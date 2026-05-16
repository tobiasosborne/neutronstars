#=
Pre-computed atmosphere spectrum grid for fast rendering.

Stores emergent intensity I_ν(cos θ_e) at a grid of (T_eff, B, θ_B) values.
For rendering, interpolates between grid points to get the spectrum at
each surface element's local conditions.

Phase 3 upgrade: replaces the modified-blackbody placeholder (f_col × B_ν).
Phase 3d (bead D1): adds the θ_B axis. Previously the grid hardcoded
θ_B = 0 (pole-on opacities) for every (T,B) cell, so every magnetic
render silently used pole-on radiative transfer regardless of which
surface element was being rendered — the dominant variable for the
X/O-mode opacity ratio was missing. The grid now has a third axis
θ_B_grid (e.g. {0, π/4, π/2}) and `lookup_spectrum` takes a per-pixel
θ_B argument. A backward-compat 2-axis overload (no θ_B) builds a
single-θ_B=0 grid and `lookup_spectrum(grid, T, B, ν, μ)` (no θ_B)
still works on those grids, returning the θ_B=0 spectrum.
=#

module AtmosphereGrid

using Printf
using ..PhysicalConstants: h, k_B, σ_SB
using ..GauntFactor: GauntTable
using ..RTAtmosphere: solve_atmosphere, AtmosphereResult
using ..MagneticAtmosphere: solve_magnetic_atmosphere, MagneticAtmosphereResult
using ..BlackbodyAtmosphere: planck_Bnu
using ..FeautrierSolver: gauss_legendre_half

export AtmosphereSpectrumGrid, AtmosphereGridProvenance, build_atmosphere_grid, lookup_spectrum

"""
    AtmosphereGridProvenance

Metadata describing how an `AtmosphereSpectrumGrid` was built. Critical
for catching stale caches when HDF5-on-disk caching ships — without
this, a grid saved with one set of solver knobs cannot be distinguished
from a grid saved with another.

`build_time` is the epoch-seconds Unix timestamp as a string (from
`string(time())`); we avoid `Dates` to keep the dependency surface
small. Convert with `Dates.unix2datetime(parse(Float64, s))` if needed.
`code_sha` is the `git rev-parse HEAD` of the source tree at build
time, or `"unknown"` if git is unavailable.
For non-magnetic grids (`B == 0` everywhere) the magnetic-only fields
`tol_flux` and `flux_damping` are recorded as `NaN`.

`theta_B_grid` records the actual θ_B sample points the grid was built
on (post-D1; previously a single scalar). A single-element vector means
a pole-on grid (back-compat path).
"""
struct AtmosphereGridProvenance
    code_sha::String             # git HEAD SHA at build time (or "unknown")
    build_time::String           # epoch seconds (string(time())); avoids Dates dep
    K::Int                       # frequency grid size
    M::Int                       # angular grid size
    N::Int                       # depth grid size
    max_iter::Int
    tol_T::Float64               # ΔT/T tolerance for Rybicki
    tol_flux::Float64            # flux tolerance (magnetic only; NaN for non-magnetic)
    flux_damping::Float64        # NaN for non-magnetic
    gaunt_path::String           # path of the loaded Gaunt table (proxy for content; sha would be better)
    theta_B_grid::Vector{Float64}# axis of B-normal angles [rad]
    converged::Array{Bool, 3}    # per-(iT, iB, iθB) convergence flag
    iters::Array{Int, 3}         # per-(iT, iB, iθB) iteration count
end

"""
Pre-computed atmosphere spectrum grid.

Stores emergent intensity I_ν(ν, cos θ_e) at a grid of (T_eff, B, θ_B) values.
The ν_grid and μ_grid are common across all models.

The `provenance` field records solver knobs, code SHA, and build time so
that on-disk caches (a future feature) can be invalidated when any of
those change.
"""
struct AtmosphereSpectrumGrid
    T_grid::Vector{Float64}              # T_eff values [K]
    B_grid::Vector{Float64}              # B field values [G] (0.0 = non-magnetic)
    θ_B_grid::Vector{Float64}            # B-normal angle values [rad] (length≥1)
    g_s::Float64                         # surface gravity (common)
    ν_grid::Vector{Float64}              # K frequencies [Hz]
    μ_grid::Vector{Float64}              # M angle cosines
    # I_cache[iT, iB, iθB] = K × M matrix of emergent intensity (non-magnetic)
    # or K × M matrix of summed X+O modes (magnetic).
    I_cache::Array{Matrix{Float64}, 3}   # (nT, nB, nθB) array of K×M matrices
    provenance::AtmosphereGridProvenance
end

"""
Return the current git HEAD SHA for the source tree, or `"unknown"` if
git is not available (e.g. when running from a tarball or in a CI
sandbox without git in `PATH`).
"""
function _git_sha()
    try
        return readchomp(`git -C $(@__DIR__) rev-parse HEAD`)
    catch
        return "unknown"
    end
end

"""
    build_atmosphere_grid(T_grid, B_grid, θ_B_grid, g_s, gaunt; K=50, M=8, N=100,
                          max_iter=80, tol_T=1e-3, tol_flux=1e-2,
                          flux_damping=0.5, gaunt_path="unknown", verbose=true)

Pre-compute atmosphere models at each (T_eff, B, θ_B) grid point.
Returns an `AtmosphereSpectrumGrid` for fast lookup during rendering.

Solver knobs (`tol_T`, `tol_flux`, `flux_damping`) and the
`gaunt_path` provenance string are recorded on the returned grid via
`AtmosphereGridProvenance`. They do NOT change defaults previously used
(non-magnetic uses its own internal `tol=1e-6`; magnetic uses `tol_T`
for the temperature convergence and the solver-default `tol_flux` /
`flux_damping` unless overridden here).

For `B == 0` cells the result is θ_B-independent, so the solver is
called only once per (T, B=0) pair and the result is broadcast across
the θ_B axis.

A backward-compat 2-axis overload `build_atmosphere_grid(T_grid, B_grid,
g_s, gaunt; ...)` is provided that auto-uses `θ_B_grid = [0.0]`.
"""
function build_atmosphere_grid(T_grid::Vector{Float64},
                                B_grid::Vector{Float64},
                                θ_B_grid::Vector{Float64},
                                g_s::Float64,
                                gaunt::GauntTable;
                                K::Int=50, M::Int=8, N::Int=100,
                                max_iter::Int=80,
                                tol_T::Float64=1e-3,
                                tol_flux::Float64=1e-2,
                                flux_damping::Float64=0.5,
                                gaunt_path::String="unknown",
                                verbose::Bool=true)
    @assert !isempty(θ_B_grid) "θ_B_grid must contain at least one angle"
    @assert all(0 .<= θ_B_grid .<= π) "θ_B angles must lie in [0, π]"
    @assert issorted(θ_B_grid) "θ_B_grid must be sorted ascending for interpolation"

    nT = length(T_grid)
    nB = length(B_grid)
    nθB = length(θ_B_grid)

    verbose && @printf("Building atmosphere grid: %d T × %d B × %d θ_B = %d models\n",
                        nT, nB, nθB, nT * nB * nθB)

    # Use the first non-magnetic model to establish common grids
    r_ref = solve_atmosphere(T_grid[1], g_s, gaunt; K=K, M=M, N=N,
                              max_iter=30, tol=1e-6, verbose=false)
    ν_grid = r_ref.ν_grid
    μ_grid = r_ref.μ_grid

    I_cache = Array{Matrix{Float64}, 3}(undef, nT, nB, nθB)
    converged_grid = falses(nT, nB, nθB)
    iters_grid = zeros(Int, nT, nB, nθB)

    # Track whether any magnetic models were built; if the grid is purely
    # non-magnetic, the magnetic-only provenance fields are recorded as NaN.
    any_magnetic = false

    total = nT * nB * nθB
    step = 0
    for (iB, B) in enumerate(B_grid)
        for (iT, T_eff) in enumerate(T_grid)
            if B == 0.0
                # Non-magnetic atmosphere: θ_B-independent. Solve once
                # and broadcast across the θ_B axis.
                step += nθB
                verbose && @printf("  [%d/%d] T_eff=%.2e K, B=0 (θ_B-indep.) ... ",
                                    step, total, T_eff)
                r = solve_atmosphere(T_eff, g_s, gaunt; K=K, M=M, N=N,
                                      max_iter=30, tol=1e-6, verbose=false)
                I_one = copy(r.I_emergent)  # K × M
                for iθB in 1:nθB
                    # Each θ_B cell gets its own copy so callers cannot
                    # accidentally mutate the cache through aliasing.
                    I_cache[iT, iB, iθB] = copy(I_one)
                    converged_grid[iT, iB, iθB] = r.converged
                    # `AtmosphereResult` does not expose an explicit iteration count;
                    # record max_iter=30 as the budget actually used here.
                    iters_grid[iT, iB, iθB] = 30
                end
                verbose && @printf("converged=%s, F/σT⁴=%.3f\n", r.converged,
                    _flux_ratio(r.I_emergent, μ_grid, ν_grid, T_eff))
            else
                # Magnetic atmosphere (sum X + O modes), one solve per θ_B.
                any_magnetic = true
                for (iθB, θ_B) in enumerate(θ_B_grid)
                    step += 1
                    verbose && @printf("  [%d/%d] T_eff=%.2e K, B=%.2e G, θ_B=%.1f° ... ",
                                        step, total, T_eff, B, rad2deg(θ_B))
                    r = solve_magnetic_atmosphere(T_eff, g_s, B, θ_B, gaunt;
                            K=K, M=M, N=N, max_iter=max_iter, tol=tol_T,
                            flux_tol=tol_flux, flux_damping=flux_damping,
                            verbose=false)
                    # Sum both modes: I_total[k, j] = I_X[k, j] + I_O[k, j]
                    I_total = r.I_emergent[:, :, 1] .+ r.I_emergent[:, :, 2]
                    I_cache[iT, iB, iθB] = I_total  # K × M
                    converged_grid[iT, iB, iθB] = r.converged
                    iters_grid[iT, iB, iθB] = r.n_iterations
                    verbose && @printf("converged=%s, iters=%d\n", r.converged, r.n_iterations)
                end
            end
        end
    end

    prov = AtmosphereGridProvenance(
        _git_sha(),
        string(time()),
        K, M, N,
        max_iter,
        tol_T,
        any_magnetic ? tol_flux : NaN,
        any_magnetic ? flux_damping : NaN,
        gaunt_path,
        copy(θ_B_grid),
        converged_grid,
        iters_grid,
    )

    return AtmosphereSpectrumGrid(copy(T_grid), copy(B_grid), copy(θ_B_grid), g_s,
                                   copy(ν_grid), copy(μ_grid), I_cache, prov)
end

"""
    build_atmosphere_grid(T_grid, B_grid, g_s, gaunt; kwargs...)

Backward-compat 2-axis overload: builds a single-θ_B grid at θ_B=0.0.
Existing callers that pass only `(T_grid, B_grid)` continue to work but
will silently get pole-on opacities for every magnetic cell. Pass an
explicit `θ_B_grid` to take advantage of the new axis.
"""
build_atmosphere_grid(T_grid::Vector{Float64},
                       B_grid::Vector{Float64},
                       g_s::Float64,
                       gaunt::GauntTable;
                       kwargs...) =
    build_atmosphere_grid(T_grid, B_grid, [0.0], g_s, gaunt; kwargs...)

"""
    lookup_spectrum(grid, T_eff, B, θ_B, ν_obs, cos_θe) → I_ν [erg/s/cm²/Hz/sr]

Look up the emergent specific intensity at given conditions by
interpolation in (T_eff, B, θ_B), and linear interpolation in angle and
nearest-neighbour in frequency (matching the original Phase 3 lookup
fidelity).

For T/B/θ_B outside the grid, uses nearest-neighbour (clamping).
When the grid has a single θ_B sample, the θ_B interpolation collapses
to that one cell (back-compat path).
"""
function lookup_spectrum(grid::AtmosphereSpectrumGrid,
                          T_eff::Float64, B::Float64, θ_B::Float64,
                          ν_obs::Vector{Float64}, cos_θe::Float64)
    nν = length(ν_obs)
    I_out = zeros(nν)

    # Find bracketing indices for T, B, θ_B
    iT, fT = _bracket_interp(grid.T_grid, T_eff)
    iB, fB = _bracket_interp(grid.B_grid, B)
    iθB, fθB = _bracket_interp(grid.θ_B_grid, θ_B)

    # Find bracketing angle index
    iμ, fμ = _bracket_interp(grid.μ_grid, cos_θe)

    nT = length(grid.T_grid)
    nB = length(grid.B_grid)
    nθB = length(grid.θ_B_grid)
    nμ = length(grid.μ_grid)

    # Trilinear interpolation in (T, B, θ_B); linear in angle; nearest in ν.
    # 8 corner models in the (T, B, θ_B) cube.
    for (wT, jT) in ((1-fT, iT), (fT, min(iT+1, nT)))
        for (wB, jB) in ((1-fB, iB), (fB, min(iB+1, nB)))
            for (wθ, jθ) in ((1-fθB, iθB), (fθB, min(iθB+1, nθB)))
                w = wT * wB * wθ
                w < 1e-10 && continue

                I_model = grid.I_cache[jT, jB, jθ]  # K × M

                # Interpolate in angle for each observed frequency
                for iν in 1:nν
                    # Find nearest frequency in grid
                    kν = _nearest(grid.ν_grid, ν_obs[iν])

                    # Angle interpolation
                    I_at_ν = (1-fμ) * I_model[kν, iμ] + fμ * I_model[kν, min(iμ+1, nμ)]
                    I_out[iν] += w * max(I_at_ν, 0.0)
                end
            end
        end
    end

    return I_out
end

"""
    lookup_spectrum(grid, T_eff, B, ν_obs, cos_θe) → I_ν

Backward-compat 4-arg overload (no θ_B). Uses θ_B = 0.0; on a
single-θ_B=0 grid (the back-compat build path) this returns exactly the
same answer the pre-D1 code would have. On a multi-θ_B grid it returns
the pole-on cell, which is almost certainly NOT what the caller wants —
prefer the 5-arg form.
"""
lookup_spectrum(grid::AtmosphereSpectrumGrid,
                T_eff::Float64, B::Float64,
                ν_obs::Vector{Float64}, cos_θe::Float64) =
    lookup_spectrum(grid, T_eff, B, 0.0, ν_obs, cos_θe)

# --- Internal helpers ---

"""
Compute F/σT⁴ from emergent specific intensity.

The hemispheric bolometric flux is
    F = ∫₀^∞ dν · 2π ∫₀¹ I_ν(μ) μ dμ
which we evaluate with the same Gauss-Legendre quadrature used by the
Feautrier solver. `gauss_legendre_half(M)` returns nodes/weights on the
upward hemisphere μ ∈ [0,1] (so the 2π · ∫₀¹ already covers the half-sphere
of emergent rays; no extra factor of 2 is needed). Frequency integration is
a simple trapezoid-style left-Riemann sum on the existing ν_grid.
"""
function _flux_ratio(I_em, μ, ν_grid, T_eff)
    _, w = gauss_legendre_half(length(μ))
    F = 0.0
    for k in 1:length(ν_grid)-1
        dν = ν_grid[k+1] - ν_grid[k]
        for j in 1:length(μ)
            F += 2π * μ[j] * w[j] * I_em[k, j] * dν
        end
    end
    return F / (σ_SB * T_eff^4)
end

"""Find bracketing index and interpolation fraction in a sorted grid."""
function _bracket_interp(grid::Vector{Float64}, x::Float64)
    n = length(grid)
    if n == 1 || x <= grid[1]
        return 1, 0.0
    elseif x >= grid[end]
        return n, 0.0
    end
    i = searchsortedlast(grid, x)
    i = clamp(i, 1, n-1)
    f = (x - grid[i]) / (grid[i+1] - grid[i])
    return i, f
end

"""Find nearest index in grid."""
function _nearest(grid::Vector{Float64}, x::Float64)::Int
    return argmin(abs.(grid .- x))
end

end # module
