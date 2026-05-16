#=
Pre-computed atmosphere spectrum grid for fast rendering.

Stores emergent intensity I_ν(cos θ_e) at a grid of (T_eff, B) values.
For rendering, interpolates between grid points to get the spectrum at
each surface element's local conditions.

Phase 3 upgrade: replaces the modified-blackbody placeholder (f_col × B_ν).
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
"""
struct AtmosphereGridProvenance
    code_sha::String           # git HEAD SHA at build time (or "unknown")
    build_time::String         # epoch seconds (string(time())); avoids Dates dep
    K::Int                     # frequency grid size
    M::Int                     # angular grid size
    N::Int                     # depth grid size
    max_iter::Int
    tol_T::Float64             # ΔT/T tolerance for Rybicki
    tol_flux::Float64          # flux tolerance (magnetic only; NaN for non-magnetic)
    flux_damping::Float64      # NaN for non-magnetic
    gaunt_path::String         # path of the loaded Gaunt table (proxy for content; sha would be better)
    theta_B::Float64           # currently grids are pole-on θ_B=0; this records that
    converged::Matrix{Bool}    # per-(iT, iB) convergence flag
    iters::Matrix{Int}         # per-(iT, iB) iteration count
end

"""
Pre-computed atmosphere spectrum grid.

Stores emergent intensity I_ν(ν, cos θ_e) at a grid of (T_eff, B) values.
The ν_grid and μ_grid are common across all models.

The `provenance` field records solver knobs, code SHA, and build time so
that on-disk caches (a future feature) can be invalidated when any of
those change.
"""
struct AtmosphereSpectrumGrid
    T_grid::Vector{Float64}           # T_eff values [K]
    B_grid::Vector{Float64}           # B field values [G] (0.0 = non-magnetic)
    g_s::Float64                      # surface gravity (common)
    ν_grid::Vector{Float64}           # K frequencies [Hz]
    μ_grid::Vector{Float64}           # M angle cosines
    # I_cache[iT, iB] = K × M matrix of emergent intensity (non-magnetic)
    # or K × M × 2 array (magnetic, two modes summed)
    I_cache::Array{Matrix{Float64}, 2}  # (nT, nB) array of K×M matrices
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
    build_atmosphere_grid(T_grid, B_grid, g_s, gaunt; K=50, M=8, N=100,
                          max_iter=80, tol_T=1e-3, tol_flux=1e-2,
                          flux_damping=0.5, theta_B=0.0,
                          gaunt_path="unknown", verbose=true)

Pre-compute atmosphere models at each (T_eff, B) grid point.
Returns an AtmosphereSpectrumGrid for fast lookup during rendering.

Solver knobs (`tol_T`, `tol_flux`, `flux_damping`, `theta_B`) and the
`gaunt_path` provenance string are recorded on the returned grid via
`AtmosphereGridProvenance`. They do NOT change defaults previously used
(non-magnetic uses its own internal `tol=1e-6`; magnetic uses `tol_T`
for the temperature convergence and the solver-default `tol_flux` /
`flux_damping` unless overridden here).
"""
function build_atmosphere_grid(T_grid::Vector{Float64},
                                B_grid::Vector{Float64},
                                g_s::Float64,
                                gaunt::GauntTable;
                                K::Int=50, M::Int=8, N::Int=100,
                                max_iter::Int=80,
                                tol_T::Float64=1e-3,
                                tol_flux::Float64=1e-2,
                                flux_damping::Float64=0.5,
                                theta_B::Float64=0.0,
                                gaunt_path::String="unknown",
                                verbose::Bool=true)
    nT = length(T_grid)
    nB = length(B_grid)

    verbose && @printf("Building atmosphere grid: %d T × %d B = %d models\n", nT, nB, nT * nB)

    # Use the first non-magnetic model to establish common grids
    r_ref = solve_atmosphere(T_grid[1], g_s, gaunt; K=K, M=M, N=N,
                              max_iter=30, tol=1e-6, verbose=false)
    ν_grid = r_ref.ν_grid
    μ_grid = r_ref.μ_grid

    I_cache = Array{Matrix{Float64}}(undef, nT, nB)
    converged_grid = falses(nT, nB)
    iters_grid = zeros(Int, nT, nB)

    # Track whether any magnetic models were built; if the grid is purely
    # non-magnetic, the magnetic-only provenance fields are recorded as NaN.
    any_magnetic = false

    for (iB, B) in enumerate(B_grid)
        for (iT, T_eff) in enumerate(T_grid)
            verbose && @printf("  [%d/%d] T_eff=%.2e K, B=%.2e G ... ",
                                (iB-1)*nT + iT, nT*nB, T_eff, B)

            if B == 0.0
                # Non-magnetic atmosphere
                r = solve_atmosphere(T_eff, g_s, gaunt; K=K, M=M, N=N,
                                      max_iter=30, tol=1e-6, verbose=false)
                I_cache[iT, iB] = copy(r.I_emergent)  # K × M
                converged_grid[iT, iB] = r.converged
                # `AtmosphereResult` does not expose an explicit iteration count;
                # record max_iter=30 as the budget actually used here.
                iters_grid[iT, iB] = 30
                verbose && @printf("converged=%s, F/σT⁴=%.3f\n", r.converged,
                    _flux_ratio(r.I_emergent, μ_grid, ν_grid, T_eff))
            else
                # Magnetic atmosphere (sum X + O modes)
                any_magnetic = true
                r = solve_magnetic_atmosphere(T_eff, g_s, B, theta_B, gaunt;
                        K=K, M=M, N=N, max_iter=max_iter, tol=tol_T,
                        flux_tol=tol_flux, flux_damping=flux_damping,
                        verbose=false)
                # Sum both modes: I_total[k, j] = I_X[k, j] + I_O[k, j]
                I_total = r.I_emergent[:, :, 1] .+ r.I_emergent[:, :, 2]
                I_cache[iT, iB] = I_total  # K × M
                converged_grid[iT, iB] = r.converged
                iters_grid[iT, iB] = r.n_iterations
                verbose && @printf("converged=%s, iters=%d\n", r.converged, r.n_iterations)
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
        theta_B,
        converged_grid,
        iters_grid,
    )

    return AtmosphereSpectrumGrid(copy(T_grid), copy(B_grid), g_s,
                                   copy(ν_grid), copy(μ_grid), I_cache, prov)
end

"""
    lookup_spectrum(grid, T_eff, B, ν_obs, cos_θe) → I_ν [erg/s/cm²/Hz/sr]

Look up the emergent specific intensity at given conditions by bilinear
interpolation in (T_eff, B) and linear interpolation in angle.

For T/B outside the grid, uses nearest-neighbour (clamping).
"""
function lookup_spectrum(grid::AtmosphereSpectrumGrid,
                          T_eff::Float64, B::Float64,
                          ν_obs::Vector{Float64}, cos_θe::Float64)
    nν = length(ν_obs)
    I_out = zeros(nν)

    # Find bracketing indices for T
    iT, fT = _bracket_interp(grid.T_grid, T_eff)
    # Find bracketing indices for B
    iB, fB = _bracket_interp(grid.B_grid, B)

    # Find bracketing angle index
    iμ, fμ = _bracket_interp(grid.μ_grid, cos_θe)

    # Bilinear interpolation in (T, B), linear in angle
    # 4 corner models
    for (wT, jT) in ((1-fT, iT), (fT, min(iT+1, length(grid.T_grid))))
        for (wB, jB) in ((1-fB, iB), (fB, min(iB+1, length(grid.B_grid))))
            w = wT * wB
            w < 1e-10 && continue

            I_model = grid.I_cache[jT, jB]  # K × M

            # Interpolate in angle for each observed frequency
            for iν in 1:nν
                # Find nearest frequency in grid
                kν = _nearest(grid.ν_grid, ν_obs[iν])

                # Angle interpolation
                I_at_ν = (1-fμ) * I_model[kν, iμ] + fμ * I_model[kν, min(iμ+1, length(grid.μ_grid))]
                I_out[iν] += w * max(I_at_ν, 0.0)
            end
        end
    end

    return I_out
end

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
    if x <= grid[1]
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
