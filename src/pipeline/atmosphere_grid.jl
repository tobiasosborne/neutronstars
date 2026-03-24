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

export AtmosphereSpectrumGrid, build_atmosphere_grid, lookup_spectrum

"""
Pre-computed atmosphere spectrum grid.

Stores emergent intensity I_ν(ν, cos θ_e) at a grid of (T_eff, B) values.
The ν_grid and μ_grid are common across all models.
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
end

"""
    build_atmosphere_grid(T_grid, B_grid, g_s, gaunt; K=50, M=8, N=100, verbose=true)

Pre-compute atmosphere models at each (T_eff, B) grid point.
Returns an AtmosphereSpectrumGrid for fast lookup during rendering.
"""
function build_atmosphere_grid(T_grid::Vector{Float64},
                                B_grid::Vector{Float64},
                                g_s::Float64,
                                gaunt::GauntTable;
                                K::Int=50, M::Int=8, N::Int=100,
                                max_iter::Int=80,
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

    for (iB, B) in enumerate(B_grid)
        for (iT, T_eff) in enumerate(T_grid)
            verbose && @printf("  [%d/%d] T_eff=%.2e K, B=%.2e G ... ",
                                (iB-1)*nT + iT, nT*nB, T_eff, B)

            if B == 0.0
                # Non-magnetic atmosphere
                r = solve_atmosphere(T_eff, g_s, gaunt; K=K, M=M, N=N,
                                      max_iter=30, tol=1e-6, verbose=false)
                I_cache[iT, iB] = copy(r.I_emergent)  # K × M
                verbose && @printf("converged=%s, F/σT⁴=%.3f\n", r.converged,
                    _flux_ratio(r.I_emergent, μ_grid, ν_grid, T_eff))
            else
                # Magnetic atmosphere (sum X + O modes)
                r = solve_magnetic_atmosphere(T_eff, g_s, B, 0.0, gaunt;
                        K=K, M=M, N=N, max_iter=max_iter, tol=1e-3, verbose=false)
                # Sum both modes: I_total[k, j] = I_X[k, j] + I_O[k, j]
                I_total = r.I_emergent[:, :, 1] .+ r.I_emergent[:, :, 2]
                I_cache[iT, iB] = I_total  # K × M
                verbose && @printf("converged=%s, iters=%d\n", r.converged, r.n_iterations)
            end
        end
    end

    return AtmosphereSpectrumGrid(copy(T_grid), copy(B_grid), g_s,
                                   copy(ν_grid), copy(μ_grid), I_cache)
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

"""
    lookup_spectrum_scalar(grid, T_eff, B, ν, cos_θe) → I_ν

Single-frequency version for per-ray evaluation.
"""
function lookup_spectrum_scalar(grid::AtmosphereSpectrumGrid,
                                 T_eff::Float64, B::Float64,
                                 ν::Float64, cos_θe::Float64)::Float64
    iT, fT = _bracket_interp(grid.T_grid, T_eff)
    iB, fB = _bracket_interp(grid.B_grid, B)
    iμ, fμ = _bracket_interp(grid.μ_grid, cos_θe)
    kν = _nearest(grid.ν_grid, ν)

    I = 0.0
    for (wT, jT) in ((1-fT, iT), (fT, min(iT+1, length(grid.T_grid))))
        for (wB, jB) in ((1-fB, iB), (fB, min(iB+1, length(grid.B_grid))))
            w = wT * wB
            w < 1e-10 && continue
            I_model = grid.I_cache[jT, jB]
            I_at_ν = (1-fμ) * I_model[kν, iμ] + fμ * I_model[kν, min(iμ+1, length(grid.μ_grid))]
            I += w * max(I_at_ν, 0.0)
        end
    end
    return I
end

# --- Internal helpers ---

"""Compute F/σT⁴ from emergent intensity."""
function _flux_ratio(I_em, μ, ν_grid, T_eff)
    F = 0.0
    for k in 1:length(ν_grid)-1
        dν = ν_grid[k+1] - ν_grid[k]
        for j in 1:length(μ)
            F += 2π * μ[j] * I_em[k, j] * (μ[j] > 0 ? 2.0 : 0.0) * dν  # approximate weight
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
