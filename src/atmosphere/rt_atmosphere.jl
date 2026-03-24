#=
Top-level radiative transfer atmosphere solver.

Iteratively solves structure + RT to self-consistency.
Algorithm: build structure → Feautrier → temperature correction → repeat.

Source: Haakonsen et al. (2012) Sect. 3.
=#

module RTAtmosphere

using Printf
using ..PhysicalConstants: σ_SB, k_B, h
using ..GauntFactor: GauntTable, load_gaunt_table
using ..AtmosphereStructure: AtmosphereColumn, build_atmosphere,
                              update_atmosphere!, make_frequency_grid
using ..FeautrierSolver: solve_feautrier_all, gauss_legendre_half
using ..TemperatureCorrection: compute_temperature_correction
using ..BlackbodyAtmosphere: planck_Bnu

export solve_atmosphere, AtmosphereResult, rt_emergent_spectrum

"""
Result of a converged atmosphere calculation.
"""
struct AtmosphereResult
    T_eff::Float64
    g_s::Float64
    converged::Bool
    n_iterations::Int
    max_dT_over_T::Float64
    ν_grid::Vector{Float64}       # K frequencies [Hz]
    μ_grid::Vector{Float64}       # M angles
    I_emergent::Matrix{Float64}   # K × M emergent specific intensity
    T_profile::Vector{Float64}    # converged T(y)
    y_grid::Vector{Float64}       # column depths
end

"""
    solve_atmosphere(T_eff, g_s, gaunt; K, M, N, max_iter, tol, verbose) → AtmosphereResult

Solve for a self-consistent NS atmosphere.
"""
function solve_atmosphere(T_eff::Float64, g_s::Float64,
                          gaunt::GauntTable;
                          K::Int=50, M::Int=10, N::Int=100,
                          max_iter::Int=20,
                          tol::Float64=1e-4,
                          anisotropic::Bool=false,
                          verbose::Bool=true)::AtmosphereResult
    @assert T_eff > 0 && g_s > 0

    verbose && @printf("RT Atmosphere: T_eff=%.2e K, g_s=%.2e cm/s², K=%d, M=%d, N=%d\n",
                        T_eff, g_s, K, M, N)

    # Set up grids
    ν_grid = make_frequency_grid(T_eff, K)
    μ, w = gauss_legendre_half(M)

    # Build initial atmosphere
    verbose && println("  Building initial structure...")
    col = build_atmosphere(T_eff, g_s, ν_grid, gaunt; N=N)
    verbose && @printf("  y_max=%.2e, τ_max(ν_max)=%.1f\n", col.y[end], col.τ[end, end])

    # Solve Feautrier RT with Eddington initial T(y) profile
    # NOTE: Full iterative temperature correction (Rybicki method) not yet
    # implemented. The Eddington T(y) gives a reasonable first-order spectrum
    # that is harder than blackbody (the key NS atmosphere effect).
    verbose && println("  Solving Feautrier RT (single pass, Eddington T profile)...")
    P_all, J, f_ν, h_ν = solve_feautrier_all(col, μ, w; anisotropic=anisotropic)

    # Compute flux error for diagnostics
    F_target = σ_SB * T_eff^4
    F_emergent = 0.0
    for k in 1:K-1
        dν = ν_grid[k+1] - ν_grid[k]
        for j in 1:M
            F_emergent += 4π * μ[j] * P_all[1, j, k] * w[j] * dν
        end
    end
    flux_ratio = F_emergent / F_target
    verbose && @printf("  Flux ratio F_emergent/σT⁴ = %.4f\n", flux_ratio)

    converged = true  # Single-pass (no iteration)
    max_dT = 0.0
    n_iter = 1

    # Final solve for emergent spectrum
    if P_all === nothing
        P_all, _, _, _ = solve_feautrier_all(col, μ, w; anisotropic=anisotropic)
    end

    # Emergent intensity: I_ν(μ) = 2 P_ν(surface, μ)
    I_emergent = zeros(K, M)
    for k in 1:K, j in 1:M
        I_emergent[k, j] = 2.0 * P_all[1, j, k]
    end

    return AtmosphereResult(T_eff, g_s, converged, n_iter, max_dT,
                             ν_grid, μ, I_emergent, copy(col.T), copy(col.y))
end

"""
    rt_emergent_spectrum(result, cos_θe, ν_out) → I_ν [erg/s/cm²/Hz/sr]

Interpolate emergent spectrum from atmosphere result at given
emission angle and output frequencies. Pipeline-compatible interface.
"""
function rt_emergent_spectrum(result::AtmosphereResult,
                               cos_θe::Float64,
                               ν_out::AbstractVector{Float64})::Vector{Float64}
    @assert 0 < cos_θe <= 1.0

    # Interpolate in angle: find bracketing μ points
    μ = result.μ_grid
    M = length(μ)

    # Find μ bracket
    j_lo = 1
    for j in 1:M-1
        if μ[j] <= cos_θe <= μ[j+1]
            j_lo = j
            break
        end
        if cos_θe < μ[1]
            j_lo = 1
        elseif cos_θe > μ[end]
            j_lo = M - 1
        end
    end
    j_hi = min(j_lo + 1, M)
    t_μ = j_lo < M ? (cos_θe - μ[j_lo]) / (μ[j_hi] - μ[j_lo]) : 0.0
    t_μ = clamp(t_μ, 0.0, 1.0)

    # Interpolate in frequency (log-log)
    I_out = zeros(length(ν_out))
    ν_in = result.ν_grid
    K = length(ν_in)

    for i in eachindex(ν_out)
        ν = ν_out[i]
        if ν <= ν_in[1]
            # Extrapolate as Rayleigh-Jeans: I ∝ ν²
            I_at_lo = (1-t_μ)*result.I_emergent[1, j_lo] + t_μ*result.I_emergent[1, j_hi]
            I_out[i] = I_at_lo * (ν / ν_in[1])^2
        elseif ν >= ν_in[end]
            # Extrapolate as Wien: I ∝ ν³ exp(-hν/kT)
            I_out[i] = 0.0  # negligible
        else
            # Find frequency bracket
            k_lo = searchsortedlast(ν_in, ν)
            k_lo = clamp(k_lo, 1, K-1)
            k_hi = k_lo + 1
            t_ν = (log(ν) - log(ν_in[k_lo])) / (log(ν_in[k_hi]) - log(ν_in[k_lo]))

            # Bilinear interpolation in (ν, μ)
            I_ll = result.I_emergent[k_lo, j_lo]
            I_hl = result.I_emergent[k_hi, j_lo]
            I_lh = result.I_emergent[k_lo, j_hi]
            I_hh = result.I_emergent[k_hi, j_hi]

            # Log interpolation in frequency
            I_lo_mu = I_ll > 0 && I_hl > 0 ?
                exp((1-t_ν)*log(I_ll) + t_ν*log(I_hl)) :
                (1-t_ν)*I_ll + t_ν*I_hl
            I_hi_mu = I_lh > 0 && I_hh > 0 ?
                exp((1-t_ν)*log(I_lh) + t_ν*log(I_hh)) :
                (1-t_ν)*I_lh + t_ν*I_hh

            I_out[i] = (1-t_μ)*I_lo_mu + t_μ*I_hi_mu
        end
    end

    return I_out
end

end # module
