#=
Top-level radiative transfer atmosphere solver.

Iteratively solves structure + RT to self-consistency:
  1. Build initial atmosphere (Eddington T profile)
  2. Solve Feautrier RT at all frequencies
  3. Compute Rybicki temperature correction ΔT
  4. Apply correction, update opacities and structure
  5. Check convergence: max|ΔT/T| < tol
  6. Repeat until converged

Source: Haakonsen et al. (2012) Sect. 3.
=#

module RTAtmosphere

using Printf
using ..PhysicalConstants: σ_SB, k_B, h, m_p
using ..GauntFactor: GauntTable, load_gaunt_table
using ..AtmosphereStructure: AtmosphereColumn, build_atmosphere,
                              update_atmosphere!, make_frequency_grid
using ..FeautrierSolver: solve_feautrier_all, gauss_legendre_half
using ..TemperatureCorrection: compute_temperature_correction
using ..BlackbodyAtmosphere: planck_Bnu
using ..HydrogenOpacity: kappa_ff, sigma_thomson, total_opacity

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

Solve for a self-consistent NS atmosphere using iterative Rybicki
temperature correction.
"""
function solve_atmosphere(T_eff::Float64, g_s::Float64,
                          gaunt::GauntTable;
                          K::Int=50, M::Int=10, N::Int=100,
                          max_iter::Int=30,
                          tol::Float64=1e-4,
                          anisotropic::Bool=true,
                          verbose::Bool=true)::AtmosphereResult
    @assert T_eff > 0 && g_s > 0

    verbose && @printf("RT Atmosphere: T_eff=%.2e K, g_s=%.2e cm/s², K=%d, M=%d, N=%d\n",
                        T_eff, g_s, K, M, N)

    # Set up grids
    ν_grid = make_frequency_grid(T_eff, K)
    μ, w = gauss_legendre_half(M)

    # Build initial atmosphere with Eddington T(y)
    verbose && println("  Building initial structure...")
    col = build_atmosphere(T_eff, g_s, ν_grid, gaunt; N=N)
    verbose && @printf("  y_max=%.2e, τ_max(ν_max)=%.1f, N=%d depths\n",
                        col.y[end], maximum(col.τ[end, :]), col.N)

    converged = false
    max_dT = Inf
    n_iter = 0

    for iter in 1:max_iter
        n_iter = iter

        # Solve Feautrier RT at all frequencies
        P_all, J, f_ν, h_ν = solve_feautrier_all(col, μ, w; anisotropic=anisotropic)

        # Compute flux for diagnostics
        F_target = σ_SB * T_eff^4
        F_bol = _bolometric_flux(P_all, μ, w, ν_grid)
        flux_ratio = F_bol / F_target

        # Compute Rybicki temperature correction
        ΔT = compute_temperature_correction(col, f_ν, h_ν, J)

        # Dampen large corrections to aid convergence
        max_dT = maximum(abs.(ΔT ./ col.T))
        if max_dT > 0.3
            damp = 0.3 / max_dT
            ΔT .*= damp
            max_dT = 0.3
        end

        verbose && @printf("  iter %2d: max|ΔT/T|=%.2e, F/σT⁴=%.4f\n",
                            iter, max_dT, flux_ratio)

        # Check convergence
        if max_dT < tol
            converged = true
            verbose && @printf("  CONVERGED at iteration %d (max|ΔT/T|=%.2e < %.2e)\n",
                                iter, max_dT, tol)
            break
        end

        # Apply correction
        for i in 1:col.N
            col.T[i] = max(col.T[i] + ΔT[i], 0.1 * T_eff)  # floor at 10% of T_eff
        end

        # Update atmospheric structure (density, opacities, optical depths)
        _update_structure!(col, gaunt)
    end

    if !converged
        verbose && @printf("  WARNING: not converged after %d iterations (max|ΔT/T|=%.2e)\n",
                            max_iter, max_dT)
    end

    # Final Feautrier solve for emergent spectrum
    P_all, J, f_ν, h_ν = solve_feautrier_all(col, μ, w; anisotropic=anisotropic)

    # Emergent intensity: I_ν(μ) = 2 P_ν(surface, μ)
    I_emergent = zeros(K, M)
    for k in 1:K, j in 1:M
        I_emergent[k, j] = 2.0 * P_all[1, j, k]
    end

    # Final flux diagnostic
    F_bol = _bolometric_flux(P_all, μ, w, ν_grid)
    verbose && @printf("  Final F/σT⁴=%.4f\n", F_bol / (σ_SB * T_eff^4))

    return AtmosphereResult(T_eff, g_s, converged, n_iter, max_dT,
                             ν_grid, μ, I_emergent, copy(col.T), copy(col.y))
end

"""
Compute bolometric emergent flux from the Feautrier solution.
"""
function _bolometric_flux(P_all, μ, w, ν_grid)
    K = length(ν_grid)
    M = length(μ)
    F = 0.0
    for k in 1:K-1
        dν = ν_grid[k+1] - ν_grid[k]
        for j in 1:M
            # F_ν = 2π ∫₀¹ I_ν(μ) μ dμ where I = 2P at surface
            F += 2π * μ[j] * 2.0 * P_all[1, j, k] * w[j] * dν
        end
    end
    return F
end

"""
Update atmospheric structure after a temperature correction.
Recomputes density, opacities, and optical depths from the new T profile.
"""
function _update_structure!(col::AtmosphereColumn, gaunt::GauntTable)
    for i in 1:col.N
        # Ideal gas for ionised H: ρ = m_p P/(2 k_B T), μ=0.5
        col.ρ[i] = m_p * col.P[i] / (2.0 * k_B * col.T[i])

        # Recompute opacities at each frequency
        for k in 1:col.K
            ν = col.ν_grid[k]
            col.κ[i, k] = kappa_ff(ν, col.T[i], col.ρ[i], gaunt)
            col.k_total[i, k] = col.κ[i, k] + col.σ_scat
            col.ρ_alb[i, k] = col.σ_scat / col.k_total[i, k]
        end
    end

    # Recompute optical depths from surface
    for k in 1:col.K
        col.τ[1, k] = 0.0
        for i in 2:col.N
            dy = col.y[i] - col.y[i-1]
            col.τ[i, k] = col.τ[i-1, k] + 0.5 * (col.k_total[i, k] + col.k_total[i-1, k]) * dy
        end
    end
end

"""
    rt_emergent_spectrum(result, cos_θe, ν_out) → I_ν [erg/s/cm²/Hz/sr]

Interpolate emergent spectrum at given emission angle and frequencies.
"""
function rt_emergent_spectrum(result::AtmosphereResult,
                               cos_θe::Float64,
                               ν_out::AbstractVector{Float64})::Vector{Float64}
    @assert 0 < cos_θe <= 1.0

    μ = result.μ_grid
    M = length(μ)

    # Find μ bracket
    j_lo = 1
    for j in 1:M-1
        if μ[j] <= cos_θe <= μ[j+1]
            j_lo = j
            break
        end
    end
    if cos_θe < μ[1]; j_lo = 1; end
    if cos_θe > μ[end]; j_lo = M - 1; end
    j_hi = min(j_lo + 1, M)
    t_μ = j_lo < M ? clamp((cos_θe - μ[j_lo]) / (μ[j_hi] - μ[j_lo]), 0.0, 1.0) : 0.0

    # Interpolate in frequency (log-log)
    ν_in = result.ν_grid
    K = length(ν_in)
    I_out = zeros(length(ν_out))

    for i in eachindex(ν_out)
        ν = ν_out[i]
        if ν <= ν_in[1]
            I_at = (1-t_μ)*result.I_emergent[1, j_lo] + t_μ*result.I_emergent[1, j_hi]
            I_out[i] = I_at * (ν / ν_in[1])^2  # Rayleigh-Jeans
        elseif ν >= ν_in[end]
            I_out[i] = 0.0
        else
            k_lo = searchsortedlast(ν_in, ν)
            k_lo = clamp(k_lo, 1, K-1)
            k_hi = k_lo + 1
            t_ν = (log(ν) - log(ν_in[k_lo])) / (log(ν_in[k_hi]) - log(ν_in[k_lo]))

            I_ll = result.I_emergent[k_lo, j_lo]
            I_hl = result.I_emergent[k_hi, j_lo]
            I_lh = result.I_emergent[k_lo, j_hi]
            I_hh = result.I_emergent[k_hi, j_hi]

            I_lo = I_ll > 0 && I_hl > 0 ? exp((1-t_ν)*log(I_ll) + t_ν*log(I_hl)) : (1-t_ν)*I_ll + t_ν*I_hl
            I_hi = I_lh > 0 && I_hh > 0 ? exp((1-t_ν)*log(I_lh) + t_ν*log(I_hh)) : (1-t_ν)*I_lh + t_ν*I_hh
            I_out[i] = (1-t_μ)*I_lo + t_μ*I_hi
        end
    end

    return I_out
end

end # module
