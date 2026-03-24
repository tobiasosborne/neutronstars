#=
Normal mode decomposition and angle-averaged opacities in a magnetic field.

Computes the polarisation-dependent opacity at each angle θ_B, then
performs the angle and frequency averaging to get Rosseland mean
opacities K_∥ and K_⊥.

Key fix vs previous version:
1. Polarization weights from full cold-plasma dielectric tensor
   (DielectricTensor module) instead of ad-hoc Eq. 26 approximation.
2. Angle-average each mode SEPARATELY (P&C 2003 Eq. 30), then combine
   modes via harmonic mean. Previous code combined modes first (wrong).

Source: Potekhin & Chabrier (2003) ApJ 585, 955, Eqs. 25-30.
=#

module MagneticModes

using ..PhysicalConstants: e_charge, m_e, m_p, c, h, ħ, k_B, σ_T
using ..MagneticFF: cyclotron_freq_e, cyclotron_freq_p,
                     sigma_total_alpha, sigma_ff_alpha, sigma_pp_alpha,
                     sigma_scat_alpha
using ..DielectricTensor: polarization_weights_full
using ..HydrogenOpacity: dBnu_dT
using ..BlackbodyAtmosphere: planck_Bnu

export mode_opacity, mode_absorption, mode_scattering, effective_opacity
export kappa_parallel_mono, kappa_perp_mono
export rosseland_magnetic, make_rosseland_frequency_grid

const m_H = m_p + m_e

"""
    mode_opacity(j, ν, θ_B, B, T, ρ) → κ_j [cm²/g]

Opacity for normal mode j (1=extraordinary, 2=ordinary) at angle θ_B.
P&C 2003 Eq. (27): κ_j = (1/m_H) Σ_α |e_{j,α}|² σ_α
"""
function mode_opacity(j::Int, ν::Float64, θ_B::Float64,
                      B::Float64, T::Float64, ρ::Float64)::Float64
    @assert j ∈ (1, 2)
    ω = 2π * ν
    n_e = ρ / m_H

    w1, w2 = polarization_weights_full(ω, B, θ_B, n_e)
    w = j == 1 ? w1 : w2

    κ = 0.0
    for (idx, α) in enumerate((-1, 0, 1))
        σ = sigma_total_alpha(α, ω, B, T, ρ)
        κ += w[idx] * σ
    end

    return κ / m_H
end

"""
    mode_absorption(j, ν, θ_B, B, T, ρ) → κ_j^abs [cm²/g]

Absorption-only opacity for normal mode j (free-free + proton-proton).
Same polarization weight decomposition as mode_opacity, but excludes scattering.
"""
function mode_absorption(j::Int, ν::Float64, θ_B::Float64,
                          B::Float64, T::Float64, ρ::Float64)::Float64
    @assert j ∈ (1, 2)
    ω = 2π * ν
    n_e = ρ / m_H

    w1, w2 = polarization_weights_full(ω, B, θ_B, n_e)
    w = j == 1 ? w1 : w2

    κ = 0.0
    for (idx, α) in enumerate((-1, 0, 1))
        σ = sigma_ff_alpha(α, ω, B, T, ρ) + sigma_pp_alpha(α, ω, B, T, ρ)
        κ += w[idx] * σ
    end

    return κ / m_H
end

"""
    mode_scattering(j, ν, θ_B, B, T, ρ) → σ_j^scat [cm²/g]

Scattering-only opacity for normal mode j (magnetic Thomson).
Same polarization weight decomposition as mode_opacity, but only scattering.
"""
function mode_scattering(j::Int, ν::Float64, θ_B::Float64,
                          B::Float64, T::Float64, ρ::Float64)::Float64
    @assert j ∈ (1, 2)
    ω = 2π * ν
    n_e = ρ / m_H

    w1, w2 = polarization_weights_full(ω, B, θ_B, n_e)
    w = j == 1 ? w1 : w2

    σ_scat = 0.0
    for (idx, α) in enumerate((-1, 0, 1))
        σ_scat += w[idx] * sigma_scat_alpha(α, ω, B)
    end

    return σ_scat / m_H
end

"""
    effective_opacity(ν, θ_B, B, T, ρ) → κ_eff [cm²/g]

Effective non-polarised opacity (harmonic mean of two modes) at a single angle.
"""
function effective_opacity(ν::Float64, θ_B::Float64,
                           B::Float64, T::Float64, ρ::Float64)::Float64
    κ1 = mode_opacity(1, ν, θ_B, B, T, ρ)
    κ2 = mode_opacity(2, ν, θ_B, B, T, ρ)

    if κ1 <= 0 || κ2 <= 0
        return max(κ1, κ2)
    end
    return 2.0 / (1.0/κ1 + 1.0/κ2)
end

"""
    kappa_parallel_perp_mono(ν, B, T, ρ; N_θ=40) → (κ_∥, κ_⊥) [cm²/g]

Monochromatic angle-averaged opacities for transport along and across B.

CORRECT implementation: angle-average each mode j separately (P&C 2003 Eq. 30),
then combine modes via harmonic mean.

  1/κ_j^∥ = (3/4) ∫₀^π cos²θ sinθ / κ_j(θ) dθ
  1/κ_j^⊥ = (3/2) ∫₀^π sin³θ / κ_j(θ) dθ
  κ_eff^∥ = 2/(1/κ_1^∥ + 1/κ_2^∥)    [harmonic mean over modes]
"""
function kappa_parallel_perp_mono(ν::Float64, B::Float64, T::Float64,
                                   ρ::Float64; N_θ::Int=40)
    dθ = π / N_θ

    # Accumulate 1/κ integrals for each mode separately
    inv_κ1_par = 0.0
    inv_κ2_par = 0.0
    inv_κ1_perp = 0.0
    inv_κ2_perp = 0.0

    for i in 1:N_θ
        θ = (i - 0.5) * dθ
        sinθ = sin(θ)
        cosθ = cos(θ)

        κ1 = mode_opacity(1, ν, θ, B, T, ρ)
        κ2 = mode_opacity(2, ν, θ, B, T, ρ)

        weight_par = cosθ^2 * sinθ * dθ
        weight_perp = sinθ^3 * dθ

        if κ1 > 0
            inv_κ1_par += weight_par / κ1
            inv_κ1_perp += weight_perp / κ1
        end
        if κ2 > 0
            inv_κ2_par += weight_par / κ2
            inv_κ2_perp += weight_perp / κ2
        end
    end

    inv_κ1_par *= 3.0 / 4.0
    inv_κ2_par *= 3.0 / 4.0
    inv_κ1_perp *= 3.0 / 2.0
    inv_κ2_perp *= 3.0 / 2.0

    # Per-mode angle-averaged opacities
    κ1_par = inv_κ1_par > 0 ? 1.0 / inv_κ1_par : 0.0
    κ2_par = inv_κ2_par > 0 ? 1.0 / inv_κ2_par : 0.0
    κ1_perp = inv_κ1_perp > 0 ? 1.0 / inv_κ1_perp : 0.0
    κ2_perp = inv_κ2_perp > 0 ? 1.0 / inv_κ2_perp : 0.0

    # Combine modes via harmonic mean for unpolarised radiation
    κ_par = harmonic_mean_2(κ1_par, κ2_par)
    κ_perp = harmonic_mean_2(κ1_perp, κ2_perp)

    return κ_par, κ_perp
end

# Backwards-compatible wrappers
function kappa_parallel_mono(ν::Float64, B::Float64, T::Float64,
                              ρ::Float64; N_θ::Int=40)::Float64
    κ_par, _ = kappa_parallel_perp_mono(ν, B, T, ρ; N_θ=N_θ)
    return κ_par
end

function kappa_perp_mono(ν::Float64, B::Float64, T::Float64,
                          ρ::Float64; N_θ::Int=40)::Float64
    _, κ_perp = kappa_parallel_perp_mono(ν, B, T, ρ; N_θ=N_θ)
    return κ_perp
end

"""
    make_rosseland_frequency_grid(B, T; N_base=200, N_resonance=20) → ν_grid

Construct a log-spaced frequency grid with extra points near cyclotron resonances.
"""
function make_rosseland_frequency_grid(B::Float64, T::Float64;
                                       N_base::Int=200, N_resonance::Int=20)
    ν_min = 0.01 * k_B * T / h
    ν_max = 200.0 * k_B * T / h
    base = [10.0^logν for logν in range(log10(ν_min), log10(ν_max), length=N_base)]

    # Cluster near electron cyclotron frequency
    ω_ce = e_charge * B / (m_e * c)
    ν_ce = ω_ce / (2π)
    if ν_min < ν_ce < ν_max
        resonance = [ν_ce * 10.0^δ for δ in range(-0.3, 0.3, length=N_resonance)]
        append!(base, resonance)
    end

    # Cluster near proton cyclotron frequency
    ν_cp = ν_ce * m_e / m_p
    if ν_min < ν_cp < ν_max
        resonance = [ν_cp * 10.0^δ for δ in range(-0.3, 0.3, length=N_resonance)]
        append!(base, resonance)
    end

    sort!(base)
    unique!(base)
    return base
end

"""
    rosseland_magnetic(B, T, ρ, ν_grid) → (K_∥, K_⊥) [cm²/g]

Rosseland mean opacities along and across B.
These are columns 13-14 of the Potekhin tables.

Uses the combined kappa_parallel_perp_mono which correctly angle-averages
each mode separately before combining them.
"""
function rosseland_magnetic(B::Float64, T::Float64, ρ::Float64,
                             ν_grid::AbstractVector{Float64})
    num_par = 0.0   # ∫ (dB/dT / κ_∥) dν
    num_perp = 0.0
    den = 0.0       # ∫ dB/dT dν

    for i in 1:length(ν_grid)-1
        ν = 0.5 * (ν_grid[i] + ν_grid[i+1])
        dν = ν_grid[i+1] - ν_grid[i]
        db = dBnu_dT(ν, T)

        κ_par, κ_perp = kappa_parallel_perp_mono(ν, B, T, ρ)

        if κ_par > 0
            num_par += db / κ_par * dν
        end
        if κ_perp > 0
            num_perp += db / κ_perp * dν
        end
        den += db * dν
    end

    K_par = den > 0 && num_par > 0 ? den / num_par : 0.0
    K_perp = den > 0 && num_perp > 0 ? den / num_perp : 0.0

    return K_par, K_perp
end

# --- Internal helpers ---

"""Harmonic mean of two positive values. Returns 0 if either is zero."""
function harmonic_mean_2(a::Float64, b::Float64)::Float64
    if a <= 0 || b <= 0
        return max(a, b)
    end
    return 2.0 / (1.0/a + 1.0/b)
end

end # module
