#=
Magnetic free-free cross-sections for fully ionised hydrogen.

Implements the three basic polarisation cross-sections (α = -1, 0, +1)
in a quantizing magnetic field, including both electron-photon and
proton-photon processes.

Source: Potekhin & Chabrier (2003) ApJ 585, 955, Eqs. 33-52.
Local: refs/potekhin_chabrier_2003_ff_opacity.pdf
=#

module MagneticFF

using SpecialFunctions: besselk
using ..PhysicalConstants: e_charge, m_e, m_p, c, h, ħ, k_B, σ_T
using ..MagneticCoulomb: coulomb_log_magnetic, coulomb_log_classical_safe

export cyclotron_freq_e, cyclotron_freq_p, beta_e
export sigma_ff_alpha, sigma_pp_alpha, sigma_scat_alpha
export sigma_total_alpha, coulomb_log_classical, coulomb_log_proton

"Electron cyclotron frequency [rad/s]. ω_ce = eB/(m_e c)."
cyclotron_freq_e(B::Float64) = e_charge * B / (m_e * c)

"Proton cyclotron frequency [rad/s]. ω_cp = eB/(m_p c)."
cyclotron_freq_p(B::Float64) = e_charge * B / (m_p * c)

"Quantization parameter β_e = ℏω_ce/(k_BT). Eq. (1)."
beta_e(B::Float64, T::Float64) = ħ * cyclotron_freq_e(B) / (k_B * T)

"""
    coulomb_log_classical(u) → Λ_cl

Classical (non-magnetic) Coulomb logarithm. P&C 2003 Eq. (43).
Λ_cl = e^{u/2} K₀(u/2)
"""
function coulomb_log_classical(u::Float64)::Float64
    u2 = u / 2.0
    if u2 > 500.0
        return 1.0  # asymptotic
    end
    return exp(u2) * besselk(0, u2)
end

"""
    coulomb_log_proton(u) → Λ_pp

Proton-proton Coulomb logarithm. P&C 2003 Eq. (48).
Analytic fit accurate to ~1.5%.
"""
function coulomb_log_proton(u::Float64)::Float64
    return 0.6 * log(22.0 * u^(-1) + 9.0 * u^(-0.3)) + 0.4 * sqrt(π * u)
end

"""
    nu_ff_alpha(α, ω, B, T, n_e) → ν_α^ff [s⁻¹]

Effective free-free absorption frequency. P&C 2003 Eq. (41).
ν_α^ff = (4/3)√(2π/(m_e k_BT)) × (n_e e⁴/(ℏω)) × (1-e^{-u}) × Λ_α^ff
"""
function nu_ff_alpha(α::Int, ω::Float64, B::Float64,
                     T::Float64, n_e::Float64)::Float64
    u = ħ * ω / (k_B * T)
    u = max(u, 1e-30)

    # Use full magnetic Coulomb logarithm when field is quantizing
    β = beta_e(B, T)
    if β > 0.1
        Λ = coulomb_log_magnetic(α, u, β)
    else
        Λ = coulomb_log_classical_safe(u)
    end

    stimulated = u < 500.0 ? (1.0 - exp(-u)) : 1.0

    pref = 4.0/3.0 * sqrt(2π / (m_e * k_B * T))
    return pref * n_e * e_charge^4 / (ħ * ω) * stimulated * Λ
end

"""
    nu_pp(ω, B, T, n_p) → ν^pp [s⁻¹]

Effective proton-proton collision frequency. P&C 2003 Eq. (47).
"""
function nu_pp(ω::Float64, T::Float64, n_p::Float64)::Float64
    u = ħ * ω / (k_B * T)
    u = max(u, 1e-30)

    Λ = coulomb_log_proton(u)
    stimulated = u < 500.0 ? (1.0 - exp(-u)) : 1.0

    pref = 256.0/3.0 * sqrt(π / (m_p * k_B * T))
    return pref * n_p * e_charge^4 / (ħ * ω) * (k_B * T / (m_p * c^2)) * stimulated * Λ
end

"""
    sigma_ff_alpha(α, ω, B, T, ρ) → σ_α^ff [cm²]

Magnetic free-free cross-section per atom for polarisation α.
P&C 2003 Eq. (51), combined electron+proton treatment.

σ_α^ff ≈ ω² / [(ω+αω_ce)²(ω-αω_cp)² + ω²ν̃_α²] × 4πe²ν_α^ff/(m_e c)
"""
function sigma_ff_alpha(α::Int, ω::Float64, B::Float64,
                        T::Float64, ρ::Float64)::Float64
    @assert α ∈ (-1, 0, 1)
    @assert ω > 0 && B >= 0 && T > 0 && ρ > 0

    if B < 1e6  # non-magnetic limit
        return _sigma_ff_nonmag(ω, T, ρ)
    end

    n_e = ρ / (m_p + m_e)
    ω_ce = cyclotron_freq_e(B)
    ω_cp = cyclotron_freq_p(B)

    ν_ff = nu_ff_alpha(α, ω, B, T, n_e)

    # Natural radiative widths (Eq. 39 and analogous for proton)
    ν_e_s = 2.0/3.0 * e_charge^2 / (m_e * c^3) * ω^2
    ν_p_s = 2.0/3.0 * e_charge^2 / (m_p * c^3) * ω^2

    # Proton-proton damping
    ν_pp_val = nu_pp(ω, T, n_e)

    # Damping frequencies
    ν_e = ν_e_s + ν_pp_val
    ν_p = ν_p_s + ν_pp_val

    # Combined damping (Eq. 52)
    ν_tilde = ν_ff
    if ω > 0
        ν_tilde += (1.0 + α * ω_ce / ω) * ν_p + (1.0 - α * ω_cp / ω) * ν_e
    end

    # Resonance denominators
    denom_e = ω + α * ω_ce   # electron cyclotron resonance
    denom_p = ω - α * ω_cp   # proton cyclotron resonance

    # Full denominator (Eq. 51)
    denom = denom_e^2 * denom_p^2 + ω^2 * ν_tilde^2

    # Cross-section
    σ = ω^2 / denom * 4π * e_charge^2 * ν_ff / (m_e * c)

    return max(σ, 0.0)
end

"""
Non-magnetic free-free cross-section (B→0 limit) per atom.
"""
function _sigma_ff_nonmag(ω::Float64, T::Float64, ρ::Float64)::Float64
    n_e = ρ / (m_p + m_e)
    u = ħ * ω / (k_B * T)
    Λ = coulomb_log_classical(max(u, 1e-30))
    stimulated = u < 500.0 ? (1.0 - exp(-u)) : 1.0

    # From Eq. (41) + (37) with ω_ce → 0:
    # σ^ff = (4πe²/(m_e c ω²)) × ν_ff
    ν_ff = 4.0/3.0 * sqrt(2π / (m_e * k_B * T)) * n_e * e_charge^4 / (ħ * ω) * stimulated * Λ
    return 4π * e_charge^2 * ν_ff / (m_e * c * ω^2)
end

"""
    sigma_pp_alpha(α, ω, B, T, ρ) → σ_α^pp [cm²]

Proton-proton collision cross-section. P&C 2003 Eq. (46).
"""
function sigma_pp_alpha(α::Int, ω::Float64, B::Float64,
                        T::Float64, ρ::Float64)::Float64
    if B < 1e6
        return 0.0  # negligible at B=0
    end

    n_p = ρ / (m_p + m_e)
    ω_cp = cyclotron_freq_p(B)

    ν_pp_val = nu_pp(ω, T, n_p)

    # Damping (Eq. 53): ν_{p,α} = (m_e/m_p) ν_α^ff + ν_p^s + ν^pp
    ν_p_s = 2.0/3.0 * e_charge^2 / (m_p * c^3) * ω^2
    ν_ff_scaled = m_e / m_p * nu_ff_alpha(α, ω, B, T, n_p)
    ν_p_alpha = ν_ff_scaled + ν_p_s + ν_pp_val

    denom = (ω - α * ω_cp)^2 + ν_p_alpha^2
    σ = 4π * e_charge^2 * ν_pp_val / (m_p * c) / denom

    return max(σ, 0.0)
end

"""
    sigma_scat_alpha(α, ω, B) → σ_α^{e,s} [cm²]

Magnetic electron scattering cross-section. P&C 2003 Eq. (33).
"""
function sigma_scat_alpha(α::Int, ω::Float64, B::Float64)::Float64
    if B < 1e6
        return σ_T / 3.0  # each polarisation gets 1/3 of Thomson
    end

    ω_ce = cyclotron_freq_e(B)
    ν_e_s = 2.0/3.0 * e_charge^2 / (m_e * c^3) * ω^2

    denom = (ω + α * ω_ce)^2 + ν_e_s^2
    return ω^2 / denom * σ_T
end

"""
    sigma_total_alpha(α, ω, B, T, ρ) → σ_α^total [cm²]

Total cross-section (absorption + scattering) for polarisation α.
For fully ionised hydrogen (x_H = 0): σ_α = σ_α^ff + σ_α^pp + σ_α^scat.
"""
function sigma_total_alpha(α::Int, ω::Float64, B::Float64,
                           T::Float64, ρ::Float64)::Float64
    return sigma_ff_alpha(α, ω, B, T, ρ) +
           sigma_pp_alpha(α, ω, B, T, ρ) +
           sigma_scat_alpha(α, ω, B)
end

end # module
