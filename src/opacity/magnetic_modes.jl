#=
Normal mode decomposition and angle-averaged opacities in a magnetic field.

Computes the polarisation-dependent opacity at each angle θ_B, then
performs the angle and frequency averaging to get Rosseland mean
opacities K_∥ and K_⊥.

Source: Potekhin & Chabrier (2003) ApJ 585, 955, Eqs. 25-30.
=#

module MagneticModes

using ..PhysicalConstants: e_charge, m_e, m_p, c, h, ħ, k_B, σ_T
using ..MagneticFF: cyclotron_freq_e, cyclotron_freq_p,
                     sigma_total_alpha, sigma_ff_alpha, sigma_pp_alpha,
                     sigma_scat_alpha
using ..HydrogenOpacity: dBnu_dT
using ..BlackbodyAtmosphere: planck_Bnu

export mode_opacity, effective_opacity
export kappa_parallel_mono, kappa_perp_mono
export rosseland_magnetic

const m_H = m_p + m_e

"""
    polarisation_weights(ω, B, θ_B, n_e) → (w1, w2)

Compute |e_{j,α}|² weights for the two normal modes (j=1: extraordinary, j=2: ordinary).

In the simplified cold-plasma limit (Eq. 26), the mode parameter q determines
the decomposition. For |q| >> 1 (far from crossings):
  Mode 1 (X): |e_{1,+1}|² ≈ sin²θ/(1+cos²θ), |e_{1,-1}|² ≈ 2cos²θ/(1+cos²θ), |e_{1,0}|² ≈ sin²θcos²θ/(...)
  Mode 2 (O): complementary weights

Returns: (w1[3], w2[3]) where w[α+2] = |e_{j,α}|² for α ∈ {-1, 0, +1}.
"""
function polarisation_weights(ω::Float64, B::Float64,
                               θ_B::Float64, n_e::Float64)
    ω_ce = cyclotron_freq_e(B)
    ω_pl2 = 4π * n_e * e_charge^2 / m_e  # plasma frequency squared

    cosθ = cos(θ_B)
    sinθ = sin(θ_B)

    # Cold plasma approximation for q (Eq. 26 simplified)
    # q ≈ (ω_ce/ω) × sinθ²/(2cosθ) for ω >> ω_pl
    if abs(cosθ) < 1e-10
        # θ_B ≈ 90°: modes are pure linear (along and across B)
        # Mode 1: α=0 (along B), Mode 2: mix of ±1
        w1 = [0.0, 1.0, 0.0]  # [α=-1, α=0, α=+1]
        w2 = [0.5, 0.0, 0.5]
        return w1, w2
    end

    q = ω_ce * sinθ^2 / (2ω * abs(cosθ))

    # For large |q|: mode 1 is nearly linearly polarised along B
    # For small |q|: modes are nearly circularly polarised
    q2 = q^2
    sq = sqrt(1.0 + q2)

    # Stokes parameter decomposition (Ginzburg 1970, §10)
    # Mode 1 (extraordinary): favours α = -1 at small θ, α = 0 at large θ
    # Mode 2 (ordinary): complementary
    f = q / (sq + 1.0)  # ∈ [0, 1) as q goes from 0 to ∞

    # Weights for mode 1 (extraordinary)
    cos2 = cosθ^2
    sin2 = sinθ^2

    # Simplified decomposition valid away from mode crossings
    w1_m1 = cos2 / (1.0 + cos2) * (1.0 - f)     # α = -1
    w1_0  = sin2 * cos2 / (1.0 + cos2) * f * 2   # α = 0
    w1_p1 = 1.0 / (1.0 + cos2) * (1.0 - f)       # α = +1

    # Normalise
    s1 = w1_m1 + w1_0 + w1_p1
    if s1 > 0
        w1_m1 /= s1; w1_0 /= s1; w1_p1 /= s1
    end

    # Mode 2 (ordinary): complementary
    w2_m1 = 1.0 / (1.0 + cos2) * f           # α = -1
    w2_0  = sin2 / (1.0 + cos2)               # α = 0
    w2_p1 = cos2 / (1.0 + cos2) * f           # α = +1

    s2 = w2_m1 + w2_0 + w2_p1
    if s2 > 0
        w2_m1 /= s2; w2_0 /= s2; w2_p1 /= s2
    end

    return [w1_m1, w1_0, w1_p1], [w2_m1, w2_0, w2_p1]
end

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

    w1, w2 = polarisation_weights(ω, B, θ_B, n_e)
    w = j == 1 ? w1 : w2

    κ = 0.0
    for (idx, α) in enumerate((-1, 0, 1))
        σ = sigma_total_alpha(α, ω, B, T, ρ)
        κ += w[idx] * σ
    end

    return κ / m_H
end

"""
    effective_opacity(ν, θ_B, B, T, ρ) → κ_eff [cm²/g]

Effective non-polarised opacity (harmonic mean of two modes).
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
    kappa_parallel_mono(ν, B, T, ρ; N_θ=20) → κ_∥(ν) [cm²/g]

Monochromatic angle-averaged opacity for transport along B.
P&C 2003 Eq. (30): 1/κ_∥ = (3/4) ∫₀^π cos²θ/κ(θ) sinθ dθ
"""
function kappa_parallel_mono(ν::Float64, B::Float64, T::Float64,
                              ρ::Float64; N_θ::Int=20)::Float64
    # Gauss-Legendre quadrature on [0, π] for the integral
    inv_κ = 0.0
    dθ = π / N_θ

    for i in 1:N_θ
        θ = (i - 0.5) * dθ
        κ_eff = effective_opacity(ν, θ, B, T, ρ)
        if κ_eff > 0
            inv_κ += cos(θ)^2 / κ_eff * sin(θ) * dθ
        end
    end

    inv_κ *= 3.0 / 4.0

    return inv_κ > 0 ? 1.0 / inv_κ : 0.0
end

"""
    kappa_perp_mono(ν, B, T, ρ; N_θ=20) → κ_⊥(ν) [cm²/g]

Monochromatic angle-averaged opacity for transport across B.
P&C 2003 Eq. (30): 1/κ_⊥ = (3/2) ∫₀^π sin³θ/κ(θ) dθ
"""
function kappa_perp_mono(ν::Float64, B::Float64, T::Float64,
                          ρ::Float64; N_θ::Int=20)::Float64
    inv_κ = 0.0
    dθ = π / N_θ

    for i in 1:N_θ
        θ = (i - 0.5) * dθ
        κ_eff = effective_opacity(ν, θ, B, T, ρ)
        if κ_eff > 0
            inv_κ += sin(θ)^3 / κ_eff * dθ
        end
    end

    inv_κ *= 3.0 / 2.0

    return inv_κ > 0 ? 1.0 / inv_κ : 0.0
end

"""
    rosseland_magnetic(B, T, ρ, ν_grid) → (K_∥, K_⊥) [cm²/g]

Rosseland mean opacities along and across B.
These are columns 13-14 of the Potekhin tables.
"""
function rosseland_magnetic(B::Float64, T::Float64, ρ::Float64,
                             ν_grid::AbstractVector{Float64})
    num_par = 0.0  # ∫ (1/κ_∥) dB/dT dν
    num_perp = 0.0
    den = 0.0      # ∫ dB/dT dν

    for i in 1:length(ν_grid)-1
        ν = 0.5 * (ν_grid[i] + ν_grid[i+1])
        dν = ν_grid[i+1] - ν_grid[i]
        db = dBnu_dT(ν, T)

        κ_par = kappa_parallel_mono(ν, B, T, ρ)
        κ_perp = kappa_perp_mono(ν, B, T, ρ)

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

end # module
