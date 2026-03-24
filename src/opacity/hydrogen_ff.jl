#=
Free-free and Thomson scattering opacities for fully ionised hydrogen.
B = 0 (non-magnetic) case.

Source: Haakonsen et al. (2012) ApJ 749:52, Eqs. 11-12.
Local: refs/haakonsen_2012_mcphac.pdf
=#

module HydrogenOpacity

using ..PhysicalConstants: e_charge, m_e, m_p, h, c, k_B, σ_T
using ..GauntFactor: GauntTable, gaunt_ff

export kappa_ff, sigma_thomson, total_opacity, scattering_albedo
export rosseland_mean, dBnu_dT

"""
    kappa_ff(ν, T, ρ, table) → κ_ν [cm²/g]

Free-free absorption opacity for fully ionised hydrogen.
Haakonsen et al. (2012) Eq. 12.

κ_ν = (4e⁶/(3m_e hc)) √(2π/(3k_B m_e)) T^{-1/2} ρ/(m_p+m_e)²
      × ν^{-3} (1 - exp(-hν/(k_BT))) g̃_ff
"""
function kappa_ff(ν::Float64, T::Float64, ρ::Float64,
                  table::GauntTable)::Float64
    @assert ν > 0 && T > 0 && ρ > 0

    gff = gaunt_ff(ν, T, table)

    # Prefactor (CGS constants)
    # 4e⁶/(3 m_e h c) × √(2π/(3 k_B m_e)) / (m_p + m_e)²
    pref = 4.0 * e_charge^6 / (3.0 * m_e * h * c) *
           sqrt(2π / (3.0 * k_B * m_e)) / (m_p + m_e)^2

    x = h * ν / (k_B * T)
    stimulated = x < 500.0 ? (1.0 - exp(-x)) : 1.0

    κ = pref * T^(-0.5) * ρ * ν^(-3) * stimulated * gff

    @assert isfinite(κ) && κ >= 0 "κ_ff non-finite: $κ at ν=$ν, T=$T, ρ=$ρ"
    return κ
end

"""
    sigma_thomson() → σ [cm²/g]

Reduced Thomson scattering opacity for fully ionised hydrogen.
Haakonsen et al. (2012) Eq. 11 with f_ion = 1:
σ_{T,t} = σ_T / (m_p + m_e)
"""
function sigma_thomson()::Float64
    return σ_T / (m_p + m_e)
end

"""
    total_opacity(ν, T, ρ, table) → k_ν [cm²/g]

Total extinction opacity: k_ν = κ_ff + σ_thomson.
"""
function total_opacity(ν::Float64, T::Float64, ρ::Float64,
                       table::GauntTable)::Float64
    return kappa_ff(ν, T, ρ, table) + sigma_thomson()
end

"""
    scattering_albedo(ν, T, ρ, table) → ρ_ν

Scattering albedo: ρ_ν = σ / (κ + σ). Haakonsen Eq. 6.
"""
function scattering_albedo(ν::Float64, T::Float64, ρ::Float64,
                           table::GauntTable)::Float64
    σ = sigma_thomson()
    k = kappa_ff(ν, T, ρ, table) + σ
    return σ / k
end

"""
    dBnu_dT(ν, T) → dB_ν/dT [erg s⁻¹ cm⁻² Hz⁻¹ sr⁻¹ K⁻¹]

Derivative of Planck function with respect to temperature.
"""
function dBnu_dT(ν::Float64, T::Float64)::Float64
    x = h * ν / (k_B * T)
    if x > 500.0
        return 0.0
    end
    ex = exp(x)
    return 2h^2 * ν^4 / (c^2 * k_B * T^2) * ex / (ex - 1)^2
end

"""
    rosseland_mean(T, ρ, ν_grid, table) → k_R [cm²/g]

Rosseland mean opacity: 1/k_R = ∫(1/k_ν)(dB_ν/dT)dν / ∫(dB_ν/dT)dν
"""
function rosseland_mean(T::Float64, ρ::Float64,
                        ν_grid::AbstractVector{Float64},
                        table::GauntTable)::Float64
    num = 0.0  # ∫ (1/k_ν) dB/dT dν
    den = 0.0  # ∫ dB/dT dν

    for i in 1:length(ν_grid)-1
        ν = 0.5 * (ν_grid[i] + ν_grid[i+1])
        dν = ν_grid[i+1] - ν_grid[i]
        k = total_opacity(ν, T, ρ, table)
        db = dBnu_dT(ν, T)
        num += db / k * dν
        den += db * dν
    end

    return den > 0 ? den / num : sigma_thomson()
end

end # module
