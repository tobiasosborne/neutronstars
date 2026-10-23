#=
Free-free and Thomson scattering opacities for fully ionised hydrogen.
B = 0 (non-magnetic) case.

Source: Haakonsen et al. (2012) ApJ 749:52, Eqs. 11-12.
Local: refs/haakonsen_2012_mcphac.pdf
=#

module HydrogenOpacity

using ..PhysicalConstants: e_charge, m_e, m_p, m_H, h, c, k_B, σ_T
using ..GauntFactor: GauntTable, gaunt_ff

export kappa_ff, sigma_thomson, total_opacity, scattering_albedo
export rosseland_mean, dBnu_dT

"""
    kappa_ff(ν, T, ρ, table) → κ_ν [cm²/g]

Free-free absorption opacity for fully ionised hydrogen.
Haakonsen et al. (2012) Eq. 12.

κ_ν = (4e⁶/(3m_e hc)) √(2π/(3k_B m_e)) T^{-1/2} ρ/(m_p+m_e)²
      × ν^{-3} (1 - exp(-hν/(k_BT))) g̃_ff

`::Real` accepted (not only `Float64`) so this kernel composes through
ForwardDiff `Dual` numbers — see also the E6 magnetic-chain refactor.
"""
function kappa_ff(ν::Real, T::Real, ρ::Real, table::GauntTable)
    gff = gaunt_ff(ν, T, table)

    # Prefactor (CGS constants) — pure Float64, no AD dependence
    pref = 4.0 * e_charge^6 / (3.0 * m_e * h * c) *
           sqrt(2π / (3.0 * k_B * m_e)) / m_H^2

    x = h * ν / (k_B * T)
    stimulated = x < 500.0 ? (1.0 - exp(-x)) : one(x)

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
    return σ_T / m_H
end

"""
    total_opacity(ν, T, ρ, table) → k_ν [cm²/g]

Total extinction opacity: k_ν = κ_ff + σ_thomson.
"""
function total_opacity(ν::Real, T::Real, ρ::Real, table::GauntTable)
    return kappa_ff(ν, T, ρ, table) + sigma_thomson()
end

"""
    scattering_albedo(ν, T, ρ, table) → ρ_ν

Scattering albedo: ρ_ν = σ / (κ + σ). Haakonsen Eq. 6.
"""
function scattering_albedo(ν::Real, T::Real, ρ::Real, table::GauntTable)
    σ = sigma_thomson()
    k = kappa_ff(ν, T, ρ, table) + σ
    return σ / k
end

"""
    dBnu_dT(ν, T) → dB_ν/dT [erg s⁻¹ cm⁻² Hz⁻¹ sr⁻¹ K⁻¹]

Derivative of Planck function with respect to temperature.
"""
function dBnu_dT(ν::Real, T::Real)
    x = h * ν / (k_B * T)
    if x > 500.0
        return zero(x)
    end
    ex = exp(x)
    return 2h^2 * ν^4 / (c^2 * k_B * T^2) * ex / (ex - 1)^2
end

"""
    rosseland_mean(T, ρ, ν_grid, table) → k_R [cm²/g]

Rosseland mean opacity: 1/k_R = ∫(1/k_ν)(dB_ν/dT)dν / ∫(dB_ν/dT)dν
"""
function rosseland_mean(T::Real, ρ::Real,
                        ν_grid::AbstractVector{<:Real},
                        table::GauntTable)
    # Promote so the accumulators carry the Dual-type if any input is Dual.
    R = promote_type(typeof(T), typeof(ρ), eltype(ν_grid))
    num = zero(R)  # ∫ (1/k_ν) dB/dT dν
    den = zero(R)  # ∫ dB/dT dν

    for i in 1:length(ν_grid)-1
        ν = 0.5 * (ν_grid[i] + ν_grid[i+1])
        dν = ν_grid[i+1] - ν_grid[i]
        k = total_opacity(ν, T, ρ, table)
        db = dBnu_dT(ν, T)
        num += db / k * dν
        den += db * dν
    end

    return den > 0 ? den / num : R(sigma_thomson())
end

end # module
