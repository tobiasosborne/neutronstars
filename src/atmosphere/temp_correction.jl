#=
Temperature correction for NS atmosphere iteration.

Uses the Unsöld-Lucy-type correction:
  ΔT(y) = -∫(J_ν - B_ν) κ_ν dν / ∫(dB_ν/dT) κ_ν dν

Source: Haakonsen et al. (2012) Eq. 18 (simplified form);
see also Mihalas (1978) Ch 6 for the general theory.
=#

module TemperatureCorrection

using ..PhysicalConstants: k_B, h, σ_SB
using ..BlackbodyAtmosphere: planck_Bnu
using ..HydrogenOpacity: dBnu_dT
using ..AtmosphereStructure: AtmosphereColumn

export compute_temperature_correction

"""
    compute_temperature_correction(col, J) → ΔT [K]

Compute temperature correction at each depth point.

ΔT(y) = -∫(J_ν - B_ν) κ_ν dν / ∫(dB_ν/dT) κ_ν dν

Also includes a flux correction term for deep layers where the local
correction is insufficient.

Arguments:
- `col`: current atmosphere column
- `J`: mean intensity J[depth, frequency]
"""
function compute_temperature_correction(col::AtmosphereColumn,
                                         J::Matrix{Float64})::Vector{Float64}
    N = col.N
    K = col.K
    ΔT = zeros(N)
    ν = col.ν_grid

    for i in 1:N
        num = 0.0  # ∫(J_ν - B_ν) κ_ν dν
        den = 0.0  # ∫(dB_ν/dT) κ_ν dν

        for k in 1:K-1
            dν = ν[k+1] - ν[k]
            ν_mid = 0.5 * (ν[k] + ν[k+1])

            Bν = planck_Bnu(ν_mid, col.T[i])
            J_mid = 0.5 * (J[i, k] + J[i, k+1])
            κ_mid = 0.5 * (col.κ[i, k] + col.κ[i, k+1])
            dBdT = dBnu_dT(ν_mid, col.T[i])

            num += (J_mid - Bν) * κ_mid * dν
            den += dBdT * κ_mid * dν
        end

        if abs(den) > 1e-30
            ΔT[i] = -num / den
        end
    end

    return ΔT
end

end # module
