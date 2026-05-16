#=
Temperature correction for NS atmosphere iteration — Rybicki method.

Computes ΔT(y) that drives the atmosphere toward radiative equilibrium
(constant integrated flux at all depths). Uses the global coupling method
of Rybicki (1971) as implemented by Haakonsen et al. (2012) Appendix A.

The Rybicki method couples all depths and frequencies through a global
matrix system, unlike the local Unsöld-Lucy correction which diverges
for scattering-dominated atmospheres.

Source: Haakonsen et al. (2012) ApJ 749:52, Appendix A, Eqs. A1-A33.
        Rybicki, G.B. (1971) JQSRT 11, 589.
=#

module TemperatureCorrection

using LinearAlgebra
using ..PhysicalConstants: k_B, h, σ_SB
using ..Tridiag: tridiag_lu_forward!, tridiag_lu_back!
using ..BlackbodyAtmosphere: planck_Bnu
using ..HydrogenOpacity: dBnu_dT
using ..AtmosphereStructure: AtmosphereColumn

export compute_temperature_correction

"""
    compute_temperature_correction(col, f_ν, h_ν, J) → ΔT [K]

Rybicki temperature correction. Given the current atmosphere structure
and Feautrier solution (Eddington factors f_ν, h_ν, and mean intensity J),
compute the temperature correction ΔT(y) at each depth.

Uses isotropic scattering approximation for the temperature correction
equation (Haakonsen Eq. A2), while the Feautrier solver uses the full
anisotropic scattering.

The method:
1. For each frequency k, build tridiagonal T_k (Eqs. A19-A21, A29-A32)
2. For each k, solve T_k⁻¹ K_k and T_k⁻¹ U_k
3. Accumulate dense W matrix and rhs from all frequencies
4. Solve (W + I) J̄ = rhs
5. ΔT = J̄ - B̄

Arguments:
- col: atmosphere column (T, y, κ, k_total, ρ_alb, ν_grid)
- f_ν: Eddington factors f[depth, freq] from Feautrier
- h_ν: flux Eddington factors h[depth, freq] (surface value used for BC)
- J: mean intensity J[depth, freq] from Feautrier
"""
function compute_temperature_correction(col::AtmosphereColumn{Tp},
                                         f_ν::AbstractMatrix{<:Real},
                                         h_ν::AbstractMatrix{<:Real},
                                         J::AbstractMatrix{<:Real}
                                         ) where {Tp<:Real}
    N = col.N
    nν = col.nν
    ν = col.ν_grid

    # Working real type: promotes column, Eddington factors, intensity.
    R = promote_type(Tp, eltype(f_ν), eltype(h_ν), eltype(J))

    # Frequency quadrature weights (trapezoid on the log-spaced grid)
    b = zeros(R, nν)
    for k in 1:nν
        if k == 1
            b[k] = 0.5 * (ν[2] - ν[1])
        elseif k == nν
            b[k] = 0.5 * (ν[nν] - ν[nν-1])
        else
            b[k] = 0.5 * (ν[k+1] - ν[k-1])
        end
    end

    # Pre-compute per-depth weighted sums for B̄ and denominator
    denom = zeros(R, N)
    B_bar = zeros(R, N)
    for i in 1:N
        for k in 1:nν
            dBdT = dBnu_dT(ν[k], col.T[i])
            κ_k = col.κ[i, k]  # absorption opacity only
            denom[i] += dBdT * κ_k * b[k]
            B_bar[i] += planck_Bnu(ν[k], col.T[i]) * κ_k * b[k]
        end
        if denom[i] > 0
            B_bar[i] /= denom[i]
        end
    end

    W = zeros(R, N, N)
    rhs = zeros(R, N)

    for k in 1:nν
        T_diag = zeros(R, N)
        T_sub = zeros(R, N-1)
        T_sup = zeros(R, N-1)
        U_k = zeros(R, N)
        K_k = zeros(R, N)

        _build_rybicki_system!(T_diag, T_sub, T_sup, U_k, K_k,
                                col, k, f_ν, h_ν, ν, denom, b, B_bar)

        # Precompute LU factorization of T_k (Thomas algorithm forward sweep)
        lu_diag, lu_sup, lu_rhs_K = tridiag_lu_forward!(T_sub, T_diag, T_sup, K_k)

        # Solve T_k⁻¹ K_k via back-substitution
        x_k = tridiag_lu_back!(lu_diag, lu_sup, lu_rhs_K)

        # V_k diagonal: (V_k)_i = κ_k b_k / denom_i
        V_k = zeros(R, N)
        for i in 1:N
            V_k[i] = denom[i] > 0 ? col.κ[i, k] * b[k] / denom[i] : zero(R)
        end

        # Accumulate rhs += V_k .* x_k
        for i in 1:N
            rhs[i] += V_k[i] * x_k[i]
        end

        # For W: need (T_k⁻¹)_{i,j} × U_{k,j} for each j where U_{k,j} ≠ 0
        for j in 1:N
            if abs(U_k[j]) < 1e-30
                continue  # skip zero columns
            end
            e_j = zeros(R, N)
            e_j[j] = one(R)
            _, _, lu_rhs_j = tridiag_lu_forward!(T_sub, T_diag, T_sup, e_j)
            z_j = tridiag_lu_back!(lu_diag, lu_sup, lu_rhs_j)

            for i in 1:N
                W[i, j] += V_k[i] * U_k[j] * z_j[i]
            end
        end
    end

    # Solve (W + I) J̄ = rhs  →  (W + I) is dense N×N
    for i in 1:N
        W[i, i] += one(R)
    end

    J_bar = W \ rhs

    # Temperature correction: ΔT_i = J̄_i - B̄_i  (from Eq. A5)
    ΔT = J_bar .- B_bar

    return ΔT
end

"""
Build the tridiagonal matrix T_k and vectors U_k, K_k for frequency index k.
Implements Haakonsen Eqs. A19-A23 with boundary conditions A29-A33.
"""
function _build_rybicki_system!(T_diag::AbstractVector, T_sub::AbstractVector,
                                 T_sup::AbstractVector,
                                 U_k::AbstractVector, K_k::AbstractVector,
                                 col, k, f_ν, h_ν, ν, denom, b, B_bar)
    N = col.N

    for i in 2:N-1
        Δ1 = (col.k_total[i, k] + col.k_total[i-1, k]) * (col.y[i] - col.y[i-1])
        Δ2 = (col.k_total[i, k] + col.k_total[i+1, k]) * (col.y[i+1] - col.y[i])
        π_sum = Δ1 + Δ2

        ρ_k = col.ρ_alb[i, k]
        f_k = f_ν[i, k]

        # Tridiagonal elements (Eqs. A19-A21)
        T_sub[i-1] = -8 * f_ν[i-1, k] / (π_sum * Δ1)
        T_sup[i] = -8 * f_ν[i+1, k] / (π_sum * Δ2)
        T_diag[i] = ((1 - ρ_k) / f_k + 8 / (π_sum * Δ1) + 8 / (π_sum * Δ2)) * f_k

        # U_k (Eq. A22): -(dB_k/dT)(1-ρ_k)
        dBdT = dBnu_dT(ν[k], col.T[i])
        U_k[i] = -dBdT * (1 - ρ_k)

        # K_k (Eq. A23)
        B_k = planck_Bnu(ν[k], col.T[i])
        B_bar_i = B_bar[i]
        K_k[i] = (B_k - dBdT * B_bar_i) * (1 - ρ_k)
    end

    # Surface boundary condition (i=1, Eqs. A28-A30)
    Δ_surf = 0.5 * (col.k_total[2, k] + col.k_total[1, k]) * (col.y[2] - col.y[1])
    h_k = h_ν[1, k]
    f_k1 = f_ν[1, k]
    f_k2 = f_ν[2, k]

    T_diag[1] = h_k + f_k1 / Δ_surf
    T_sup[1] = -f_k2 / Δ_surf

    # Surface BC for U_k and K_k: both zero
    U_k[1] = zero(eltype(U_k))
    K_k[1] = zero(eltype(K_k))

    # Bottom boundary condition (i=N): J_ν = B_ν
    T_diag[N] = one(eltype(T_diag))
    if N > 1
        T_sub[N-1] = zero(eltype(T_sub))
    end
    U_k[N] = zero(eltype(U_k))
    K_k[N] = planck_Bnu(ν[k], col.T[N])
end

end # module
