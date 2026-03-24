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
function compute_temperature_correction(col::AtmosphereColumn,
                                         f_ν::Matrix{Float64},
                                         h_ν::Matrix{Float64},
                                         J::Matrix{Float64})::Vector{Float64}
    N = col.N
    K = col.K
    ν = col.ν_grid

    # Frequency quadrature weights (trapezoid on the log-spaced grid)
    b = zeros(K)
    for k in 1:K
        if k == 1
            b[k] = 0.5 * (ν[2] - ν[1])
        elseif k == K
            b[k] = 0.5 * (ν[K] - ν[K-1])
        else
            b[k] = 0.5 * (ν[k+1] - ν[k-1])
        end
    end

    # Pre-compute per-depth weighted sums for B̄ and denominator
    # denom_i = Σ_k (dB_k/dT)(y_i) × κ_k(y_i) × b_k
    # B̄_i = Σ_k B_k(T₀(y_i)) × κ_k(y_i) × b_k / denom_i
    denom = zeros(N)
    B_bar = zeros(N)
    for i in 1:N
        for k in 1:K
            dBdT = dBnu_dT(ν[k], col.T[i])
            κ_k = col.κ[i, k]  # absorption opacity only
            denom[i] += dBdT * κ_k * b[k]
            B_bar[i] += planck_Bnu(ν[k], col.T[i]) * κ_k * b[k]
        end
        if denom[i] > 0
            B_bar[i] /= denom[i]
        end
    end

    # Accumulate W matrix (N×N dense) and rhs (N-vector) over frequencies
    #
    # From Eq. A24: [I + Σ_k V_k T_k⁻¹ diag(U_k)] J̄ = Σ_k V_k T_k⁻¹ K_k
    #
    # V_k is diagonal (N×N), T_k⁻¹ is full (N×N from tridiagonal), diag(U_k) is diagonal.
    # W[i,j] += V_{k,i} × (T_k⁻¹)_{i,j} × U_{k,j}
    # rhs[i] += V_{k,i} × (T_k⁻¹ K_k)_i
    #
    # Computing (T_k⁻¹)_{i,j} requires solving T_k z = e_j for each column j.
    # We precompute the LU factors and do N back-substitutions per frequency.
    W = zeros(N, N)
    rhs = zeros(N)

    for k in 1:K
        T_diag = zeros(N)
        T_sub = zeros(N-1)
        T_sup = zeros(N-1)
        U_k = zeros(N)
        K_k = zeros(N)

        _build_rybicki_system!(T_diag, T_sub, T_sup, U_k, K_k,
                                col, k, f_ν, h_ν, ν, denom, b)

        # Precompute LU factorization of T_k (Thomas algorithm forward sweep)
        lu_diag, lu_sup, lu_rhs_K = _tridiag_lu_forward(T_sub, T_diag, T_sup, K_k)

        # Solve T_k⁻¹ K_k via back-substitution
        x_k = _tridiag_lu_back(lu_diag, lu_sup, lu_rhs_K)

        # V_k diagonal: (V_k)_i = κ_k b_k / denom_i
        V_k = zeros(N)
        for i in 1:N
            V_k[i] = denom[i] > 0 ? col.κ[i, k] * b[k] / denom[i] : 0.0
        end

        # Accumulate rhs += V_k .* x_k
        for i in 1:N
            rhs[i] += V_k[i] * x_k[i]
        end

        # For W: need (T_k⁻¹)_{i,j} × U_{k,j} for each j where U_{k,j} ≠ 0
        # Solve T_k z_j = e_j for each j, then W[i,j] += V_k[i] * U_k[j] * z_j[i]
        for j in 1:N
            if abs(U_k[j]) < 1e-30
                continue  # skip zero columns
            end
            # Solve T_k × z_j = e_j (unit vector)
            e_j = zeros(N)
            e_j[j] = 1.0
            _, _, lu_rhs_j = _tridiag_lu_forward(T_sub, T_diag, T_sup, e_j)
            z_j = _tridiag_lu_back(lu_diag, lu_sup, lu_rhs_j)

            # W[i,j] += V_k[i] * U_k[j] * z_j[i]
            for i in 1:N
                W[i, j] += V_k[i] * U_k[j] * z_j[i]
            end
        end
    end

    # Solve (W + I) J̄ = rhs  →  (W + I) is dense N×N
    for i in 1:N
        W[i, i] += 1.0
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
function _build_rybicki_system!(T_diag, T_sub, T_sup, U_k, K_k,
                                 col, k, f_ν, h_ν, ν, denom, b)
    N = col.N

    for i in 2:N-1
        # Optical depth differences (Eqs. A12-A14)
        # NOTE: Paper defines Δ WITHOUT the 1/2 factor:
        #   Δ₁ = [k(y_i) + k(y_{i-1})] × (y_i - y_{i-1})
        Δ1 = (col.k_total[i, k] + col.k_total[i-1, k]) * (col.y[i] - col.y[i-1])
        Δ2 = (col.k_total[i, k] + col.k_total[i+1, k]) * (col.y[i+1] - col.y[i])
        π_sum = Δ1 + Δ2

        ρ_k = col.ρ_alb[i, k]
        f_k = f_ν[i, k]

        # Tridiagonal elements (Eqs. A19-A21)
        # Note: T_k acts on (f_k J_k), not J_k directly
        T_sub[i-1] = -8.0 * f_ν[i-1, k] / (π_sum * Δ1)
        T_sup[i] = -8.0 * f_ν[i+1, k] / (π_sum * Δ2)
        T_diag[i] = ((1.0 - ρ_k) / f_k + 8.0 / (π_sum * Δ1) + 8.0 / (π_sum * Δ2)) * f_k

        # Hmm wait — re-reading Eqs. A19-A21 more carefully:
        # The equation is for J_k (not f_k J_k).
        # (T_k)_{i,i-1} = -8 f_k(y_{i-1}) / (π Δ₁)
        # (T_k)_{i,i+1} = -8 f_k(y_{i+1}) / (π Δ₂)
        # (T_k)_{i,i} = [(1-ρ_k)/f_k(y_i) + 8/(πΔ₁) + 8/(πΔ₂)] f_k(y_i)
        # These act on J_k directly.

        # U_k (Eq. A22): -(dB_k/dT)(1-ρ_k)
        dBdT = dBnu_dT(ν[k], col.T[i])
        U_k[i] = -dBdT * (1.0 - ρ_k)

        # K_k (Eq. A23)
        B_k = planck_Bnu(ν[k], col.T[i])
        B_bar_i = denom[i] > 0 ? sum(planck_Bnu(ν[kp], col.T[i]) * col.κ[i, kp] * b[kp] for kp in 1:col.K) / denom[i] : B_k
        K_k[i] = (B_k - dBdT * B_bar_i) * (1.0 - ρ_k)
    end

    # Surface boundary condition (i=1, Eqs. A28-A30)
    # Discretized surface BC: (f₂J₂ - f₁J₁)/Δ_surf = h_k J₁
    # → (h_k + f₁/Δ) J₁ - (f₂/Δ) J₂ = 0  [pure radiation constraint]
    # Surface Δ uses the 0.5 mean-opacity convention (matching McPHAC CalcTt.c):
    #   Δ_surf = 0.5 × (k₁ + k₂) × (y₂ - y₁)
    # This differs from the interior Δ which omits the 0.5 factor (Eq. A12).
    Δ_surf = 0.5 * (col.k_total[2, k] + col.k_total[1, k]) * (col.y[2] - col.y[1])
    h_k = h_ν[1, k]
    f_k1 = f_ν[1, k]
    f_k2 = f_ν[2, k]

    T_diag[1] = h_k + f_k1 / Δ_surf
    T_sup[1] = -f_k2 / Δ_surf

    # Surface BC for U_k and K_k: both zero (McPHAC CalcUt.c line 10, CalcKt.c line 15).
    # The surface temperature correction comes through the global W coupling,
    # not through direct local thermal terms.
    U_k[1] = 0.0
    K_k[1] = 0.0

    # Bottom boundary condition (i=N, Eqs. A31-A33): J_ν = B_ν
    T_diag[N] = 1.0
    if N > 1
        T_sub[N-1] = 0.0
    end
    U_k[N] = 0.0  # temperature correction doesn't affect bottom BC
    K_k[N] = planck_Bnu(ν[k], col.T[N])
end

"""
LU forward sweep for tridiagonal system. Returns modified diagonal, super-diagonal, and RHS.
The sub-diagonal multipliers are applied to both the matrix and the RHS.
"""
function _tridiag_lu_forward(sub::Vector{Float64}, diag::Vector{Float64},
                              sup::Vector{Float64}, d::Vector{Float64})
    N = length(diag)
    lu_diag = copy(diag)
    lu_sup = copy(sup)
    lu_rhs = copy(d)

    for i in 2:N
        if abs(lu_diag[i-1]) < 1e-30
            continue
        end
        m = sub[i-1] / lu_diag[i-1]
        lu_diag[i] -= m * lu_sup[i-1]
        lu_rhs[i] -= m * lu_rhs[i-1]
    end

    return lu_diag, lu_sup, lu_rhs
end

"""
Back-substitution after LU forward sweep.
"""
function _tridiag_lu_back(lu_diag::Vector{Float64}, lu_sup::Vector{Float64},
                           lu_rhs::Vector{Float64})::Vector{Float64}
    N = length(lu_diag)
    x = zeros(N)
    if abs(lu_diag[N]) > 1e-30
        x[N] = lu_rhs[N] / lu_diag[N]
    end
    for i in N-1:-1:1
        if abs(lu_diag[i]) > 1e-30
            x[i] = (lu_rhs[i] - lu_sup[i] * x[i+1]) / lu_diag[i]
        end
    end
    return x
end

"""
Solve tridiagonal system Ax = d using Thomas algorithm (convenience wrapper).
"""
function _tridiag_solve(sub::Vector{Float64}, diag::Vector{Float64},
                         sup::Vector{Float64}, d::Vector{Float64})::Vector{Float64}
    lu_diag, lu_sup, lu_rhs = _tridiag_lu_forward(sub, diag, sup, d)
    return _tridiag_lu_back(lu_diag, lu_sup, lu_rhs)
end

end # module
