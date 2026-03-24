#=
Feautrier radiative transfer solver.

Solves the second-order RT equation at a single frequency:
  μ_j² d²P_j/dτ² = Σ_{j'} R_{j,j'} P_{j'} - T_j

via block-tridiagonal elimination (Auer 1976).

Source: Haakonsen et al. (2012) Sect. 3.3, Eqs. 13-17.
=#

module FeautrierSolver

using LinearAlgebra
using ..PhysicalConstants: k_B, h
using ..BlackbodyAtmosphere: planck_Bnu
using ..AtmosphereStructure: AtmosphereColumn

export solve_feautrier_all, gauss_legendre_half

"""
    gauss_legendre_half(M) → (μ, w)

Gauss-Legendre quadrature abscissae and weights on [0, 1].
Returns M points (outward hemisphere only).
"""
function gauss_legendre_half(M::Int)
    # Transform from [-1,1] to [0,1]: x = (t+1)/2, w' = w/2
    μ_full, w_full = _gauss_legendre(2M)
    μ = Float64[]
    w = Float64[]
    for i in eachindex(μ_full)
        if μ_full[i] > 0
            push!(μ, μ_full[i])
            push!(w, w_full[i])
        end
    end
    return μ[1:M], w[1:M]
end

"""
    solve_feautrier_all(col, μ, w; anisotropic=true) → P_all, J, f_ν, h_ν

Solve Feautrier RT at all frequencies. Returns:
- P_all[i, j, k]: Feautrier mean at depth i, angle j, frequency k
- J[i, k]: mean intensity
- f_ν[i, k]: Eddington factor
- h_ν[i, k]: flux Eddington factor
"""
function solve_feautrier_all(col::AtmosphereColumn,
                             μ::Vector{Float64},
                             w::Vector{Float64};
                             anisotropic::Bool=true)
    N = col.N
    M = length(μ)
    K = col.K

    P_all = zeros(N, M, K)
    J = zeros(N, K)
    f_ν = zeros(N, K)
    h_ν = zeros(N, K)

    for k in 1:K
        # Find depth range where τ < τ_max (80)
        N_eff = N
        for i in 2:N
            if col.τ[i, k] > 80.0
                N_eff = i
                break
            end
        end
        N_eff = max(N_eff, 3)

        P_k = _solve_single_frequency(col, k, N_eff, μ, w, anisotropic)

        # Store results (pad beyond N_eff with diffusion value)
        Bν_bottom = planck_Bnu(col.ν_grid[k], col.T[N_eff])
        for j in 1:M
            for i in 1:N_eff
                P_all[i, j, k] = P_k[i, j]
            end
            for i in N_eff+1:N
                P_all[i, j, k] = planck_Bnu(col.ν_grid[k], col.T[i])
            end
        end

        # Compute J, f, h at this frequency
        for i in 1:N
            J[i, k] = sum(P_all[i, j, k] * w[j] for j in 1:M)
            if J[i, k] > 0
                f_ν[i, k] = sum(μ[j]^2 * P_all[i, j, k] * w[j] for j in 1:M) / J[i, k]
                h_ν[i, k] = sum(μ[j] * P_all[i, j, k] * w[j] for j in 1:M) / J[i, k]
            else
                f_ν[i, k] = 1.0 / 3.0
                h_ν[i, k] = 0.0
            end
        end
    end

    return P_all, J, f_ν, h_ν
end

"""
Solve Feautrier at a single frequency index k, for depth points 1:N_eff.
Returns P[i, j] array.
"""
function _solve_single_frequency(col::AtmosphereColumn, k::Int,
                                  N_eff::Int,
                                  μ::Vector{Float64},
                                  w::Vector{Float64},
                                  anisotropic::Bool)
    M = length(μ)
    ν = col.ν_grid[k]

    # Build dtau between adjacent depth points (in units of total opacity × dy)
    dtau = zeros(N_eff)
    for i in 2:N_eff
        dy = col.y[i] - col.y[i-1]
        dtau[i] = 0.5 * (col.k_total[i, k] + col.k_total[i-1, k]) * dy
    end

    # Block tridiagonal system: A_i P_{i-1} + B_i P_i + C_i P_{i+1} = Q_i
    # Forward sweep: store modified B̃ and Q̃
    B_tilde = [zeros(M, M) for _ in 1:N_eff]
    Q_tilde = [zeros(M) for _ in 1:N_eff]
    C_store = [zeros(M, M) for _ in 1:N_eff]

    for i in 1:N_eff
        A = zeros(M, M)
        B = zeros(M, M)
        C = zeros(M, M)
        Q = zeros(M)

        Bν = planck_Bnu(ν, col.T[i])
        ρ_alb = col.ρ_alb[i, k]

        if i == 1
            # Surface boundary: I_- = 0 → special form
            # From Auer (1976): P = ½I+ at surface
            dt_next = dtau[2]
            for j in 1:M
                δ = dt_next / μ[j]
                C[j, j] = -1.0 / δ
                B[j, j] = 1.0 + 1.0 / δ
            end
        elseif i == N_eff
            # Bottom boundary: diffusion approximation P = B_ν
            for j in 1:M
                B[j, j] = 1.0
            end
            Q .= Bν
            # No scattering terms needed (isotropic at depth)
            B_tilde[i] = B
            Q_tilde[i] = Q
            if i > 1
                # Eliminate A
                A_mat = zeros(M, M)
                dt_prev = dtau[i]
                for j in 1:M
                    δ = dt_prev / μ[j]
                    A_mat[j, j] = -1.0 / (δ * dt_prev / μ[j])
                end
                # Actually at bottom BC: P = B_ν regardless, so just set it
            end
            B_tilde[i] = B
            Q_tilde[i] = Q
            C_store[i] = C
            # Skip elimination for bottom BC
            if i > 1
                Binv_prev = B_tilde[i-1] \ I(M)
                B_tilde[i] = B  # Keep as identity
                Q_tilde[i] = Q  # Keep as B_ν
            end
            continue
        else
            # Interior point
            dt_prev = dtau[i]
            dt_next = dtau[min(i+1, N_eff)]
            dt_avg = 0.5 * (dt_prev + dt_next)
            for j in 1:M
                coeff_A = 1.0 / (dt_prev * dt_avg) * μ[j]^2
                coeff_C = 1.0 / (dt_next * dt_avg) * μ[j]^2
                A[j, j] = -coeff_A
                C[j, j] = -coeff_C
                B[j, j] = coeff_A + coeff_C + 1.0
            end
        end

        # Add scattering redistribution to B (off-diagonal)
        # Haakonsen Eq. 4,7: S contains 2/(κ+σ) × Σ P × dσ/dμ' × w
        # dσ/dμ' = σ_T × (3/16)(3 + 3μ²μ'² - μ² - μ'²) (aniso, Eq. 16)
        # dσ/dμ' = σ_T × 1/2 (isotropic, Eq. 17 with g=1/2)
        # Factor in R matrix: 2ρ × dσ/(σ×dμ') × w
        if i < N_eff  # not at bottom BC
            for j in 1:M, jp in 1:M
                if anisotropic
                    # Haakonsen Eq. 16: cross = (3/16)(3+3μ²μ'²-μ²-μ'²)
                    cross = 3.0/16.0 * (3.0 + 3.0*μ[j]^2*μ[jp]^2 - μ[j]^2 - μ[jp]^2)
                else
                    # Isotropic: dσ/dμ' = σ/2, so cross = 1/2
                    cross = 0.5
                end
                B[j, jp] -= 2.0 * ρ_alb * cross * w[jp]
            end

            # Thermal source
            Q .= (1.0 - ρ_alb) * Bν
        end

        C_store[i] = copy(C)

        # Forward elimination
        if i == 1
            B_tilde[i] = copy(B)
            Q_tilde[i] = copy(Q)
        else
            # B̃_i = B_i - A_i × B̃_{i-1}⁻¹ × C_{i-1}
            tmp = B_tilde[i-1] \ C_store[i-1]
            B_tilde[i] = B - A * tmp
            Q_tilde[i] = Q - A * (B_tilde[i-1] \ Q_tilde[i-1])
        end
    end

    # Back-substitution
    P = zeros(N_eff, M)
    P[N_eff, :] = B_tilde[N_eff] \ Q_tilde[N_eff]

    for i in N_eff-1:-1:1
        P[i, :] = B_tilde[i] \ (Q_tilde[i] - C_store[i] * P[i+1, :])
    end

    # Sanity: P should be non-negative (specific intensity average)
    for i in 1:N_eff, j in 1:M
        P[i, j] = max(P[i, j], 0.0)
    end

    return P
end

"""
Gauss-Legendre quadrature on [-1, 1] via eigenvalue method.
"""
function _gauss_legendre(n::Int)
    β = [i / sqrt(4i^2 - 1) for i in 1:n-1]
    T = SymTridiagonal(zeros(n), β)
    vals, vecs = eigen(T)
    w = 2.0 .* vecs[1, :].^2
    return vals, w
end

end # module
