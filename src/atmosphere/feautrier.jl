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
using ..SolverDefaults: TAU_DIFFUSION

export solve_feautrier_all, solve_feautrier_all_adaptive, gauss_legendre_half, _log_interp

"""
    gauss_legendre_half(M) → (μ, w)

Gauss-Legendre quadrature abscissae and weights on [0, 1].
Returns M points (outward hemisphere only).

For M ∈ {4, 6, 8, 10, 12, 16} a precomputed lookup (`_GL_HALF_TABLE`) is
returned without an `eigen` call. Other M fall back to `_compute_gl_half`.
The cached vectors are read-only; callers must not mutate them.
"""
function gauss_legendre_half(M::Int)
    if haskey(_GL_HALF_TABLE, M)
        return _GL_HALF_TABLE[M]
    end
    return _compute_gl_half(M)
end

"""
    _compute_gl_half(M) → (μ, w)

Eigenvalue-based computation of Gauss-Legendre nodes/weights on [0, 1].
Used to bootstrap `_GL_HALF_TABLE` and as the fallback for non-tabulated M.
"""
function _compute_gl_half(M::Int)
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

# Precomputed Gauss-Legendre nodes/weights on [0,1] for M values used in
# the codebase (rt_atmosphere.jl, magnetic_atmosphere.jl). Values generated
# from `_compute_gl_half` printed at 17 significant digits — bit-identical
# to the eigen path when read back as Float64 literals.
const _GL_HALF_TABLE = Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}(
    4 => (Float64[0.18343464249565056, 0.52553240991632921, 0.79666647741362717, 0.96028985649753618],
          Float64[0.36268378337836232, 0.3137066458778876, 0.22238103445337404, 0.10122853629037609]),
    6 => (Float64[0.1252334085114698, 0.36783149899818013, 0.58731795428661759, 0.76990267419430447, 0.90411725637047469, 0.98156063424671924],
          Float64[0.24914704581340252, 0.23349253653835614, 0.20316742672306526, 0.16007832854334608, 0.10693932599531807, 0.047175336386511696]),
    8 => (Float64[0.095012509837637538, 0.28160355077925947, 0.45801677765722726, 0.61787624440264421, 0.75540440835500311, 0.86563120238783187, 0.94457502307323271, 0.98940093499164994],
          Float64[0.1894506104550675, 0.18260341504492311, 0.16915651939500195, 0.14959598881657718, 0.12462897125553475, 0.095158511682492911, 0.062253523938648039, 0.027152459411754308]),
    10 => (Float64[0.076526521133497782, 0.2277858511416454, 0.37370608871541966, 0.51086700195082713, 0.63605368072651525, 0.74633190646015113, 0.83911697182221889, 0.91223442825132606, 0.96397192727791381, 0.993128599185095],
           Float64[0.15275338713072575, 0.14917298647260305, 0.14209610931838021, 0.13168863844917597, 0.1181945319615197, 0.10193011981724084, 0.083276741576705129, 0.062672048334109318, 0.04060142980038714, 0.017614007139152146]),
    12 => (Float64[0.06405689286260563, 0.19111886747361662, 0.31504267969616384, 0.4337935076260454, 0.54542147138883967, 0.64809365193697532, 0.74012419157855436, 0.82000198597390295, 0.88641552700440096, 0.93827455200273269, 0.97472855597130958, 0.99518721999702142],
           Float64[0.12793819534675185, 0.12583745634682769, 0.12167047292780363, 0.11550566805372615, 0.10744427011596527, 0.097618652104114634, 0.086190161531953122, 0.073346481411080675, 0.05929858491543713, 0.044277438817419704, 0.028531388628933532, 0.012341229799987145]),
    16 => (Float64[0.048307665687739032, 0.14447196158279674, 0.23928736225213731, 0.33186860228212778, 0.4213512761306355, 0.50689990893222947, 0.5877157572407623, 0.66304426693021512, 0.7321821187402896, 0.79448379596794239, 0.84936761373256997, 0.89632115576605209, 0.93490607593773967, 0.96476225558750639, 0.98561151154526838, 0.99726386184948157],
           Float64[0.096540088514727174, 0.095638720079274875, 0.09384439908080483, 0.091173878695764515, 0.087652093004403867, 0.083311924226947665, 0.078193895787070089, 0.072345794108848338, 0.065822222776361586, 0.058684093478535822, 0.050998059262376189, 0.042835898022226787, 0.034273862913021425, 0.025392065309262017, 0.016274394730905473, 0.0070186100094700617]),
)

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
    nν = col.nν

    P_all = zeros(N, M, nν)
    J = zeros(N, nν)
    f_ν = zeros(N, nν)
    h_ν = zeros(N, nν)

    for k in 1:nν
        # Find depth range where τ < τ_max (TAU_DIFFUSION)
        N_eff = N
        for i in 2:N
            if col.τ[i, k] > TAU_DIFFUSION
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
    solve_feautrier_all_adaptive(col, μ, w; N_ν=200, anisotropic=true) → P_all, J, f_ν, h_ν

Frequency-adaptive Feautrier solver (McPHAC GetColumnsNu.c approach).
For each frequency k:
  1. Find y_max where τ_k = TAU_DIFFUSION (=80)
  2. Create N_k = min(N_ν, i_max) log-spaced depth points from y_min to y_max
  3. Interpolate T, k_total, ρ_alb to the per-frequency grid
  4. Solve Feautrier on the per-frequency grid
  5. Interpolate J, f_ν, h_ν, P_all back to the common grid

Status (2026-05-16): On by default in `solve_atmosphere` after a parameter
sweep at T_eff=1e6 K, g=2e14 cm/s² against an N=400 non-adaptive reference:

  N        non-adaptive maxΔI/I    adaptive maxΔI/I    ratio
  50       0.1461                  0.1371              0.94×
  100      0.0462                  0.0377              0.82×
  150      0.0200                  0.0116              0.58×
  200      0.0119                  0.0073              0.61×

Both solvers agree to <1% at N=400 (0.00712 max bulk residual), confirming
they converge to the same limit. The earlier "~7% → <1%" docstring
claim was over-stated, but adaptive does reduce the discretization-error
prefactor by 0.6–0.9× at production N. Pass `adaptive=false` to opt out.
"""
function solve_feautrier_all_adaptive(col::AtmosphereColumn,
                                       μ::Vector{Float64},
                                       w::Vector{Float64};
                                       N_ν::Int=200,
                                       anisotropic::Bool=true)
    N = col.N
    M = length(μ)
    nν = col.nν

    P_all = zeros(N, M, nν)
    J = zeros(N, nν)
    f_ν = zeros(N, nν)
    h_ν = zeros(N, nν)

    for k in 1:nν
        ν = col.ν_grid[k]

        # Find y_max where τ_k reaches TAU_DIFFUSION
        i_max = N
        for i in 2:N
            if col.τ[i, k] > TAU_DIFFUSION
                i_max = i
                break
            end
        end
        i_max = max(i_max, 3)
        y_min = col.y[1]
        y_max = col.y[i_max]

        # Create log-spaced per-frequency grid
        N_k = min(N_ν, i_max)  # don't use more points than needed
        logy_min = log10(max(y_min, 1e-30))
        logy_max = log10(y_max)
        y_nu = [10.0^logy for logy in range(logy_min, logy_max, length=N_k)]

        # Interpolate T, k_total, ρ_alb to per-frequency grid
        T_nu = _log_interp(col.y[1:i_max], col.T[1:i_max], y_nu)
        k_total_nu = _log_interp(col.y[1:i_max],
                                   [col.k_total[i, k] for i in 1:i_max], y_nu)
        ρ_alb_nu = _log_interp(col.y[1:i_max],
                                 [col.ρ_alb[i, k] for i in 1:i_max], y_nu)

        # Build dtau on the per-frequency grid
        dtau_nu = zeros(N_k)
        for i in 2:N_k
            dy = y_nu[i] - y_nu[i-1]
            dtau_nu[i] = 0.5 * (k_total_nu[i] + k_total_nu[i-1]) * dy
        end

        # Solve Feautrier on the per-frequency grid
        P_k = _solve_on_grid(N_k, M, dtau_nu, μ, w, ν, T_nu, ρ_alb_nu, anisotropic)

        # Compute J, f, h on the per-frequency grid
        J_nu = zeros(N_k)
        f_nu = zeros(N_k)
        h_nu = zeros(N_k)
        for i in 1:N_k
            J_nu[i] = sum(P_k[i, j] * w[j] for j in 1:M)
            if J_nu[i] > 0
                f_nu[i] = sum(μ[j]^2 * P_k[i, j] * w[j] for j in 1:M) / J_nu[i]
                h_nu[i] = sum(μ[j] * P_k[i, j] * w[j] for j in 1:M) / J_nu[i]
            else
                f_nu[i] = 1.0 / 3.0
                h_nu[i] = 0.0
            end
        end

        # Interpolate J, f, h back to the common grid
        J_common = _log_interp(y_nu, J_nu, col.y[1:i_max])
        f_common = _log_interp(y_nu, f_nu, col.y[1:i_max])
        h_common = _log_interp(y_nu, h_nu, col.y[1:i_max])

        for i in 1:i_max
            J[i, k] = J_common[i]
            f_ν[i, k] = clamp(f_common[i], 0.0, 1.0)
            h_ν[i, k] = h_common[i]
        end
        # Deep layers: diffusion limit
        for i in i_max+1:N
            J[i, k] = planck_Bnu(ν, col.T[i])
            f_ν[i, k] = 1.0 / 3.0
            h_ν[i, k] = 0.0
        end

        # Also store P_all for emergent intensity (interpolate surface values)
        for j in 1:M
            P_surface = _log_interp(y_nu, [P_k[ii, j] for ii in 1:N_k], col.y[1:i_max])
            for i in 1:i_max
                P_all[i, j, k] = P_surface[i]
            end
            for i in i_max+1:N
                P_all[i, j, k] = planck_Bnu(ν, col.T[i])
            end
        end
    end

    return P_all, J, f_ν, h_ν
end

"""Solve Feautrier on an arbitrary depth grid (resampled version)."""
function _solve_on_grid(N_eff, M, dtau, μ, w, ν, T, ρ_alb, anisotropic)
    B_tilde = [zeros(M, M) for _ in 1:N_eff]
    Q_tilde = [zeros(M) for _ in 1:N_eff]
    C_store = [zeros(M, M) for _ in 1:N_eff]

    for i in 1:N_eff
        A = zeros(M, M)
        B = zeros(M, M)
        C = zeros(M, M)
        Q = zeros(M)

        Bν = planck_Bnu(ν, T[i])
        ρ = ρ_alb[i]

        if i == 1
            # Surface: pure radiation BC
            dt_next = dtau[2]
            for j in 1:M
                δ = dt_next / μ[j]
                C[j, j] = -1.0 / δ
                B[j, j] = 1.0 + 1.0 / δ
            end
        elseif i == N_eff
            # Bottom: LTE
            for j in 1:M; B[j, j] = 1.0; end
            Q .= Bν
            B_tilde[i] = B; Q_tilde[i] = Q; C_store[i] = C
            continue
        else
            # Interior
            dt_prev = dtau[i]
            dt_next = dtau[min(i+1, N_eff)]
            dt_avg = 0.5 * (dt_prev + dt_next)
            for j in 1:M
                coeff_A = μ[j]^2 / (dt_prev * dt_avg)
                coeff_C = μ[j]^2 / (dt_next * dt_avg)
                A[j, j] = -coeff_A
                C[j, j] = -coeff_C
                B[j, j] = coeff_A + coeff_C + 1.0
            end
        end

        # Scattering + source (not at surface)
        if i > 1 && i < N_eff
            if anisotropic
                for j in 1:M, jp in 1:M
                    B[j, jp] -= 2.0 * ρ * 0.5 * w[jp]
                end
            else
                for j in 1:M
                    B[j, j] -= ρ
                end
            end
            Q .= (1.0 - ρ) * Bν
        end

        C_store[i] = copy(C)

        if i == 1
            B_tilde[i] = copy(B)
            Q_tilde[i] = copy(Q)
        else
            tmp = B_tilde[i-1] \ C_store[i-1]
            B_tilde[i] = B - A * tmp
            Q_tilde[i] = Q - A * (B_tilde[i-1] \ Q_tilde[i-1])
        end
    end

    P = zeros(N_eff, M)
    P[N_eff, :] = B_tilde[N_eff] \ Q_tilde[N_eff]
    for i in N_eff-1:-1:1
        P[i, :] = B_tilde[i] \ (Q_tilde[i] - C_store[i] * P[i+1, :])
    end
    for i in 1:N_eff, j in 1:M
        P[i, j] = max(P[i, j], 0.0)
    end
    return P
end

"""Log-linear interpolation: log-space in x, linear in y."""
function _log_interp(x_old::AbstractVector{Float64}, y_old::AbstractVector{Float64},
                      x_new::AbstractVector{Float64};
                      min_x::Float64 = eps(Float64))
    @assert all(x -> x >= 0, x_old) "_log_interp requires non-negative x_old"
    n_old = length(x_old)
    n_new = length(x_new)
    y_new = zeros(n_new)
    logx_old = log10.(max.(x_old, min_x))

    for i in 1:n_new
        logx = log10(max(x_new[i], min_x))
        # Find bracketing index
        j = searchsortedlast(logx_old, logx)
        j = clamp(j, 1, n_old - 1)
        if j >= n_old
            y_new[i] = y_old[end]
        else
            f = (logx - logx_old[j]) / max(logx_old[j+1] - logx_old[j], min_x)
            f = clamp(f, 0.0, 1.0)
            y_new[i] = (1-f) * y_old[j] + f * y_old[j+1]
        end
    end
    return y_new
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
            B_tilde[i] = B
            Q_tilde[i] = Q
            C_store[i] = C
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

        # Add scattering redistribution and thermal source
        # Surface (i=1): pure radiation BC — no scattering, no thermal source.
        #   McPHAC: K[0]=0, U[0]=0 (CalcK.c line 16, CalcU.c line 9).
        #   Outgoing radiation comes through tridiagonal coupling to interior.
        # Interior (2 ≤ i < N_eff): Q = (1-ρ) B_ν, with scattering in B matrix.
        if i > 1 && i < N_eff
            for j in 1:M, jp in 1:M
                if anisotropic
                    cross = 3.0/16.0 * (3.0 + 3.0*μ[j]^2*μ[jp]^2 - μ[j]^2 - μ[jp]^2)
                else
                    cross = 0.5
                end
                B[j, jp] -= 2.0 * ρ_alb * cross * w[jp]
            end
            Q .= (1.0 - ρ_alb) * Bν
        end
        # i=1: Q=0, B unchanged (pure radiation BC)

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
