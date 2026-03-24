#=
Magnetic atmosphere solver with two-mode radiative transfer.

In a magnetized NS atmosphere, radiation propagates in two normal modes
(extraordinary/X and ordinary/O) with distinct opacities. Each mode
obeys its own Feautrier equation; they couple through the temperature
correction (both modes contribute to the total radiative flux).

Source: Suleimanov, Potekhin & Werner (2009) A&A 500, 891.
        Haakonsen et al. (2012) ApJ 749:52 (McPHAC), for B=0 limit.
=#

module MagneticAtmosphere

using Printf
using LinearAlgebra
using ..PhysicalConstants: σ_SB, k_B, h, m_p, m_e
using ..GauntFactor: GauntTable
using ..AtmosphereStructure: AtmosphereStructure, make_frequency_grid
using ..BlackbodyAtmosphere: planck_Bnu
using ..HydrogenOpacity: kappa_ff, sigma_thomson, dBnu_dT
using ..MagneticModes: mode_opacity, mode_absorption, mode_scattering
using ..FeautrierSolver: gauss_legendre_half

export solve_magnetic_atmosphere, MagneticAtmosphereResult

const m_H = m_p + m_e

"""
Result of a converged magnetic atmosphere calculation.
"""
struct MagneticAtmosphereResult
    T_eff::Float64
    g_s::Float64
    B::Float64              # magnetic field strength [G]
    θ_B::Float64            # angle between B and surface normal [rad]
    converged::Bool
    n_iterations::Int
    max_dT_over_T::Float64
    ν_grid::Vector{Float64}
    μ_grid::Vector{Float64}
    I_emergent::Array{Float64, 3}   # K × M × 2 (freq, angle, mode)
    T_profile::Vector{Float64}
    y_grid::Vector{Float64}
end

"""
    solve_magnetic_atmosphere(T_eff, g_s, B, θ_B, gaunt; ...) → MagneticAtmosphereResult

Solve a plane-parallel magnetized NS atmosphere with two-mode RT.
Each mode j (1=X, 2=O) has its own opacity κ_j(ν,T,ρ,B,θ_B).
Temperature correction enforces Σ_j ∫ κ_j(J_j - B_ν/2) dν = 0.

For B=0, this should recover the non-magnetic result (both modes identical).
"""
function solve_magnetic_atmosphere(T_eff::Float64, g_s::Float64,
                                    B::Float64, θ_B::Float64,
                                    gaunt::GauntTable;
                                    K::Int=50, M::Int=8, N::Int=200,
                                    max_iter::Int=30,
                                    tol::Float64=1e-4,
                                    verbose::Bool=true)::MagneticAtmosphereResult
    @assert T_eff > 0 && g_s > 0 && B >= 0
    @assert 0 <= θ_B <= π

    verbose && @printf("Magnetic RT: T_eff=%.2e K, g_s=%.2e, B=%.2e G, θ_B=%.1f°\n",
                        T_eff, g_s, B, rad2deg(θ_B))

    ν_grid = make_frequency_grid(T_eff, K)
    μ, w = gauss_legendre_half(M)

    # Build initial atmosphere column (Eddington T profile)
    y, T, ρ, P = _build_initial_column(T_eff, g_s, N, ν_grid, gaunt)

    verbose && @printf("  N=%d, y_max=%.2e, K=%d frequencies\n", N, y[end], K)

    # Compute opacities for both modes at all (depth, frequency) points
    # κ_j[i, k]: absorption opacity for mode j at depth i, frequency k
    # σ_j[i, k]: scattering opacity (mode-dependent for B>0)
    # For B>0: use mode_opacity(j, ν, θ_B, B, T, ρ)
    # For B=0: both modes have κ_ff + σ_T (recover non-magnetic)

    κ = zeros(N, K, 2)      # absorption opacity per mode
    k_total = zeros(N, K, 2) # total opacity per mode
    ρ_alb = zeros(N, K, 2)  # scattering albedo per mode
    τ = zeros(N, K, 2)      # optical depth per mode

    _compute_magnetic_opacities!(κ, k_total, ρ_alb, τ, y, T, ρ, ν_grid, B, θ_B, gaunt)

    converged = false
    max_dT = Inf
    n_iter = 0

    # Storage for Feautrier solutions (per mode)
    P_all = zeros(N, M, K, 2)  # Feautrier P[depth, angle, freq, mode]
    J = zeros(N, K, 2)          # mean intensity per mode
    f_ν = zeros(N, K, 2)        # Eddington factor per mode
    h_ν = zeros(N, K, 2)        # flux Eddington factor per mode

    for iter in 1:max_iter
        n_iter = iter

        # Solve Feautrier RT for each mode independently
        for j in 1:2
            _solve_feautrier_mode!(view(P_all, :, :, :, j),
                                   view(J, :, :, j),
                                   view(f_ν, :, :, j),
                                   view(h_ν, :, :, j),
                                   y, T, ν_grid, μ, w,
                                   view(k_total, :, :, j),
                                   view(ρ_alb, :, :, j))
        end

        # Compute flux diagnostic
        F_bol = _bolometric_flux_2mode(P_all, μ, w, ν_grid)
        flux_ratio = F_bol / (σ_SB * T_eff^4)

        # Rybicki temperature correction (radiative equilibrium)
        ΔT = _rybicki_two_mode(N, K, ν_grid, y, T, κ, k_total, ρ_alb, f_ν, h_ν, J)

        # Uniform flux correction: gently scale ALL temperatures toward
        # the correct emergent flux.  ΔT_flux/T = -α(F/F_target - 1)/4
        # α is small (0.1) so it doesn't destabilize Rybicki.
        if abs(flux_ratio - 1.0) > 0.02
            α_flux = 0.1
            for i in 1:N
                ΔT[i] += -α_flux * T[i] * (flux_ratio - 1.0) / (4.0 * flux_ratio)
            end
        end

        max_dT = maximum(abs.(ΔT ./ T))
        if max_dT > 0.3
            ΔT .*= 0.3 / max_dT
            max_dT = 0.3
        end

        verbose && @printf("  iter %2d: max|ΔT/T|=%.2e, F/σT⁴=%.4f\n",
                            iter, max_dT, flux_ratio)

        if max_dT < tol && abs(flux_ratio - 1.0) < 0.05
            converged = true
            verbose && @printf("  CONVERGED at iteration %d (flux error %.1f%%)\n",
                                iter, abs(flux_ratio - 1.0) * 100)
            break
        end
        # Also converge if flux is good and ΔT oscillations are small
        if max_dT < 10*tol && abs(flux_ratio - 1.0) < 0.03
            converged = true
            verbose && @printf("  CONVERGED at iteration %d (flux error %.1f%%, ΔT/T=%.1e)\n",
                                iter, abs(flux_ratio - 1.0) * 100, max_dT)
            break
        end

        # Apply correction
        for i in 1:N
            T[i] = max(T[i] + ΔT[i], 0.1 * T_eff)
        end

        # Update density and opacities
        for i in 1:N
            ρ[i] = m_p * P[i] / (2.0 * k_B * T[i])
        end
        _compute_magnetic_opacities!(κ, k_total, ρ_alb, τ, y, T, ρ, ν_grid, B, θ_B, gaunt)
    end

    # Final Feautrier solve
    for j in 1:2
        _solve_feautrier_mode!(view(P_all, :, :, :, j),
                               view(J, :, :, j),
                               view(f_ν, :, :, j),
                               view(h_ν, :, :, j),
                               y, T, ν_grid, μ, w,
                               view(k_total, :, :, j),
                               view(ρ_alb, :, :, j))
    end

    # Emergent intensity: I_j(ν, μ) = 2P_j(surface, μ)
    I_emergent = zeros(K, M, 2)
    for k in 1:K, jj in 1:M, mode in 1:2
        I_emergent[k, jj, mode] = 2.0 * P_all[1, jj, k, mode]
    end

    verbose && @printf("  Final F/σT⁴=%.4f\n",
                        _bolometric_flux_2mode(P_all, μ, w, ν_grid) / (σ_SB * T_eff^4))

    return MagneticAtmosphereResult(T_eff, g_s, B, θ_B, converged, n_iter, max_dT,
                                     ν_grid, μ, I_emergent, copy(T), copy(y))
end

# --- Internal functions ---

"""Build initial atmosphere column using the non-magnetic infrastructure."""
function _build_initial_column(T_eff, g_s, N, ν_grid, gaunt)
    # Reuse the non-magnetic build_atmosphere for the initial structure
    # (proper Eddington T profile with actual Rosseland mean)
    col = AtmosphereStructure.build_atmosphere(T_eff, g_s, ν_grid, gaunt; N=N)
    return copy(col.y), copy(col.T), copy(col.ρ), copy(col.P)
end

"""Compute magnetic opacities for both modes at all (depth, freq) points."""
function _compute_magnetic_opacities!(κ, k_total, ρ_alb, τ, y, T, ρ, ν_grid, B, θ_B, gaunt)
    N = length(y)
    K = length(ν_grid)

    for i in 1:N, k in 1:K
        ν = ν_grid[k]
        if B > 0
            # Magnetic: mode-dependent opacities, properly separated
            for j in 1:2
                κ_abs = mode_absorption(j, ν, θ_B, B, T[i], ρ[i])
                σ_scat = mode_scattering(j, ν, θ_B, B, T[i], ρ[i])
                κ[i, k, j] = max(κ_abs, 1e-30)
                k_total[i, k, j] = κ[i, k, j] + σ_scat
                ρ_alb[i, k, j] = σ_scat / k_total[i, k, j]
            end
        else
            # B=0: both modes identical (recover non-magnetic)
            κ_ff = kappa_ff(ν, T[i], ρ[i], gaunt)
            σ_s = sigma_thomson()
            for j in 1:2
                κ[i, k, j] = κ_ff
                k_total[i, k, j] = κ_ff + σ_s
                ρ_alb[i, k, j] = σ_s / k_total[i, k, j]
            end
        end
    end

    # Compute optical depth from surface for each mode
    for j in 1:2, k in 1:K
        τ[1, k, j] = 0.0
        for i in 2:N
            dy = y[i] - y[i-1]
            τ[i, k, j] = τ[i-1, k, j] + 0.5 * (k_total[i, k, j] + k_total[i-1, k, j]) * dy
        end
    end
end

"""Solve Feautrier RT for a single polarization mode."""
function _solve_feautrier_mode!(P_out, J_out, f_out, h_out,
                                 y, T, ν_grid, μ, w,
                                 k_tot_mode, ρ_alb_mode)
    N = length(y)
    M = length(μ)
    K = length(ν_grid)

    for k in 1:K
        ν = ν_grid[k]

        # Build dtau for this mode
        dtau = zeros(N)
        for i in 2:N
            dy = y[i] - y[i-1]
            dtau[i] = 0.5 * (k_tot_mode[i, k] + k_tot_mode[i-1, k]) * dy
        end

        # Find effective depth (τ < 80)
        N_eff = N
        τ_cum = 0.0
        for i in 2:N
            τ_cum += dtau[i]
            if τ_cum > 80.0
                N_eff = i
                break
            end
        end
        N_eff = max(N_eff, 3)

        # Solve Feautrier for this frequency and mode
        P_k = _feautrier_single(N_eff, M, dtau, μ, w, ν, T,
                                 view(ρ_alb_mode, :, k))

        # Store results (pad beyond N_eff with B_ν/2 per mode)
        for j in 1:M
            for i in 1:N_eff
                P_out[i, j, k] = P_k[i, j]
            end
            for i in N_eff+1:N
                P_out[i, j, k] = 0.5 * planck_Bnu(ν, T[i])
            end
        end

        # Compute J, f, h
        for i in 1:N
            J_out[i, k] = sum(P_out[i, j, k] * w[j] for j in 1:M)
            if J_out[i, k] > 0
                f_out[i, k] = sum(μ[j]^2 * P_out[i, j, k] * w[j] for j in 1:M) / J_out[i, k]
                h_out[i, k] = sum(μ[j] * P_out[i, j, k] * w[j] for j in 1:M) / J_out[i, k]
            else
                f_out[i, k] = 1.0/3.0
                h_out[i, k] = 0.0
            end
        end
    end
end

"""Solve Feautrier at single frequency for ONE polarization mode.
Per-mode Planck source is B_ν/2 (unpolarized emission split equally)."""
function _feautrier_single(N_eff, M, dtau, μ, w, ν, T, ρ_alb)
    B_tilde = [zeros(M, M) for _ in 1:N_eff]
    Q_tilde = [zeros(M) for _ in 1:N_eff]
    C_store = [zeros(M, M) for _ in 1:N_eff]

    for i in 1:N_eff
        A = zeros(M, M)
        B = zeros(M, M)
        C = zeros(M, M)
        Q = zeros(M)

        Bν_half = 0.5 * planck_Bnu(ν, T[i])  # per-mode blackbody = B_ν/2
        ρ = ρ_alb[i]

        if i == 1
            # Surface: pure radiation BC (no source, no scattering)
            dt_next = dtau[2]
            for j in 1:M
                δ = dt_next / μ[j]
                C[j, j] = -1.0 / δ
                B[j, j] = 1.0 + 1.0 / δ
            end
        elseif i == N_eff
            # Bottom: P_j = B_ν/2 (LTE per mode)
            for j in 1:M
                B[j, j] = 1.0
            end
            Q .= Bν_half
            B_tilde[i] = B
            Q_tilde[i] = Q
            C_store[i] = C
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

        # Scattering and source (interior only, not surface)
        if i > 1 && i < N_eff
            for j in 1:M, jp in 1:M
                B[j, jp] -= 2.0 * ρ * 0.5 * w[jp]
            end
            Q .= (1.0 - ρ) * Bν_half  # per-mode source = (1-ρ) B_ν/2
        end

        C_store[i] = copy(C)

        # Forward elimination
        if i == 1
            B_tilde[i] = copy(B)
            Q_tilde[i] = copy(Q)
        else
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
    for i in 1:N_eff, j in 1:M
        P[i, j] = max(P[i, j], 0.0)
    end
    return P
end

"""
Rybicki temperature correction for two-mode atmosphere.
Each mode j contributes to the correction through its own T_k, U_k, K_k matrices.
The V_k weights sum over both modes.

ΔT enforces: Σ_j ∫ κ_j (J_j - B_ν/2) dν = 0 at each depth.
"""
function _rybicki_two_mode(N, K, ν, y, T, κ, k_total, ρ_alb, f_ν, h_ν, J)
    # Frequency quadrature weights
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

    # Denominator and B̄ — sum over BOTH modes with per-mode Planck = B_ν/2
    # denom_i = Σ_j Σ_k (dB_ν/dT)/2 × κ_j × b_k
    # B̄_i = Σ_j Σ_k (B_ν/2) × κ_j × b_k / denom_i
    denom = zeros(N)
    B_bar = zeros(N)
    for i in 1:N
        for k in 1:K
            dBdT_half = 0.5 * dBnu_dT(ν[k], T[i])
            Bν_half = 0.5 * planck_Bnu(ν[k], T[i])
            for j in 1:2
                denom[i] += dBdT_half * κ[i, k, j] * b[k]
                B_bar[i] += Bν_half * κ[i, k, j] * b[k]
            end
        end
        if denom[i] > 0
            B_bar[i] /= denom[i]
        end
    end

    # Accumulate W and rhs over both modes and all frequencies
    W = zeros(N, N)
    rhs = zeros(N)

    for j in 1:2
        for k in 1:K
            T_diag = zeros(N)
            T_sub = zeros(N-1)
            T_sup = zeros(N-1)
            U_k = zeros(N)
            K_k = zeros(N)

            _build_rybicki_system_magnetic!(T_diag, T_sub, T_sup, U_k, K_k,
                                            N, k, j, ν, y, T, κ, k_total, ρ_alb,
                                            f_ν, h_ν, denom, b, B_bar)

            # LU solve and accumulate
            lu_diag, lu_sup, lu_rhs_K = _tridiag_lu_forward(T_sub, T_diag, T_sup, K_k)
            x_k = _tridiag_lu_back(lu_diag, lu_sup, lu_rhs_K)

            # V_k: weight for this mode and frequency
            V_k = zeros(N)
            for i in 1:N
                V_k[i] = denom[i] > 0 ? κ[i, k, j] * b[k] / denom[i] : 0.0
            end

            # Accumulate rhs
            for i in 1:N
                rhs[i] += V_k[i] * x_k[i]
            end

            # Accumulate W
            for jcol in 1:N
                abs(U_k[jcol]) < 1e-30 && continue
                e_j = zeros(N); e_j[jcol] = 1.0
                _, _, lu_rhs_j = _tridiag_lu_forward(T_sub, T_diag, T_sup, e_j)
                z_j = _tridiag_lu_back(lu_diag, lu_sup, lu_rhs_j)
                for i in 1:N
                    W[i, jcol] += V_k[i] * U_k[jcol] * z_j[i]
                end
            end
        end
    end

    # Solve (W + I) J̄ = rhs
    for i in 1:N
        W[i, i] += 1.0
    end
    J_bar = W \ rhs

    return J_bar .- B_bar
end

"""Build Rybicki system for mode j, frequency k.
Per-mode Planck source is B_ν/2. B̄ is the GLOBAL average (both modes)."""
function _build_rybicki_system_magnetic!(T_diag, T_sub, T_sup, U_k, K_k,
                                          N, k, j, ν, y, T, κ, k_total, ρ_alb,
                                          f_ν, h_ν, denom, b, B_bar)
    for i in 2:N-1
        Δ1 = (k_total[i, k, j] + k_total[i-1, k, j]) * (y[i] - y[i-1])
        Δ2 = (k_total[i, k, j] + k_total[i+1, k, j]) * (y[i+1] - y[i])
        π_sum = Δ1 + Δ2

        ρ_k = ρ_alb[i, k, j]
        f_k = f_ν[i, k, j]

        T_sub[i-1] = -8.0 * f_ν[i-1, k, j] / (π_sum * Δ1)
        T_sup[i] = -8.0 * f_ν[i+1, k, j] / (π_sum * Δ2)
        T_diag[i] = ((1.0 - ρ_k) / f_k + 8.0 / (π_sum * Δ1) + 8.0 / (π_sum * Δ2)) * f_k

        # Per-mode: dB_j/dT = (dB_ν/dT)/2, B_j = B_ν/2
        dBdT_half = 0.5 * dBnu_dT(ν[k], T[i])
        U_k[i] = -dBdT_half * (1.0 - ρ_k)

        Bν_half = 0.5 * planck_Bnu(ν[k], T[i])
        # K_k uses the GLOBAL B̄ (precomputed from both modes), not per-mode
        K_k[i] = (Bν_half - dBdT_half * B_bar[i]) * (1.0 - ρ_k)
    end

    # Surface: pure radiation BC
    Δ_surf = 0.5 * (k_total[2, k, j] + k_total[1, k, j]) * (y[2] - y[1])
    T_diag[1] = h_ν[1, k, j] + f_ν[1, k, j] / Δ_surf
    T_sup[1] = -f_ν[2, k, j] / Δ_surf
    U_k[1] = 0.0
    K_k[1] = 0.0

    # Bottom: J_j = B_j = B_ν/2 (per-mode LTE)
    T_diag[N] = 1.0
    if N > 1; T_sub[N-1] = 0.0; end
    U_k[N] = 0.0
    K_k[N] = 0.5 * planck_Bnu(ν[k], T[N])
end

"""Bolometric flux from two-mode Feautrier solution."""
function _bolometric_flux_2mode(P_all, μ, w, ν_grid)
    K = length(ν_grid)
    M = length(μ)
    F = 0.0
    for k in 1:K-1
        dν = ν_grid[k+1] - ν_grid[k]
        for j in 1:M, mode in 1:2
            F += 2π * μ[j] * 2.0 * P_all[1, j, k, mode] * w[j] * dν
        end
    end
    return F
end

# --- Tridiagonal solver (copied from TemperatureCorrection for self-containedness) ---

function _tridiag_lu_forward(sub, diag, sup, d)
    N = length(diag)
    lu_diag = copy(diag); lu_sup = copy(sup); lu_rhs = copy(d)
    for i in 2:N
        abs(lu_diag[i-1]) < 1e-30 && continue
        m = sub[i-1] / lu_diag[i-1]
        lu_diag[i] -= m * lu_sup[i-1]
        lu_rhs[i] -= m * lu_rhs[i-1]
    end
    return lu_diag, lu_sup, lu_rhs
end

function _tridiag_lu_back(lu_diag, lu_sup, lu_rhs)
    N = length(lu_diag)
    x = zeros(N)
    abs(lu_diag[N]) > 1e-30 && (x[N] = lu_rhs[N] / lu_diag[N])
    for i in N-1:-1:1
        abs(lu_diag[i]) > 1e-30 && (x[i] = (lu_rhs[i] - lu_sup[i] * x[i+1]) / lu_diag[i])
    end
    return x
end

end # module
