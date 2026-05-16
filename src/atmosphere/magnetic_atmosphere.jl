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
using ..PhysicalConstants: σ_SB, k_B, h, m_p, m_e, m_H
using ..Tridiag: tridiag_lu_forward!, tridiag_lu_back!
using ..GauntFactor: GauntTable
using ..AtmosphereStructure: AtmosphereStructure, make_frequency_grid, density_from_PT
using ..BlackbodyAtmosphere: planck_Bnu
using ..HydrogenOpacity: kappa_ff, sigma_thomson, dBnu_dT
using ..MagneticModes: mode_absorption, mode_scattering, mode_opacity_split, effective_opacity
using ..FeautrierSolver: gauss_legendre_half
using ..DielectricTensor: vacuum_resonance_density, vacuum_resonance_pjump
using ..SolverDefaults: TAU_DIFFUSION, TAU_EFFECTIVE_MIN, T_SURFACE_FRAC_MCPHAC,
                         Y_MAX_SEMIINFINITE, FLUX_DAMPING_DEFAULT,
                         FLUX_SCALE_CLAMP_LO, FLUX_SCALE_CLAMP_HI,
                         DELTA_T_OVER_T_CAP, OPACITY_FLOOR

export solve_magnetic_atmosphere, MagneticAtmosphereResult, magnetic_emergent_spectrum

"""
Result of a converged magnetic atmosphere calculation.

Bead E15: `column` snapshots the full converged magnetic column as a
NamedTuple `(y, T, ρ, P, κ, k_total, ρ_alb, τ)`. The magnetic solver
operates on bare arrays (no dedicated struct), so a NamedTuple is the
lowest-friction container. The mode axis is the last axis of κ, k_total,
ρ_alb, τ (shape `N × nν × 2`). `T_profile` and `y_grid` are retained for
backward compatibility and are equivalent to `column.T` and `column.y`.
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
    I_emergent::Array{Float64, 3}   # nν × M × 2 (freq, angle, mode)
    T_profile::Vector{Float64}      # == column.T (E15: retained for compat)
    y_grid::Vector{Float64}         # == column.y (E15: retained for compat)
    column::NamedTuple              # full converged column (E15)
    # Vacuum-resonance mode-conversion diagnostics (SPW09 Eq. 16-17).
    # `P_jump[k]` is the non-conversion probability at frequency `k`;
    # `NaN` means no resonance layer was found inside the depth grid
    # (no swap applied). `pjump_enabled` records whether the post-
    # processing was active.
    pjump_enabled::Bool
    P_jump::Vector{Float64}         # length nν, NaN where no resonance in grid
end

"""
    solve_magnetic_atmosphere(T_eff, g_s, B, θ_B, gaunt; ...) → MagneticAtmosphereResult
    solve_magnetic_atmosphere(T_eff, g_s, B, θ_B; ...) → MagneticAtmosphereResult

Solve a plane-parallel magnetized NS atmosphere with two-mode RT.
Each mode j (1=X, 2=O) has its own opacity κ_j(ν,T,ρ,B,θ_B).
Temperature correction enforces Σ_j ∫ κ_j(J_j - B_ν/2) dν = 0.

For B=0, this should recover the non-magnetic result (both modes identical).

Bead D10: `gaunt` is now `Union{GauntTable, Nothing}`; a Gaunt table is only
required when `B == 0` (the non-magnetic limit uses `kappa_ff`, which needs
the table). For pure-magnetic runs (`B > 0`) the mode opacities use their
own Coulomb logarithm and never touch the table, so callers may either pass
`nothing` for the 5-arg form or use the 4-arg overload, which forwards
`gaunt = nothing`. The B=0 branch in `_compute_magnetic_opacities!` errors
explicitly if `gaunt === nothing`.
"""
function solve_magnetic_atmosphere(T_eff::Float64, g_s::Float64,
                                    B::Float64, θ_B::Float64,
                                    gaunt::Union{GauntTable, Nothing};
                                    nν::Int=50, M::Int=8, N::Int=200,
                                    max_iter::Int=30,
                                    tol::Float64=1e-4,
                                    flux_tol::Float64=1e-2,
                                    flux_damping::Float64=FLUX_DAMPING_DEFAULT,
                                    apply_pjump::Bool=false,
                                    verbose::Bool=true)::MagneticAtmosphereResult
    @assert T_eff > 0 && g_s > 0 && B >= 0
    @assert 0 <= θ_B <= π

    verbose && @printf("Magnetic RT: T_eff=%.2e K, g_s=%.2e, B=%.2e G, θ_B=%.1f°\n",
                        T_eff, g_s, B, rad2deg(θ_B))

    ν_grid = make_frequency_grid(T_eff, nν)
    μ, w = gauss_legendre_half(M)

    # Build initial atmosphere column (Eddington T profile)
    y, T, ρ, P = _build_initial_column(T_eff, g_s, N, ν_grid, gaunt, B, θ_B)

    verbose && @printf("  N=%d, y_max=%.2e, nν=%d frequencies\n", N, y[end], nν)

    # Compute opacities for both modes at all (depth, frequency) points
    # κ_j[i, k]: absorption opacity for mode j at depth i, frequency k
    # σ_j[i, k]: scattering opacity (mode-dependent for B>0)
    # For B>0: split true absorption from scattering extinction.
    # For B=0: both modes have κ_ff + σ_T (recover non-magnetic)

    κ = zeros(N, nν, 2)      # absorption opacity per mode
    k_total = zeros(N, nν, 2) # total opacity per mode
    ρ_alb = zeros(N, nν, 2)  # scattering albedo per mode
    τ = zeros(N, nν, 2)      # optical depth per mode

    _compute_magnetic_opacities!(κ, k_total, ρ_alb, τ, y, T, ρ, ν_grid, B, θ_B, gaunt)

    converged = false
    max_dT = Inf
    n_iter = 0
    flux_clamp_warned = false  # emit clamp warning at most once per call

    # Storage for Feautrier solutions (per mode)
    P_all = zeros(N, M, nν, 2)  # Feautrier P[depth, angle, freq, mode]
    J = zeros(N, nν, 2)          # mean intensity per mode
    f_ν = zeros(N, nν, 2)        # Eddington factor per mode
    h_ν = zeros(N, nν, 2)        # flux Eddington factor per mode

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

        # Compute two-mode Rybicki temperature correction
        ΔT = _rybicki_two_mode(N, nν, ν_grid, y, T, κ, k_total, ρ_alb, f_ν, h_ν, J)

        # Ad-hoc grey-body flux scaling (NOT the Suleimanov-Potekhin-Werner 2009
        # Avrett-Krook correction, despite the name we've been using). Motivation:
        # T_eff ∝ F^{1/4}, so a global flux deficit suggests a uniform temperature
        # bump. Heuristic — clamped to [FLUX_SCALE_CLAMP_LO, FLUX_SCALE_CLAMP_HI]
        # to keep iterations from overshooting; works for B<~10^13 G but may
        # explain the unresolved θ_B=45° mismatch noted in notes/HANDOFF.md
        # "Immediate Next Task". Implementing actual SPW09 §2 (depth-resolved
        # Avrett-Krook from their Eq. 19-22 + Kurucz 1970 surface correction)
        # is deferred — see notes/decisions.md D10.
        F_bol = _bolometric_flux_2mode(P_all, μ, w, ν_grid)
        flux_ratio = F_bol / (σ_SB * T_eff^4)
        if B > 0 && isfinite(flux_ratio) && flux_ratio > 0
            raw_flux_scale = flux_ratio^(-0.25)
            raw_scale = 1.0 + flux_damping * (raw_flux_scale - 1.0)
            if !flux_clamp_warned && !(0.5 < raw_scale < 2.0)
                @warn "magnetic atmosphere flux correction clamped" iter raw_flux_scale flux_ratio raw_scale flux_damping
                flux_clamp_warned = true
            end
            flux_scale = clamp(raw_scale, FLUX_SCALE_CLAMP_LO, FLUX_SCALE_CLAMP_HI)
            ΔT .+= (flux_scale - 1.0) .* T
        end

        max_dT = maximum(abs.(ΔT ./ T))
        if max_dT > DELTA_T_OVER_T_CAP
            ΔT .*= DELTA_T_OVER_T_CAP / max_dT
            max_dT = DELTA_T_OVER_T_CAP
        end

        verbose && @printf("  iter %2d: max|ΔT/T|=%.2e, F/σT⁴=%.4f\n",
                            iter, max_dT, flux_ratio)

        if max_dT < tol && abs(flux_ratio - 1.0) < flux_tol
            converged = true
            verbose && @printf("  CONVERGED at iteration %d\n", iter)
            break
        end

        # Apply correction
        for i in 1:N
            T[i] = max(T[i] + ΔT[i], 0.1 * T_eff)
        end

        # Update density and opacities
        for i in 1:N
            ρ[i] = density_from_PT(P[i], T[i])
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
    I_emergent = zeros(nν, M, 2)
    for k in 1:nν, jj in 1:M, mode in 1:2
        I_emergent[k, jj, mode] = 2.0 * P_all[1, jj, k, mode]
    end

    # Vacuum-resonance mode-conversion post-processing per SPW09 Eqs. (16)-(17)
    # / van Adelsberg & Lai (2006) Eqs. (4), (33)-(34). The Feautrier solve
    # above treated the two modes as decoupled along the entire column. SPW09's
    # actual prescription (their §2 around Eq. 16) is to compute that
    # independent solve and then apply the X<->O mixing at the resonance
    # density:
    #     I_X^obs = P_jump · I_X + (1 - P_jump) · I_O
    #     I_O^obs = (1 - P_jump) · I_X + P_jump · I_O
    # `P_jump[k]` is set to NaN for frequencies whose resonance density falls
    # outside the depth grid (no resonance crossing in the modelled column ⇒
    # no swap applied). Photons emitted from layers above the resonance are
    # not represented in this post-processing approximation (would require
    # depth-resolved swap during the RT integration; deferred — see TODO).
    P_jump_vec = fill(NaN, nν)
    if apply_pjump && B > 0
        I_temp = similar(I_emergent)
        I_temp .= I_emergent
        for k in 1:nν
            ω = 2π * ν_grid[k]
            ρ_V = vacuum_resonance_density(ω, B, θ_B)
            # Find the resonance index: ρ(y) is monotonically increasing with
            # y in our log-spaced column. We require ρ[1] < ρ_V < ρ[N_eff_safe]
            # so the surface (ρ[1]) is below the resonance and the resonance
            # lies inside the integration domain.
            if !(ρ_V > ρ[1] && ρ_V < ρ[end])
                continue   # no resonance crossing in column ⇒ P_jump stays NaN
            end
            i_V = searchsortedfirst(ρ, ρ_V)
            i_V = clamp(i_V, 2, N)
            # Linear gradient dρ/dy across the bracketing pair.
            dρ_dy = (ρ[i_V] - ρ[i_V-1]) / (y[i_V] - y[i_V-1])
            if !(dρ_dy > 0 && isfinite(dρ_dy))
                continue
            end
            P = vacuum_resonance_pjump(ω, B, θ_B, T[i_V], dρ_dy)
            P_jump_vec[k] = P
            for jj in 1:M
                I_X = I_temp[k, jj, 1]
                I_O = I_temp[k, jj, 2]
                I_emergent[k, jj, 1] = P * I_X + (1 - P) * I_O
                I_emergent[k, jj, 2] = (1 - P) * I_X + P * I_O
            end
        end
        if verbose
            n_active = count(isfinite, P_jump_vec)
            if n_active > 0
                P_finite = filter(isfinite, P_jump_vec)
                @printf("  P_jump: applied at %d/%d freqs, range [%.3f, %.3f], mean=%.3f\n",
                        n_active, nν, minimum(P_finite), maximum(P_finite),
                        sum(P_finite)/length(P_finite))
            else
                @printf("  P_jump: no frequency had a vacuum resonance inside the column\n")
            end
        end
    end

    verbose && @printf("  Final F/σT⁴=%.4f\n",
                        _bolometric_flux_2mode(P_all, μ, w, ν_grid) / (σ_SB * T_eff^4))

    # Bead E15: snapshot the full converged column (bare arrays, no struct
    # in the magnetic solver). copy.() detaches from the mutable scratch
    # arrays that the iteration loop has been overwriting in place.
    column = (y       = copy(y),
              T       = copy(T),
              ρ       = copy(ρ),
              P       = copy(P),
              κ       = copy(κ),
              k_total = copy(k_total),
              ρ_alb   = copy(ρ_alb),
              τ       = copy(τ))

    return MagneticAtmosphereResult(T_eff, g_s, B, θ_B, converged, n_iter, max_dT,
                                     ν_grid, μ, I_emergent,
                                     column.T, column.y, column,
                                     apply_pjump, P_jump_vec)
end

# Bead D10: 4-arg convenience overload for pure-magnetic runs (B > 0). Forwards
# `gaunt = nothing` so the B>0 branch (which doesn't touch the Gaunt table) can
# run without a 75 kB free-free table being loaded from disk. Julia method
# dispatch picks this overload when the caller supplies four positional args
# and the 5-arg method when `gaunt` is supplied explicitly.
solve_magnetic_atmosphere(T_eff::Float64, g_s::Float64,
                          B::Float64, θ_B::Float64; kwargs...) =
    solve_magnetic_atmosphere(T_eff, g_s, B, θ_B, nothing; kwargs...)

"""
    magnetic_emergent_spectrum(result::MagneticAtmosphereResult,
                                ν_obs::Vector{Float64},
                                cos_θe::Float64;
                                polarized::Bool=false) → Vector{Float64}
                                                          or Matrix{Float64}

Interpolate the magnetic atmosphere's emergent intensity onto observation
frequencies and a single angle (cos θe = μ).

When `polarized=false` (default): returns a `length(ν_obs)`-length vector
of mode-summed intensity `I_X + I_O`. This is what intensity-only
renderers want.

When `polarized=true`: returns a `length(ν_obs) × 2` matrix (column 1 = X,
column 2 = O). Use this for future polarization-aware downstream code
(e.g. Kerr ray tracing with parallel-transported basis).

For both cases, frequency interpolation is log-log linear (matching
`RTAtmosphere.rt_emergent_spectrum`); angle interpolation is linear in μ.
Below the solver's `ν_grid[1]`, intensity is extended via Rayleigh-Jeans
`(ν/ν₁)²` scaling (matches the non-magnetic spectrum). Above
`ν_grid[end]`, the intensity is clamped to zero (Wien tail decays
faster than any local extrapolation can capture). In μ, values outside
`[μ_grid[1], μ_grid[end]]` are clamped to the boundary value (no
extrapolation).

This is the magnetic counterpart of `RTAtmosphere.rt_emergent_spectrum`,
which only handles the non-magnetic `nν × M` case (see
`src/atmosphere/rt_atmosphere.jl`). The nν × M × 2 magnetic case used to
be re-implemented inline in every script
(see `scripts/render_rxj1856_visible_magnetic_atmosphere.jl`); promoting
it here makes it testable and reusable.
"""
function magnetic_emergent_spectrum(result::MagneticAtmosphereResult,
                                     ν_obs::AbstractVector{Float64},
                                     cos_θe::Float64;
                                     polarized::Bool=false)
    @assert 0 < cos_θe <= 1.0

    μ = result.μ_grid
    Mμ = length(μ)

    # μ bracket (mirrors rt_emergent_spectrum)
    j_lo = 1
    for j in 1:Mμ-1
        if μ[j] <= cos_θe <= μ[j+1]
            j_lo = j
            break
        end
    end
    if cos_θe < μ[1]; j_lo = 1; end
    if cos_θe > μ[end]; j_lo = Mμ - 1; end
    j_hi = min(j_lo + 1, Mμ)
    t_μ = j_lo < Mμ ? clamp((cos_θe - μ[j_lo]) / (μ[j_hi] - μ[j_lo]), 0.0, 1.0) : 0.0

    ν_in = result.ν_grid
    nν = length(ν_in)
    N_out = length(ν_obs)

    # Per-mode storage (always compute both, sum only if !polarized)
    out_per_mode = zeros(N_out, 2)

    for mode in 1:2
        for i in 1:N_out
            ν = ν_obs[i]
            if ν <= ν_in[1]
                I_lo_b = result.I_emergent[1, j_lo, mode]
                I_hi_b = result.I_emergent[1, j_hi, mode]
                I_at = (1.0 - t_μ) * I_lo_b + t_μ * I_hi_b
                out_per_mode[i, mode] = I_at * (ν / ν_in[1])^2  # Rayleigh-Jeans
            elseif ν >= ν_in[end]
                out_per_mode[i, mode] = 0.0
            else
                k_lo = clamp(searchsortedlast(ν_in, ν), 1, nν - 1)
                k_hi = k_lo + 1
                t_ν = (log(ν) - log(ν_in[k_lo])) / (log(ν_in[k_hi]) - log(ν_in[k_lo]))

                I_ll = result.I_emergent[k_lo, j_lo, mode]
                I_hl = result.I_emergent[k_hi, j_lo, mode]
                I_lh = result.I_emergent[k_lo, j_hi, mode]
                I_hh = result.I_emergent[k_hi, j_hi, mode]

                I_lo = (I_ll > 0 && I_hl > 0) ?
                       exp((1.0 - t_ν) * log(I_ll) + t_ν * log(I_hl)) :
                       (1.0 - t_ν) * I_ll + t_ν * I_hl
                I_hi = (I_lh > 0 && I_hh > 0) ?
                       exp((1.0 - t_ν) * log(I_lh) + t_ν * log(I_hh)) :
                       (1.0 - t_ν) * I_lh + t_ν * I_hh
                out_per_mode[i, mode] = (1.0 - t_μ) * I_lo + t_μ * I_hi
            end
        end
    end

    if polarized
        return out_per_mode
    else
        return out_per_mode[:, 1] .+ out_per_mode[:, 2]
    end
end

# --- Internal functions ---

"""Build initial atmosphere column.

For B = 0 we reuse the non-magnetic `build_atmosphere` (which needs the Gaunt
table for the Rosseland-mean Eddington T profile). For B > 0 the magnetic
mode opacities don't use the Gaunt table at all, so `gaunt` may be `nothing`:
we build the log-spaced y grid and hydrostatic P = g_s * y directly and use
`_magnetic_eddington_temperature` for the initial T (it doesn't touch the
Gaunt table either). Bead D10.
"""
function _build_initial_column(T_eff, g_s, N, ν_grid, gaunt, B=0.0, θ_B=0.0)
    y_max = 1e2
    max_y = Y_MAX_SEMIINFINITE  # Suleimanov+ 2009 semi-infinite atmosphere depth scale.

    if B <= 0
        # Non-magnetic: defer to build_atmosphere (requires a Gaunt table).
        gaunt === nothing && error(
            "Non-magnetic (B=0) initial column requires a Gaunt table; " *
            "pass `gaunt = load_gaunt_table(\"refs/code/McPHAC/gffgu.dat\")`.")
        col = AtmosphereStructure.build_atmosphere(T_eff, g_s, ν_grid, gaunt;
                                                   N=N, y_max=y_max)
        return copy(col.y), copy(col.T), copy(col.ρ), copy(col.P)
    end

    # B > 0 path. Grow the column until every mode/frequency reaches the
    # diffusion-depth cutoff (Feautrier needs τ ≳ TAU_DIFFUSION at the bottom).
    # When a Gaunt table is available we reuse `build_atmosphere` (which has
    # its own non-magnetic τ ≥ TAU_DIFFUSION auto-grow at the highest frequency); the
    # resulting (y, P) seeds the magnetic τ-growth loop with a deeper column,
    # which empirically helps convergence. When `gaunt === nothing` (pure-
    # magnetic call) we construct (y, P) directly from log-spacing and
    # hydrostatic equilibrium — no Gaunt table is needed because we discard
    # `build_atmosphere`'s T profile and use `_magnetic_eddington_temperature`
    # instead.
    while true
        if gaunt === nothing
            y_min = 1e-9
            y_loc = [10.0^logy for logy in range(log10(y_min), log10(y_max), length=N)]
            P_loc = g_s .* y_loc
        else
            col = AtmosphereStructure.build_atmosphere(T_eff, g_s, ν_grid, gaunt;
                                                       N=N, y_max=y_max)
            y_loc = copy(col.y)
            P_loc = copy(col.P)
        end

        T_mag = _magnetic_eddington_temperature(y_loc, T_eff, g_s, ν_grid, B, θ_B)
        ρ_mag = [density_from_PT(P_loc[i], T_mag[i]) for i in 1:N]

        nν = length(ν_grid)
        κ = zeros(N, nν, 2)
        k_total = zeros(N, nν, 2)
        ρ_alb = zeros(N, nν, 2)
        τ = zeros(N, nν, 2)
        _compute_magnetic_opacities!(κ, k_total, ρ_alb, τ,
                                     y_loc, T_mag, ρ_mag, ν_grid, B, θ_B, gaunt)

        τ_min = minimum(τ[end, :, :])
        if τ_min >= TAU_DIFFUSION || y_loc[end] >= max_y
            return y_loc, T_mag, ρ_mag, P_loc
        end

        scale = (TAU_DIFFUSION / max(τ_min, 1.0))^1.2
        y_max = min(max(y_loc[end] * max(scale, 1.5), y_loc[end] * 1.5), max_y)
    end
end

"""Magnetic grey atmosphere initial guess at the local field-normal angle."""
function _magnetic_eddington_temperature(y, T_eff, g_s, ν_grid, B, θ_B)
    N = length(y)
    T = zeros(N)
    # Initial guess for the surface temperature. The T_SURFACE_FRAC_MCPHAC
    # value (0.265) is the McPHAC non-magnetic surface limit (H12 §3.1);
    # justified as a starting point for the iteration, not as a magnetic-
    # atmosphere boundary value. See notes/decisions.md D9 for the rationale
    # (and why we don't yet use the grey-atmosphere limit (1/2)^{1/4} T_eff
    # = 0.841 T_eff).
    T[1] = T_SURFACE_FRAC_MCPHAC * T_eff

    P1 = g_s * y[1]
    ρ1 = density_from_PT(P1, T[1])
    k_R_prev = _rosseland_directional_magnetic(T[1], ρ1, ν_grid, B, θ_B)
    τ_R = 0.0

    for i in 2:N
        dy = y[i] - y[i-1]
        P_i = g_s * y[i]

        # Fixed-point iteration because the Rosseland opacity depends on T and ρ.
        T_i = T[i-1]
        k_R_i = k_R_prev
        for _ in 1:4
            τ_trial = τ_R + 0.5 * (k_R_prev + k_R_i) * dy
            T_i = (0.75 * T_eff^4 * (τ_trial + 2.0/3.0))^0.25
            ρ_i = density_from_PT(P_i, T_i)
            k_R_i = _rosseland_directional_magnetic(T_i, ρ_i, ν_grid, B, θ_B)
        end

        τ_R += 0.5 * (k_R_prev + k_R_i) * dy
        T[i] = (0.75 * T_eff^4 * (τ_R + 2.0/3.0))^0.25
        k_R_prev = k_R_i
    end

    return T
end

function _rosseland_directional_magnetic(T, ρ, ν_grid, B, θ_B)
    numerator = 0.0
    denominator = 0.0

    for k in 1:length(ν_grid)-1
        ν = 0.5 * (ν_grid[k] + ν_grid[k+1])
        dν = ν_grid[k+1] - ν_grid[k]
        weight = dBnu_dT(ν, T)
        κ_eff = effective_opacity(ν, θ_B, B, T, ρ)
        if κ_eff > 0
            numerator += weight / κ_eff * dν
            denominator += weight * dν
        end
    end

    return denominator > 0 && numerator > 0 ? denominator / numerator : 0.0
end

"""Compute magnetic opacities for both modes at all (depth, freq) points."""
function _compute_magnetic_opacities!(κ, k_total, ρ_alb, τ, y, T, ρ, ν_grid, B, θ_B, gaunt)
    N = length(y)
    nν = length(ν_grid)

    for i in 1:N, k in 1:nν
        ν = ν_grid[k]
        if B > 0
            # Magnetic: mode-dependent opacities
            for j in 1:2
                κ_abs_raw, σ_scat_raw = mode_opacity_split(j, ν, θ_B, B, T[i], ρ[i])
                κ_abs = max(κ_abs_raw, OPACITY_FLOOR)
                σ_scat = max(σ_scat_raw, 0.0)
                total = max(κ_abs + σ_scat, OPACITY_FLOOR)

                κ[i, k, j] = κ_abs
                k_total[i, k, j] = total
                ρ_alb[i, k, j] = σ_scat / total
            end
        else
            # B=0: both modes identical (recover non-magnetic).
            # Bead D10: the Gaunt table is only consulted on this branch
            # (kappa_ff needs it). Pure-magnetic callers can pass
            # `gaunt = nothing`, but a B=0 run with no table is an error.
            gaunt === nothing && error(
                "Non-magnetic (B=0) path requires a Gaunt table; " *
                "pass `gaunt = load_gaunt_table(\"refs/code/McPHAC/gffgu.dat\")` " *
                "or use a B>0 magnetic configuration.")
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
    for j in 1:2, k in 1:nν
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
    nν = length(ν_grid)

    # Preallocate scratch (reused across frequencies). Sized to N (max N_eff).
    B_tilde = [zeros(M, M) for _ in 1:N]
    Q_tilde = [zeros(M)    for _ in 1:N]
    C_store = [zeros(M, M) for _ in 1:N]
    A_work  = zeros(M, M)
    B_work  = zeros(M, M)
    C_work  = zeros(M, M)
    Q_work  = zeros(M)
    workmat = zeros(M, M)   # factorization scratch for B_tilde[i-1]
    tmp_mat = zeros(M, M)   # B_tilde[i-1] \ C_store[i-1]
    tmp_vec = zeros(M)      # B_tilde[i-1] \ Q_tilde[i-1]
    P_buf   = zeros(N, M)   # Feautrier P solution per frequency
    p_col   = zeros(M)      # column buffer for back-substitution
    rhs_buf = zeros(M)      # rhs buffer for back-substitution
    dtau    = zeros(N)

    for k in 1:nν
        ν = ν_grid[k]

        # Build dtau for this mode (reuse buffer)
        dtau[1] = 0.0
        for i in 2:N
            dy = y[i] - y[i-1]
            dtau[i] = 0.5 * (k_tot_mode[i, k] + k_tot_mode[i-1, k]) * dy
        end

        # Find effective depth.  In scattering-dominated layers, τ_total=80
        # can still be above the thermalization depth, so require effective
        # optical depth ∫sqrt(κ_abs k_total)dy to be safely large as well.
        N_eff = N
        τ_cum = 0.0
        τ_eff_cum = 0.0
        for i in 2:N
            τ_cum += dtau[i]
            κ_abs_prev = k_tot_mode[i-1, k] * (1.0 - ρ_alb_mode[i-1, k])
            κ_abs_curr = k_tot_mode[i, k] * (1.0 - ρ_alb_mode[i, k])
            κ_abs_mid = max(0.5 * (κ_abs_prev + κ_abs_curr), 0.0)
            k_tot_mid = max(0.5 * (k_tot_mode[i-1, k] + k_tot_mode[i, k]), 0.0)
            dy = y[i] - y[i-1]
            τ_eff_cum += sqrt(κ_abs_mid * k_tot_mid) * dy
            if τ_cum > TAU_DIFFUSION && τ_eff_cum > TAU_EFFECTIVE_MIN
                N_eff = i
                break
            end
        end
        N_eff = max(N_eff, 3)

        # Solve Feautrier for this frequency and mode (writes into P_buf)
        _feautrier_single!(P_buf, N_eff, M, dtau, μ, w, ν, T,
                           view(ρ_alb_mode, :, k),
                           B_tilde, Q_tilde, C_store,
                           A_work, B_work, C_work, Q_work,
                           workmat, tmp_mat, tmp_vec, p_col, rhs_buf)

        # Store results (pad beyond N_eff with B_ν/2 per mode)
        for j in 1:M
            for i in 1:N_eff
                P_out[i, j, k] = P_buf[i, j]
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
Per-mode Planck source is B_ν/2 (unpolarized emission split equally).

Writes solution into `P_out` (first `N_eff × M` entries). All other arrays
are caller-allocated scratch reused across frequencies."""
function _feautrier_single!(P_out, N_eff, M, dtau, μ, w, ν, T, ρ_alb,
                            B_tilde, Q_tilde, C_store,
                            A, B, C, Q,
                            workmat, tmp_mat, tmp_vec, p_col, rhs_buf)
    for i in 1:N_eff
        # Reset per-step working matrices/vectors
        fill!(A, 0.0)
        fill!(B, 0.0)
        fill!(C, 0.0)
        fill!(Q, 0.0)

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
            B_tilde[i] .= B
            Q_tilde[i] .= Q
            C_store[i] .= C  # C is zero here
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
        # Use Hummer (1962) dipole phase function for Thomson scattering, matching
        # the non-magnetic Feautrier (src/atmosphere/feautrier.jl:393-399).  This
        # ensures the B=0 limit recovers the non-magnetic result exactly: at B=0
        # both modes have ρ_alb = σ_T/(κ_ff+σ_T) and the per-mode kernel sums to
        # the unpolarized dipole kernel.  For B>0 this is still a local-per-mode
        # approximation (inter-mode redistribution = SPW09 Eq. 8 is the full
        # physics), but it is the *same* approximation as the non-magnetic case,
        # which is the consistency requirement.
        if i > 1 && i < N_eff
            for j in 1:M, jp in 1:M
                cross = 3.0/16.0 * (3.0 + 3.0*μ[j]^2*μ[jp]^2 - μ[j]^2 - μ[jp]^2)
                B[j, jp] -= 2.0 * ρ * cross * w[jp]
            end
            Q .= (1.0 - ρ) * Bν_half  # per-mode source = (1-ρ) B_ν/2
        end

        C_store[i] .= C

        # Forward elimination (write directly into preallocated B_tilde[i], Q_tilde[i])
        if i == 1
            B_tilde[i] .= B
            Q_tilde[i] .= Q
        else
            # Factor B_tilde[i-1] once and reuse for both C and Q solves.
            workmat .= B_tilde[i-1]
            F = lu!(workmat)
            tmp_mat .= C_store[i-1]
            ldiv!(F, tmp_mat)             # tmp_mat = B_tilde[i-1] \ C_store[i-1]
            tmp_vec .= Q_tilde[i-1]
            ldiv!(F, tmp_vec)             # tmp_vec = B_tilde[i-1] \ Q_tilde[i-1]
            # B_tilde[i] = B - A * tmp_mat
            mul!(B_tilde[i], A, tmp_mat, -1.0, 0.0)
            B_tilde[i] .+= B
            # Q_tilde[i] = Q - A * tmp_vec
            mul!(Q_tilde[i], A, tmp_vec, -1.0, 0.0)
            Q_tilde[i] .+= Q
        end
    end

    # Back-substitution (in-place into P_out[i, :])
    # Solve B_tilde[N_eff] * p_col = Q_tilde[N_eff]
    workmat .= B_tilde[N_eff]
    p_col .= Q_tilde[N_eff]
    F = lu!(workmat)
    ldiv!(F, p_col)
    for j in 1:M
        P_out[N_eff, j] = p_col[j]
    end
    for i in N_eff-1:-1:1
        # rhs = Q_tilde[i] - C_store[i] * P_out[i+1, :]
        for j in 1:M
            p_col[j] = P_out[i+1, j]
        end
        rhs_buf .= Q_tilde[i]
        mul!(rhs_buf, C_store[i], p_col, -1.0, 1.0)
        # solve B_tilde[i] * p_col = rhs_buf
        workmat .= B_tilde[i]
        p_col .= rhs_buf
        F = lu!(workmat)
        ldiv!(F, p_col)
        for j in 1:M
            P_out[i, j] = p_col[j]
        end
    end
    for i in 1:N_eff, j in 1:M
        P_out[i, j] = max(P_out[i, j], 0.0)
    end
    return P_out
end

"""
Rybicki temperature correction for two-mode atmosphere.
Each mode j contributes to the correction through its own T_k, U_k, K_k matrices.
The V_k weights sum over both modes.

ΔT enforces: Σ_j ∫ κ_j (J_j - B_ν/2) dν = 0 at each depth.
"""
function _rybicki_two_mode(N, nν, ν, y, T, κ, k_total, ρ_alb, f_ν, h_ν, J)
    # Frequency quadrature weights
    b = zeros(nν)
    for k in 1:nν
        if k == 1
            b[k] = 0.5 * (ν[2] - ν[1])
        elseif k == nν
            b[k] = 0.5 * (ν[nν] - ν[nν-1])
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
        for k in 1:nν
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
        for k in 1:nν
            T_diag = zeros(N)
            T_sub = zeros(N-1)
            T_sup = zeros(N-1)
            U_k = zeros(N)
            K_k = zeros(N)

            _build_rybicki_system_magnetic!(T_diag, T_sub, T_sup, U_k, K_k,
                                            N, k, j, ν, y, T, κ, k_total, ρ_alb,
                                            f_ν, h_ν, denom, b, B_bar)

            # LU solve and accumulate
            lu_diag, lu_sup, lu_rhs_K = tridiag_lu_forward!(T_sub, T_diag, T_sup, K_k)
            x_k = tridiag_lu_back!(lu_diag, lu_sup, lu_rhs_K)

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
                _, _, lu_rhs_j = tridiag_lu_forward!(T_sub, T_diag, T_sup, e_j)
                z_j = tridiag_lu_back!(lu_diag, lu_sup, lu_rhs_j)
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
    nν = length(ν_grid)
    M = length(μ)
    F = 0.0
    for k in 1:nν-1
        dν = ν_grid[k+1] - ν_grid[k]
        for j in 1:M, mode in 1:2
            F += 2π * μ[j] * 2.0 * P_all[1, j, k, mode] * w[j] * dν
        end
    end
    return F
end

end # module
