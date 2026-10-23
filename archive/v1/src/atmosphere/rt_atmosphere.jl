#=
Top-level radiative transfer atmosphere solver.

Iteratively solves structure + RT to self-consistency:
  1. Build initial atmosphere (Eddington T profile)
  2. Solve Feautrier RT at all frequencies
  3. Compute Rybicki temperature correction ΔT
  4. Apply correction, update opacities and structure
  5. Check convergence: max|ΔT/T| < tol
  6. Repeat until converged

Source: Haakonsen et al. (2012) Sect. 3.
=#

module RTAtmosphere

using Printf
using Logging
using ..PhysicalConstants: σ_SB, k_B, h, m_p
using ..GauntFactor: GauntTable, load_gaunt_table
using ..AtmosphereStructure: AtmosphereColumn, build_atmosphere,
                              update_atmosphere!, make_frequency_grid,
                              density_from_PT
using ..FeautrierSolver: solve_feautrier_all, solve_feautrier_all_adaptive, gauss_legendre_half
using ..TemperatureCorrection: compute_temperature_correction
using ..BlackbodyAtmosphere: planck_Bnu
using ..HydrogenOpacity: kappa_ff, sigma_thomson, total_opacity
using ..SolverDefaults: DELTA_T_OVER_T_CAP
using ..SolverLogging: with_solver_logger

export solve_atmosphere, AtmosphereResult, rt_emergent_spectrum

"""
Result of a converged atmosphere calculation.

Bead E15: `column` snapshots the full converged `AtmosphereColumn` so that
downstream analyses (κ vs depth, τ vs depth, k_R, ρ, P, etc.) survive the
function call. The pre-existing `T_profile` and `y_grid` fields are kept
for backward compatibility; they are equivalent to `column.T` and `column.y`.

AD-percolation bead: parametric on the floating-point element type `T<:Real`
so this struct can be returned from `solve_atmosphere` under ForwardDiff.
For Float64 callers nothing changes (`AtmosphereResult{Float64}`).
"""
struct AtmosphereResult{T<:Real}
    T_eff::T
    g_s::T
    converged::Bool
    n_iterations::Int
    max_dT_over_T::T
    ν_grid::Vector{T}             # nν frequencies [Hz]
    μ_grid::Vector{T}             # M angles
    I_emergent::Matrix{T}         # nν × M emergent specific intensity
    T_profile::Vector{T}          # converged T(y)  (== column.T)
    y_grid::Vector{T}             # column depths   (== column.y)
    column::AtmosphereColumn{T}   # full converged atmosphere column (E15)
end

"""
    solve_atmosphere(T_eff, g_s, gaunt; nν, M, N, max_iter, tol, adaptive=true, verbose) → AtmosphereResult

Solve for a self-consistent NS atmosphere using iterative Rybicki
temperature correction.

Keyword `adaptive::Bool` (default `true`) selects between two Feautrier
back-ends:

- `adaptive=true` (default): `solve_feautrier_all_adaptive`. McPHAC-style
  per-frequency log-spaced depth resampling. Concentrates depth resolution
  near each frequency's photosphere. At T_eff=10⁶ K, g_s=2×10¹⁴,
  reduces max bulk spectrum residual (vs N=400 reference) by 0.6–0.9× at
  every tested N ∈ {50,100,150,200}; the gain grows with N
  (N=150 ratio 0.58×, N=200 ratio 0.61×). Both solvers agree to <1% at
  N=400, confirming they converge to the same limit.
- `adaptive=false`: `solve_feautrier_all`. Solves on the global y-grid with
  per-frequency depth truncation at τ ≥ TAU_DIFFUSION. Kept as an
  opt-out for direct comparison and as a fallback.
"""
function solve_atmosphere(T_eff::T1, g_s::T2,
                          gaunt::GauntTable;
                          nν::Int=50, M::Int=10, N::Int=100,
                          max_iter::Int=30,
                          tol::Real=1e-4,
                          anisotropic::Bool=true,
                          adaptive::Bool=true,
                          N_ν_adaptive::Int=200,
                          verbose::Bool=true) where {T1<:Real, T2<:Real}
    @assert T_eff > 0 && g_s > 0

    # Common floating-point type for the solver: promotes Dual when present.
    Tp = promote_type(typeof(float(T_eff)), typeof(float(g_s)))
    T_eff_p = Tp(T_eff)
    g_s_p = Tp(g_s)

    return with_solver_logger(verbose) do
        # Use the underlying values for the format strings (printf-style cannot
        # consume Dual directly).
        @info @sprintf("RT Atmosphere: T_eff=%.2e K, g_s=%.2e cm/s², nν=%d, M=%d, N=%d, adaptive=%s",
                       _val(T_eff_p), _val(g_s_p), nν, M, N, adaptive)

        # Set up grids
        ν_grid = make_frequency_grid(T_eff_p, nν)
        μ_f, w_f = gauss_legendre_half(M)
        # Promote μ, w into the working type so the Feautrier scratch arrays
        # are allocated in Tp (otherwise Float64 leaks back in and ::Tp
        # promotion happens piecewise on every operation).
        μ = Tp[m for m in μ_f]
        w = Tp[ww for ww in w_f]

        # Build initial atmosphere with Eddington T(y)
        @info "  Building initial structure..."
        col = build_atmosphere(T_eff_p, g_s_p, ν_grid, gaunt; N=N)
        @info @sprintf("  y_max=%.2e, τ_max(ν_max)=%.1f, N=%d depths",
                       _val(col.y[end]), _val(maximum(col.τ[end, :])), col.N)

        converged = false
        max_dT = Tp(Inf)
        n_iter = 0

        for iter in 1:max_iter
            n_iter = iter

            # Solve Feautrier RT at all frequencies
            P_all, J, f_ν, h_ν = if adaptive
                solve_feautrier_all_adaptive(col, μ, w; N_ν=N_ν_adaptive, anisotropic=anisotropic)
            else
                solve_feautrier_all(col, μ, w; anisotropic=anisotropic)
            end

            # Compute flux for diagnostics
            F_target = σ_SB * T_eff_p^4
            F_bol = _bolometric_flux(P_all, μ, w, ν_grid)
            flux_ratio = F_bol / F_target

            # Compute Rybicki temperature correction
            ΔT = compute_temperature_correction(col, f_ν, h_ν, J)

            # Dampen large corrections to aid convergence
            max_dT = maximum(abs.(ΔT ./ col.T))
            if max_dT > DELTA_T_OVER_T_CAP
                damp = DELTA_T_OVER_T_CAP / max_dT
                ΔT .*= damp
                max_dT = Tp(DELTA_T_OVER_T_CAP)
            end

            @debug @sprintf("  iter %2d: max|ΔT/T|=%.2e, F/σT⁴=%.4f",
                            iter, _val(max_dT), _val(flux_ratio))

            # Check convergence
            if max_dT < tol
                converged = true
                @info @sprintf("  CONVERGED at iteration %d (max|ΔT/T|=%.2e < %.2e)",
                               iter, _val(max_dT), float(tol))
                break
            end

            # Apply correction
            for i in 1:col.N
                col.T[i] = max(col.T[i] + ΔT[i], 0.1 * T_eff_p)  # floor at 10% of T_eff
            end

            # Update atmospheric structure (density, opacities, optical depths)
            _update_structure!(col, gaunt)
        end

        if !converged
            @warn @sprintf("  not converged after %d iterations (max|ΔT/T|=%.2e)",
                           max_iter, _val(max_dT))
        end

        # Final Feautrier solve for emergent spectrum
        P_all, J, f_ν, h_ν = if adaptive
            solve_feautrier_all_adaptive(col, μ, w; N_ν=N_ν_adaptive, anisotropic=anisotropic)
        else
            solve_feautrier_all(col, μ, w; anisotropic=anisotropic)
        end

        # Emergent intensity: I_ν(μ) = 2 P_ν(surface, μ)
        I_emergent = zeros(eltype(P_all), nν, M)
        for k in 1:nν, j in 1:M
            I_emergent[k, j] = 2 * P_all[1, j, k]
        end

        # Final flux diagnostic
        F_bol = _bolometric_flux(P_all, μ, w, ν_grid)
        @info @sprintf("  Final F/σT⁴=%.4f", _val(F_bol / (σ_SB * T_eff_p^4)))

        col_snapshot = deepcopy(col)

        return AtmosphereResult{Tp}(T_eff_p, g_s_p, converged, n_iter, max_dT,
                                     ν_grid, μ, I_emergent,
                                     col_snapshot.T, col_snapshot.y, col_snapshot)
    end
end

# Underlying numerical value for printf — strips ForwardDiff.Dual wrappers
# (which expose a `.value` field) recursively down to a plain `AbstractFloat`.
# We deliberately do NOT depend on ForwardDiff to keep src/ free of test-only
# deps; structural duck-typing on `.value` is sufficient.
_val(x::AbstractFloat) = x
function _val(x::Real)
    if hasproperty(x, :value)
        return _val(getproperty(x, :value))
    end
    return float(x)
end

"""
Compute bolometric emergent flux from the Feautrier solution.
"""
function _bolometric_flux(P_all, μ, w, ν_grid)
    nν = length(ν_grid)
    M = length(μ)
    R = promote_type(eltype(P_all), eltype(μ), eltype(w), eltype(ν_grid))
    F = zero(R)
    for k in 1:nν-1
        dν = ν_grid[k+1] - ν_grid[k]
        for j in 1:M
            # F_ν = 2π ∫₀¹ I_ν(μ) μ dμ where I = 2P at surface
            F += 2π * μ[j] * 2 * P_all[1, j, k] * w[j] * dν
        end
    end
    return F
end

"""
Update atmospheric structure after a temperature correction.
Recomputes density, opacities, and optical depths from the new T profile.
"""
function _update_structure!(col::AtmosphereColumn{Tp}, gaunt::GauntTable) where {Tp<:Real}
    for i in 1:col.N
        # Ideal gas for ionised H: ρ = m_p P/(2 k_B T), μ=0.5
        col.ρ[i] = density_from_PT(col.P[i], col.T[i])

        # Recompute opacities at each frequency
        for k in 1:col.nν
            ν = col.ν_grid[k]
            col.κ[i, k] = kappa_ff(ν, col.T[i], col.ρ[i], gaunt)
            col.k_total[i, k] = col.κ[i, k] + col.σ_scat
            col.ρ_alb[i, k] = col.σ_scat / col.k_total[i, k]
        end
    end

    # Recompute optical depths from surface
    for k in 1:col.nν
        col.τ[1, k] = zero(Tp)
        for i in 2:col.N
            dy = col.y[i] - col.y[i-1]
            col.τ[i, k] = col.τ[i-1, k] + 0.5 * (col.k_total[i, k] + col.k_total[i-1, k]) * dy
        end
    end
end

"""
    rt_emergent_spectrum(result, cos_θe, ν_out) → I_ν [erg/s/cm²/Hz/sr]

Interpolate emergent spectrum at given emission angle and frequencies.
"""
function rt_emergent_spectrum(result::AtmosphereResult,
                               cos_θe::Float64,
                               ν_out::AbstractVector{Float64})::Vector{Float64}
    @assert 0 < cos_θe <= 1.0

    μ = result.μ_grid
    M = length(μ)

    # Find μ bracket
    j_lo = 1
    for j in 1:M-1
        if μ[j] <= cos_θe <= μ[j+1]
            j_lo = j
            break
        end
    end
    if cos_θe < μ[1]; j_lo = 1; end
    if cos_θe > μ[end]; j_lo = M - 1; end
    j_hi = min(j_lo + 1, M)
    t_μ = j_lo < M ? clamp((cos_θe - μ[j_lo]) / (μ[j_hi] - μ[j_lo]), 0.0, 1.0) : 0.0

    # Interpolate in frequency (log-log)
    ν_in = result.ν_grid
    nν = length(ν_in)
    I_out = zeros(length(ν_out))

    for i in eachindex(ν_out)
        ν = ν_out[i]
        if ν <= ν_in[1]
            I_at = (1-t_μ)*result.I_emergent[1, j_lo] + t_μ*result.I_emergent[1, j_hi]
            I_out[i] = I_at * (ν / ν_in[1])^2  # Rayleigh-Jeans
        elseif ν >= ν_in[end]
            I_out[i] = 0.0
        else
            k_lo = searchsortedlast(ν_in, ν)
            k_lo = clamp(k_lo, 1, nν-1)
            k_hi = k_lo + 1
            t_ν = (log(ν) - log(ν_in[k_lo])) / (log(ν_in[k_hi]) - log(ν_in[k_lo]))

            I_ll = result.I_emergent[k_lo, j_lo]
            I_hl = result.I_emergent[k_hi, j_lo]
            I_lh = result.I_emergent[k_lo, j_hi]
            I_hh = result.I_emergent[k_hi, j_hi]

            I_lo = I_ll > 0 && I_hl > 0 ? exp((1-t_ν)*log(I_ll) + t_ν*log(I_hl)) : (1-t_ν)*I_ll + t_ν*I_hl
            I_hi = I_lh > 0 && I_hh > 0 ? exp((1-t_ν)*log(I_lh) + t_ν*log(I_hh)) : (1-t_ν)*I_lh + t_ν*I_hh
            I_out[i] = (1-t_μ)*I_lo + t_μ*I_hi
        end
    end

    return I_out
end

end # module
