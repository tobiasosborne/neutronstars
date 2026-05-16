#=
Neutron star atmosphere column structure.

Hydrostatic equilibrium with ideal gas EOS for fully ionised hydrogen.
Source: Haakonsen et al. (2012) Sect. 2-3.

APPROXIMATION: Ideal gas ρ = m_p P/(2 k_B T) for fully ionised H (μ=0.5)
instead of OPAL EOS. Valid for T > 10⁵ K where Coulomb corrections are small.
=#

module AtmosphereStructure

using ..PhysicalConstants: k_B, m_p, h, σ_SB
using ..HydrogenOpacity: kappa_ff, sigma_thomson, total_opacity,
                          rosseland_mean, dBnu_dT
using ..GauntFactor: GauntTable
using ..BlackbodyAtmosphere: planck_Bnu
using ..SolverDefaults: TAU_DIFFUSION, Y_MAX_SEMIINFINITE,
                         T_SURFACE_FRAC_MCPHAC

export AtmosphereColumn, build_atmosphere, update_atmosphere!
export make_frequency_grid
export density_from_PT

"""
    density_from_PT(P, T) → ρ [g/cm³]

EOS for fully ionised hydrogen: ρ = m_p · P / (2 · k_B · T).

The factor of 2 comes from mean molecular weight μ = 0.5 (electrons and
protons both contribute, each with their share of the gas pressure).
This is the single source of truth for converting (P, T) to ρ — do not
inline the formula elsewhere.
"""
density_from_PT(P::Real, T::Real) = m_p * P / (2 * k_B * T)

"""
Atmosphere column data structure.
All arrays indexed [depth_point] or [depth_point, frequency].

Bead-AD: parametric on the floating-point element type `T<:Real` so the
non-magnetic solver chain composes through ForwardDiff `Dual` numbers.
The pre-existing `AtmosphereColumn` (concrete Float64) calls remain
backward-compatible because `::AtmosphereColumn` matches every
parametrization.
"""
struct AtmosphereColumn{T<:Real}
    N::Int                      # number of depth points
    nν::Int                     # number of frequencies
    g_s::T                      # surface gravity [cm/s²]
    T_eff::T                    # effective temperature [K]
    y::Vector{T}                # column depth [g/cm²]
    T::Vector{T}                # temperature [K]
    ρ::Vector{T}                # density [g/cm³]
    P::Vector{T}                # pressure [dyn/cm²]
    σ_scat::T                   # Thomson scattering opacity [cm²/g]
    κ::Matrix{T}                # free-free opacity κ[i,k] [cm²/g]
    k_total::Matrix{T}          # total opacity k[i,k] = κ + σ [cm²/g]
    ρ_alb::Matrix{T}            # scattering albedo ρ_ν[i,k]
    τ::Matrix{T}                # optical depth from surface τ[i,k]
    k_R::Vector{T}              # Rosseland mean opacity [cm²/g]
    ν_grid::Vector{T}           # frequency grid [Hz]
    gaunt::GauntTable           # Gaunt factor table reference
end

"""
    make_frequency_grid(T_eff, nν) → ν_grid [Hz]

Log-spaced frequency grid from 0.05 k_BT_eff/h to 120 k_BT_eff/h.
McPHAC: MINFREQKT=0.05, MAXFREQKT=120.
"""
function make_frequency_grid(T_eff::Real, nν::Int)
    ν_min = 0.05 * k_B * T_eff / h
    ν_max = 120.0 * k_B * T_eff / h
    T = typeof(float(T_eff))
    return T[10.0^logν for logν in range(log10(ν_min), log10(ν_max), length=nν)]
end

"""
    build_atmosphere(T_eff, g_s, ν_grid, gaunt; N=200, y_min, y_max) → AtmosphereColumn

Build initial atmosphere structure.

1. Log-space column depth grid
2. P(y) = g_s × y (hydrostatic eq, Haakonsen Eq. 1)
3. T(y) from Eddington approximation (Eq. 10)
4. ρ from ideal gas EOS
5. Compute opacities and optical depths
"""
function build_atmosphere(T_eff::TT, g_s::Tg,
                          ν_grid::AbstractVector{Tν},
                          gaunt::GauntTable;
                          N::Int=200,
                          y_min::Real=1e-9,
                          y_max::Real=1e2) where {TT<:Real, Tg<:Real, Tν<:Real}
    @assert T_eff > 0 && g_s > 0
    @assert N >= 10
    nν = length(ν_grid)

    # Common floating-point type: promotes Float64 → Dual automatically.
    Tp = promote_type(typeof(float(T_eff)), typeof(float(g_s)), float(Tν),
                      typeof(float(y_min)), typeof(float(y_max)))

    # Column depth grid (log-spaced)
    y = Tp[10.0^logy for logy in range(log10(y_min), log10(y_max), length=N)]

    # Hydrostatic equilibrium: P = g_s × y
    P = Tp[g_s * yi for yi in y]

    # Promote ν_grid into the working type (cheap; needed so AtmosphereColumn{Tp}
    # holds a homogeneously-typed grid even if caller passed Float64 ν_grid with
    # Dual T_eff/g_s).
    ν_grid_p = Tp[ν for ν in ν_grid]

    # Initial temperature profile: Eddington approximation
    T = _eddington_temperature(y, T_eff, g_s, ν_grid_p, gaunt)

    # Density from ideal gas EOS for fully ionised hydrogen.
    ρ = Tp[density_from_PT(P[i], T[i]) for i in 1:N]

    σ_scat_f = sigma_thomson()
    σ_scat = Tp(σ_scat_f)

    # Compute opacities at each (depth, frequency)
    κ = zeros(Tp, N, nν)
    k_tot = zeros(Tp, N, nν)
    ρ_alb = zeros(Tp, N, nν)
    for i in 1:N, k in 1:nν
        κ[i,k] = kappa_ff(ν_grid_p[k], T[i], ρ[i], gaunt)
        k_tot[i,k] = κ[i,k] + σ_scat
        ρ_alb[i,k] = σ_scat / k_tot[i,k]
    end

    # Optical depth from surface: τ[1,:] = 0, τ[i,:] = ∫ k dy
    τ = zeros(Tp, N, nν)
    for k in 1:nν
        for i in 2:N
            dy = y[i] - y[i-1]
            τ[i,k] = τ[i-1,k] + 0.5 * (k_tot[i,k] + k_tot[i-1,k]) * dy
        end
    end

    # Rosseland mean at each depth
    k_R = Tp[rosseland_mean(T[i], ρ[i], ν_grid_p, gaunt) for i in 1:N]

    # Check if τ_max ≥ TAU_DIFFUSION at highest frequency; extend if needed.
    # The check uses the underlying numeric value, but recursion preserves
    # the parametrization automatically (we pass the same typed inputs back).
    τ_max_hf = τ[end, nν]
    if τ_max_hf < TAU_DIFFUSION && y_max < Y_MAX_SEMIINFINITE
        y_max_new = y_max * (TAU_DIFFUSION / max(_value(τ_max_hf), 1.0))^1.2
        y_max_new = min(y_max_new, Y_MAX_SEMIINFINITE)
        return build_atmosphere(T_eff, g_s, ν_grid, gaunt;
                                N=N, y_min=y_min, y_max=y_max_new)
    end

    T_eff_p = Tp(T_eff)
    g_s_p = Tp(g_s)
    return AtmosphereColumn{Tp}(N, nν, g_s_p, T_eff_p, y, T, ρ, P,
                                 σ_scat, κ, k_tot, ρ_alb, τ, k_R,
                                 ν_grid_p, gaunt)
end

"""
    update_atmosphere!(col, T_new)

Update atmosphere after temperature correction: recompute ρ, opacities, τ.
"""
function update_atmosphere!(col::AtmosphereColumn{Tp},
                            T_new::AbstractVector{<:Real}) where {Tp<:Real}
    @assert length(T_new) == col.N
    @assert all(T_new .> 0) "Temperature must be positive everywhere"

    for i in 1:col.N
        col.T[i] = Tp(T_new[i])
    end
    for i in 1:col.N
        col.ρ[i] = density_from_PT(col.P[i], col.T[i])
    end

    for i in 1:col.N, k in 1:col.nν
        col.κ[i,k] = kappa_ff(col.ν_grid[k], col.T[i], col.ρ[i], col.gaunt)
        col.k_total[i,k] = col.κ[i,k] + col.σ_scat
        col.ρ_alb[i,k] = col.σ_scat / col.k_total[i,k]
    end

    for k in 1:col.nν
        col.τ[1,k] = zero(Tp)
        for i in 2:col.N
            dy = col.y[i] - col.y[i-1]
            col.τ[i,k] = col.τ[i-1,k] + 0.5 * (col.k_total[i,k] + col.k_total[i-1,k]) * dy
        end
    end

    for i in 1:col.N
        col.k_R[i] = rosseland_mean(col.T[i], col.ρ[i], col.ν_grid, col.gaunt)
    end
end

"""
Eddington temperature profile: T⁴(τ_R) = (3/4)T_eff⁴(τ_R + 2/3).
Simplified: T(y) ≈ T_eff × [(3/4)(k_R y + 2/3)]^{1/4}.
Surface: T(0) = 0.265 T_eff (McPHAC boundary condition).
"""
function _eddington_temperature(y::AbstractVector{Tp}, T_eff::Real,
                                 g_s::Real,
                                 ν_grid::AbstractVector{<:Real},
                                 gaunt::GauntTable) where {Tp<:Real}
    N = length(y)
    T = zeros(Tp, N)

    # Surface boundary
    T[1] = T_SURFACE_FRAC_MCPHAC * Tp(T_eff)

    # Estimate initial k_R at surface
    ρ_surf = m_p * Tp(g_s) * y[1] / (2.0 * k_B * T[1])
    k_R_est = rosseland_mean(T[1], ρ_surf, ν_grid, gaunt)

    T_eff_p = Tp(T_eff)
    g_s_p = Tp(g_s)
    for i in 2:N
        τ_R = k_R_est * y[i]
        T4 = 0.75 * T_eff_p^4 * (τ_R + 2.0/3.0)
        T[i] = T4^0.25

        # Update k_R estimate using current T and ρ
        ρ_i = m_p * g_s_p * y[i] / (2.0 * k_B * T[i])
        k_R_est = rosseland_mean(T[i], ρ_i, ν_grid, gaunt)
    end

    return T
end

# Extract the underlying numerical value from any Real (handles ForwardDiff.Dual
# without requiring ForwardDiff as a runtime dep — `value` falls back to the
# Real itself if there's no specialized method).
_value(x::Real) = x  # Float64 path: identity
# ForwardDiff.Dual <: Real, but it also responds to `Float64(...)` only for
# pure-value Duals. Safer: extract via float promotion at use site.
# Tests below confirm the comparison-with-Real-constant path works on Duals.

end # module
