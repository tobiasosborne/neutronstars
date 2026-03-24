#=
Neutron star atmosphere column structure.

Hydrostatic equilibrium with ideal gas EOS for fully ionised hydrogen.
Source: Haakonsen et al. (2012) Sect. 2-3.

APPROXIMATION: Ideal gas ρ = m_p P/(k_B T) instead of OPAL EOS.
Valid for T > 10⁵ K where Coulomb corrections are small.
=#

module AtmosphereStructure

using ..PhysicalConstants: k_B, m_p, h, σ_SB
using ..HydrogenOpacity: kappa_ff, sigma_thomson, total_opacity,
                          rosseland_mean, dBnu_dT
using ..GauntFactor: GauntTable
using ..BlackbodyAtmosphere: planck_Bnu

export AtmosphereColumn, build_atmosphere, update_atmosphere!
export make_frequency_grid

"""
Atmosphere column data structure.
All arrays indexed [depth_point] or [depth_point, frequency].
"""
mutable struct AtmosphereColumn
    N::Int                      # number of depth points
    K::Int                      # number of frequencies
    g_s::Float64                # surface gravity [cm/s²]
    T_eff::Float64              # effective temperature [K]
    y::Vector{Float64}          # column depth [g/cm²]
    T::Vector{Float64}          # temperature [K]
    ρ::Vector{Float64}          # density [g/cm³]
    P::Vector{Float64}          # pressure [dyn/cm²]
    σ_scat::Float64             # Thomson scattering opacity [cm²/g]
    κ::Matrix{Float64}          # free-free opacity κ[i,k] [cm²/g]
    k_total::Matrix{Float64}    # total opacity k[i,k] = κ + σ [cm²/g]
    ρ_alb::Matrix{Float64}      # scattering albedo ρ_ν[i,k]
    τ::Matrix{Float64}          # optical depth from surface τ[i,k]
    k_R::Vector{Float64}        # Rosseland mean opacity [cm²/g]
    ν_grid::Vector{Float64}     # frequency grid [Hz]
    gaunt::GauntTable           # Gaunt factor table reference
end

"""
    make_frequency_grid(T_eff, K) → ν_grid [Hz]

Log-spaced frequency grid from 0.05 k_BT_eff/h to 120 k_BT_eff/h.
McPHAC: MINFREQKT=0.05, MAXFREQKT=120.
"""
function make_frequency_grid(T_eff::Float64, K::Int)::Vector{Float64}
    ν_min = 0.05 * k_B * T_eff / h
    ν_max = 120.0 * k_B * T_eff / h
    return [10.0^logν for logν in range(log10(ν_min), log10(ν_max), length=K)]
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
function build_atmosphere(T_eff::Float64, g_s::Float64,
                          ν_grid::Vector{Float64},
                          gaunt::GauntTable;
                          N::Int=200,
                          y_min::Float64=1e-9,
                          y_max::Float64=1e2)::AtmosphereColumn
    @assert T_eff > 0 && g_s > 0
    @assert N >= 10
    K = length(ν_grid)

    # Column depth grid (log-spaced)
    y = [10.0^logy for logy in range(log10(y_min), log10(y_max), length=N)]

    # Hydrostatic equilibrium: P = g_s × y
    P = g_s .* y

    # Initial temperature profile: Eddington approximation
    # T(y) integrated from surface using dT/dy = (3/16) k_R (T/T_eff)³ T_eff
    # Start with boundary condition T(0) = 0.265 T_eff (McPHAC)
    T = _eddington_temperature(y, T_eff, g_s, ν_grid, gaunt)

    # Density from ideal gas EOS: P = n k_B T = (ρ/m_p) k_B T
    ρ = [m_p * P[i] / (k_B * T[i]) for i in 1:N]

    σ_scat = sigma_thomson()

    # Compute opacities at each (depth, frequency)
    κ = zeros(N, K)
    k_tot = zeros(N, K)
    ρ_alb = zeros(N, K)
    for i in 1:N, k in 1:K
        κ[i,k] = kappa_ff(ν_grid[k], T[i], ρ[i], gaunt)
        k_tot[i,k] = κ[i,k] + σ_scat
        ρ_alb[i,k] = σ_scat / k_tot[i,k]
    end

    # Optical depth from surface: τ[1,:] = 0, τ[i,:] = ∫ k dy
    τ = zeros(N, K)
    for k in 1:K
        for i in 2:N
            dy = y[i] - y[i-1]
            τ[i,k] = τ[i-1,k] + 0.5 * (k_tot[i,k] + k_tot[i-1,k]) * dy
        end
    end

    # Rosseland mean at each depth
    k_R = [rosseland_mean(T[i], ρ[i], ν_grid, gaunt) for i in 1:N]

    # Check if τ_max ≥ 80 at highest frequency; extend if needed
    τ_max_hf = τ[end, K]
    if τ_max_hf < 80.0 && y_max < 1e5
        y_max_new = y_max * (80.0 / max(τ_max_hf, 1.0))^1.2
        y_max_new = min(y_max_new, 1e5)
        return build_atmosphere(T_eff, g_s, ν_grid, gaunt;
                                N=N, y_min=y_min, y_max=y_max_new)
    end

    return AtmosphereColumn(N, K, g_s, T_eff, y, T, ρ, P,
                             σ_scat, κ, k_tot, ρ_alb, τ, k_R,
                             ν_grid, gaunt)
end

"""
    update_atmosphere!(col, T_new)

Update atmosphere after temperature correction: recompute ρ, opacities, τ.
"""
function update_atmosphere!(col::AtmosphereColumn, T_new::Vector{Float64})
    @assert length(T_new) == col.N
    @assert all(T_new .> 0) "Temperature must be positive everywhere"

    col.T .= T_new
    for i in 1:col.N
        col.ρ[i] = m_p * col.P[i] / (k_B * col.T[i])
    end

    for i in 1:col.N, k in 1:col.K
        col.κ[i,k] = kappa_ff(col.ν_grid[k], col.T[i], col.ρ[i], col.gaunt)
        col.k_total[i,k] = col.κ[i,k] + col.σ_scat
        col.ρ_alb[i,k] = col.σ_scat / col.k_total[i,k]
    end

    for k in 1:col.K
        col.τ[1,k] = 0.0
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
function _eddington_temperature(y::Vector{Float64}, T_eff::Float64,
                                 g_s::Float64, ν_grid::Vector{Float64},
                                 gaunt::GauntTable)::Vector{Float64}
    N = length(y)
    T = zeros(N)

    # Surface boundary
    T[1] = 0.265 * T_eff

    # Estimate initial k_R at surface
    ρ_surf = m_p * g_s * y[1] / (k_B * T[1])
    k_R_est = rosseland_mean(T[1], ρ_surf, ν_grid, gaunt)

    # Integrate inward using Eddington: T⁴ = (3/4)T_eff⁴(τ_R + 2/3)
    # τ_R ≈ k_R × y (crude but sufficient for initial guess)
    for i in 2:N
        τ_R = k_R_est * y[i]
        T4 = 0.75 * T_eff^4 * (τ_R + 2.0/3.0)
        T[i] = T4^0.25

        # Update k_R estimate using current T and ρ
        ρ_i = m_p * g_s * y[i] / (k_B * T[i])
        k_R_est = rosseland_mean(T[i], ρ_i, ν_grid, gaunt)
    end

    return T
end

end # module
