#=
Tolman-Oppenheimer-Volkoff equation solver.

Integrates the TOV equations from the centre of a neutron star outward
to find M(R) for a given central density and EOS.

Source: Shapiro & Teukolsky (1983) Ch 5-6; Potekhin+ (2013) Sect. 6.
The TOV equations are:

    dP/dr = -G(ε + P)(m + 4πr³P/c²) / [r(rc² - 2Gm)]
    dm/dr = 4πr²ε/c²

where ε = ρc² is the total energy density, P is pressure, m is enclosed
gravitational mass, and r is the circumferential radius.

Integration: 4th-order Runge-Kutta with adaptive step, from centre outward.
Boundary: at r→0, expand to avoid 0/0 singularity.
Termination: P drops to surface pressure P_surface.
=#

module TOVSolver

using ..PhysicalConstants: G, c, M_sun
using ..BSkEOS: BSkParams, pressure_of_density, density_of_pressure,
                energy_density_of_density

export solve_tov, TOVResult

"""
Result of a TOV integration.
"""
struct TOVResult
    M::Float64         # Total gravitational mass [g]
    R::Float64         # Circumferential radius [cm]
    M_solar::Float64   # Mass in solar masses
    R_km::Float64      # Radius in km
    g_surface::Float64 # Surface gravity [cm s⁻²]
    compactness::Float64 # u = 2GM/(Rc²)
    ρ_c::Float64       # Central density [g cm⁻³]
    # Radial profiles
    r_profile::Vector{Float64}  # r [cm]
    m_profile::Vector{Float64}  # m(r) [g]
    P_profile::Vector{Float64}  # P(r) [dyn cm⁻²]
    ρ_profile::Vector{Float64}  # ρ(r) [g cm⁻³]
end

"""
    solve_tov(ρ_c, params::BSkParams; kwargs...) → TOVResult

Solve the TOV equations for central density ρ_c [g cm⁻³].

# Keyword arguments
- `P_surface`: surface pressure [dyn cm⁻²] (default: P(10⁶ g/cm³))
- `dr_init`: initial step size [cm] (default: 100 cm)
- `dr_max`: maximum step size [cm] (default: 1000 cm = 10 m)
- `rtol`: relative tolerance for adaptive stepping (default: 1e-6)
"""
function solve_tov(ρ_c::Float64, params::BSkParams;
                   P_surface::Float64=pressure_of_density(1e6, params),
                   dr_init::Float64=100.0,
                   dr_max::Float64=1000.0,
                   rtol::Float64=1e-6)::TOVResult

    @assert ρ_c > 1e10 "Central density too low: $ρ_c (need > 10¹⁰ g/cm³)"
    @assert ρ_c < 1e16 "Central density too high: $ρ_c"

    P_c = pressure_of_density(ρ_c, params)
    ε_c = energy_density_of_density(ρ_c)

    # === Central boundary expansion ===
    # At r=0: m(r) ≈ (4π/3)ε_c r³/c², P(r) ≈ P_c - 2πG(ε_c+P_c)(ε_c+3P_c)r²/(3c⁴)
    r0 = dr_init
    m0 = 4π/3 * ε_c / c^2 * r0^3
    dP_coeff = -2π * G * (ε_c + P_c) * (ε_c + 3P_c) / (3 * c^4)
    P0 = P_c + dP_coeff * r0^2

    if P0 <= P_surface || P0 <= 0
        error("Central expansion already below surface: P_c=$P_c, P0=$P0")
    end

    # Storage for profiles
    r_prof = Float64[0.0, r0]
    m_prof = Float64[0.0, m0]
    P_prof = Float64[P_c, P0]
    ρ_prof = Float64[ρ_c, density_of_pressure(P0, params)]

    r = r0
    m = m0
    P = P0

    # === RK4 integration ===
    dr = dr_init
    max_steps = 10_000_000
    for step in 1:max_steps
        # Adaptive step: reduce dr near surface where P changes fast
        dr = min(dr_max, dr_init * max(1.0, (P / P_surface)^0.3))
        dr = min(dr, 0.01 * r)  # Don't jump more than 1% of current radius

        # RK4 step
        k1_m, k1_P = tov_rhs(r, m, P, params)
        r_half = r + 0.5*dr
        m_half1 = m + 0.5*dr*k1_m
        P_half1 = P + 0.5*dr*k1_P

        if P_half1 <= 0; P_half1 = P_surface * 0.1; end

        k2_m, k2_P = tov_rhs(r_half, m_half1, P_half1, params)
        m_half2 = m + 0.5*dr*k2_m
        P_half2 = P + 0.5*dr*k2_P

        if P_half2 <= 0; P_half2 = P_surface * 0.1; end

        k3_m, k3_P = tov_rhs(r_half, m_half2, P_half2, params)
        r_end = r + dr
        m_end1 = m + dr*k3_m
        P_end1 = P + dr*k3_P

        if P_end1 <= 0; P_end1 = P_surface * 0.1; end

        k4_m, k4_P = tov_rhs(r_end, m_end1, P_end1, params)

        m_new = m + dr/6 * (k1_m + 2k2_m + 2k3_m + k4_m)
        P_new = P + dr/6 * (k1_P + 2k2_P + 2k3_P + k4_P)

        r += dr

        # Check for surface
        if P_new <= P_surface || P_new <= 0
            # Interpolate to find exact surface
            # Linear interpolation: P(r) crosses P_surface between (r-dr) and r
            frac = (P - P_surface) / (P - P_new)
            r_surf = (r - dr) + frac * dr
            m_surf = m + frac * (m_new - m)

            push!(r_prof, r_surf)
            push!(m_prof, m_surf)
            push!(P_prof, P_surface)
            push!(ρ_prof, density_of_pressure(P_surface, params))

            return _build_result(r_surf, m_surf, ρ_c, r_prof, m_prof, P_prof, ρ_prof)
        end

        m = m_new
        P = P_new

        # Store profile (subsample to keep memory bounded)
        if step % 100 == 0 || P < 10 * P_surface
            push!(r_prof, r)
            push!(m_prof, m)
            push!(P_prof, P)
            push!(ρ_prof, density_of_pressure(P, params))
        end
    end

    error("TOV integration did not reach surface after $max_steps steps")
end

"""
TOV right-hand side: returns (dm/dr, dP/dr).
"""
function tov_rhs(r::Float64, m::Float64, P::Float64,
                 params::BSkParams)
    ρ = density_of_pressure(P, params)
    ε = energy_density_of_density(ρ)

    # dm/dr = 4πr²ε/c²
    dm_dr = 4π * r^2 * ε / c^2

    # dP/dr = -G(ε + P)(m + 4πr³P/c²) / [r(rc² - 2Gm)]
    numerator = -G * (ε + P) * (m + 4π * r^3 * P / c^2)
    denominator = r * (r * c^2 - 2G * m)

    @assert denominator > 0 "Horizon reached at r=$(r/1e5) km! (2Gm/rc² ≥ 1)"

    dP_dr = numerator / denominator

    return dm_dr, dP_dr
end

function _build_result(R, M, ρ_c, r_prof, m_prof, P_prof, ρ_prof)
    R_km = R / 1e5
    M_solar = M / M_sun
    u = 2G * M / (R * c^2)  # compactness
    g_surf = G * M / (R^2 * sqrt(1 - u))  # surface gravity (redshifted)

    @assert isfinite(M_solar) "Mass diverged"
    @assert isfinite(R_km) "Radius diverged"
    @assert M_solar > 0 "Mass must be positive"
    @assert R_km > 0 "Radius must be positive"
    @assert u < 1.0 "Compactness u=$u ≥ 1 (black hole!)"

    return TOVResult(
        M, R, M_solar, R_km, g_surf, u, ρ_c,
        r_prof, m_prof, P_prof, ρ_prof
    )
end

"""
    mass_radius_curve(params::BSkParams; n_points=50) → (M_solar[], R_km[])

Compute the M(R) relation by varying central density.
"""
function mass_radius_curve(params::BSkParams;
                           log_ρc_min::Float64=14.0,
                           log_ρc_max::Float64=15.8,
                           n_points::Int=50)
    M_arr = Float64[]
    R_arr = Float64[]
    ρc_arr = Float64[]

    for i in 1:n_points
        log_ρc = log_ρc_min + (i-1) * (log_ρc_max - log_ρc_min) / (n_points - 1)
        ρ_c = 10.0^log_ρc

        try
            result = solve_tov(ρ_c, params)
            push!(M_arr, result.M_solar)
            push!(R_arr, result.R_km)
            push!(ρc_arr, ρ_c)
        catch e
            @warn "TOV failed at ρ_c=$ρ_c: $e"
        end
    end

    return M_arr, R_arr, ρc_arr
end

end # module
