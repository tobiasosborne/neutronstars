#=
Analytical representations of unified BSk equations of state.
Source: Potekhin et al., A&A 560, A48 (2013), Eq. (3), Table 2.

P(ρ) fit valid for 10⁶ ≲ ρ ≲ 10¹⁶ g cm⁻³.
Typical fit error of P: ≈1% for ξ ≳ 6, max 3.7% at ξ = 9.51.
=#

module BSkEOS

using ..PhysicalConstants: c, m_u

export pressure_of_density, energy_density_of_density, pseudo_enthalpy
export BSk19_params, BSk20_params, BSk21_params

"""
    BSk EOS parameters a₁…a₂₃ for Eq. (3) of Potekhin+ (2013).
"""
struct BSkParams
    a::NTuple{23, Float64}
    name::String
end

# Table 2: BSk21 parameters (preferred — stiffest, M_max ≈ 2.27 M_☉)
const BSk21_params = BSkParams((
    4.857,   # a1
    6.981,   # a2
    0.00706, # a3
    0.19351, # a4
    4.085,   # a5
    12.065,  # a6
    10.521,  # a7
    1.5905,  # a8
    4.104,   # a9
    -28.726, # a10
    2.0845,  # a11
    4.89,    # a12
    14.302,  # a13
    22.881,  # a14
    -1.7690, # a15
    0.989,   # a16
    15.313,  # a17
    0.091,   # a18
    4.68,    # a19
    11.65,   # a20
    -0.086,  # a21
    10.0,    # a22
    14.15,   # a23
), "BSk21")

# Table 2: BSk20 parameters
const BSk20_params = BSkParams((
    4.078,   # a1
    7.587,   # a2
    0.00839, # a3
    0.21695, # a4
    3.614,   # a5
    11.942,  # a6
    13.751,  # a7
    1.3373,  # a8
    3.606,   # a9
    -22.996, # a10
    1.6229,  # a11
    4.88,    # a12
    14.274,  # a13
    23.560,  # a14
    -1.5564, # a15
    2.095,   # a16
    15.294,  # a17
    0.084,   # a18
    6.36,    # a19
    11.67,   # a20
    -0.042,  # a21
    14.8,    # a22
    14.18,   # a23
), "BSk20")

# Table 2: BSk19 parameters
const BSk19_params = BSkParams((
    3.916,   # a1
    7.701,   # a2
    0.00858, # a3
    0.22114, # a4
    3.269,   # a5
    11.964,  # a6
    13.349,  # a7
    1.3683,  # a8
    3.254,   # a9
    -12.953, # a10
    0.9237,  # a11
    6.20,    # a12
    14.383,  # a13
    16.693,  # a14
    -1.0514, # a15
    2.486,   # a16
    15.362,  # a17
    0.085,   # a18
    6.23,    # a19
    11.68,   # a20
    -0.029,  # a21
    20.1,    # a22
    14.19,   # a23
), "BSk19")

"""
    pressure_of_density(ρ, params::BSkParams) → P

Pressure P [dyn cm⁻²] as a function of gravitational mass density ρ [g cm⁻³].
Potekhin et al. (2013) Eq. (3).

Variables: ξ = log₁₀(ρ/g cm⁻³), ζ = log₁₀(P/dyn cm⁻²).
"""
function pressure_of_density(ρ::Float64, params::BSkParams)::Float64
    @assert ρ > 0 "Density must be positive, got $ρ"
    @assert ρ < 1e17 "Density exceeds NS range, got $ρ"

    a = params.a
    ξ = log10(ρ)

    # Potekhin+ (2013) Eq. (3)
    ζ = (a[1] + a[2]*ξ + a[3]*ξ^3) / (1 + a[4]*ξ) *
        (exp(a[5]*(ξ - a[6])) + 1)^(-1) +
        (a[7] + a[8]*ξ) * (exp(a[9]*(a[6] - ξ)) + 1)^(-1) +
        (a[10] + a[11]*ξ) * (exp(a[12]*(a[13] - ξ)) + 1)^(-1) +
        (a[14] + a[15]*ξ) * (exp(a[16]*(a[17] - ξ)) + 1)^(-1) +
        a[18] / (1 + (a[19]*(ξ - a[20]))^2) +
        a[21] / (1 + (a[22]*(ξ - a[23]))^2)

    P = 10.0^ζ

    @assert isfinite(P) "P($ρ) returned non-finite: ζ=$ζ"
    @assert P > 0 "P($ρ) must be positive, got $P"
    return P
end

"""
    energy_density_of_density(ρ) → ε

Total energy density ε = ρc² [erg cm⁻³] (rest mass + internal energy).
For the TOV equation in the form using ε and P.

Note: The BSk fit gives P(ρ) where ρ is the gravitational mass density,
which already includes internal energy contributions. For the TOV equation,
ε = ρc² is the correct energy density to use with this parametrisation.
See Potekhin+ (2013) Sect. 3.1, Eq. (1).
"""
function energy_density_of_density(ρ::Float64)::Float64
    @assert ρ > 0 "Density must be positive"
    return ρ * c^2
end

"""
    pseudo_enthalpy(P, ρ, params::BSkParams) → H

Pseudo-enthalpy H = ∫₀ᴾ dP'/(ρ(P')c² + P').
Potekhin+ (2013) Eq. (2). Computed numerically via trapezoidal rule.
"""
function pseudo_enthalpy(P::Float64, params::BSkParams;
                         n_points::Int=1000)::Float64
    @assert P > 0 "Pressure must be positive"

    log_P_min = log10(pressure_of_density(1e6, params))
    log_P_max = log10(P)

    if log_P_max <= log_P_min
        return 0.0
    end

    # Integrate using trapezoidal rule in log P space
    dlog_P = (log_P_max - log_P_min) / n_points
    H = 0.0
    for i in 0:n_points
        log_Pi = log_P_min + i * dlog_P
        Pi = 10.0^log_Pi
        # Need ρ(P) — invert P(ρ) via bisection
        ρi = density_of_pressure(Pi, params)
        εi = energy_density_of_density(ρi)
        integrand = 1.0 / (εi + Pi)
        weight = (i == 0 || i == n_points) ? 0.5 : 1.0
        H += weight * integrand * Pi * log(10.0) * dlog_P
    end

    @assert isfinite(H) "Pseudo-enthalpy diverged at P=$P"
    return H
end

"""
    density_of_pressure(P, params::BSkParams) → ρ

Invert P(ρ) via bisection. Returns ρ [g cm⁻³].
"""
function density_of_pressure(P::Float64, params::BSkParams;
                             rtol::Float64=1e-10)::Float64
    @assert P > 0 "Pressure must be positive"

    # Bracket: ρ ∈ [1e4, 1e16]
    ρ_lo = 1e4
    ρ_hi = 1e16

    P_lo = pressure_of_density(ρ_lo, params)
    P_hi = pressure_of_density(ρ_hi, params)

    # Extend bracket if needed
    while P < P_lo && ρ_lo > 1e0
        ρ_lo /= 10.0
        P_lo = pressure_of_density(ρ_lo, params)
    end
    while P > P_hi && ρ_hi < 1e17
        ρ_hi *= 10.0
        P_hi = pressure_of_density(ρ_hi, params)
    end

    @assert P >= P_lo && P <= P_hi "P=$P outside bracket [$P_lo, $P_hi]"

    # Bisection in log space
    log_lo = log10(ρ_lo)
    log_hi = log10(ρ_hi)

    for _ in 1:200
        log_mid = 0.5 * (log_lo + log_hi)
        ρ_mid = 10.0^log_mid
        P_mid = pressure_of_density(ρ_mid, params)

        if P_mid < P
            log_lo = log_mid
        else
            log_hi = log_mid
        end

        if (log_hi - log_lo) / abs(log_mid) < rtol
            return 10.0^(0.5 * (log_lo + log_hi))
        end
    end

    error("Bisection did not converge for P=$P after 200 iterations")
end

end # module
