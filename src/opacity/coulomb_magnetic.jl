#=
Magnetic Coulomb logarithm for free-free absorption in quantizing B fields.

Implements the Born-approximation Coulomb logarithm Λ_α^ff that accounts
for Landau quantization of electron orbits. This replaces the classical
(non-magnetic) Λ_cl when β_e = ℏω_ce/(k_BT) >> 1.

Source: Potekhin & Chabrier (2003) ApJ 585, 955, Eqs. (44a-e).
Local: refs/potekhin_chabrier_2003_ff_opacity.pdf
=#

module MagneticCoulomb

using SpecialFunctions: besselk
using QuadGK: quadgk
using ..PhysicalConstants: ħ, k_B, e_charge, m_e, c

export coulomb_log_magnetic, coulomb_log_classical_safe

"""
    coulomb_log_magnetic(α, u, β_e) → Λ_α^ff

Magnetic Coulomb logarithm. P&C 2003 Eq. (44a).
Λ_α^ff = (3/4) e^{u/2} Σ_n ∫₀^∞ Q_n^α(β_e, u, y) dy

Arguments:
- α ∈ {-1, 0, +1}: polarisation index
- u = ℏω/(k_BT): dimensionless photon energy
- β_e = ℏω_ce/(k_BT): quantization parameter
"""
function coulomb_log_magnetic(α::Int, u::Float64, β_e::Float64)::Float64
    @assert α ∈ (-1, 0, 1)
    @assert u > 0 && β_e > 0

    # If field is non-quantizing, fall back to classical
    if β_e < 0.1
        return coulomb_log_classical_safe(u)
    end

    prefactor = 0.75 * exp(min(u/2, 300.0))

    # Sum over Landau levels n. For β_e >> 1, n=0 dominates.
    # Include n from -N_max to +N_max.
    N_max = max(2, ceil(Int, 3.0 / β_e + 2))  # few terms needed when β_e large
    N_max = min(N_max, 50)  # cap for safety

    total = 0.0
    for n in -N_max:N_max
        # Integrate Q_n^α over y ∈ [0, ∞)
        integrand(y) = _Q_n_alpha(n, α, β_e, u, y)

        # Adaptive quadrature; integrand decays for large y
        val, _ = quadgk(integrand, 0.0, 20.0; rtol=1e-4, maxevals=500)
        total += val
    end

    Λ = prefactor * total
    return max(Λ, 0.01)  # floor to avoid log(0)
end

"""
Q_n^α(β_e, u, y) — integrand for the magnetic Coulomb logarithm.
P&C 2003 Eq. (44b-e).
"""
function _Q_n_alpha(n::Int, α::Int, β_e::Float64, u::Float64, y::Float64)::Float64
    # Eq. (44d): auxiliary variables
    exp_neg_β = exp(-min(β_e, 500.0))
    θ = (1.0 + exp_neg_β) / max(1.0 - exp_neg_β, 1e-30)

    ζ = sqrt(1.0 + 2θ*y + y^2)

    # Eq. (44e): x_n
    x_n = abs(u - n * β_e) / sqrt(0.25 + y / β_e)

    # Suppress contributions where Bessel function argument is very large
    if x_n > 300.0
        return 0.0
    end

    # Boltzmann suppression for |n| > 0 Landau levels
    # Factor: [(y+θ+ζ) sinh(β_e/2)]^{-|n|}
    abs_n = abs(n)
    if abs_n > 0
        sinh_half = 0.5 * (exp(min(β_e/2, 300.0)) - exp(-min(β_e/2, 300.0)))
        base = (y + θ + ζ) * sinh_half
        if base <= 0 || log(base) * abs_n > 300.0
            return 0.0
        end
        landau_factor = base^(-abs_n)
    else
        landau_factor = 1.0
    end

    # Eq. (44c): A_n^α
    if α == 0
        # A_n^0 = x_n K₁(x_n) / (y + β_e/4)
        K1_val = x_n > 0 ? besselk(1, x_n) : 1.0  # K₁(x) ~ 1/x as x→0
        A = x_n * K1_val / (y + β_e / 4.0)
    else
        # A_n^{±1} = (y + θ + |n|ζ) / ζ² × K₀(x_n)
        K0_val = x_n > 1e-10 ? besselk(0, x_n) : 10.0  # K₀(x) ~ -ln(x) as x→0
        A = (y + θ + abs_n * ζ) / ζ^2 * K0_val
    end

    # Eq. (44b): Q_n^α = y / [ζ × ...] × A_n^α
    Q = y / ζ * landau_factor * A

    return isfinite(Q) ? max(Q, 0.0) : 0.0
end

"""
    coulomb_log_classical_safe(u) → Λ_cl

Classical Coulomb logarithm with safe handling of edge cases.
Λ_cl = e^{u/2} K₀(u/2). P&C 2003 Eq. (43).
"""
function coulomb_log_classical_safe(u::Float64)::Float64
    u2 = max(u / 2.0, 1e-30)
    if u2 > 300.0
        # Asymptotic: K₀(x) ~ √(π/(2x)) e^{-x}, so e^x K₀(x) ~ √(π/(2x))
        return sqrt(π / (2u2))
    end
    return exp(u2) * besselk(0, u2)
end

end # module
