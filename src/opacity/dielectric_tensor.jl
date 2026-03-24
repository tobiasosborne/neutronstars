#=
Full cold-plasma dielectric tensor and normal-mode polarization vectors.

Computes the mode parameter q from the Stix parameters (S, D, P)
and derives the cyclic polarization weights |e_{j,α}|² from K_j.

The q-based approach is numerically stable at all densities because q
doesn't depend on density in the dominant term (it's set by ω_ce/ω).
The n²-based approach suffers catastrophic cancellation at low density
where S ≈ n² ≈ 1.

Source: Ginzburg (1970) §10, Potekhin & Chabrier (2003) ApJ 585, 955,
        Eq. 25-27. Shafranov (1967) for polarization vectors.
=#

module DielectricTensor

using ..PhysicalConstants: e_charge, m_e, m_p, c, ħ

export stix_parameters, polarization_weights_full

"""
    stix_parameters(ω, B, n_e) → (S, D, P)

Cold-plasma Stix parameters including both electron and proton contributions.
"""
function stix_parameters(ω::Float64, B::Float64, n_e::Float64)
    @assert ω > 0 && B > 0 && n_e > 0

    ω_ce = e_charge * B / (m_e * c)
    ω_cp = e_charge * B / (m_p * c)
    ω_pe² = 4π * n_e * e_charge^2 / m_e
    ω_pi² = 4π * n_e * e_charge^2 / m_p

    ω² = ω * ω
    de = ω² - ω_ce^2
    dp = ω² - ω_cp^2

    S = 1.0 - ω_pe² / de - ω_pi² / dp
    D = ω_ce * ω_pe² / (ω * de) - ω_cp * ω_pi² / (ω * dp)
    P = 1.0 - ω_pe² / ω² - ω_pi² / ω²

    return S, D, P
end

"""
    polarization_weights_full(ω, B, θ_B, n_e) → (w1, w2)

Compute |e_{j,α}|² for the two normal modes using the q-based approach.

Mode parameter q (P&C 2003 Eq. 25, for cold plasma without vacuum):
  q = (P - S) sin²θ / (2D cosθ)

Transverse polarization ratios (Shafranov 1967):
  K₁ = -q + √(1 + q²)   (extraordinary)
  K₂ = -q - √(1 + q²)   (ordinary)

Longitudinal component (from wave equation with n² = S - K_j D):
  K_{z,j} = -K_j (S - K_j D) sinθ cosθ / (P - (S - K_j D) sin²θ)

Returns w1[3], w2[3] where w[1]=|e_{j,-1}|², w[2]=|e_{j,0}|², w[3]=|e_{j,+1}|².
"""
function polarization_weights_full(ω::Float64, B::Float64,
                                    θ_B::Float64, n_e::Float64)
    S, D, P = stix_parameters(ω, B, n_e)

    cosθ = cos(θ_B)
    sinθ = sin(θ_B)

    # θ ≈ 0: quasi-longitudinal, modes are circularly polarised
    if abs(sinθ) < 1e-8
        w1 = [0.0, 0.0, 1.0]   # extraordinary: right-circular (α = +1)
        w2 = [1.0, 0.0, 0.0]   # ordinary: left-circular (α = -1)
        return w1, w2
    end

    # θ ≈ π/2: quasi-transverse
    if abs(cosθ) < 1e-8
        w1 = [0.0, 1.0, 0.0]   # extraordinary: along B (α = 0)
        w2 = [0.5, 0.0, 0.5]   # ordinary: circular mix (α = ±1)
        return w1, w2
    end

    # Compute q from Eq. 25 (full formula, no approximation)
    # q = (P - S) sin²θ / (2D cosθ)
    if abs(D) < 1e-30
        # No magnetic splitting → linear polarisation
        w1 = [0.0, 1.0, 0.0]
        w2 = [0.5, 0.0, 0.5]
        return w1, w2
    end

    q = (P - S) * sinθ^2 / (2.0 * D * cosθ)

    # K values from quadratic K² + 2qK - 1 = 0
    sq = sqrt(1.0 + q^2)
    K1 = -q + sq     # extraordinary
    K2 = -q - sq     # ordinary (= -1/K1)

    # Longitudinal components K_{z,j}
    # n² = S - K_j D, then K_z = -K(S-KD)sinθcosθ / (P - (S-KD)sin²θ)
    w1 = compute_weights_from_K(K1, S, D, P, sinθ, cosθ)
    w2 = compute_weights_from_K(K2, S, D, P, sinθ, cosθ)

    return w1, w2
end

"""
Compute |e_{j,α}|² from the transverse polarization ratio K_j.
"""
function compute_weights_from_K(K::Float64, S::Float64, D::Float64,
                                 P::Float64, sinθ::Float64, cosθ::Float64)
    # Refractive index for this mode: n² = S - K*D
    n² = S - K * D

    # Longitudinal component
    denom_z = P - n² * sinθ^2
    if abs(denom_z) < 1e-30 * (1.0 + abs(P) + abs(n²))
        # E_z dominates: mode is linearly polarised along B
        return [0.0, 1.0, 0.0]
    end

    Kz = -n² * sinθ * cosθ * K / denom_z

    # Normalisation factor
    norm = 1.0 + K^2 + Kz^2

    # Cyclic components
    wp1 = (K + 1.0)^2 / (2.0 * norm)   # |e_{+1}|²
    wm1 = (K - 1.0)^2 / (2.0 * norm)   # |e_{-1}|²
    w0  = Kz^2 / norm                    # |e_0|²

    # Safety clamp and renormalize
    wp1 = max(0.0, wp1)
    wm1 = max(0.0, wm1)
    w0  = max(0.0, w0)
    s = wm1 + w0 + wp1
    if s > 0
        wm1 /= s; w0 /= s; wp1 /= s
    end

    return [wm1, w0, wp1]
end

end # module
