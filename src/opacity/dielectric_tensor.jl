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

using ..PhysicalConstants: e_charge, m_e, m_p, c, ħ, α_fine, B_Q, keV

export stix_parameters, vacuum_coefficients
export cold_polarization_parameter, vacuum_polarization_parameter
export vacuum_mode_parameters
export polarization_weights_cold, polarization_weights_vacuum, polarization_weights_full
export vacuum_resonance_density, vacuum_resonance_pjump

"""
    stix_parameters(ω, B, n_e) → (S, D, P)

Cold-plasma Stix parameters including both electron and proton contributions.

Sign convention: this uses the Ginzburg (1970) / Lai-Ho convention where
`D ≡ g` is the off-diagonal element of the Hermitian dielectric tensor in
van Adelsberg & Lai (2006) Eq. (10), `[ε_pl] = [[ε, ig, 0], [-ig, ε, 0],
[0, 0, η]]`. Identification with our Stix notation: S = ε, D = g, P = η.

For a single electron species (q = -e), g = +ω_ce ω_pe² / (ω(ω² - ω_ce²)),
so D_here = +ω_ce ω_pe² / (ω · de). This is *opposite in sign* to the
textbook Stix (1992) definition D_Stix = (R-L)/2, which carries an extra
minus from the electron charge. The two conventions differ only by a global
sign for D, which cancels in every observable (q ∝ 1/D appears squared via
K₁K₂ = -1, and |e_α|² is built from K = -q ± √(1+q²)). See van Adelsberg &
Lai (2006) Eqs. (10, 19) and P&C 2003 Eq. (25) for the underlying tensor
definition.
"""
function stix_parameters(ω::Float64, B::Float64, n_e::Float64)::NTuple{3, Float64}
    ω_ce = e_charge * B / (m_e * c)
    ω_cp = e_charge * B / (m_p * c)
    ω_pe² = 4π * n_e * e_charge^2 / m_e
    ω_pi² = 4π * n_e * e_charge^2 / m_p

    ω² = ω * ω
    de = ω² - ω_ce^2
    dp = ω² - ω_cp^2

    S = 1.0 - ω_pe² / de - ω_pi² / dp
    # D = g = ε_xy / i in the Ginzburg/Lai tensor (Eq. 10 of vAL2006).
    # Electron contribution: +ω_ce ω_pe²/(ω·de); ion contribution carries
    # opposite sign because protons have opposite charge from electrons.
    D = ω_ce * ω_pe² / (ω * de) - ω_cp * ω_pi² / (ω * dp)
    P = 1.0 - ω_pe² / ω² - ω_pi² / ω²

    return S, D, P
end

"""
    polarization_weights_cold(ω, B, θ_B, n_e) → (w1, w2)

Compute |e_{j,α}|² for the two normal modes of a magnetised cold plasma
(no vacuum-polarisation correction), using the q-based approach of
Ginzburg (1970) §10 / Potekhin & Chabrier (2003) Eqs. 25-27.

Each of `w1`, `w2` is an `NTuple{3, Float64}` (allocation-free).

Mode parameter q (P&C 2003 Eq. 25, cold plasma, no vacuum terms):
  q = (P - S) sin²θ / (2D cosθ)

Transverse polarization ratios (Shafranov 1967; quadratic K² + 2qK - 1 = 0):
  K₁ = -q + √(1 + q²)   (extraordinary)
  K₂ = -q - √(1 + q²)   (ordinary)

Longitudinal component (from the wave equation with n² = S - K_j D):
  K_{z,j} = -K_j (S - K_j D) sinθ cosθ / (P - (S - K_j D) sin²θ)

Returns w1[3], w2[3] where w[1]=|e_{j,-1}|², w[2]=|e_{j,0}|², w[3]=|e_{j,+1}|².

For magnetised neutron-star atmospheres where B ≳ 10⁶ G, prefer
`polarization_weights_vacuum`, which adds the QED vacuum-polarisation
correction (Potekhin, Lai & Chabrier 2004).
"""
function polarization_weights_cold(ω::Float64, B::Float64,
                                   θ_B::Float64, n_e::Float64)::NTuple{2, NTuple{3, Float64}}
    S, D, P = stix_parameters(ω, B, n_e)

    cosθ = cos(θ_B)
    sinθ = sin(θ_B)

    # θ ≈ 0: quasi-longitudinal, modes are circularly polarised
    if abs(sinθ) < 1e-8
        w1 = (0.0, 0.0, 1.0)   # extraordinary: right-circular (α = +1)
        w2 = (1.0, 0.0, 0.0)   # ordinary: left-circular (α = -1)
        return w1, w2
    end

    # θ ≈ π/2: quasi-transverse
    if abs(cosθ) < 1e-8
        w1 = (0.0, 1.0, 0.0)   # extraordinary: along B (α = 0)
        w2 = (0.5, 0.0, 0.5)   # ordinary: circular mix (α = ±1)
        return w1, w2
    end

    # Compute q from Eq. 25 (full formula, no approximation)
    # q = (P - S) sin²θ / (2D cosθ)
    if abs(D) < 1e-30
        # No magnetic splitting → linear polarisation
        w1 = (0.0, 1.0, 0.0)
        w2 = (0.5, 0.0, 0.5)
        return w1, w2
    end

    q = cold_polarization_parameter(ω, B, θ_B, n_e)

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
    polarization_weights_full(ω, B, θ_B, n_e) → (w1, w2)

DEPRECATED compatibility wrapper. Equivalent to
`polarization_weights_cold(ω, B, θ_B, n_e)`. Kept for backward compatibility
with callers that predate the vacuum-polarisation refactor (commit 4dc9a5d).
New code should call either `polarization_weights_cold` (no vacuum effects,
pure P&C 2003) or `polarization_weights_vacuum` (vacuum-polarised, Potekhin,
Lai & Chabrier 2004) explicitly.
"""
function polarization_weights_full(ω::Float64, B::Float64,
                                    θ_B::Float64, n_e::Float64)::NTuple{2, NTuple{3, Float64}}
    return polarization_weights_cold(ω, B, θ_B, n_e)
end

"""
    polarization_weights_vacuum(ω, B, θ_B, n_e) → (w1, w2)

Compute |e_{j,α}|² including the QED vacuum-polarization correction to the
real dielectric/permeability tensors. Vacuum coefficients (a, q, m) follow
Potekhin, Lai & Chabrier (2004) Appendix Eqs. (A7)-(A9); the B-frame cyclic
mode-vector construction follows van Adelsberg & Lai (2006) Eqs. (16)-(18)
(equivalent to PLC 2004 Eqs. (20)-(22) / Ho & Lai 2003).

Falls back to `polarization_weights_cold` for B < 10⁶ G, near the analytic
limits (sinθ → 0, cosθ → 0, D → 0), and when the vacuum-shifted permittivities
εp = S + a or ηp = P + a + q approach zero.

Each of `w1`, `w2` is an `NTuple{3, Float64}` (allocation-free).
"""
function polarization_weights_vacuum(ω::Float64, B::Float64,
                                     θ_B::Float64, n_e::Float64)::NTuple{2, NTuple{3, Float64}}
    if B < 1e6
        return polarization_weights_cold(ω, B, θ_B, n_e)
    end

    S, D, P = stix_parameters(ω, B, n_e)

    cosθ = cos(θ_B)
    sinθ = sin(θ_B)

    # Keep the analytic limiting weights where the β expression is singular.
    if abs(sinθ) < 1e-8 || abs(cosθ) < 1e-8 || abs(D) < 1e-30
        return polarization_weights_cold(ω, B, θ_B, n_e)
    end

    a_hat, q_vac, m_vac = vacuum_coefficients(B)
    εp = S + a_hat
    ηp = P + a_hat + q_vac

    if abs(εp) < 1e-30 || abs(ηp) < 1e-30
        return polarization_weights_cold(ω, B, θ_B, n_e)
    end

    β, r, K1, K2, Kz1, Kz2 = vacuum_mode_parameters(ω, B, θ_B, n_e)

    w1 = compute_weights_from_K_Kz(K1, Kz1, sinθ, cosθ)
    w2 = compute_weights_from_K_Kz(K2, Kz2, sinθ, cosθ)

    return w1, w2
end

"""
    vacuum_coefficients(B) → (a_hat, q, m)

Vacuum polarization coefficients from Potekhin, Lai & Chabrier (2004)
Appendix Eqs. (A7)-(A9).  These fits recover the Euler-Heisenberg weak-field
limits and remain accurate for arbitrary B in the low-frequency limit.

Verified verbatim against PLC2004 Appendix A on 2026-05-16: each constant
(0.25487, 0.75, 1.2, 1.33, 0.56, 3.75, 2.7) and prefactor (-2α/9π, 7α/45π,
-α/3π) matches the published Eqs. (A7), (A8), (A9). Max relative error vs
the exact Heyl-Hernquist expression: 1.1% (A7), 2.3% (A8), 4.2% (A9).
"""
function vacuum_coefficients(B::Float64)::NTuple{3, Float64}
    b = B / B_Q
    if b == 0.0
        return 0.0, 0.0, 0.0
    end

    # PLC2004 Eq. (A7): â ≈ -(2α/9π) ln[1 + (b²/5)(1 + 0.25487 b^(3/4))/(1 + 0.75 b^(5/4))]
    a_hat = -2.0 * α_fine / (9π) *
            log(1.0 + (b^2 / 5.0) *
                (1.0 + 0.25487 * b^(3/4)) /
                (1.0 + 0.75 * b^(5/4)))
    # PLC2004 Eq. (A8): q ≈ (7α/45π) b² (1 + 1.2 b)/(1 + 1.33 b + 0.56 b²)
    q = 7.0 * α_fine / (45π) * b^2 *
        (1.0 + 1.2 * b) /
        (1.0 + 1.33 * b + 0.56 * b^2)
    # PLC2004 Eq. (A9): m ≈ -(α/3π) b² / (3.75 + 2.7 b^(5/4) + b²)
    m = -α_fine / (3π) * b^2 /
        (3.75 + 2.7 * b^(5/4) + b^2)

    return a_hat, q, m
end

"""
Cold-plasma polarization parameter q used in the P&C 2003 mode quadratic.

Contract: the expression q = (P - S) sin²θ / (2 D cosθ) is singular when
`D → 0` (unmagnetised limit) or `cosθ → 0` (strictly transverse propagation).
In those regimes the mode quadratic degenerates and callers must dispatch on
the analytic limit themselves before calling this function — `polarization_weights_cold`
does so at lines ~97/104/112 and `polarization_weights_vacuum` at line ~176.

Rather than returning a signed `Inf` (which was correct-by-accident: downstream
`isfinite` guards happened to catch it), this routine now throws a
`DomainError` so misuse is loud rather than silent.
"""
function cold_polarization_parameter(ω::Float64, B::Float64,
                                     θ_B::Float64, n_e::Float64)::Float64
    S, D, P = stix_parameters(ω, B, n_e)
    cosθ = cos(θ_B)
    sinθ = sin(θ_B)
    if abs(D) < 1e-30 || abs(cosθ) < 1e-30
        throw(DomainError((D, cosθ),
            "cold_polarization_parameter is singular when |D| < 1e-30 or " *
            "|cosθ| < 1e-30; caller must handle the analytic limit " *
            "(see polarization_weights_cold / polarization_weights_vacuum)."))
    end
    return (P - S) * sinθ^2 / (2.0 * D * cosθ)
end

"""
    vacuum_polarization_parameter(ω, B, θ_B, n_e) → β

Mode polarization parameter including vacuum dielectric and magnetic
permeability corrections.  This is the β-like parameter in Potekhin, Lai &
Chabrier (2004), written with the sign convention of the local cold-plasma
quadratic.
"""
function vacuum_polarization_parameter(ω::Float64, B::Float64,
                                       θ_B::Float64, n_e::Float64)::Float64
    S, D, P = stix_parameters(ω, B, n_e)
    a_hat, q_vac, m_vac = vacuum_coefficients(B)

    εp = S + a_hat
    ηp = P + a_hat + q_vac
    cosθ = cos(θ_B)
    sinθ = sin(θ_B)

    if abs(D) < 1e-30 || abs(cosθ) < 1e-30 || abs(ηp) < 1e-30
        return cold_polarization_parameter(ω, B, θ_B, n_e)
    end

    numerator = ηp - εp + D^2 / εp + ηp * m_vac / (1.0 + a_hat)
    return numerator / (2.0 * D) * (εp / ηp) * sinθ^2 / cosθ
end

"""
    vacuum_mode_parameters(ω, B, θ_B, n_e) → (β, r, K1, K2, Kz1, Kz2)

Vacuum-corrected normal-mode parameters from Potekhin, Lai & Chabrier
(2004).  The `K1/K2` branch labels follow their X/O convention,
`K_j = β[1 + (-1)^j sqrt(1 + r/β²)]`, evaluated in a numerically stable
form away from and at β=0.
"""
function vacuum_mode_parameters(ω::Float64, B::Float64,
                                θ_B::Float64, n_e::Float64)::NTuple{6, Float64}
    S, D, P = stix_parameters(ω, B, n_e)
    a_hat, q_vac, m_vac = vacuum_coefficients(B)
    sinθ = sin(θ_B)
    cosθ = cos(θ_B)

    β = vacuum_polarization_parameter(ω, B, θ_B, n_e)
    r = 1.0 + m_vac / (1.0 + a_hat) * sinθ^2
    r = max(r, 0.0)
    sq = sqrt(β^2 + r)

    if β > 0
        K1 = β - sq
        K2 = β + sq
    elseif β < 0
        K1 = β + sq
        K2 = β - sq
    else
        K1 = -sqrt(r)
        K2 = sqrt(r)
    end

    Kz1 = _vacuum_Kz(K1, S, D, P, a_hat, q_vac, sinθ, cosθ)
    Kz2 = _vacuum_Kz(K2, S, D, P, a_hat, q_vac, sinθ, cosθ)

    return β, r, K1, K2, Kz1, Kz2
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
        return (0.0, 1.0, 0.0)
    end

    Kz = -n² * sinθ * cosθ * K / denom_z

    return compute_weights_from_K_Kz(K, Kz, 0.0, 1.0)
end

function _vacuum_Kz(K::Float64, S::Float64, D::Float64, P::Float64,
                    a_hat::Float64, q_vac::Float64,
                    sinθ::Float64, cosθ::Float64)::Float64
    εp = S + a_hat
    ηp = P + a_hat + q_vac
    denom = εp * sinθ^2 + ηp * cosθ^2
    if abs(denom) < 1e-30 * (1.0 + abs(εp) + abs(ηp))
        return copysign(Inf, -((εp - ηp) * K * cosθ + D) * sinθ)
    end
    return -((εp - ηp) * K * cosθ + D) / denom * sinθ
end

function compute_weights_from_K_Kz(K::Float64, Kz::Float64,
                                   sinθ::Float64, cosθ::Float64)::NTuple{3, Float64}
    # Cyclic mode-vector weights in the B-frame, from van Adelsberg & Lai
    # (2006) Eqs. (22)-(23). The k-frame eigenvector (Eq. 17) is
    #   e_j = (1/√N)(i K_j, 1, i K_{z,j})
    # with N = 1 + |K_j|² + |K_{z,j}|². Rotating about y by -θ_kB so the
    # new z'-axis lies along B gives B-frame components
    #   e_x' ∝ K cosθ + Kz sinθ   (the "transverse" combination below)
    #   e_y' ∝ 1
    #   e_z' ∝ -K sinθ + Kz cosθ  (the "longitudinal" combination, sign
    #                              flipped here — only its square matters)
    # The cyclic projections e_± = (1/√2)(e_x' ± i e_y') then yield
    # |e_±|² = (transverse ± 1)² / (2N), |e_0|² = longitudinal² / N,
    # which matches vAL2006 (22)-(23) verbatim for real K, Kz.
    if !isfinite(K) || !isfinite(Kz)
        return (0.0, 1.0, 0.0)
    end

    norm = 1.0 + K^2 + Kz^2
    transverse = K * cosθ + Kz * sinθ
    longitudinal = K * sinθ - Kz * cosθ

    wp1 = (transverse + 1.0)^2 / (2.0 * norm)  # |e_{+1}|², vAL2006 Eq. (22), upper sign
    wm1 = (transverse - 1.0)^2 / (2.0 * norm)  # |e_{-1}|², vAL2006 Eq. (22), lower sign
    w0  = longitudinal^2 / norm                # |e_0|²,    vAL2006 Eq. (23)

    # Safety clamp and renormalize
    wp1 = max(0.0, wp1)
    wm1 = max(0.0, wm1)
    w0  = max(0.0, w0)
    s = wm1 + w0 + wp1
    if s > 0
        wm1 /= s; w0 /= s; wp1 /= s
    end

    return (wm1, w0, wp1)
end

# ---------------------------------------------------------------------------
# Vacuum-resonance mode conversion (P_jump)
#
# Source: van Adelsberg & Lai (2006) MNRAS 373, 1495 — Eqs. (2), (3), (4).
#         Suleimanov, Potekhin & Werner (2009) A&A 500, 891 — Eqs. (13), (16),
#         (17) (which reproduce the vAL2006 formulae verbatim, with the angle
#         labelled θ_kB ≡ θ_B in our codebase).
# ---------------------------------------------------------------------------

"""
    vacuum_resonance_density(ω, B, θ_B; Y_e=1.0, f_B=1.0) → ρ_V [g/cm³]

Density at which the plasma and vacuum contributions to the dielectric
tensor cancel, i.e. the vacuum-resonance crossing point.

This is van Adelsberg & Lai (2006) Eq. (2):

```
ρ_V = 0.96 · Y_e⁻¹ · E_1² · B_14² · f_B⁻²   [g cm⁻³]
```

where `E_1 = E/(1 keV)`, `B_14 = B/(10¹⁴ G)`, and `f_B ≈ 1` is a slowly
varying function of `B` (vAL2006 references Lai & Ho 2003 Eq. 2.41 for the
exact expression). We take `f_B = 1` by default, an approximation that is
accurate to ≲ 10% across `B = 10¹² – 10¹⁵ G` and dominates none of the
subsequent P_jump physics. SPW09 Eq. (13) uses the same expression
through their `E_V(ρ)` relation.

The resonance density does **not** depend on `θ_B`; `θ_B` is accepted in
the signature for symmetry with `vacuum_resonance_pjump`.

Inputs are CGS: `ω` [rad s⁻¹], `B` [G]. Returns `ρ_V` [g/cm³] ≥ 0.
"""
function vacuum_resonance_density(ω::Float64, B::Float64, θ_B::Float64;
                                  Y_e::Float64=1.0, f_B::Float64=1.0)::Float64
    @assert ω > 0 "ω must be positive, got $ω"
    @assert B >= 0 "B must be non-negative, got $B"
    @assert Y_e > 0 "Y_e must be positive, got $Y_e"
    @assert f_B > 0 "f_B must be positive, got $f_B"

    E_keV = ħ * ω / keV       # photon energy in keV
    B_14 = B / 1.0e14
    ρ_V = 0.96 / Y_e * E_keV^2 * B_14^2 / f_B^2
    @assert isfinite(ρ_V) "vacuum_resonance_density returned non-finite for ω=$ω B=$B"
    return ρ_V
end

"""
    vacuum_resonance_pjump(ω, B, θ_B, T, dρ_dy; Z=1.0, A=1.0, f_B=1.0) → P_jump ∈ [0, 1]

Probability that a photon undergoes mode conversion (X ↔ O) when it
crosses the vacuum resonance at density `ρ_V`. The plane-parallel
column-depth gradient is supplied as `dρ_dy` [(g cm⁻³)/(g cm⁻²)] = [cm⁻¹]
at the resonance layer; the density scale length along the ray (normal
incidence, `cos θ_e = 1` for the column build) is then `H_ρ = 1/|dρ/dy|`
in cm (because `dz = dy/ρ` and `ρ d ln ρ/dz = dρ/dz = ρ · dρ/dy`, so
`H_ρ ≡ |d ln ρ/dz|⁻¹ = 1/|dρ/dy|`).

This is van Adelsberg & Lai (2006) Eqs. (3)-(4):

```
E_ad   = 2.52 · [ f_B · tan θ_kB · |1 - (E_Bi/E)²| ]^(2/3) · (1 cm / H_ρ)^(1/3)   [keV]
P_jump = exp[ -π/2 · (E / E_ad)³ ]
```

with the ion-cyclotron energy `E_Bi = 0.63 · (Z/A) · B_14` keV. SPW09 use
the same formulae (Eqs. 16-17). `Z = A = 1` for fully ionised hydrogen.

Limits:
  * `dρ_dy → 0` (very gentle gradient) ⇒ `H_ρ → ∞` ⇒ `E_ad → 0` ⇒
    `P_jump → 0` (fully adiabatic, modes preserved).
  * `dρ_dy → ∞` (discontinuous jump)   ⇒ `H_ρ → 0` ⇒ `E_ad → ∞` ⇒
    `P_jump → 1` (fully non-adiabatic, modes swap).
  * `θ_B → 0` ⇒ `tan θ_B → 0` ⇒ `E_ad → 0` ⇒ `P_jump → 0` (only photons
    propagating exactly along B see no mode mixing).

Returns 0 for non-physical inputs (`dρ_dy ≤ 0`, NaN gradient, `θ_B = 0`
or `π`) — caller treats this as "no resonance crossing in this layer".
"""
function vacuum_resonance_pjump(ω::Float64, B::Float64, θ_B::Float64,
                                T::Float64, dρ_dy::Float64;
                                Z::Float64=1.0, A::Float64=1.0,
                                f_B::Float64=1.0)::Float64
    @assert ω > 0 "ω must be positive, got $ω"
    @assert B >= 0 "B must be non-negative, got $B"
    @assert T > 0 "T must be positive, got $T"
    @assert isfinite(dρ_dy) "dρ_dy must be finite, got $dρ_dy"

    # Degenerate angle: tan θ_kB = 0 (along-B) or undefined at π/2; modes
    # don't mix at exactly θ_B = 0 (vAL2006 Eq. 3 prefactor vanishes).
    if abs(sin(θ_B)) < 1e-12 || abs(cos(θ_B)) < 1e-12
        return abs(cos(θ_B)) < 1e-12 ? 1.0 : 0.0
        # cos θ_B = 0 (strictly transverse): tan θ_B → ∞ ⇒ E_ad → ∞ ⇒ P_jump → 1.
        # sin θ_B = 0 (along B):              tan θ_B = 0 ⇒ E_ad → 0   ⇒ P_jump → 0.
    end

    # No physical resonance crossing in this layer (caller upstream should
    # have already screened on y_V being inside the grid; this is defence
    # in depth).
    if !(dρ_dy > 0) || !isfinite(dρ_dy)
        return 0.0
    end

    H_ρ = 1.0 / dρ_dy           # cm; derivation in the docstring.
    @assert H_ρ > 0 && isfinite(H_ρ) "H_ρ must be positive finite, got $H_ρ"

    E_keV = ħ * ω / keV
    B_14 = B / 1.0e14
    E_Bi = 0.63 * (Z / A) * B_14          # ion cyclotron energy, keV
    cyc_factor = abs(1.0 - (E_Bi / E_keV)^2)

    tanθ = tan(θ_B)
    # vAL2006 Eq. 3: argument of (·)^(2/3) is always ≥ 0; floor at 0 in
    # case of round-off at the cyclotron itself (E = E_Bi exactly).
    arg = f_B * abs(tanθ) * cyc_factor
    if arg <= 0
        # E_ad → 0 ⇒ (E/E_ad)³ → ∞ ⇒ P_jump → 0 (perfectly adiabatic).
        return 0.0
    end

    E_ad = 2.52 * arg^(2.0/3.0) * (1.0 / H_ρ)^(1.0/3.0)   # keV
    @assert isfinite(E_ad) && E_ad >= 0 "E_ad must be finite & non-negative, got $E_ad"

    if E_ad <= 0
        return 0.0
    end

    # vAL2006 Eq. 4 / SPW09 Eq. 17.
    exponent = -0.5 * π * (E_keV / E_ad)^3
    P = exp(exponent)
    # exp() can underflow to 0 cleanly; clamp tiny round-off above 1.
    P = clamp(P, 0.0, 1.0)
    @assert isfinite(P) "P_jump returned non-finite, got $P (E_ad=$E_ad)"
    return P
end

end # module
