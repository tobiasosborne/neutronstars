# v1 dependency audit (rule 7: own the pipeline)

Date: 2026-10-23. Scope: production code paths of `archive/v1/src/**` and
`archive/v1/scripts/**` (tests and `verification/` ignored). The two
production entry points are `scripts/render_rxj1856_visible_atmosphere.jl`
(B = 0 H atmosphere) and `scripts/render_rxj1856_visible_magnetic_atmosphere.jl`
(B = 1e12 G), plus the library renderers `render_neutron_star` /
`render_spectral_cube` in `src/pipeline/render.jl`. All `file:line`
references are relative to `archive/v1/` unless they start with `refs/`.

Provenance classes follow `truth/SCHEMA.md` (measured, convention,
derivable, fit, table, code). One extra label, **placeholder**, marks
"typical" numbers that are none of these; it is not a valid `truth/` class and
every placeholder must be replaced or dropped.

Classification rule used: an item is *derivable* if v1 evaluates a formula
that follows from stated physics (even if under an approximation such as the
Born approximation); it is *fit* if the runtime number is a published
coefficient that approximates someone else's calculation; *table* / *code* if
the runtime number is someone else's numerical output.

## Part A — what v1 consumed at runtime

| # | file:line | what it is | class | physical calculation it stands in for / notes |
|---|-----------|------------|-------|-----------------------------------------------|
| A1 | `src/constants.jl:10-62`, `src/opacity/gaunt_ff.jl:20` | CODATA constants c, G, h, ħ, k_B, m_e, m_p, m_u, e, σ_T, α, r_e, eV, σ_SB, Ry | measured | Irreducible. Defects: `e_charge = 4.80320451e-10` (l.41) is the pre-2018 value (CODATA 2018 gives e·c/10 = 4.803204712570263e-10 esu, Δ = 5e-8); σ_SB is typed in, not derived from h, k, c. |
| A2 | `src/constants.jl:54` | M_sun = 1.98892e33 g, docstring "IAU 2015 nominal" | measured | Docstring is false: IAU 2015 B3 fixes GM_sun, giving M_sun = GM_sun/G ≈ 1.9884e33 g (Δ ≈ 2.6e-4). Cancels in both RX J1856 scripts (M is built from M_sun and divided back). |
| A3 | `src/constants.jl:56`, `src/utils/units.jl:93-103` | pc, au, mas | convention | IAU exact definitions. Irreducible. |
| A4 | `scripts/render_rxj1856_visible_atmosphere.jl:54`, `scripts/render_rxj1856_visible_magnetic_atmosphere.jl:54`, `src/opacity/gaunt_ff.jl:36-124` | **McPHAC copy of Sutherland (1998) `gffgu.dat`** (`refs/code/McPHAC/gffgu.dat`, 3321 rows), bilinear interpolation, silent clamp to u, γ² ∈ [1e-4, 1e4] (l.101-103) | table | Thermally averaged non-relativistic free-free Gaunt factor ḡ_ff(u, γ²) from exact Coulomb (Sommerfeld) wavefunctions. Column labels of the published table are swapped (`refs/haakonsen_2012_mcphac.txt` footnote 2). In the magnetic script the table only seeds the initial (y, P) column (`src/atmosphere/magnetic_atmosphere.jl:469-512`); B > 0 opacities do not read it. |
| A5 | `src/opacity/hydrogen_ff.jl:165-179` | Kramers free-free prefactor, Haakonsen+2012 Eq. 12 | derivable | Classical bremsstrahlung × ḡ_ff. Owned once A4 is owned. |
| A6 | `src/opacity/hydrogen_ff.jl:188-190`, `src/opacity/magnetic_ff.jl:318-334` | Thomson σ_T/m_H; magnetic electron scattering PC03 Eq. 33 | derivable | Classical damped-oscillator scattering. Note: the local PC03 is ApJ 585, 955; `refs/INDEX.md` mislabels it as A&A 399, 1007. |
| A7 | `src/atmosphere/atm_structure.jl:160`, used again in `src/atmosphere/magnetic_atmosphere.jl:509` | Ideal fully ionised H gas ρ = m_p P/(2 k_B T), also for B > 0 | derivable | Approximation: no Coulomb non-ideality, no Landau quantisation of the electron gas, no partial ionisation. Mixes m_p here with m_H in the opacities (`magnetic_ff.jl:234`). |
| A8 | `src/atmosphere/blackbody.jl:29-40`, `src/atmosphere/feautrier.jl:63-75` | Planck function; hard-coded Gauss–Legendre half-range nodes and weights | derivable | GL nodes are already recomputed by eigenvalues at `feautrier.jl:506`; the hard-coded table is redundant. |
| A9 | `src/atmosphere/feautrier.jl:463`, `src/atmosphere/magnetic_atmosphere.jl:779,1022`, `src/opacity/magnetic_modes.jl:53-58` | Azimuth-averaged Rayleigh phase function 3/16(3+3μ²μ'²−μ²−μ'²); per-mode scattering with no X↔O redistribution | derivable | Correct for B = 0 Thomson. For B > 0 it is an approximation: magnetic scattering is mode- and angle-dependent and couples the modes. |
| A10 | `src/atmosphere/atm_structure.jl:197-202` | RT frequency grid 0.05–120 k_B T_eff/h (McPHAC MINFREQKT, MAXFREQKT) | convention | Numerical choice. For RX J1856 (T_eff ≈ 5.2e5 K) it starts at ≈ 550 nm, which caused the v1 visible-band failure (HANDOFF). |
| A11 | `src/atmosphere/rt_atmosphere.jl:275-292`, `src/atmosphere/magnetic_atmosphere.jl:411-429`, `scripts/render_rxj1856_visible_atmosphere.jl:79` | Below the grid: Rayleigh–Jeans ν² extrapolation; above the grid: 0; μ clamped to [μ_1, μ_M] and cos α clamped to ≥ 1e-6 | derivable | Silent extrapolation in place of solving the RT at the needed ν and μ. Fix: put those ν in the RT grid and do a formal solution at arbitrary μ. |
| A12 | `src/atmosphere/solver_defaults.jl:25`, used in `atm_structure.jl:338` | T(0) = 0.265 T_eff, "empirical (McPHAC)" | code | Initial-guess seed taken from McPHAC behaviour. It should not affect a converged model; that has to be tested. |
| A13 | `src/atmosphere/solver_defaults.jl:19-43`, `rt_atmosphere.jl:159`, `magnetic_atmosphere.jl:199`, `coulomb_magnetic.jl:35,45-46,64`, `magnetic_ff.jl:230,293,319`, `dielectric_tensor.jl:167` | Solver numerics: τ_diff = 80, τ_eff,min = 10, y_max = 1e5 g cm⁻², damping 0.5, flux clamps 0.9/1.1, \|ΔT/T\| ≤ 0.3, opacity floor 1e-30, T floor 0.1 T_eff, B < 1e6 G → non-magnetic, β_e < 0.1 → classical, Landau sum capped at 50, Λ floor 0.01 | convention | Numerical controls, not physics. Several are silent clamps (rule 5): T floor, Λ floor, flux clamp. |
| A14 | `src/atmosphere/magnetic_atmosphere.jl:160-180` | "Ad-hoc grey-body flux scaling" (code comment), not SPW09 | convention | Iteration heuristic. It shapes the answer whenever the model is not fully converged. The production B = 1e12 G run uses tol = 7e-4 (script l.59). |
| A15 | `src/opacity/coulomb_magnetic.jl:128-135` | Classical Coulomb log Λ_cl = e^{u/2} K₀(u/2), PC03 Eq. 43 | derivable | Born Maxwell average of PC03 Eq. 42. |
| A16 | `src/opacity/coulomb_magnetic.jl:31-120` | Magnetic Coulomb log Λ_α^ff, PC03 Eq. 44a–e, evaluated by quadgk | derivable | v1 computes this itself; it is not a fit. It is Born-only, though. v1 omits the non-Born rescaling by ḡ_ff(Hummer)/ḡ_Born and the finite-proton-mass quantum correction factor that PC03 applies to Eqs. 41/51 (PC03 §4.3.2, local). As a result the B → 0 limit of the magnetic path (Born Λ_cl) differs from the B = 0 path (Sutherland table). |
| A17 | `src/opacity/magnetic_ff.jl:158-311` | Magnetic ff, pp and damping, PC03 Eqs. 37–41, 46–47, 51–53 | derivable | Classical cold-plasma plus PC03 Appendix A/B results, written as closed forms. |
| A18 | `src/opacity/magnetic_ff.jl:173-175` | **Λ_pp ≈ 0.6 ln(22/u + 9u^-0.3) + 0.4√(πu)**, PC03 Eq. 48 | fit | Fit (max error 1.5%) to the Maxwell average of the proton–proton quadrupole Born cross-section, PC03 App. A Eqs. A4–A6. |
| A19 | `src/opacity/dielectric_tensor.jl:45-389`, `src/opacity/magnetic_modes.jl:37-89` | Cold-plasma Stix tensor (e and p), mode parameter q/β, K_j, K_z,j, cyclic weights \|e_j,α\|² | derivable | PC03 Eq. 25–27, vAL06 Eq. 10, 16–23. Weights are renormalised after a clamp (l.380-386). |
| A20 | `src/opacity/dielectric_tensor.jl:209-229` | **Vacuum-polarisation coefficients â, q, m**, PLCH04 Eqs. A7–A9 | fit | Fit (max error 1.1/2.3/4.2%) to the Heyl–Hernquist one-loop Euler–Heisenberg result (PLCH04 Eqs. A2–A6, local). The docstring's error figures are copied from the paper; v1 never ran a comparison. |
| A21 | `src/opacity/dielectric_tensor.jl:424-436` | **ρ_V = 0.96 Y_e⁻¹ E₁² B₁₄² f_B⁻²**, f_B := 1 | fit | Rounded closed form of the vacuum-resonance condition from Lai & Ho. f_B (Ho & Lai 2003 Eq. 2.41) is replaced by 1. Optional path: `apply_pjump = false` by default (`magnetic_atmosphere.jl:100`). |
| A22 | `src/opacity/dielectric_tensor.jl:482-539` | **E_ad = 2.52(f_B tanθ\|1−(E_Bi/E)²\|)^{2/3}(1 cm/H_ρ)^{1/3}**, E_Bi = 0.63 B₁₄ keV, P_jump = exp[−(π/2)(E/E_ad)³] | fit | Rounded Landau–Zener adiabaticity result (Lai & Ho 2002). Same optional path as A21. |
| A23 | `src/atmosphere/magnetic_atmosphere.jl:280-312` | Depth weighting of the X↔O swap by exp(−τ_V/μ) ("Eddington–Barbier") | derivable | v1's own post-processing approximation of mode conversion inside the RT. Optional (as A21). |
| A24 | `src/eos/bsk_eos.jl:26-104` | BSk19/20/21 coefficients a₁…a₂₃, Potekhin+2013 Eq. 3 / Table 2 | fit | **Not on any render path.** `render.jl:15` imports `solve_tov` but never calls it. |
| A25 | `src/eos/tov.jl` | TOV integrator | derivable | Not on any render path. |
| A26 | `src/eos/potekhin_table_reader.jl` | Reader for Potekhin magnetised H EOS and opacity tables | table | **Not on a production path.** It is not `include`d by `src/NeutronStar.jl` (header l.4-8), so it is verification-only. |
| A27 | `src/constants.jl:70` | ρ_nuc = 2.8e14 g cm⁻³ | placeholder | Unused anywhere in `src/`. |
| A28 | `src/surface/dipole.jl:76-79` | \|B\| = (B_p/2)√(1+3cos²θ_B) | derivable | Flat-space centred dipole. The GR correction to the dipole field at r ≈ 3 r_g is ignored. |
| A29 | `src/surface/dipole.jl:57-66` | **T(θ_B) = T_eq + (T_p − T_eq) cos²θ_B** | fit | An unsourced two-parameter ansatz; the v1 docstring itself says it is not Greenstein–Hartke. It stands in for T_s(θ) from anisotropic heat conduction through a magnetised envelope. |
| A30 | `src/pipeline/render.jl:351-369` | Magnetic colatitude θ_B passed to the atmosphere as "angle between B and the surface normal" | derivable | Computed wrongly. For a dipole the field–normal angle Θ obeys tanΘ = ½ tanθ_B (follows from the B_r, B_θ in CLAUDE.md v1 Module 4). |
| A31 | `src/geodesics/schwarzschild.jl:1-16, 65-200`; `b_max = 1.05 b_crit` l.149 | Beloborodov (2002) cosine relation, parallel rays from infinity | derivable | Approximation (~1%) to the exact bending integral. The header claims "Also implements exact integral … Eq. (9)"; **that claim is false**, because no such function exists in the file. PRD R1 requires exact geodesics from a finite moving camera. |
| A32 | `scripts/*:36-48`, `src/pipeline/render.jl:131-138, 320-371` | 1+z = (1−u)^{-1/2}, I_obs = I_em/(1+z)³, M from (z, R∞), g_s = GM(1+z)/R² | derivable | Correct. |
| A33 | `src/pipeline/render.jl:381` | pixel_scale = 2R/(d N), "approximate" | derivable | Wrong. The image half-width is b_max = 1.05 R/√(1−u), not R. |
| A34 | `scripts/render_rxj1856_visible_atmosphere.jl:36-41` and magnetic script l.31-42 | **RX J1856: T∞ = 4.3e5 K, R∞ = 17 km, z_g = 0.22** (Ho+2007) | code | These are outputs of Ho's thin magnetic H-atmosphere code fitted to Chandra LETG data **at d = 140 pc** (`refs/ho_2007_rxj1856_thin_atm.txt` l.39-41, Table 1), for B = 4e12 G. v1 did not rescale to the measured 123 pc (`walter_2010`), so R∞ would become ≈ 14.9 km. They stand in for our own fit of our own surface model to the measured spectrum and photometry. |
| A35 | magnetic script l.47-48 | B = 1e12 G, θ_B = π/4 | placeholder | Chosen "because it converges" (comment l.44-46). Ho+2007 gives (3–4)e12 G; spin-down gives 1.5e13 G (`obs.rxj1856-dipole-field`). |
| A36 | `src/pipeline/render.jl:34-42` | NSParams example: 1.4 M_☉, 12 km, 1e12 G, T_p = 1.5e6 K, T_eq = 3e5 K, α = 0.3, 100 pc | placeholder | "Typical NS" numbers with no source. |
| A37 | `src/atmosphere/blackbody.jl:54-65`, `render.jl:134` (p = 0 hard-coded) | Colour correction f_col and limb-darkening exponent p ("Suleimanov+ for f_col") | fit | Defaults f_col = 1, p = 0, so the legacy path is in effect a pure blackbody. f_col is a fit to atmosphere-model output; p is unsourced. |
| A38 | `src/data/cvrl_cie1931_2deg.csv`, loaded at `scripts/*:70` | CIE 1931 2° CMFs, 360–830 nm at 1 nm | convention | Definitional standard table. Irreducible. |
| A39 | `src/colorimetry/cie_srgb.jl:135-137` | XYZ→linear-sRGB matrix (3.2406255, −1.537208, …) "Exact IEC" | convention | The values match neither the local source (`refs/stokes_1996_srgb_equations.txt`: 3.2410, −1.5374, …) nor the 4-digit IEC print. They are the 7-digit inverse of the primaries/D65 matrix. Derive them from the primaries and white point instead (convention, local in `truth/entries/colour*.toml`). |
| A40 | `src/colorimetry/cie_srgb.jl:147-150` | sRGB transfer 12.92 / 0.0031308 / 1.055 / 2.4 | convention | Irreducible. |
| A41 | `src/colorimetry/cie_srgb.jl:163-166`; key = 0.36 at `scripts/*:97/90`, `render.jl:164,421` | Reinhard global operator, key 0.36 (paper default 0.18), L_avg = log-average | convention | Display choice. PRD R6 requires it to be documented and reversible. Per-channel clamp at l.198-200 caused hue banding. |
| A42 | `scripts/*:65`; `render.jl:114-116,180-192` | Visible grid 360–830 nm × 180, cube grid 0.01 eV–100 keV, false-colour bands 0.1/0.5/2/10 keV | convention | Output sampling choices. |
| A43 | `src/pipeline/atmosphere_grid.jl:329, 490-527` | Grid lookup clamps T, B, θ_B to the nearest edge outside the grid | convention | Silent extrapolation (rule 5). |

Counts (43 rows): measured 2, convention 10, derivable 17, fit 7, table 2,
code 2, placeholder 3. Fit/table/code items that actually ran in the
production RX J1856 renders: A4 (Gaunt table), A12 (T(0) seed), A18 (Λ_pp,
magnetic script), A20 (vacuum coefficients, magnetic script), A34 (Ho+2007
parameters). A21–A23 only run with `apply_pjump = true`. A24 and A26 never
ran in production. A29 and A37 ran only in the library renderer
`render_spectral_cube` / `render_neutron_star`.

## Part B — route to ownership

Policy: DROP only when the flyby in PRD.md does not consume the quantity.
Every other item gets a layer-by-layer RECOMPUTE chain. Each layer gives its
local source (or "missing", with the paper to fetch), a rough size in lines of
code (LoC), and the verification number it can be checked against. arXiv IDs
marked "(confirm)" are from memory and must be checked when the paper is
fetched.

| item | route | dependency chain (layer: physics — source — size — verification) |
|------|-------|-------------------------------------------------------------------|
| A4 Gaunt table | RECOMPUTE | **L0** ḡ_ff(u, γ²) = Maxwell average over initial electron energy of g_ff(η_i, η_f) — usage local (`haakonsen_2012` Eq. 12), averaging integral: missing, fetch van Hoof et al. 2014, MNRAS 444, 420, arXiv:1407.1083 — 60 LoC — verify against `refs/code/McPHAC/gffgu.dat` (e.g. γ² = 1e-4, u = 1e-4 → 5.524282; becomes a verification fixture, with the column swap noted).<br>**L1** exact non-relativistic dipole g_ff from Sommerfeld's formula, hypergeometric ₂F₁ with complex parameters (η = Z/ka₀); needs recurrences or arbitrary precision at large η — missing: Karzas & Latter 1961, ApJS 6, 167 and van Hoof 2014 §2 — 200 LoC — verify the Born limit g → (√3/π) ln[(p_f+p_i)/(p_f−p_i)], local in PC03 Eq. 42, and the thermal Born limit ḡ_Born = (√3/π)Λ_cl (PC03 Eq. 43, local).<br>**L2** Coulomb continuum wavefunctions and dipole matrix element — missing: Bethe & Salpeter 1957 (book); Born limit local (PC03 Eq. 42) — 0 extra LoC if L1 uses the closed form; 150 LoC for an independent numerical check by ODE integration of the radial Coulomb functions. The same code also supplies the non-Born ratio ḡ_ff/ḡ_Born that PC03 applies to Λ_α^ff, which replaces the Hummer 1988 (ApJ 327, 477) Padé fit. |
| A12 T(0) = 0.265 T_eff | RECOMPUTE | **L0** grey Eddington profile T⁴ = ¾T_eff⁴(τ + ⅔), T(0) = 2^{-1/4} T_eff — local (`haakonsen_2012` Eq. 10 context; `atm_structure.jl:326` already states it) — 5 LoC — verify that the converged spectrum is unchanged (to tolerance) when the seed changes from 0.265 to 0.841 T_eff. Exact Hopf value: missing, Mihalas 1978 (book, manual download). |
| A16 (omitted corrections to Λ_α^ff) | RECOMPUTE | **L0** non-Born rescaling Λ_α → Λ_α ḡ_ff/ḡ_Born — local (PC03 p. 965 text) — 10 LoC once A4 is owned.<br>**L1** finite-proton-mass Λ^e_α, PC03 Eqs. B19–B22 (thermal average over Landau states n, s and k), replacing the unstated "correction factor of order unity" for Eqs. 41/51 — local — 200 LoC — verify against PC03 Fig. 6 curves (local; values read off the figure, so tolerance should be loose).<br>**L2** Coulomb matrix elements ṽ_{ns,n's'} between Landau states (PC03 Eqs. B6–B8) — local, but they cite Potekhin 1994, J. Phys. B 27, 1073 Eqs. A3–A10 (missing; fetch) — 150 LoC — verify the α = 0 case reproduces Eq. 44 (PC03 states it "well reproduces"; local).<br>**L3** Landau wavefunctions Φ_ns — local (`ruder_1994`) — 50 LoC — verify normalisation and orthogonality analytically. |
| A18 Λ_pp fit | RECOMPUTE | **L0** Maxwell average of σ_pp = 2σ*_pp − σ_× (PC03 Eqs. A4–A6, with F_pp(p_i) as stated after A6) — local — 60 LoC (1-D quadrature) — verify against PC03 Eq. 48 to ≤ 1.5% (local) and the asymptotes PC03 states.<br>**L1** quadrupole expansion of the two-proton current j_eff (Eqs. A1–A3) — local — analytic, 0 LoC.<br>Note: PC03 Eq. 54 (local) shows pp collisions are negligible in outer layers, so this is ownership debt, not accuracy debt. |
| A20 PLCH04 A7–A9 | RECOMPUTE | **L0** closed form â, â+q, m from X(ξ) = 2 lnΓ(ξ/2) − … (PLCH04 Eqs. A2–A6) — local — 60 LoC (lgamma, digamma, one quadrature ∫X dξ) — verify the weak-field limit A1 (â = −2α_f b²/45π, q = 7α_f b²/45π, m = −4α_f b²/45π, local) and the A7–A9 maxima (1.1% at b = 0.07, 2.3% at b = 0.4, 4.2% at b = 0.3, local). The constant 0.0329199 in A2 must be re-derived (it should follow from requiring â → 0 as b → 0), not copied.<br>**L1** one-loop Euler–Heisenberg effective Lagrangian (proper-time integral) → χ_vac, μ tensors — missing: Heyl & Hernquist 1997, J. Phys. A 30, 6485 (as cited by PLCH04) and Euler & Heisenberg 1936, Z. Phys. 98, 714 — 100 LoC (direct numerical proper-time integral as an independent check of L0) — verify against Kohri & Yamada 2002 Fig. 2 (missing; fetch) and L0.<br>Validity: ħω ≪ m_e c² (local, PLCH04 App.), which is satisfied for optical to soft X-ray. |
| A21 ρ_V, f_B | RECOMPUTE | **L0** resonance density: the root in ρ where the vacuum-corrected mode parameter (our A19 plus A20 tensor) makes the modes circular / the plasma and vacuum terms cancel — local statement (vAL06 Eq. 1–2 and text; PLCH04 mode equations) — 40 LoC (root find on our own tensor) — verify against vAL06 Eq. 2 with f_B → 1 for b ≪ 1.<br>**L1** f_B as a function of B — missing: Ho & Lai 2003, MNRAS 338, 233, Eq. 2.41, arXiv:astro-ph/0201380 (confirm) — only needed as a verification number, because L0 makes f_B implicit. |
| A22 E_ad, P_jump, E_Bi | RECOMPUTE | **L0** E_Bi = ħeB/(m_p c) — derivable from CODATA — 1 LoC — verify vAL06's 0.63 B₁₄ keV (local).<br>**L1** Landau–Zener jump probability P = exp(−πγ_ad/2), with γ_ad computed from our own tensor (mode splitting at the crossing and the derivative of the mode parameter along the ray, with H_ρ from our atmosphere column) — statement local (vAL06 Eq. 3–4, SPW09 Eq. 16–17); derivation missing: Lai & Ho 2002, ApJ 566, 373, arXiv:astro-ph/0108127 (confirm) — 80 LoC — verify vAL06's 2.52 coefficient in its stated limit.<br>**L2** independent check: integrate the two-mode amplitude evolution equations through the resonance — missing: Lai & Ho 2003 / Ho & Lai 2003 mode-evolution equations — 150 LoC — verify that L1 and L2 agree.<br>**L3** replace A23 (exp(−τ_V/μ) swap) with mode conversion inside the two-mode RT at the resonance depth, as vAL06 §3 does — local (`van_adelsberg_lai_2006`) — 150 LoC — verify vAL06 published spectra and figures (local). |
| A24 BSk fit (+ A25 TOV) | DROP | The flyby uses only the Schwarzschild exterior plus surface radiance. M and R enter only as (M, R) or (z, R∞), and both come from the spectral/photometric fit in A34 or from the uncertainty ensemble (PRD M5). g_s = GM(1+z)/R² follows from them. Interior structure is never consumed. **Caveat:** if B is ever taken from spin-down with our own R (A35), the moment of inertia I enters, and I needs interior structure. Then the chain is TOV (local `potekhin_2013`, `eos.tov-*` entries) ← EOS of dense matter, and BSk would have to be recomputed from the BSk energy-density functional (missing: Goriely, Chamel & Pearson 2013; very deep). The recommended way to keep it dropped is to treat B as a fitted parameter, not to derive it from spin-down. |
| A26 Potekhin tables | DROP (as runtime) | Never runtime. Keep it as a verification fixture under `verification/fixtures/` for the magnetic EOS and opacity work, consistent with the old rule. |
| A29 T(θ) ansatz | RECOMPUTE | **L0** local T_s at each surface patch from (T_b, B, Θ_Bn, g): 1-D plane-parallel envelope with constant flux, d lnT/d lnP = (3Pκ/16g)(T_s⁴/T⁴) (PPP2015 Eq. 38, local), with the anisotropic flux projection (GH 1983 Eqs. 6–10, local, OCR poor; PPP2015 §4.4, local); top boundary = our atmosphere (F = σT_s⁴); bottom boundary ρ_b = 1e10 g cm⁻³ at T_b — 150 LoC — verify PPP2015 Eq. 37 (Gudmundsson fit T_b = 1.288e8 (T_s6⁴/g₁₄)^0.455 K, local) for B = 0 Fe, PPP2015 Fig. 11 and the App. B analytic T_b–T_s relation (local).<br>**L1** conductive opacity κ_cond = 16σT³/(3ρκ_e) combined with radiative opacity (PPP2015 Eq. 26, local) — 20 LoC.<br>**L2** magnetised electron conduction tensor κ_∥, κ_⊥, κ_∧ from energy-averaged τ_ij(ε) (PPP2015 Eqs. 31–36, local) — 120 LoC — verify the B → 0 limit and PPP2015 Fig. 9 (local).<br>**L3** electron collision frequencies: e-i for non-degenerate and degenerate electrons (PPP2015 Eqs. A.1, A.3, local), e-e (A.2, local); magnetic Coulomb logarithms Λ_∥, Λ_⊥ — missing: Potekhin 1999, A&A 351, 787, arXiv:astro-ph/9909100 (confirm); e-phonon in the crystal — missing: Baiko & Yakovlev 1995, Astron. Lett. 21, 702; ion structure factor S(q) of the liquid OCP — compute from our own OCP MC (L5) — 300 LoC — verify PPP2015 Fig. 4 and 9 (local).<br>**L4** envelope EOS: magnetised partially degenerate electron gas (PC01, `potekhin_chabrier_2001_eos`, local) — 300 LoC — verify PC01 published values (local).<br>**L5** ion Coulomb free energy and S(q): PC01 uses the Potekhin & Chabrier 2000 fit (PRE 62, 8554, missing); recompute with OCP Monte Carlo (Ewald + Metropolis) — missing: Slattery, Doolen & DeWitt 1982, PRA 26, 2255 (verification table); Frenkel & Smit (book) — 200 LoC — verify Debye–Hückel at Γ → 0 (analytic) and SDD82 Table I.<br>**L6** radiative opacity at the envelope top: our A4/A16–A18 for H. For an Fe or C envelope, bound-free opacities are needed (missing: Opacity Project; deep). A light-element (H/He) envelope avoids this, but its composition is a parameter.<br>**L7** T_b, B_p, obliquity and envelope composition: free parameters, constrained by measured flux and pulse amplitude (`obs.rxj1856-pulse-amplitude`, local) as PRD R4 demands. Cooling theory is not needed.<br>Total ≈ 1200 LoC. |
| A34 Ho+2007 T∞, R∞, z | RECOMPUTE | **L0** our own maximum-likelihood fit of (T∞ or T_b, R∞, z, N_H, B, envelope parameters) to measured data: HST optical/UV photometry (`kaplan_2011`, `obs.*` entries, local), distance 123 pc (`walter_2010`, local), Chandra LETG spectrum (paper local as `burwitz_2003`; **the raw spectrum and responses are missing**, fetch from the Chandra archive using the ObsIDs in Burwitz 2003 / Ho 2007) — 300 LoC — verify that our blackbody fit reproduces Burwitz's R∞/d (`obs.rxj1856-xray-bb-r-over-d`, local), and that with Ho-like model assumptions we recover Ho+2007 Table 1 (T∞ = 4.34(2)e5 K, R∞ = 17 d₁₄₀ km, z = 0.22(2), local).<br>**L1** forward model = our atmosphere grid (A4, A16–A22) plus exact GR (A31), already owned in the chain above.<br>**L2** X-ray ISM photoabsorption σ_ISM(E) — missing: Wilms, Allen & McCray 2000, ApJ 542, 914, arXiv:astro-ph/0008425 (confirm), which is itself built on Verner et al. 1996 (ApJ 465, 487) fits to photoionisation cross-sections. Recompute route: atomic photoionisation of H, He, C, N, O, Ne, Fe from Hartree–Fock or R-matrix continuum (deep; source missing: Opacity Project / Verner 1996) — 200 LoC with abundances — verify Morrison & McCammon 1983 (missing) and Wilms 2000 tables.<br>**L3** optical extinction curve (needed for R5) — see Part C, M2.<br>**L4** instrument responses (ARF/RMF, HST throughputs): measured, irreducible. |
| A35 B = 1e12 G | RECOMPUTE | Make B a fitted parameter of A34, or the spin-down value (`obs.rxj1856-dipole-field`; B = 3.2e19 √(P Ṗ) is a convention that assumes vacuum dipole spin-down with assumed I and R, see the A24 caveat). |
| A27 ρ_nuc, A36 NSParams defaults | DROP | Not consumed by the flyby. Every number comes from RX J1856 data, the A34 fit, or the labelled ensemble. |
| A37 f_col, p | DROP | Milestone 1 uses a pure blackbody (f_col = 1 and p = 0 by definition, a convention). Milestone 3 gets the angle- and frequency-dependent intensity directly from the RT solution. A colour-correction factor is then a derived diagnostic, never an input. |
| A31 Beloborodov | RECOMPUTE (exact) | **L0** exact Schwarzschild null geodesic ODE or bending integral from a finite-r moving camera — local (`pechenick_1983`, `beloborodov_2002` Eq. 9; `gr.bending-integral`, `gr.null-geodesic-*` entries) — 200 LoC — verify the flat-space limit (visible fraction exactly ½), Pechenick Fig. 1, and Beloborodov's stated ≈1% accuracy of the cosine relation (local). |
| A11, A13, A14, A43 numerics | RECOMPUTE (method) | Not physics inputs. Replace each silent clamp or extrapolation with either (a) computing the requested ν or μ, or (b) a loud error (rule 5). A10's grid must contain the full camera band × (1+z) × the Doppler range of the flyby. |
| A30, A33 bugs | RECOMPUTE | Field–normal angle from the dipole components (tanΘ = ½ tanθ_B in flat space; use GR dipole components when A28 is upgraded — source missing: Wasserman & Shapiro 1983, ApJ 265, 1036). Pixel scale from the camera model. |
| A39 sRGB matrix | RECOMPUTE (convention) | Invert the primaries/D65 matrix built from the chromaticities in `colour.srgb-*-primary` and `colour.srgb-d65-white-*` (local). Verify against Stokes 1996 to its printed 4 decimals (local). |

Measured and convention items (A1–A3, A38, A40–A42) are irreducible. The
derivable items A5–A9, A15, A17, A19, A28 and A32 need only re-implementation
from their local equations; they carry no ownership debt.

## Part C — ownership debt per milestone

| milestone | blocking fit/table/code items (rank = order to own) | non-v1 debt that appears here |
|-----------|-----------------------------------------------------|-------------------------------|
| **M1** blackbody + exact GR + real sky | 1. **A34** (T∞, R∞, z from Ho's code). A blackbody disc still needs a T and an R. Unblock by declaring them ensemble inputs labelled `code` (or by using the measured d and the X-ray blackbody R∞/d as the stated input) until A34 is recomputed. No opacity, EOS or atmosphere fit is consumed. A31 (exact geodesics) is derivable and already sourced locally. | Gaia DR3 photometry → spectral radiance per star needs a stellar SED model per (BP−RP, G). If taken from a published model grid (e.g. synthetic spectra), that is a **new table dependency**. Owning it means a blackbody or own-atmosphere fit to Gaia XP spectra, which are measured. The diffuse Milky Way map is measured (licence open, PRD). |
| **M2** Earth-distance photometry (R5) | 1. **A34** again: the model's T, R, z must come from our own fit, otherwise R5 is circular, because Ho's numbers were themselves fitted to these data. | (a) HST filter throughputs: measured, missing (STScI synphot files). (b) STMAG zero point: convention, local (`obs.stmag-definition`). (c) **Interstellar extinction**: Kaplan 2011 uses A_V = N_H/1.79e21 (Predehl & Schmitt 1995, A&A 293, 889, a fit to X-ray halo data) and the Cardelli, Clayton & Mathis 1989 (ApJ 345, 245) curve, a fit to measured extinction curves (`obs.approx-kaplan2011-extinction`, local). Own route: fit A_V to our own model with the extinction curve as a measured-curve dataset, or recompute grain extinction by Mie theory with a dust model (missing: Draine 2003 ARA&A; deep). A_V ≲ 0.1 mag (local), so it is small but must be owned or stated. |
| **M3** surface models | In order: 1. **A4** Gaunt factors (blocks the B = 0 H atmosphere, the first model; 300 LoC; local fixture exists). 2. **A12** seed (5 LoC). 3. **A20** vacuum coefficients (60 LoC, all local). 4. **A18** Λ_pp (60 LoC, local). 5. **A16 corrections** (non-Born plus finite proton mass; 400 LoC; Potekhin 1994 missing). 6. **A21/A22** vacuum resonance and mode conversion (120–270 LoC; Lai & Ho 2002/2003 missing). Items 3–6 matter only for the magnetic H model. 7. **A29** surface temperature map via the envelope chain (≈ 1200 LoC; Potekhin 1999, Baiko & Yakovlev 1995, Slattery+1982 missing), unless the map is parametrised and fitted to the pulse profile. 8. **A34** final fit with the owned models (R5). | The condensed-surface model (PRD scope c) will bring **new fits**: `van_adelsberg_2005` and `potekhin_2012` (both local) give fitting formulae for condensed-surface reflectivity and emissivity. Own route: the Fresnel-type reflectivity of the magnetised condensed medium from its dielectric tensor (vAL05 derivation, local) plus the condensed-matter EOS for surface density (missing: Medin & Lai 2006/2007). Audit those papers before M3 starts. |

M4 and M5 add no new physics inputs beyond M3. The display model (A41, R6)
is convention.

## Summary

Of the 43 numerical inputs v1 consumed or carried on production paths, 29 are
already acceptable under rule 7: 2 measured, 10 conventions and 17 derivable
from local equations. A few of these carry defects, including a pre-2018
electron charge, a mislabelled M_sun, an sRGB matrix that matches no local
source, a misused dipole angle, a wrong pixel scale and an "exact geodesic"
claim that is false. Seven are fits, two are tables, two are other people's
code output, and three are unsourced placeholders.

Only five of the fit/table/code items actually ran in the production RX J1856
renders: the Sutherland/McPHAC Gaunt table, the McPHAC T(0) seed, the PC03
Λ_pp fit, the PLCH04 vacuum-polarisation fits, and Ho+2007's T∞/R∞/z.

The TOV solver, the BSk fit and the Potekhin tables were never on a render
path. They can be dropped, because the flyby needs only the exterior metric
and surface radiance. That stays true as long as B is fitted rather than
derived from spin-down with a computed moment of inertia.

Every remaining fit has a recompute chain whose top layers are stated in local
references:
- PC03 Appendices A and B give Λ_pp and the finite-proton-mass Λ^e_α.
- PLCH04 Eqs. A2–A6 give the exact Heyl–Hernquist vacuum coefficients.
- vAL06 gives the resonance and Landau–Zener structure.
- PPP2015 and PC01 give envelope transport and EOS.

The missing sources are concentrated in a short fetch list: van Hoof 2014 and
Karzas & Latter 1961 for Gaunt factors; Potekhin 1994; Lai & Ho 2002 and Ho &
Lai 2003; Heyl & Hernquist 1997; Potekhin 1999; Baiko & Yakovlev 1995;
Slattery, Doolen & DeWitt 1982; Wilms+2000 and Verner+1996; the Chandra LETG
data; and the HST throughputs.

Milestone 1 carries essentially no physics-fit debt. Its one issue is that the
RX J1856 T, R and z are Ho's code output at d = 140 pc, and they must be
labelled as ensemble inputs.

Milestone 2's debt is that same parameter set, which becomes circular under
R5, plus the interstellar extinction law.

Milestone 3 carries the bulk of the debt. In order it is: Gaunt factors (about
300 LoC, the first blocker), the vacuum-polarisation and mode-conversion fits
(about 400 LoC), the magnetic free-free corrections (about 500 LoC), the
surface temperature map through the envelope chain (about 1200 LoC including
the OCP Monte Carlo), and the condensed-surface fits still to be audited.
