# Design Decisions

## D1: Start with Phase 2 tracer bullet using blackbody placeholder
**Date:** 2026-03-24
**Decision:** Implement end-to-end pipeline with modified blackbody atmosphere before self-consistent RT.
**Rationale:** Working pipeline with known approximations beats perfect opacity module with no pipeline. The blackbody placeholder (f_col correction from Suleimanov+ 2009) is a standard observational approximation.

## D2: CGS units throughout physics code
**Date:** 2026-03-24
**Decision:** All physics functions use CGS (erg, cm, g, s, Gauss).
**Rationale:** Standard in astrophysics literature. All reference papers use CGS. Avoids unit conversion errors when comparing against published values.

## D3: Schwarzschild first, Kerr deferred to Phase 4
**Date:** 2026-03-24
**Decision:** Ray tracing starts with Schwarzschild metric (no rotation).
**Rationale:** Most NS spin periods >> R/c ≈ 40 μs, so rotation effects are small. Schwarzschild has exact elliptic integral solutions (Pechenick+ 1983). Kerr adds significant complexity (Dexter & Agol 2009).

## D4: Fully ionised H atmosphere only in Phase 2-3
**Date:** 2026-03-24
**Decision:** Defer partial ionisation (Potekhin & Chabrier 2004) to Phase 4.
**Rationale:** Partial ionisation matters for T_eff < 5×10⁵ K and B > 10¹² G. The canonical tracer bullet target (T_pole = 10⁶ K) is fully ionised.

## D5: Julia native — no runtime dependency on external Fortran/C
**Date:** 2026-03-24
**Decision:** All physics implemented in Julia. External codes (McPHAC, Potekhin Fortran) used only for verification.
**Rationale:** Reproducibility and transparency. Every equation must be traceable in our code.

## D6: Preserve cold-plasma mode API and add explicit vacuum mode weights
**Date:** 2026-05-14
**Decision:** Keep `polarization_weights_full` as the cold-plasma compatibility wrapper and add `polarization_weights_vacuum` for QED vacuum-polarized local mode vectors. Magnetic opacities use the vacuum-aware weights.
**Rationale:** Existing callers and tests that mean "P&C 2003 cold plasma" keep their behavior, while production magnetic opacities use the physically updated local mode vectors. This avoids silently changing every low-level dielectric caller's convention.

## D7: Use figure digitization as validation, not as the sole ground truth
**Date:** 2026-05-14
**Decision:** Use digitized Suleimanov Fig. 2 curves as regression and smoke-validation artifacts, with published equations and reference implementations remaining the primary ground truth.
**Rationale:** Fig. 2 digitization is useful for catching dex-scale errors, but branch labels, annotation pixels, resonance markers, and anti-aliased line strokes make it unsuitable for percent-level validation without stronger curve extraction.

## D8: Surface temperature — phenomenological cos² ansatz, not Greenstein-Hartke
**Date:** 2026-05-16
**Context:** Code review (reviews/03_code.md finding B13) flagged that
`surface_temperature` in `src/surface/dipole.jl` uses
`T(θ_B) = T_eq + (T_pole − T_eq) cos²(θ_B)` but cited Greenstein &
Hartke (1983) Eq. 1. GH 1983 actually derive `T(θ) ∝ |cos θ|^{1/4}`
for a pure dipole with anisotropic thermal conductivity, and their
Eq. 1 is the dipole *field magnitude*, not the temperature.

**Decision:** Keep the cos² form for now and fix the labelling
(docstring + module header) to stop mis-attributing it. Document the
true GH form as a deferred physics option in `notes/approximations.md`
(see A11).

**Rationale:** The cos² ansatz has two free parameters (T_pole, T_eq)
and smoothly interpolates between them, which is what every existing
render and the Tier-1 tests in `test/runtests.jl` assume. The true GH
|cos θ|^{1/4} profile has only one parameter and goes to zero at the
equator; adopting it would (a) break the existing API and tests, and
(b) be physically wrong for a NS atmosphere with finite electron
thermal conductivity, where T does not vanish at the magnetic equator.
A future deliberate physics decision (with self-consistent thermal
transport, e.g. Potekhin, Pons & Page 2015) should supersede both.

**Scope:** Documentation/labelling only — function body and tests
unchanged.

## D9: Magnetic atmosphere surface T initial guess
**Date:** 2026-05-16
**Context:** Code review (reviews/03_code.md B1) noted that
`_magnetic_eddington_temperature` in
`src/atmosphere/magnetic_atmosphere.jl` initialises the surface
temperature with `T[1] = 0.265 * T_eff`, the empirical McPHAC
non-magnetic surface limit (H12 §3.1). For a magnetic atmosphere this
has no derivation — the magnetic Rosseland mean K_⊥ differs from k_R
so the surface limit also differs.

**Decision:** Keep 0.265 as the iteration starting value. The Rybicki
iteration converges (B=10¹², θ_B=π/4 in ~60 iter; tested smoke
passes flux F/σT⁴ within rtol=0.05). The alternative analytic grey
limit T[1] = (1/2)^{1/4} T_eff = 0.841 T_eff is unproven to converge
faster and would invalidate the existing convergence baseline.

**Rationale:** This is an initial-guess heuristic, not a physical
boundary condition. Future work (HANDOFF.md "Immediate Next Task")
on the θ_B=45° Suleimanov Fig 2 mismatch may revisit this — if the
mismatch turns out to be initial-guess-sensitive, swap to the grey
limit and re-benchmark.

**Scope:** Documentation/labelling only — no code change.

## D10: Magnetic atmosphere flux correction — ad-hoc grey scaling, not SPW09 Avrett-Krook
**Date:** 2026-05-16
**Context:** reviews/03_code.md B3 (Convergent finding #1) noted that
`_iterate_magnetic_rybicki` in `src/atmosphere/magnetic_atmosphere.jl`
applies ΔT += (flux_scale - 1) × T with
flux_scale = clamp(1 + flux_damping × (flux_ratio^{-1/4} − 1), 0.9, 1.1).
This was labelled "Suleimanov-style" in HANDOFF and commit history. SPW09
§2 actually uses depth-resolved Avrett-Krook (their Eq. 19-22) + a Kurucz
1970 surface correction — neither of which the code implements.

**Decision:** Keep the ad-hoc grey scaling for now and relabel it honestly
(comments + HANDOFF + CLAUDE). DO NOT implement actual Avrett-Krook in
this code-review pass. Real Avrett-Krook is a substantive physics
deliverable (~100-200 LoC, needs per-depth ε_H(m) integral, surface
correction tuning, validation against Suleimanov Fig 2 at θ_B=45°)
deserving its own session.

**Rationale:** The heuristic works for B=10^12 G (converges, flux passes
rtol=0.05). The unresolved θ_B=45° Fig 2 mismatch is the most likely
place where the heuristic limitation matters — and that's already the
project's immediate next task (HANDOFF.md). Doing both at once would
conflate physics improvement with code-review hygiene.

**Scope:** Comments, docstrings, HANDOFF.md, CLAUDE.md — no code
behaviour change.
