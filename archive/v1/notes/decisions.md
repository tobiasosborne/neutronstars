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

## D11: Suleimanov Fig. 2 verification harness must use vacuum-polarized weights
**Date:** 2026-05-16
**Context:** `verification/plot_suleimanov_fig2_clean_validation.py`
reported ~2.7-3.0 dex residuals on the upper free-free branches at
`θ_B=45°, B=10¹⁴ G` while `θ_B=5°` and `θ_B=90°` were fine.
Reviewer D9 inspected the vacuum-polarization / mode-vector physics
in `src/opacity/dielectric_tensor.jl` + `src/opacity/magnetic_modes.jl`
and concluded "physics matches the papers". The HANDOFF
notes therefore speculated that the residual was caused by the ad-hoc
grey flux correction (D10).

That diagnosis was wrong. The flux correction lives inside the
magnetic atmosphere solver and only affects the converged T profile;
it cannot influence opacity at fixed `(T, ρ, B, θ_B)`. The actual bug
was in the verification harness itself:
`verification/compute_suleimanov_fig2_opacities.jl` imported
`polarization_weights_full`, which since commit 4dc9a5d (D6) is a
deprecated compatibility wrapper around `polarization_weights_cold`
(no vacuum polarization). The script was therefore generating
cold-plasma opacities and comparing them against SPW09 Fig. 2, which
explicitly uses vacuum-polarized normal modes:

> "Opacities are calculated in the same way as in the paper by
> van Adelsberg & Lai (2006) ... The vacuum polarization effect is
> taken into account following the same work."  (SPW09 §2, p. 3,
> lines following Eq. 10; also Fig. 2 caption: "The proton cyclotron
> and vacuum resonances are also shown.")

`src/opacity/magnetic_modes.jl` already uses
`polarization_weights_vacuum`, so the production opacity pipeline was
correct — only the verification harness was wrong.

**Decision:** Replace `polarization_weights_full` with
`polarization_weights_vacuum` in
`verification/compute_suleimanov_fig2_opacities.jl`. With the fix in
place, the θ_B=45° best-fit residuals drop from 2.66-5.58 dex to
0.024-0.089 dex; every digitized branch at θ_B ∈ {5°, 45°, 90°} now
matches its best computed mode to <0.4 dex median absolute error.
Add θ_B=45° to `verification/check_suleimanov_fig2_metrics.py`
(previously only θ_B=5° and 90° were gated).

**Rationale:** The verification harness must use the same mode-vector
convention as the production code and as the reference paper. D9's
"physics matches" was a correct statement about the source code; what
was missing was a runtime check that the *verification script* was
exercising that code path.

**Scope:** `verification/compute_suleimanov_fig2_opacities.jl` (one-line
import + one-line call site), `verification/check_suleimanov_fig2_metrics.py`
(two new gated rows for θ_B=45°). No production-code change. Test suite
still passes (119/119).

## D12: Depth-resolved Avrett-Krook flux correction — blocked on missing sources
**Date:** 2026-05-16
**Context:** D10 deferred the depth-resolved Avrett-Krook flux correction
to a dedicated session. This session is that attempt. The plan was to
replace the uniform grey scaling

```
flux_scale = clamp(1 + flux_damping*(flux_ratio^{-1/4} - 1),
                   FLUX_SCALE_CLAMP_LO, FLUX_SCALE_CLAMP_HI)
ΔT_total[i] = ΔT_lambda[i] + (flux_scale - 1)*T[i]    # global, depth-flat
```

in `_iterate_magnetic_rybicki` (now inlined in `solve_magnetic_atmosphere`)
with a **per-depth** flux-deficit correction `ΔT_flux[i]` driven by the
local relative flux error `ε_H(τ_i) = (H_actual(τ_i) - H_target)/H_target`,
applied additively on top of the existing Rybicki Λ-correction. SPW09 §2
names this scheme (Avrett-Krook) and credits it to **Kurucz (1970), SAO
Special Report 309** for the precise formula. SPW09 prints only the
signature (their unnumbered formula around p. 5 line 881-884):

> "The second procedure is the Avrett-Krook flux correction, which uses
> the relative flux error ε_H(m), and it is performed in the deep layers."

and cites Kurucz 1970 for "a detailed description".

**Decision:** STOP. Revert any code change and document. The Avrett-Krook
formula is **not available from any source in `refs/`**, and no usable
proxy exists either:

1. **Kurucz 1970, SAO Special Report 309** — not in `refs/`. The classical
   reference for the Avrett-Krook correction as adapted to stellar
   atmospheres. SPW09 only cites it; the formula does not appear in
   SPW09 itself.
2. **Mihalas 1978 "Stellar Atmospheres" Ch 7 (eq. 7-46)** — not in
   `refs/`. The textbook formulation. The task brief claimed Haakonsen
   2012 Appendix A reproduces it; **this is incorrect** — Haakonsen
   Appendix A is the pure Rybicki Λ-correction (eqs A1-A33), which is
   *already implemented* in `src/atmosphere/temp_correction.jl` and
   `_rybicki_two_mode` in `magnetic_atmosphere.jl`. No flux-derivative
   pass appears anywhere in Haakonsen.
3. **Avrett 1971 / Avrett-Krook 1963 ApJ 137 874** — not in `refs/`.
4. **McPHAC source (`refs/code/McPHAC/DeltaT.c`, `CalcJtbar.c`)** —
   inspected. McPHAC uses ONLY the Rybicki Λ-correction
   `ΔT = ∫κ(J-B)dν / ∫κ(dB/dT)dν` and an optional fractional damping
   (`DTFRAC`). It has NO Avrett-Krook flux step. (Haakonsen 2012 §3.4
   confirms: "The temperature correction scheme used is the same as
   [Zavlin 1996]" — Λ-only.) So McPHAC is not a valid cross-code
   reference for the AK formula either.

Implementing a "depth-resolved flux correction" without one of these
sources would either (a) be reverse-engineered from memory (violates
Sacred Principle #1), or (b) be a heuristic of comparable epistemic
status to the current grey scaling (no net improvement, more code).

**Baseline verification (before any code change):**
```
solve_magnetic_atmosphere(1e6, 2e14, 1e12, π/4, gaunt;
                          nν=24, M=4, N=50, max_iter=20, tol=5e-4,
                          flux_tol=1e-2, flux_damping=0.5, verbose=true)
→ CONVERGED at iteration 17, Final F/σT⁴ = 0.9998
```
The current grey-scaling solver hits 4σ inside the 5% test tolerance
within the orchestrator's 20-iteration budget. The case for replacing
it on physics grounds exists (deep-layer accuracy at high B), but the
case for replacing it without the source formula does not.

**Next-session prerequisites (any one suffices to proceed):**

1. **Acquire Kurucz 1970 SAO Special Report 309.** Out of print
   (1970 SAO bulletin); not on ADS as PDF; sometimes available via
   inter-library loan from SAO/CfA. The Avrett-Krook section is ~10
   pages.
2. **Acquire Mihalas 1978 "Stellar Atmospheres" 2nd ed., Ch 7.4-7.5.**
   The "ΔT_flux" formulation is eq. 7-46 (per the orchestrator brief;
   not independently verified). ISBN 0-7167-0359-9. ~5 pages.
3. **Acquire a stellar atmosphere code that implements AK and is
   open-source.** Kurucz's ATLAS9/12 (Fortran, public) is the canonical
   implementation. Reading `atmos.for` or `synthe.for` to extract the
   AK update is feasible. ATLAS source: kurucz.harvard.edu.
4. **Decide the project will rely on the Rybicki Λ-correction alone**
   (matching McPHAC), demote SPW09's three-pass scheme to "future
   work", and just remove the grey scaling. This is a smaller change
   but needs its own verification (convergence at θ_B=45°, B=1e14 G,
   without the grey scaling crutch).

**Scope of this attempt:** No production-code change. This entry only.
Total session: ~30 min reading (SPW09 §2, Haakonsen Appendix A,
McPHAC DeltaT.c), ~5 min baseline smoke, 0 min coding. The orchestrator's
STOP criterion "physics judgment call you can't justify from the local
PDFs" was triggered on the source-availability check before any
implementation began.

**Estimated work for a future session, assuming prerequisite #2 (Mihalas
Ch 7.4-7.5):** ~150-200 LoC for a per-depth ε_H(τ) integrator + AK
update + Kurucz surface correction, plus a θ_B=45°, B=1e14 G
verification (which is what the current grey scaling has been suspected
of breaking but never proven to — see D11 for why the original
verification was misdiagnosed).
