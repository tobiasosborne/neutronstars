# Review Synthesis — 2026-05-16

Four parallel critical reviews ran on commit `318de1a` immediately after the
master/origin merge. Source files:

| # | Reviewer | File | Lines |
|---|----------|------|-------|
| 1 | Architecture | `reviews/01_architecture.md` | 213 |
| 2 | Test coverage | `reviews/02_tests.md` | 268 |
| 3 | Code smells + derivations | `reviews/03_code.md` | 331 |
| 4 | Repo tidiness | `reviews/04_tidiness.md` | 600 |

This document consolidates findings that appear in multiple reviews, ranks
them, and proposes an order of attack. It does not replace the individual
reviews — the per-finding detail (file:line, code quotes, paper references,
suggested fixes) is in the originals.

---

## Headline assessment

**The codebase is in better shape than the user fears, but the gap between
HANDOFF.md claims and what is actually verified is larger than HANDOFF.md
implies.** The physics modules are well-decomposed, the recently-merged
remote work (vacuum polarization, Eq. 44e fix, Suleimanov Fig 2 harness) is
genuinely good, and several pieces of discipline are unusually high
quality (per-mode `B_ν/2` documented everywhere, mode-separated angle
averaging, captured intermediate regression values).

The mess concentrates in three places:

1. **The pipeline layer (`AtmosphereGrid`, `Renderer`) does not respect the
   physics layer.** It silently drops `θ_B`, collapses polarization modes,
   uses nearest-neighbour interpolation, prints a numerically-wrong flux
   diagnostic, and has zero tests. Phase 4 (Kerr ray tracing) cannot land
   cleanly on this foundation.
2. **Test coverage is bimodal.** Either pinned to last-digit precision
   (change-detectors, not physics checks) or accepts 25× looser tolerance
   than the solver itself uses. The headline physics claims in HANDOFF.md
   ("matches McPHAC within 1.2%", "B=0 limit recovers non-magnetic to
   <0.1%", "F/σT⁴ = 0.99", TOV M_max, visible_fraction values) have **zero
   automated regression backing**.
3. **The working tree is a mess** — two stale duplicate CLAUDE.md files,
   `output/` half-tracked / half-untracked (120 MB of frames untracked,
   2.3 MB of redundant PPMs committed), `refs/code/McPHAC` is a 109 MB
   nested git repo masquerading as vendored source (dirty,
   not-a-submodule, has compiled `.o` files), `Manifest.toml` committed
   for a library, orphaned `proof/` directory from a 2026-03-23 one-shot
   tool run.

---

## Convergent findings (multiple reviewers flagged the same thing)

These are the issues where the reviewers, working independently, converged.

### Convergent #1 — The flux correction is mislabeled and over-clamped
- Code review **B3**: "ΔT *= 1 + flux_damping × (flux_ratio^(-1/4) − 1) is NOT from SPW09." SPW09 uses Avrett-Krook (per-depth ε_H(m)), not grey scaling.
- Code review **A13**: the `[0.9, 1.1]` clamp hides divergence — a runaway iteration would silently cap at 10% correction even when 40% is needed.
- Architecture **#21**: `flux_damping=0.5` is one of ~22 scattered magic numbers with no central knob.
- Architecture **#15**: the `_flux_ratio` diagnostic printed during grid build is itself numerically wrong (missing Gauss-Legendre weight).
- Test review **#9**, **#10**, **#11**: regression hazards — dropping the flux correction entirely, or setting `flux_damping = 0`, both silently pass the current flux test (`rtol=0.25`).

**Severity: HIGH.** This is likely connected to the unresolved `θ_B = 45°` Suleimanov Fig 2 mismatch flagged in HANDOFF.md as the next task.

### Convergent #2 — `update_atmosphere!` is exported with a wrong EOS factor
- Architecture **#2**: drops factor of 2 from μ=0.5 ionised H. Function is exported, documented, but dead — iteration uses private `_update_structure!` instead.
- Code review **A9**: also flags that the EOS formula `ρ = m_p P / (2 k_B T)` is inlined five times across the codebase with no `density_from_PT` helper.

**Severity: HIGH (latent).** Any user following the public API gets 2× density. Single fix: extract `density_from_PT(P, T)` to `AtmosphereStructure`, replace all five copies, delete the bad public function.

### Convergent #3 — The renderer silently uses pole-on physics
- Architecture **#1** (critical): `AtmosphereGrid` is built with `θ_B = 0.0`; `lookup_spectrum` doesn't take `θ_B`; the renderer computes per-pixel `θ_B` and discards it. **Every magnetic render in `output/` is pole-on opacities applied to a dipole T-map.**
- Code review acknowledges in closing notes (intentionally not duplicating).
- Test review **#16**: zero tests of the render pipeline; the bug is invisible to the test suite.

**Severity: CRITICAL.** The pulsed-fraction shape in the README hero animations is decoupled from the actual radiative transfer. Fix is structural — `θ_B` needs to become a grid axis.

### Convergent #4 — Non-magnetic Rybicki is O(N·K²) by accident
- Code review **A5**: `temp_correction.jl:73-85` precomputes `B_bar[i]` per depth (~O(N·K)), then `:197` *recomputes the same sum inside the inner k-loop*. The precomputed array is unused. The magnetic version (`magnetic_atmosphere.jl:574`) correctly uses the precomputed array.
- Test review: non-magnetic Rybicki has zero test coverage, so the inefficiency (and any future regression to a correctness bug there) is silent.

**Severity: HIGH (performance only; not a correctness bug today).** A 4–8× slowdown of the non-magnetic Rybicki step that would never be discovered through running the code, only through reading it. Two-line fix.

### Convergent #5 — Magnetic atmosphere does NOT recover non-magnetic at B=0
- Code review **B4**: the magnetic Feautrier uses an isotropic scattering kernel, the non-magnetic uses Hummer (1962) dipole. At B→0 these don't agree. The HANDOFF claim "B=0 recovers non-magnetic to <0.1%" must be true only because the scattering contribution is small at the tested conditions, not because the kernels are consistent.
- Test review **#1** (priority A): the B=0 consistency check is not in `runtests.jl`. The two solvers could silently diverge with no warning.
- Architecture: praises the single-codepath-with-branch design but doesn't verify it.

**Severity: MEDIUM.** Real physical inconsistency. The fix is to use the same kernel in both paths weighted by `ρ_alb`.

### Convergent #6 — Pipeline layer cannot scale to Phase 4
- Architecture **#1, #3, #4, #7, #9, #11**: AtmosphereGrid has no provenance metadata, no `θ_B`, no polarization axis, no concurrency, inconsistent array shapes, fragile interpolation.
- Code review **A3, A12, A18**: nearest-neighbour ν lookup with O(K)-per-call allocation, 6.5 GB garbage per 512² render.
- Test review **#16**: zero tests of the entire render pipeline.

**Severity: HIGH (long-term).** When Kerr ray tracing arrives (HANDOFF Phase 4), polarization must flow along null geodesics — there must be a polarization axis. When HDF5 caching arrives (HANDOFF backlog #3), there must be provenance. Most of these need to be fixed before either of those features lands.

### Convergent #7 — No CI, no gated verification
- Test review (entire section "Verification harness assessment"): the Suleimanov Fig 2 harness is the only thing with an exit code, but it's not run by `runtests.jl`. `eos_checks.jl` prints "PASS"/"FAIL" with no exit code — a `FAIL` line will silently scroll past.
- Architecture **#25**: tests depend on `refs/code/McPHAC/gffgu.dat` — a file inside a dirty nested repo — so clean-CI is blocked.
- Tidiness: the McPHAC dependency is itself unstable (nested repo, no `.gitmodules`).

**Severity: HIGH.** No mechanism prevents a future agent from silently breaking the 1.2% McPHAC agreement.

---

## Independent high-severity findings (one reviewer only, but real)

These didn't appear in multiple reviews, but each is concrete and load-bearing.

| ID | Reviewer | Severity | Summary |
|----|----------|----------|---------|
| Code **A1** | code | HIGH | `polarization_weights_*` returns heap-allocated `Vector{Float64}` of length 3 per call — ~600k allocations per solve. Fix: `NTuple{3,Float64}` or `SVector{3,Float64}`. |
| Code **A2** | code | HIGH | `polarization_weights_vacuum` called twice per opacity sample because `mode_absorption` and `mode_scattering` are separate entry points. ~4× speedup available via combined `mode_opacity_split`. |
| Code **A4** | code | HIGH | `_feautrier_single` allocates `Vector{Matrix{Float64}}` per frequency per mode + redundant `copy(B)` / `copy(Q)`. Preallocate once per `_solve_feautrier_mode!` call. |
| Code **A6, A7** | code | MEDIUM | Bessel-function fallbacks for tiny `x_n` are numerically wrong: `K0_val=10.0` (should be ≈23 at `x=1e-10`); `K1_val=1.0` makes `x·K_1(x)→0` instead of →1. Only triggers near cyclotron resonance — but resonance is exactly where accuracy matters. |
| Code **A23** | code | MEDIUM | Two implementations of `coulomb_log_classical` disagree at high `u`: `magnetic_ff.jl:39` returns 1.0; `coulomb_magnetic.jl:118` returns the correct `√(π/u)` asymptote. Delete the bad version. |
| Code **B11** | code | MEDIUM | `sigma_scat_alpha` non-magnetic fallback returns `σ_T/3` per polarisation. Under current `Σ|e_α|²=1` weighting, this is 3× too low. Dead code path today, but wrong if exercised. |
| Code **B13** | code | MEDIUM | `surface_temperature` uses `T_eq + (T_pole − T_eq) cos²θ` but docstring claims Greenstein-Hartke 1983, whose actual prediction is `T ∝ |cos θ|^{1/4}` — a much shallower distribution. Either rename the function or implement the actual GH formula. |
| Arch **#5, #6** | architecture | MEDIUM | Three nearly-identical interpolation routines: `rt_emergent_spectrum` (good, log-log), `lookup_spectrum` (nearest-neighbour), `lookup_spectrum_scalar` (dead duplicate), plus a fourth inlined in `scripts/render_rxj1856_visible_magnetic_atmosphere.jl`. |
| Arch **#17** | architecture | LOW | `src/eos/potekhin_table_reader.jl` is a 235-line orphaned module — never included. Delete or wire into `verification/`. |
| Test **Tier 1** | tests | HIGH | Six one-line tests would catch the most embarrassing potential regressions: TOV 1.4 M☉ → R≈12.6 km; visible_fraction(u=1/3)=0.75; solar BB → sRGB (255,252,245); dipole pole/equator limits. Total effort: ~15 minutes. |

---

## Needs paper verification (flagged honestly, can't be resolved here)

| ID | Reviewer | What | Why unverified |
|----|----------|------|----------------|
| Code **B2** | code | `compute_weights_from_K_Kz` mode-vector sign convention | `refs/ho_lai_2001_vacuum_pol.pdf` is mislabeled per HANDOFF — actually an unrelated quasar paper. Need to re-download van Adelsberg & Lai (2006) from arXiv. **This is at the heart of the unresolved θ_B=45° Fig 2 mismatch.** |
| Code **B7** | code | `vacuum_coefficients` constants `0.25487, 0.75, 1.2, 1.33, 0.56, 3.75, 2.7` | `refs/potekhin_chabrier_ho_2014_opacities.pdf` is mislabeled per HANDOFF — actually Boev & Kovalev 2014 (exciton BEC). Need to re-download PLC2004 from arXiv. |
| Code **B12** | code | Stix `D` has opposite sign to textbook convention | May be intentional Ginzburg/P&C convention; need to confirm against P&C 2003 Eqs. 22-25. |

Reviewer 3 explicitly refused to guess on these — flagging is the correct posture.

---

## Proposed order of attack

### This week — 1 day total

The fewest changes that have the highest signal:

1. **Tidiness quick wins (1 hour)** — drop-in `.gitignore` is in
   `reviews/04_tidiness.md`. `git rm --cached Manifest.toml`. Delete
   `refs/potekhin_tables/hmagnet.tar.gz` (duplicate) and
   `refs/potekhin_tables/eos_potekhin.f` (0 bytes). Recovers
   ~120 MB of working-tree noise.
2. **Delete `NEUTRON_STAR_PIPELINE.md`, rename
   `NEUTRON_STAR_CLAUDE_MD_v2.md → CLAUDE.md` (10 min).** Cherry-pick the
   Phase-2 tracer-bullet recipe from PIPELINE before deleting.
3. **Tier 1 tests (15 min)** — six one-liners listed in
   `reviews/02_tests.md`. Catches TOV / dipole / colorimetry regressions
   that are currently invisible.
4. **Fix Convergent #4 (5 min)** — `temp_correction.jl:197`: replace the
   inline B̄ sum with the precomputed `B_bar[i]`. 4–8× speedup of the
   non-magnetic Rybicki step.
5. **Tighten the magnetic flux test tolerance (1 min)** —
   `runtests.jl:116`: change `rtol=0.25` to `rtol=0.05`. If it fails
   today, that's itself a finding worth knowing.
6. **Fix the `_flux_ratio` diagnostic (10 min)** —
   `atmosphere_grid.jl:168-177`: missing Gauss-Legendre weight. Either
   delete the function and call the canonical `_bolometric_flux` from
   `RTAtmosphere`, or fix the formula. Currently lying to the user
   during grid build.

### This sprint — 1–2 weeks

The structural fixes that unblock Phase 4:

7. **Fix Convergent #2 — `update_atmosphere!` EOS bug.** Extract
   `density_from_PT(P, T)` to `AtmosphereStructure`. Replace all 5 inline
   copies. Delete or fix the broken exported function.
8. **Fix Convergent #5 — B=0 limit consistency.** Use the same scattering
   kernel in both `solve_magnetic_atmosphere` (B=0 path) and
   `solve_atmosphere`. Add the B=0 consistency test (Tier 2 in
   `reviews/02_tests.md`).
9. **Wire `verification/check_suleimanov_fig2_metrics.py` into
   `runtests.jl`** as a `@testset` that calls the script and asserts
   exit code 0. ~30 lines. Converts the Suleimanov gate from manual to
   automated.
10. **Add CI** — `.github/workflows/ci.yml` that runs `Pkg.test()` plus
    `python3 verification/check_suleimanov_fig2_metrics.py`. Add
    `[extras]` + `[targets] test = ["Test"]` to `Project.toml`.
11. **Add provenance metadata to `AtmosphereSpectrumGrid`** (Architecture
    #3) — K, M, N, tolerances, gaunt-table SHA, code SHA, θ_B,
    convergence flags per cell. Required before HDF5 caching ships.
12. **Convert `refs/code/McPHAC` to a real submodule** (Tidiness R1).
    ~30 min. Fixes the persistent `git status` dirt and the CI-blocker
    in Architecture #25.

### Before Phase 4 — 2–4 weeks

The pipeline-layer rewrite:

13. **Fix Convergent #3 — renderer θ_B drop.** Add `θ_B` as a grid axis
    (initial set `{0, π/4, π/2}`). Route through `lookup_spectrum` and
    `render_spectral_cube`. Until that lands, rename the current artifact
    to `AtmosphereSpectrumGridPoleOn` so future agents can't be misled.
14. **Fix Convergent #6 — pipeline architectural rewrite.** Introduce a
    `PhysicalSurfaceState → ObservedSpectrum` interface parameterised
    over (T_eff, B, θ_B, polarization). Carry polarization axis through
    grid + renderer; replace `_nearest` with `searchsortedlast` + log-log
    interp; preallocate Feautrier scratch space.
15. **Fix Convergent #1 — implement the actual Avrett-Krook flux
    correction** (or honestly rename the heuristic). Re-test against
    Suleimanov Fig 2 to see if the θ_B=45° mismatch closes.
16. **Re-download mislabeled PDFs from arXiv** (van Adelsberg & Lai 2006;
    Potekhin, Lai & Chabrier 2004) — required to verify Code review B2,
    B7, B12 and to close the θ_B=45° investigation.
17. **Hot-loop allocation fixes** (Code A1, A2, A3, A4). Expected
    ~5–10× speedup of magnetic atmosphere solves and renders.

### Optional / nice-to-have

- Convert `verbose::Bool` to `@info`/`@debug` (Code A15 / Arch #20).
- Centralise magic numbers into a `SolverDefaults` module (Arch #21).
- Structured error types (Arch #18).
- `Float64` → `T<:Real` for AD compatibility (Code A19).
- `proof/` directory cleanup (Tidiness R4).
- HANDOFF.md per-session pattern (Tidiness R2).

---

## What's done well (from all four reviews, deduplicated)

- Module decomposition is clean and matches the physics pipeline. New
  readers can find any piece quickly.
- Physical constants are isolated, sourced (CODATA), re-exported once.
- Per-mode `B_ν/2` discipline is consistently documented at every use
  site and reflected in code.
- The recently-merged remote work (Eq. 44e fix, vacuum polarization
  mode vectors, Suleimanov Fig 2 digitisation, magnetic-Eddington init,
  thermalisation-depth guard) is genuinely good work — clean, sourced,
  regression-tested where possible.
- `vacuum_mode_parameters` is properly decomposed into testable pieces.
- `MagneticAtmosphereResult` captures full provenance (T_eff, g_s, B,
  θ_B, iter count, max ΔT/T, grids, profiles, intensities).
- `notes/{decisions,approximations,failures}.md` is a healthy
  append-only pattern.
- `verification/VERIFICATION_LOG.md` and `verification/check_suleimanov_fig2_metrics.py`
  are the right shapes for what they do — they just need to be gated.
- `src/` directory tree is clean, snake_case throughout, matches the
  declared module structure.
- LICENSE, Project.toml, and README are all in order.
- The Suleimanov Fig 2 digitisation pipeline (Python + Julia round-trip
  with metric thresholds) is the project's best-engineered verification.
- The adaptive Feautrier solver is packaged (not yet wired in) — real
  "found a better method" work rather than a placeholder.

---

## Closing note

This is a research code that has been built in confident, rapid sessions
with clear physical motivation and a healthy notes/decisions discipline.
Most of the issues found above are not failures of effort or design taste
— they are predictable consequences of the codebase passing through a
stage where verification has not kept up with code growth. The
recently-merged Suleimanov Fig 2 harness shows the project knows how to
verify properly; the gap is that automated regression coverage hasn't
been built around any of the older Phase 2/3a/3c work.

The single highest-leverage change is probably **adding CI + Tier 1
tests + wiring the Suleimanov gate into `runtests.jl`** — about 100
lines of new test/CI code that would convert most of HANDOFF.md's
English-prose claims into enforced contracts.
