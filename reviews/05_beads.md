# Bead Registry — 2026-05-16

70 beads created from `reviews/00_synthesis.md` and the four underlying
reviews. Each bead is a task in the harness (`TaskList`); the table below
maps bead → task ID → source finding(s).

Phases run **serially** (no parallel Julia). Orchestrator is the main
Claude session; each bead is delegated to one Opus subagent at a time,
which verifies its own change.

## Phase A — Tidiness (10 beads, no Julia)

| Bead | Task # | Subject | Source |
|------|--------|---------|--------|
| A1  | 12 | drop-in `.gitignore` | Tidy quick-wins 3+4+10 |
| A2  | 13 | untrack `Manifest.toml` | Tidy quick-win 1 |
| A3  | 14 | delete redundant potekhin files | Tidy quick-win 5 |
| A4  | 15 | delete untracked frame dirs | Tidy quick-win 4 |
| A5  | 16 | rename `→ CLAUDE.md`, delete `PIPELINE` | Tidy quick-wins 7+8 |
| A6  | 17 | move `section_ns_structure.tex` to `docs/sections/` | Tidy quick-win 6 |
| A7  | 18 | add `.gitattributes` + `.editorconfig` | Tidy quick-win 10 |
| A8  | 19 | move `HANDOFF.md` → `notes/` | Tidy quick-win 9 |
| A9  | 20 | hero assets to `docs/figures/`, clean `output/` | Tidy reorg |
| A10 | 21 | McPHAC nested repo → real submodule | Tidy R1 |

## Phase B — Quick safety nets and bug fixes (15 beads, mostly Julia)

| Bead | Task # | Subject | Source |
|------|--------|---------|--------|
| B1  | 22 | Tier 1 tests (TOV/dipole/visible_fraction/colorimetry) | Test Tier 1 |
| B2  | 23 | tighten magnetic flux test `rtol=0.25 → 0.05` | Test Tier 2 #9 |
| B3  | 24 | fix Rybicki O(N·K²) recomputation | Code A5 |
| B4  | 25 | fix `_flux_ratio` diagnostic missing weight | Arch #15 |
| B5  | 26 | fix Bessel `K0_val/K1_val` fallbacks | Code A6, A7 |
| B6  | 27 | delete bad `coulomb_log_classical` | Code A23 |
| B7  | 28 | `sigma_scat_alpha` B<1e6 fallback `σ_T/3 → σ_T` | Code B11 |
| B8  | 29 | fix Greenstein-Hartke surface T claim | Code B13 |
| B9  | 30 | delete `lookup_spectrum_scalar` dead duplicate | Arch #24 |
| B10 | 31 | resolve `PotekhinTableReader` orphan | Arch #17 |
| B11 | 32 | use `PhysicalConstants.pc` literal in `render.jl` | Arch #14 partial |
| B12 | 33 | fix `polarization_weights_*` docstrings | Code A16 |
| B13 | 34 | clean Landau-level sum guard | Code A8 |
| B14 | 35 | magnetic Eddington `T(0)` hardcode `0.265` | Code B1 |
| B15 | 36 | clean `_solve_single_frequency` dead branches | Code A25 |

## Phase C — Structural fixes (7 beads, this sprint)

| Bead | Task # | Subject | Source |
|------|--------|---------|--------|
| C1  | 37 | extract `density_from_PT`, fix `update_atmosphere!` EOS bug | **Convergent #2** (Arch #2 + Code A9) |
| C2  | 38 | B=0 consistency: magnetic Feautrier dipole kernel | **Convergent #5** (Code B4) |
| C3  | 39 | wire Suleimanov gate into `runtests.jl` | review 02 |
| C4  | 40 | `[extras]/[targets]` in `Project.toml` | Tier 4 #20 |
| C5  | 41 | `.github/workflows/ci.yml` | Tier 4 #20 |
| C6  | 42 | `AtmosphereSpectrumGrid` provenance metadata | Arch #3 |
| C7  | 43 | Tier 2 tests (non-mag flux, TOV M_max sweep, Newtonian) | Test Tier 2 |

## Phase D — Pipeline rewrite (11 beads, before Phase 4)

| Bead | Task # | Subject | Source |
|------|--------|---------|--------|
| D1  | 44 | add `θ_B` grid axis (renderer drops θ_B today) | **Convergent #3** (Arch #1) |
| D2  | 45 | polarization axis preservation through grid + renderer | Arch #7 |
| D3  | 46 | `searchsortedlast` + log-log interp replacing `_nearest` | Code A3/A12/A18, Arch #5 |
| D4  | 47 | preallocate `_feautrier_single` scratch | Code A4 |
| D5  | 48 | `polarization_weights_*` return `NTuple/SVector` | Code A1 |
| D6  | 49 | `mode_opacity_split` combine abs+scat call | Code A2 |
| D7  | 50 | Avrett-Krook OR rename heuristic | **Convergent #1** (Code B3, A13) |
| D8  | 51 | re-download mislabeled PDFs (vAL2006, PLC2004) | HANDOFF caveats |
| D9  | 52 | verify B2/B7/B12 derivations against papers | Code B2, B7, B12 |
| D10 | 53 | make Gaunt table optional kwarg | Arch #13 |
| D11 | 54 | promote `magnetic_emergent_spectrum` to module | Arch #6 |

## Phase E — Architecture polish (27 beads, optional but valuable)

| Bead | Task # | Subject | Source |
|------|--------|---------|--------|
| E1  | 55 | `m_H` to `PhysicalConstants` | Code A10 |
| E2  | 56 | `flux_scale` clamp warning | Code A13 |
| E3  | 57 | move hot-loop `@assert` to outer validation | Code A14 |
| E4  | 58 | `verbose::Bool` → `@info/@debug` | Code A15, Arch #20 |
| E5  | 59 | rename `K` (freq count) → `nν` | Code A17 |
| E6  | 60 | parametric `Float64 → T<:Real` | Code A19 |
| E7  | 61 | `AtmosphereColumn` immutable | Code A21 |
| E8  | 62 | `_log_interp` clamp as kwarg | Code A24 |
| E9  | 63 | structured error types | Arch #18 |
| E10 | 64 | units handling (newtypes or Unitful) | Arch #19 |
| E11 | 65 | `SolverDefaults` module / `SolverParams` struct | Arch #21 |
| E12 | 66 | split `NSParams` → `StellarParams` + `RenderConfig` | Arch #23 |
| E13 | 67 | CMF data via Julia Artifacts | Arch #14 full |
| E14 | 68 | `NTuple` return type annotations | Code A11 |
| E15 | 69 | snapshot full `AtmosphereColumn` in result | Arch #10 |
| E16 | 70 | `Threads.@threads` grid build | Arch #11 |
| E17 | 71 | `cold_polarization_parameter` limit handling | Code B6 |
| E18 | 72 | `quadgk` error capture + warn | Code B8 |
| E19 | 73 | small-M Gauss-Legendre table or FastGaussQuadrature | Code A20 |
| E20 | 74 | centralize `n_e` definition | Code A22, Arch #22 |
| E21 | 75 | per-session HANDOFF pattern | Tidy R2 |
| E22 | 76 | reconcile zero-tables principle | Tidy R3 |
| E23 | 77 | fate of `proof/` directory | Tidy R4 |
| E24 | 78 | `scripts/render_all_hero_assets.jl` | Tidy R5 |
| E25 | 79 | factor tridiagonal solver to shared `Linalg` | Arch #8 |
| E26 | 80 | document module include order | Arch #16 |
| E27 | 81 | Tier 3 tests (sweeps + edges) | Test Tier 3 |

## Skipped (out of scope or absorbed elsewhere)

- **Test Tier 5** (performance benchmarks) — flaky on CI; iteration-count
  assertions covered indirectly by B2 (tighter flux rtol) and C7
  (TOV sweep timing implicit).
- **Arch #25** (tests depend on McPHAC checkout) — resolved structurally
  by A10 (real submodule) + C4/C5 (CI documents the dep).
- **Code A18** (argmin antipattern) — absorbed by D3 (replace _nearest).

## Orchestration rules

- **Serial execution only.** No more than one Julia process at a time.
- **One subagent per bead.** Subagent is briefed with the specific
  finding, file:line, suggested fix, and verification step.
- **Verify before marking done.** Subagent runs `using NeutronStar` and
  test suite for Julia-touching changes.
- **Check in with user** at phase boundaries (after A, B, C, D, E).
- **Commit policy**: one bead = one commit (unless trivial groups like
  pure deletions). Each commit message cites the bead ID and source
  finding.
