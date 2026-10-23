# HANDOFF — 2026-10-23 (afternoon session, end)

## State

Reboot day, two sessions. Committed and pushed (`06a70b4` + this session's close).
`truth/` is a blessed, provenanced ground-truth DB; **no pipeline code exists yet**.
Orchestration used Opus/Sonnet/Haiku subagents (no Fable) and codex GPT-6.1 sol
xhigh reviews run as background bash; recipe in `truth/harvest/HARVEST_BRIEF.md`.

- `truth/entries/*.toml` — **493 entries: 395 blessed, 97 rejected, 1 proposed**
  (`obs.mellinger-photometric-calibration`, locator fixed after r3, needs owner bless).
  New files today: gr_geodesic, camera, const_iau, radiance_display, surface_field,
  sky, sky2, display_psf. Every blessed entry has a byte-identical anchor and a codex
  accept (or applied fix for `kind = pipeline`); `reject_resolved = true` marks
  superseded r1 rejects.
- Spine `pipeline.toml`: **24 blessed / 6 rejected**. Rejected: sky-radiance (data
  files not on disk), display-exposure-tone (no glare model), envelope-tb-ts,
  atmosphere-solve, polarisation-modes, surface-radiance-condensed (M3 harvests).
- Tree: 177 green / 316 yellow / 0 red. All yellow = no code yet.
- `truth/NOTATION_DECISIONS.md`: 471 symbols; today's additions and collision
  decisions in "Additions 2026-10-23". e[esu]: code computes e = sqrt(alpha_f hbar c)
  in CGS; the 10 c e_SI bridge is a 2.7e-10 cross-check only.
- `refs/`: 62 keys, 190 checksums OK, INDEX.md regenerated. `.gitattributes` pins
  `refs/** -text`; `refs/tex/` trees are ignored (regenerable from arXiv ids).
- `truth/tools/red_tests.py` is now honest (was vacuous: missing PRD.md crashed the
  baseline). `codex_review.sh` prompt states .tex anchors are valid; do not pass
  apostrophes in its EXTRA argument (wrapper exits 127 after writing the review).
- `truth/analysis/`: six codex **gpt-6-astra xhigh** analytic arguments on the
  research-grade paths (own A34 fit, envelope T_b–T_s, exact Gaunt, condensed
  surface, vacuum resonance, extinction own route); briefs in `briefs/`, runner
  `truth/tools/codex_analysis.sh`. Launched ~13:55; check `ls truth/analysis/*.md`.

## Owner decisions needed

1. Language for Milestone-1 code (PRD recommends Julia; v1 microphysics is Julia).
2. Camera projection: pinhole vs equidistant (`gr.camera-projection-convention`).
3. Diffuse sky map: Mellinger is not starless and has no licence on disk; options
   in `pipeline.sky-radiance` notes (mask+inpaint at Gaia positions / other map /
   Gaia total minus catalogued flux).
4. Three blocked notation fixes (alias merge on `gr.notation-geometry`; two I_nu
   call-form rows on `rt.notation-mcphac`), see "r1 fix blocked" in the dictionary.
5. Jackson *Classical Electrodynamics* units appendix (manual) for
   `const.esu-from-coulomb` (rejected r2); optional given decision above.

## Next small step

1. Read the six `truth/analysis/*.md` when they land; commit them; fold their
   "missing sources" lists into `truth/GAPS.md` §2–3.
2. Sky data fetch (checksummed, `refs/data/`): Gaia DR3 rows around RX J1856,
   EDR3 passbands (VizieR), XP sampled spectra, CALSPEC `alpha_lyr_mod_002.fits`,
   Hipparcos rows; then re-propose `pipeline.sky-radiance`.
3. Milestone-1 code, red-green, in the chosen language: Schwarzschild null
   geodesic ODE with the flat-space test and Pechenick/Beloborodov benchmarks
   (`gr.dphi-max-r16`, `gr.circle180-radius-r16`, 227.18° grazing check), turning
   `gr.*` nodes green via `impl`/`verified_by`.

Gate before any commit: `python3 -I truth/check.py --render` and
`python3 -I truth/tools/red_tests.py` both clean.
