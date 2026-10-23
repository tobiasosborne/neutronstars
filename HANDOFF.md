# HANDOFF — 2026-10-23 (end of day)

## State

Reboot day. v1 is archived under `archive/v1/`; the new project has rules
(`CLAUDE.md`, 7 rules), a PRD (`PRD.md` v0.1: free-fall flyby of RX J1856
against the real sky), and a **blessed, provenanced ground-truth database**
under `truth/`. No pipeline code exists yet. Nothing is committed since the
restart; the working tree holds the whole day (125 renames + new files).

## truth/ at a glance

- `truth/entries/*.toml` — 387 entries: **299 blessed, 88 rejected, 0 proposed**.
  Every blessed entry has a byte-identical anchor in `refs/`, a provenance
  class, an ownership flag, and a codex (GPT-6.1 sol, xhigh) accept verdict
  on record (`truth/reviews/*.md`). Rejected rows are kept as negative results.
- `truth/NOTATION_DECISIONS.md` — **canonical notation dictionary**, 370 rows,
  CGS conventions. Read it before writing any code or entry.
- `truth/entries/pipeline.toml` — the DAG spine: 30 nodes rooted at
  `pipeline.rendered-frame`; 22 blessed, 8 rejected for lack of supporting
  entries (camera worldline, tetrad/projection, sky radiance, exposure/tone,
  envelope T_b-T_s, atmosphere solve, polarisation modes, condensed radiance).
- `truth/TREE.md` — recomputation tree: 146 green (measured/convention
  leaves), 241 yellow (stubs: nothing recomputed yet), 0 red.
- `truth/GAPS.md` — ranked work order for the next harvests (120 rows).
- `truth/DEPENDENCY_AUDIT.md` — what v1 consumed at runtime and the
  recompute chain for every fit/table (43 items).
- `truth/check.py` (+ `truth/tools/red_tests.py`, 7 mutations) validates
  everything; `--render` writes NOTATION.md, INDEX.md, TREE.md.

## refs/

51 sources, all checksummed (`refs/CHECKSUMS`, `refs/INDEX.md`). 30 have
arXiv TeX in `refs/<key>.tex` (+ trees under `refs/tex/`). Four inherited
mis-filed PDFs were replaced. 17 recomputation-chain papers added today.

## Next small step

1. Commit (user's call). Then `python3 -I truth/check.py --render` and
   `python3 -I truth/tools/red_tests.py` must both be clean before any commit.
2. Harvest the Milestone-1 gaps in `truth/GAPS.md` §1 (exterior line
   element, timelike geodesics, static tetrad + boost, camera projection,
   Planck function, IAU constants, Gaia passbands) with the same
   harvest -> codex r1 -> r2 workflow (`truth/tools/codex_review.sh`).
3. Then Milestone 1 code, red-green: Schwarzschild null geodesic ODE and
   its flat-space test, with `impl`/`verified_by` turning nodes green.

## Known loose ends

- Notation review applied (240 accept / 124 fix / 6 reject). 27 fix rows are
  marked "r1 fix pending" in `truth/NOTATION_DECISIONS.md` ("Review r1
  outcome" section): they need a different source/anchor/locator, i.e. a
  harvest, not a notation edit. The reviewer's reject of `e` (SI vs CGS) does
  not reproduce (3.3e-11 vs claimed 2.7e-10); row kept as disputed.
- Non-notation entries still use their source spellings; the dictionary's
  "Pipeline-level spellings" table maps them. Code must use dictionary symbols.
- Spine applier flagged repairs outside its lane: remove edge
  `opacity.approx-vacuum-low-energy -> opacity.dielectric-tensor-polarizability`;
  re-harvest the rejected `opacity.ff-*` and `rt.emergent-intensity-mcphac`.
- `eos.bsk-mmax-central-density`: anchor says 2.27 M_sun, TeX says 2.16.
- Slattery 1982 is on disk (user fetched it); CODATA 2018 and 4 scans have no TeX.
