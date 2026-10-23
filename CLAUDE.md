# Neutron Star Flyby — rules for agents

Read `PRD.md`, then `HANDOFF.md`, then the newest shard in `worklog/`.
Discover everything else as you need it. Nothing here is a tutorial.

## Rules

1. **Red-green TDD.** Write the failing test first, make it pass, refactor.
   No production code without a failing test that motivated it.

2. **Mutation testing and fuzzing where they bite.** Use them on numerical
   kernels, parsers and coordinate transforms, where a wrong sign or an
   off-by-one survives unit tests. Skip them where they only slow the loop.

3. **Ground truth is on disk, quoted byte-for-byte.** Every physical equation,
   constant, or factual claim in code or docs cites a local file under
   `refs/` and quotes it verbatim. The quote must be byte-identical to the
   local copy, so `grep -F` finds it. PDFs are quoted from their sibling
   `refs/<key>.txt` (made with `pdftotext -layout`). Observational inputs
   (Gaia, HST photometry, CIE tables) live under `refs/data/` with checksums
   in `refs/CHECKSUMS`. Nothing is cited from memory. "Well known" is not a
   source.

4. **Handoff and worklog, sharded.** `HANDOFF.md` is the current state,
   under 60 lines, rewritten each session, never appended. `worklog/YYYY-MM-DD.md`
   is one append-only shard per session. Update both before you stop.

5. **Fail loud.** No silent fallbacks, clamps, extrapolations, or NaN
   propagation. Assert on entry, assert on exit, raise with the offending
   values in the message.

6. **Verification prints numbers.** Expected value with source, computed
   value, error, tolerance. "Test passed" on its own is not evidence.

7. **Own the pipeline.** Every number the renderer consumes is either a
   measurement (Gaia, HST, CODATA), a convention (sRGB, CIE), or computed
   from first principles by code in this repo. Published fits, interpolation
   tables and other people's code outputs are verification targets, never
   runtime inputs. Every `truth/` entry declares its `provenance` class, and a
   `fit` or `table` entry must name the physics it will be recomputed from, or
   the reason the pipeline no longer needs it. Recomputation is a tree explored
   breadth-first: a stub (published fit) may stand in for a subtree so results
   come quickly, but a node is green only when every leaf below it is
   recomputed here. `truth/check.py --render` colours the tree.

## Layout (grows by discovery)

- `PRD.md` — what we are building and how we will know it is done.
- `HANDOFF.md`, `worklog/` — state and history.
- `refs/` — local sources, extracted text, data, checksums.
- `archive/v1/` — the deprecated predecessor. Use it for ground truth and
  salvage only. Never import from it at runtime.

`AGENTS.md` is a symlink to this file.
