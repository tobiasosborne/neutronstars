# Repo Tidiness Review

Reviewer scope: repo hygiene only (gitignore, working-tree noise, file
organisation, naming, doc overlap, vendored assets, build artifacts).
Code, tests, and architecture are out of scope.

Repo: `/home/tobiasosborne/Projects/neutronstars`, branch `master`,
git history 14 commits, working tree 50 MB `.git` + 320 MB working set
(of which 131 MB `output/` and 161 MB `refs/`).

---

## Executive summary — top three messes

1. **Two CLAUDE.md files, neither named CLAUDE.md.** `NEUTRON_STAR_PIPELINE.md`
   (534 lines, committed once, never touched since) and
   `NEUTRON_STAR_CLAUDE_MD_v2.md` (467 lines, committed once, never touched
   since) both literally begin with `# CLAUDE.md — Project Neutron Star`,
   both restate the same "Sacred Principles", both diverge on detail
   (different project-tree sketches, different module ordering). A new
   contributor has no signal which is canonical. README.md (67 lines)
   is the only file Claude actually onboards from, and it points at
   neither. There is no `CLAUDE.md` at the root.

2. **`output/` is half-tracked, half-untracked, with no policy.** 14
   render artefacts (PPM/PNG/GIF/DAT, ~9 MB) are committed; three
   directories of frame stills (`gif_frames/`, `gif_frames_rxj/`,
   `rotation_frames/`, ~120 MB total, 600+ files) are untracked and
   show up dirty on every `git status`. The `.gitignore` rules
   `output/*.ppm` and `output/*.png` are dead — they were either
   superseded by `git add -f`-style commits or simply never applied
   to the existing tracked files. PPMs (uncompressed bitmaps) should
   never be in git; the `rxj1856_xray.ppm` alone is 2.3 MB and is
   100% reproducible from the PNG next to it.

3. **`refs/code/McPHAC` is a nested git repo masquerading as vendored
   source, and it is dirty.** It is not a submodule (no `.gitmodules`,
   no gitlink in the index — it is checked in to our index as a
   subtree with its own `.git/` directory, which is why `git status`
   shows ` m refs/code/McPHAC`). 109 MB. Already contains compiled
   `.o` object files from an upstream build, plus three untracked
   `OUT_*` output directories from the McPHAC run, plus an unstaged
   `Makefile` change. We have no provenance ledger for what we
   modified vs upstream, and no script to refresh it. Pinning this
   as a proper submodule (or a tarball + checksum) is the only sane
   path.

The rest of the repo follows: a `proof/` ledger that is orphaned from
anything else in the project, a 9.5 MB `refs/potekhin_tables/` of
tabulated Fortran-format `.dat` files that contradict the project's own
"zero external table dependencies" principle, LaTeX build artefacts
in `docs/`, and three different naming conventions for top-level files.

---

## Working-tree inventory

`git status --short`:

| Path | Status | Recommendation |
|---|---|---|
| `refs/code/McPHAC` | ` m` (modified, nested repo dirty) | **Convert to submodule.** See "Vendored / large assets". Track upstream SHA + a thin patches/ dir if we need local mods. Either way the inner `.git/` should not be a child of ours. |
| `.claude/` | `??` untracked | **Ignore.** Tool-internal Lean-helpers from another project (`.claude/docs/lean4/`, `.claude/tools/lean4/`); nothing in this repo references them. Add `.claude/` to `.gitignore`. |
| `docs/physics_report.aux` | `??` | **Ignore.** LaTeX build artefact. |
| `docs/physics_report.log` | `??` | **Ignore.** LaTeX build artefact. |
| `docs/physics_report.out` | `??` | **Ignore.** LaTeX hyperref artefact. |
| `docs/physics_report.toc` | `??` | **Ignore.** LaTeX toc artefact. |
| `output/gif_frames/`        | `??` | **Ignore + delete.** 17 MB, ~50 PPM+PNG frame pairs for `ns_spectral_sweep.gif`. The GIF itself is committed; the frames are scratch. |
| `output/gif_frames_rxj/`    | `??` | **Ignore + delete.** 67 MB, ~50 PPM+PNG frame pairs for `rxj1856_spectral_sweep.gif`. Same story. |
| `output/rotation_frames/`   | `??` | **Ignore + delete.** 36 MB, frame pairs for `rxj1856_rotation.gif`. |

Already-tracked things that should not be tracked at all (need
`git rm --cached` + ignore rules):

| Path | Size | Why it's wrong | Recommendation |
|---|---|---|---|
| `Manifest.toml` | 5.2 KB | Pulls in transient lock state of `julia_version = "1.12.5"`. For a library (`Project.toml` declares `name = "NeutronStar"`), Julia community convention is to NOT commit Manifest. | `git rm --cached Manifest.toml`, add to `.gitignore`. |
| `output/*.ppm` (6 files) | 7.7 MB | Uncompressed bitmaps; the matching `.png` is committed. | Delete from history (or at minimum `git rm`); keep PNG. |
| `output/*.png` (8 files) | ~225 KB | Build artefacts of `render.jl`. Smaller, but still pure outputs. | Keep README hero images; move them to `docs/figures/` or `assets/` and `git rm` the rest. |
| `output/*.gif` (3 files)  | 2.2 MB | Hero animations linked from README. | Keep but move to `docs/figures/` so README can still link, and `output/` becomes purely scratch. |
| `output/*.dat` (3 files)  | 14 KB | Spectrum dumps from one-off CLI calls. | Delete — reproducible from `scripts/`. |
| `output/*.txt` (2 files)  | 1 KB | One-line diagnostic dumps from `render_rxj1856_visible_*.jl`. | Delete or move next to the script. |
| `refs/potekhin_tables/*.dat` (50 files) | 9.5 MB | Production-input tables, contradicting "zero external table dependencies" (NEUTRON_STAR_CLAUDE_MD_v2.md §5). Currently committed as ground truth despite that. | See "Vendored / large assets". |
| `refs/potekhin_tables/hmagnet.tar.gz` | 2.1 MB | Tarball containing the same `hmag*.dat` files already extracted next to it. | `git rm`. Pure duplication. |
| `refs/potekhin_tables/eos_potekhin.f` | 0 bytes | Empty file. | `git rm`. |
| `proof/externals/*.json` (19 files) | 6 KB | Provenance hashes for the literature refs (e.g. `heinke_2006 → ae57a6ce…`). Useful in principle but `refs/CHECKSUMS` already does this in one file. | Delete `proof/externals/`; CHECKSUMS is the canonical source. |
| `proof/ledger/000001.json … 000053.json` | 25 KB total | Append-only "proof initialized / external_loaded / approximation_added" events from a tool that was used once on 2026-03-23 and has not been touched since. Not referenced anywhere in src/, scripts/, tests, or HANDOFF. | Delete `proof/ledger/` (and possibly all of `proof/`). |

Things that should be **committed** but currently aren't:

- Nothing. The dirty working tree is 100% artefacts.

---

## `.gitignore` proposal

The current `.gitignore` is 6 lines and doesn't actually work (the
ignored `output/*.png` and `output/*.ppm` patterns were obviously
overridden at some point — those files are tracked). Drop-in
replacement:

```gitignore
# ============================================================
# Julia
# ============================================================
# Manifest is per-environment; commit Project.toml only.
# (Convention for libraries; see https://pkgdocs.julialang.org)
Manifest.toml
*.cov
*.jl.cov
*.jl.*.cov
*.jl.mem
deps/deps.jl
deps/build.log
deps/downloads/
deps/usr/
docs/build/
docs/site/

# ============================================================
# Render output (regeneratable from scripts/)
# ============================================================
output/
# Hero animations linked from README live in docs/figures/, not here.

# ============================================================
# LaTeX build artefacts (docs/physics_report.*)
# ============================================================
*.aux
*.bbl
*.blg
*.fdb_latexmk
*.fls
*.log
*.out
*.synctex.gz
*.toc
*.nav
*.snm
*.vrb
# Keep the compiled PDF: docs/physics_report.pdf is checked in deliberately.

# ============================================================
# Editors / OS
# ============================================================
.DS_Store
Thumbs.db
*~
.#*
\#*\#
*.swp
*.swo

# ============================================================
# Tool-local state (per-project Claude config, IDE state)
# ============================================================
.claude/
.idea/
.vscode/

# ============================================================
# Large vendored references that should not have been committed
# ============================================================
# Currently tracked; will be moved out in cleanup. See reviews/04_tidiness.md.
# refs/potekhin_tables/*.dat
# refs/code/McPHAC/
```

The two commented-out rules at the bottom are aspirational — they
need a separate `git rm --cached` step before they take effect.

A `.gitattributes` should also be added to mark binary files explicitly
(`*.ppm binary`, `*.png binary`, `*.gif binary`, `*.pdf binary`), and a
one-line `.editorconfig` for `end_of_line = lf`, `insert_final_newline =
true`. Neither exists today.

---

## File reorganisation proposal

Current top level:

```
.claude/                       (untracked, ignore)
.git/
.gitignore
HANDOFF.md                     (session-scoped, rewritten every commit)
LICENSE
Manifest.toml                  (should be ignored)
NEUTRON_STAR_CLAUDE_MD_v2.md   (467 lines, never updated)
NEUTRON_STAR_PIPELINE.md       (534 lines, never updated)
Project.toml
README.md
docs/                          (LaTeX + build artefacts)
notes/
output/                        (mixed tracked/untracked, mostly junk)
proof/                         (orphaned tooling)
refs/                          (papers + tables + nested git repo)
reviews/                       (new, this review)
scripts/
src/
test/
verification/
```

Proposed:

```
.gitignore
.gitattributes                 (NEW)
.editorconfig                  (NEW)
CLAUDE.md                      (RENAMED from NEUTRON_STAR_CLAUDE_MD_v2.md,
                                with NEUTRON_STAR_PIPELINE.md merged in
                                or deleted — see "Documentation overlap")
LICENSE
Project.toml
README.md
docs/
    physics_report.tex
    physics_report.pdf
    sections/                  (move all section_*.tex here)
        opacities.tex
        radiative_transfer.tex
        ray_tracing.tex
        ns_structure.tex       (moved out of proof/)
    figures/                   (NEW; move output/*.gif and the canonical
                                README hero PNGs here)
        ns_spectral_sweep.gif
        rxj1856_rotation.gif
        rxj1856_spectral_sweep.gif
        rxj1856_xray.png
        rxj1856_true.png
notes/
    HANDOFF.md                 (MOVED from root; this is a note, not a
                                spec)
    decisions.md
    approximations.md
    failures.md
    sessions/                  (NEW; old HANDOFFs roll over here, see
                                "Bigger refactors")
output/                        (GITIGNORED; .gitkeep only)
refs/
    INDEX.md
    CHECKSUMS
    *.pdf                      (literature)
    code/                      (vendored code as proper submodule,
                                  see "Vendored / large assets")
    tables/                    (RENAMED from potekhin_tables/, kept
                                  only if we explicitly decide to break
                                  the zero-tables principle for
                                  verification fixtures)
reviews/                       (this dir)
scripts/
src/
test/
verification/
    data/
    figures/
    *.jl, *.py
```

Specific moves:

- `proof/section_ns_structure.tex` → `docs/sections/ns_structure.tex`
  (the LaTeX root `docs/physics_report.tex` already does
  `\input{../proof/section_ns_structure}`, which is a code smell — it
  reaches outside `docs/` for one section but not for the other three).
- `docs/section_*.tex` → `docs/sections/*.tex` for uniformity.
- `HANDOFF.md` → `notes/HANDOFF.md` (or replace with `notes/sessions/`
  per-session entries, see Bigger refactors).
- `NEUTRON_STAR_CLAUDE_MD_v2.md` → `CLAUDE.md` after merging the
  non-duplicate parts of `NEUTRON_STAR_PIPELINE.md`.
- `output/rxj1856_rotation.gif`, `output/rxj1856_spectral_sweep.gif`,
  `output/ns_spectral_sweep.gif` → `docs/figures/` (so they survive
  ignoring `output/`).
- All other `output/*` → delete.
- `refs/potekhin_tables/hmagnet.tar.gz` → delete.
- `refs/potekhin_tables/eos_potekhin.f` (0 bytes) → delete.
- `proof/` → delete (or, if the JSON ledger has value, move to
  `notes/proof_ledger/` with a README explaining what tool consumes
  it; right now nothing does).

---

## Documentation overlap

There are seven prose/spec documents in this repo and the boundaries
between them are unclear.

| File | Stated purpose | Observed content | Update cadence | Recommendation |
|---|---|---|---|---|
| `README.md` (67 lines) | GitHub landing page | Goal statement, pipeline summary, RX J1856 hero images, current status, quick start, refs, license. The only doc a new reader actually starts from. | Updated alongside releases (~3 commits touched it). | **Keep.** Add a one-line pointer to `CLAUDE.md` for design intent and to `notes/HANDOFF.md` for current session state. |
| `NEUTRON_STAR_CLAUDE_MD_v2.md` (467 lines) | Effectively a CLAUDE.md (first line: `# CLAUDE.md — Project Neutron Star`). Mission, principles, project structure, pipeline modules, phasing. | Same as PIPELINE but a different revision; sketches a "graph/" subdir that does not exist; module breakdown is the canonical one. | **Never updated after the one commit that introduced it** (`8926587`, 2026-03-24). | **Rename to `CLAUDE.md`.** This is the live spec. Merge any non-duplicate content from PIPELINE.md (Phase 2 tracer-bullet recipe, Module 2c bound-free chain). |
| `NEUTRON_STAR_PIPELINE.md` (534 lines) | Also literally starts `# CLAUDE.md — Project Neutron Star: Physically Exact Visualisation from First Principles`. | Earlier draft of the same document with more end-to-end recipe detail (Phase 2 step-by-step, tracer bullet integration test). | **Never updated since initial commit** (`ff575fb`). | **Delete after extracting the few unique sections** into CLAUDE.md. There is no defensible reason to have two CLAUDE.mds, both stale, both at root, neither named CLAUDE.md. |
| `HANDOFF.md` (375 lines) | Session-to-next-session handoff. | Reads like a project status report dated 2026-05-16: working/not-working list, files changed, full reproduction commands. | **Rewritten every session** (8 commits explicitly mention it; history shows ±300 lines per session). | **Move to `notes/HANDOFF.md` and split.** See Bigger refactors — the every-session rewrite pattern destroys history. Append-only `notes/sessions/2026-05-16.md` + a thin stable `HANDOFF.md` that always points at the latest is the standard pattern. |
| `notes/decisions.md` | Per-decision ADR-style entries (D1, D2, …, D7). | Exactly that. 2.8 KB, 7 entries from 2026-03 and 2026-05. | Append-only. | **Keep as-is.** Healthy. |
| `notes/approximations.md` | Master list of named approximations with validity, source, impact. | Exactly that. Cross-references af-graph node IDs (`af node: 1.3.1`) but the graph itself isn't checked in anywhere. | Append-only. | **Keep.** Add a note that "af node" IDs are dangling — the af tool is not part of this repo. Or remove the IDs. |
| `notes/failures.md` | Failure post-mortems (F1–F4). | Exactly that. Useful. | Append-only. | **Keep.** |
| `docs/physics_report.tex` + sections | 30-page equation-by-equation code review | LaTeX doc with `section_opacities`, `section_radiative_transfer`, `section_ray_tracing`, and `section_ns_structure` (the last lives in `proof/`, all others in `docs/`). | Single commit, will likely be updated rarely. | **Keep.** Move all sections into `docs/sections/` (currently `proof/section_ns_structure.tex` is referenced via `\input{../proof/...}`, which is bad). Add `.aux/.log/.out/.toc` to gitignore (currently dirty). |
| `proof/section_ns_structure.tex` | TOV+EOS section of the physics report | Same as the other sections, just in the wrong directory. | One commit. | **Move into `docs/sections/`.** |
| `proof/meta.json`, `proof/ledger/`, `proof/externals/` | "Proof" semantic graph from an external tool (af) that was run once on 2026-03-23. | 72 tiny JSON files. Not consumed by anything in this repo. The "external" hashes duplicate `refs/CHECKSUMS`. | Frozen since 2026-03-23. | **Delete the directory.** If you want to keep the audit trail, archive to `notes/historical/proof_2026-03-23/` and never look at it again. |
| `refs/INDEX.md` | Bibkey → file → purpose table for the PDFs in `refs/`. | Exactly that. Useful. | Updated alongside new papers. | **Keep.** Cross-check entries against actual files (haven't audited line-by-line but list looks complete). |
| `verification/VERIFICATION_LOG.md` | Append-only verification results per module. | Healthy. EOS, ray tracer, colorimetry, magnetic opacity, end-to-end. | Append-only. | **Keep.** |

Net effect of the documentation cleanup:

- 2 top-level docs (README + CLAUDE) instead of 4.
- All physics-report `.tex` in `docs/sections/`.
- All notes (decisions, approximations, failures, HANDOFF, verification log)
  under `notes/` and `verification/`.

---

## Vendored / large assets

### `refs/code/McPHAC/` — 109 MB nested git repo

- Not a submodule. There is no `.gitmodules`, and `git submodule status`
  errors. It is just a directory with its own `.git/` inside. That is
  the source of the persistent ` m refs/code/McPHAC` in `git status`:
  git is reporting that the nested checkout has uncommitted local
  changes (and indeed it does — `Makefile` is modified, `OUT/.dummy` is
  deleted, three `OUT_T*/` output dirs are untracked).
- The upstream is `https://github.com/McPHAC/McPHAC.git`, pinned at
  `ad6df20 Cleaned up indentation and a few comments.`
- The checked-in copy includes pre-compiled `.o` object files
  (`CalcAF.o`, `CalcBF.o`, … ~20+ of them). These were built on
  someone's machine and contaminate the index.
- The README's quick-start example reads `refs/code/McPHAC/gffgu.dat`,
  so we DO depend on at least one file inside this tree.

Recommendations, in order of effort:

1. (Now) Add an explicit note in `refs/INDEX.md` describing how this
   directory was acquired and the pinned SHA.
2. (Soon) Convert to a proper submodule: `git rm -r --cached
   refs/code/McPHAC`, delete the inner `.git`, `git submodule add
   https://github.com/McPHAC/McPHAC.git refs/code/McPHAC`, pin to
   `ad6df20`. Now `git status` is clean and provenance is explicit.
3. (Cleanup) Add a `Makefile` patch under `refs/code/patches/` if
   the local Makefile change is meaningful; otherwise revert it.

### `refs/potekhin_tables/` — 9.5 MB of `.dat` files + tarball

- 50 files: `hmag1[1-3]_[0-9].dat` (Rosseland), `hmm1[0-2]_[0-9].dat`
  (monochromatic), `hmn1[3-5]_[0-9].dat`. Plus `hmagnet.tar.gz` which
  contains the same files re-archived (pure duplication).
- Also `eos_potekhin.f` (0 bytes, garbage).
- These tables are the kind of thing the project's own CLAUDE.md
  explicitly disclaims as production inputs ("zero external table
  dependencies"). The repo currently both swears off them and ships
  them.
- `verification/potekhin_table_comparison.jl` reads them, so they are
  legitimately needed as VERIFICATION fixtures.

Recommendations:

1. Delete the tarball (`refs/potekhin_tables/hmagnet.tar.gz`) and the
   zero-byte `.f`.
2. Rename `refs/potekhin_tables/` → `verification/fixtures/potekhin_tables/`
   to align directory with usage. They are not "references"; they are
   verification ground truth.
3. Add a `verification/fixtures/README.md` explaining what each table
   is, where it came from, and that they are NOT to be used as
   production inputs.
4. Consider Git LFS if these grow further. 9.5 MB across 50 files is
   borderline acceptable in plain git; another doubling and it
   wouldn't be.

### `verification/data/suleimanov_2009_fig2/` — 1.2 MB CSV + JSON

- `digitized_points.csv` (610 KB, 6846 lines) is figure-digitised
  points from a published PDF.
- `computed_opacities.csv` (530 KB, 4201 lines) is the output of
  `compute_suleimanov_fig2_opacities.jl` for direct comparison.
- Three small `*_metrics.json` files.

Assessment: this is the correct call. The CSVs are the actual
verification artefact and are slow to regenerate (the .jl script
involves opacity calls at many sampling points). 1.2 MB across 5
files is fine in git. **Keep.**

One nit: `computed_opacities.csv` is downstream of code that changes
fairly often. Worth adding a header comment recording the git SHA it
was generated from, so a stale CSV is easy to spot.

### `verification/figures/suleimanov_2009_fig2/` — 1.3 MB PNGs

- 5 PNGs, biggest is 400 KB.
- These are figures the user looks at to confirm validation status,
  referenced from `notes/failures.md` (F3) and `verification/VERIFICATION_LOG.md`.

Assessment: **Keep.** Acceptable size, real documentary value.

### `output/` — 131 MB total, mostly junk

Already covered above. Summary of action: ignore the whole directory;
promote 3 hero GIFs + 2 hero PNGs to `docs/figures/`; delete
everything else from history (or at minimum stop tracking).

Note that the .git directory is currently 50 MB. After `git rm`ing
the PPMs and Manifest, history will still contain them; if size
matters, a `git filter-repo` pass on `output/*.ppm` is the only way
to get the .git down. Probably not worth it unless someone is
constantly cloning fresh.

---

## Quick-win cleanups (each <5 min)

1. `git rm --cached Manifest.toml` and add to `.gitignore`. Manifest
   is per-Julia-version and per-environment; committing it means
   every contributor on a different Julia version generates noise.
2. `git rm` the LaTeX build artefacts: `docs/physics_report.{aux,log,out,toc}`
   are already untracked but are regenerated on every `latexmk` and
   will start appearing in `git status` again. Add `*.aux`, `*.log`,
   `*.out`, `*.toc`, `*.fls`, `*.fdb_latexmk` to `.gitignore`.
3. Add `.claude/` to `.gitignore`. It's tool-local Lean-helper state
   leaked from another project and has no business in this tree.
4. Add `output/gif_frames/`, `output/gif_frames_rxj/`,
   `output/rotation_frames/` (or simply `output/`) to `.gitignore`
   and `rm -rf` them. 120 MB of scratch frames recovered.
5. `git rm refs/potekhin_tables/hmagnet.tar.gz` (duplicate of
   already-extracted files) and `refs/potekhin_tables/eos_potekhin.f`
   (0 bytes). 2.1 MB recovered, one weird file gone.
6. Move `proof/section_ns_structure.tex` → `docs/sections/ns_structure.tex`
   and update the `\input{}` in `docs/physics_report.tex`. Stops the
   cross-directory `\input{../proof/...}` smell.
7. Rename `NEUTRON_STAR_CLAUDE_MD_v2.md` → `CLAUDE.md`. Even before
   merging in PIPELINE.md content, this immediately gives Claude
   (and humans) the file they expect to find.
8. Delete `NEUTRON_STAR_PIPELINE.md` after copy-pasting the few
   uniquely-useful Phase-2 paragraphs (the tracer-bullet recipe in
   §2.7) into CLAUDE.md.
9. Move `HANDOFF.md` → `notes/HANDOFF.md` and add a `notes/README.md`
   that lists what each notes file is for.
10. Create `.gitattributes` with `*.ppm binary`, `*.png binary`,
    `*.gif binary`, `*.pdf binary`, `*.dat binary`. Prevents
    accidental line-ending munging.

---

## Bigger refactors

### R1: Convert McPHAC to a real submodule

Effort: ~30 min, plus one round of testing the build path
`refs/code/McPHAC/gffgu.dat` still resolves correctly.

Steps:
- `git rm -r --cached refs/code/McPHAC`
- `rm -rf refs/code/McPHAC` (and back it up first, in case the
  Makefile/OUT changes have any value)
- `git submodule add https://github.com/McPHAC/McPHAC.git refs/code/McPHAC`
- `cd refs/code/McPHAC && git checkout ad6df20`
- `git add .gitmodules refs/code/McPHAC` (now records the gitlink)
- Document the submodule init step in README's quick-start.

Payoff: `git status` becomes clean; provenance becomes explicit;
clone size drops by ~109 MB; nobody is ever again confused about
whether modifying `refs/code/McPHAC/Makefile` is "real" work.

### R2: Stop rewriting HANDOFF.md every session

Current pattern: every session, HANDOFF.md is overwritten with
±300-line diffs. Git history of HANDOFF.md shows 8 commits, all of
them substantive rewrites. The result is that "the handoff" is
whatever the current HEAD has, and the previous session state is
buried in `git log -p`.

Proposed pattern:

- `notes/sessions/2026-05-16.md`, `notes/sessions/2026-05-14.md`,
  etc. — append-only, one file per session. Each is the snapshot of
  "what changed in this session, what works, what's broken".
- `HANDOFF.md` at root (or `notes/HANDOFF.md`) becomes a thin
  ~50-line file: "current status", "next actions", and a link to
  the most recent session file for detail.
- Or simpler: `HANDOFF.md` IS the latest session file, and at the
  start of each session you `git mv HANDOFF.md notes/sessions/<date>.md`
  and `cp` a template to a fresh HANDOFF.md.

Payoff: history of "what we knew, when" becomes browsable. The
375-line megadoc that gets re-written wholesale stops being the
unit of change.

### R3: Reconcile "zero table dependencies" with `refs/potekhin_tables/`

CLAUDE.md §5 declares zero external tables in production. In
practice, `refs/potekhin_tables/` ships 50 production-format `.dat`
files and `verification/potekhin_table_comparison.jl` reads them.

Either:

- (Honest) Move `refs/potekhin_tables/` → `verification/fixtures/potekhin_tables/`
  and update CLAUDE.md §5 to say "Published tables are used only as
  verification fixtures, stored under `verification/fixtures/`".
- (Strict) Replace the table comparison with one that recomputes
  from Potekhin's published formulae and only stores a handful of
  reference points (~10 KB) for regression.

The first is the realistic call. Either way, the current state
(directory and principle disagree) is the worst of both worlds.

### R4: Decide the fate of `proof/`

The `proof/` directory contains:
- `meta.json` (one line)
- `externals/` (19 JSON hashes of refs PDFs)
- `ledger/` (53 append-only event JSONs from 2026-03-23)
- `section_ns_structure.tex` (one section of the physics report, in
  the wrong place)

`refs/CHECKSUMS` already does the externals job. Nothing in `src/`,
`scripts/`, or `test/` reads anything in `proof/`. The `af` tool that
generated the ledger is not in this repo, and a search for "af node"
in `notes/approximations.md` shows the references go nowhere.

Recommendation: move the one `.tex` file into `docs/sections/`, then
delete `proof/` entirely. If for archival reasons the ledger feels
load-bearing, move it to `notes/historical/proof_ledger_2026-03/` and
write a one-line README explaining that this is frozen audit history
from a tool no longer in use.

### R5: Pin a render-from-scratch script

Almost every artefact in `output/` is reproducible, but you have to
know the right `scripts/render_*` or `julia -e '…'` incantation. As
part of moving `output/` to gitignored:

- Add `scripts/render_all_hero_assets.jl` that produces
  `docs/figures/{ns_spectral_sweep,rxj1856_rotation,rxj1856_spectral_sweep}.gif`
  plus the README hero PNGs, end-to-end.
- README says: "Hero animations: `julia --project=. scripts/render_all_hero_assets.jl`".

Now `output/` can be deleted with zero anxiety because anyone can
regenerate the official assets. This is what makes "`output/` is
gitignored" sustainable.

---

## Things that are fine (not all news is bad)

- `notes/{decisions,approximations,failures}.md` is a healthy
  append-only documentation pattern. Good as-is.
- `verification/VERIFICATION_LOG.md` is the right way to record
  test-vs-published-value comparisons. Good as-is.
- `refs/INDEX.md` + `refs/CHECKSUMS` for the literature PDFs is
  well-organised.
- `src/` directory tree (atmosphere/, opacity/, eos/, geodesics/,
  surface/, colorimetry/, pipeline/) is clean, snake_case throughout,
  and matches the module structure declared in CLAUDE.md.
- `Project.toml` is minimal and focused (5 deps).
- LICENSE is present and matches what README claims (GPL-3.0).

---

## Naming consistency notes

- `src/` uses `snake_case.jl` throughout — consistent.
- Top-level docs use a mix of `SCREAMING_SNAKE.md` (HANDOFF, LICENSE,
  NEUTRON_STAR_CLAUDE_MD_v2, NEUTRON_STAR_PIPELINE) and `PascalCase.md`
  (README) and `Pascal.toml` (Project, Manifest). The `_v2` suffix on
  `NEUTRON_STAR_CLAUDE_MD_v2.md` implies a `_v1` lurking somewhere
  (there isn't one — it was renamed in-place, leaving a dangling
  version marker). Cleanup: rename to `CLAUDE.md`, drop the version.
- Output files have a partial convention: `ns_*` for generic, `rxj1856_*`
  for the specific target. Internally consistent. After moving to
  `docs/figures/`, consider a subdir per target (`ns_canonical/`,
  `rxj1856/`).
- `scripts/render_rxj1856_visible_atmosphere.jl` and
  `scripts/render_rxj1856_visible_magnetic_atmosphere.jl` are 4.5 KB
  and 6.0 KB scripts that overlap heavily by name. Worth checking if
  they share enough that one of them is a drop-in superset — but
  that's code review, out of scope here.

---

## Final tally

- **Tracked files that should be deleted or gitignored:** 60+
  (Manifest.toml, all output/*.ppm/png/dat/txt except hero
  animations, refs/potekhin_tables/hmagnet.tar.gz,
  refs/potekhin_tables/eos_potekhin.f, all of proof/externals/ and
  proof/ledger/).
- **Untracked dirt that needs ignoring:** ~120 MB across
  `output/gif_frames/`, `output/gif_frames_rxj/`,
  `output/rotation_frames/`, plus `.claude/` and LaTeX build files.
- **Top-level docs that should collapse from 4 to 2:** README +
  CLAUDE.md; HANDOFF moves into `notes/`, PIPELINE deleted, the
  `_v2` suffix dropped.
- **Nested git repo masquerading as vendored source:** 1, needs
  conversion to a proper submodule.
- **Stale tool ledger orphaned from the codebase:** 1, can be deleted.

The repo isn't unsalvageable — the code tree (`src/`, `test/`,
`verification/`, `notes/decisions+approximations+failures`) is
actually well-structured. The mess is concentrated at the top level
and in `output/` + `refs/code/`. The quick wins above would fix
~80% of the visible damage in under an hour.
