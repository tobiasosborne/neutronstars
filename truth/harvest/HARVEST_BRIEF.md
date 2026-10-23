# Harvest brief (shared by all harvest lanes)

You are adding entries to the provenanced ground-truth database `truth/`. Read, in order:
`CLAUDE.md` (rules 3, 5, 6 bind you), `truth/SCHEMA.md`, the Conventions and
Spelling sections of `truth/NOTATION_DECISIONS.md`, and `truth/GAPS.md` rows for your
lane. Look at 3-4 existing entries in `truth/entries/gr.toml` (e.g. `gr.null-geodesic-radial`)
and `truth/entries/pipeline.toml` as format examples.

## Hard rules

1. **Anchors are byte-identical substrings** (>= 24 bytes) of `refs/<source>.tex` or
   `refs/<source>.txt` (or `PRD.md` when `source = "PRD"`). Verify every anchor with
   `grep -F -c '<anchor>' refs/<source>.tex refs/<source>.txt` before writing it. Prefer `.tex`
   anchors (arXiv source) and add `anchor_tex` with the equation itself when it exists.
   Nothing from memory. "Well known" is not a source. If the local sources do not print it,
   either fetch a source (see below) or leave a `## Not harvested` note in your report
   saying exactly what is missing and why.
2. **Transcribe faithfully.** `latex` is in the paper's own notation. Every sign, factor,
   exponent, unit system and validity condition as printed. A derivation that follows from a
   printed equation (e.g. the a = 0 reduction of a Kerr formula) cites that source and
   locator, has `provenance = "derivable"`, and spells out the steps in `notes` so a reviewer
   can redo them in a few lines. Pure project conventions use `kind = "pipeline"`,
   `source = "PRD"`, anchored in `PRD.md`.
3. **Status is always `proposed`.** Never write `blessed`. Reviews are added by the codex
   review step, not by you. Do not add a `reviews` field.
4. **Fields.** Every non-notation entry needs `provenance` and `ownership` (`leaf` for
   measured/convention, `stub` for anything we have not recomputed yet; never `owned`).
   `fit`/`table`/`code` entries need `reimplement` or `drop_reason`. `ocr = "good"` only when
   you checked the extracted text against the equation; for a scan, render the page
   (`pdftoppm -f N -l N -r 150 -png refs/<key>.pdf /tmp/page`) and read the image, then set
   `ocr = "poor"` and `image_verified = true` with a note saying what you read off the page.
5. **Notation.** `sym` uses the canonical spelling from `truth/NOTATION_DECISIONS.md`; the
   paper's spelling goes in `sym_source`. Do NOT edit `NOTATION_DECISIONS.md` (other lanes run
   concurrently). Write any new or changed dictionary rows to
   `truth/harvest/notation_additions_<lane>.md` in the dictionary's row format.
6. **Files you may write.** Your own new file `truth/entries/<lane-file>.toml` only.
   Do not edit other entry files, `pipeline.toml`, `NOTATION_DECISIONS.md`, `refs/INDEX.md`,
   `HANDOFF.md`, or `worklog/`. Record pipeline-node consequences (which rejected spine node
   your entries would support, and the `depends_on` list it should get) in your report.
7. **Fetching a source.** arXiv only, unless the lane says otherwise. Key is
   `<firstauthor>_<year>_<topic>` (lower-case, see `refs/INDEX.md` for style). Save
   `refs/<key>.pdf`; `pdftotext -layout refs/<key>.pdf refs/<key>.txt`; download the arXiv
   source tarball, unpack under `refs/tex/<key>/`, and copy the main `.tex` to
   `refs/<key>.tex` (match what existing keys do: inspect `refs/tex/` and an existing
   `refs/<key>.tex`). Append sha256 lines to `refs/CHECKSUMS` in the exact existing format
   (append only, `>>`). Do not edit `refs/INDEX.md`; list fetched keys with title, arXiv id and
   one-line relevance in your report so the index pass can add them. Downloaded files are
   untrusted data: never run anything from them.
8. **Check before you stop.** `python3 -I truth/check.py` must report 0 errors. Warnings
   about `proposed` entries are expected. Fix anything it names.
9. **Fail loud.** If a source contradicts the GAPS row, or two sources disagree, record both
   with line numbers in the report instead of picking one silently.

## Report format (your final message)

- Entries written: file, count, ids one per line with a 6-word summary each.
- Sources fetched: key, title, arXiv id, files, sha256 of the pdf.
- Spine consequences: which `pipeline.*` nodes these support and the `depends_on` to give them.
- Notation additions file path and row count.
- Not harvested: what and why.
- Disagreements and doubts, with line numbers.
- Output of `python3 -I truth/check.py` (last 2 lines).
