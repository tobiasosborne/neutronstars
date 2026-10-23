# truth/ — provenanced ground-truth database

Canonical store for every definition, notation, equation, constant, result,
approximation and observational datum the reboot relies on. Entries live in
`truth/entries/<domain>.toml` as `[[entry]]` tables. `truth/check.py` validates
them and renders `truth/NOTATION.md` and `truth/INDEX.md` from blessed entries.

Run: `python3 -I truth/check.py` (exit 0 = clean). `--render` also writes the
rendered files.

## Entry fields

| field | req | meaning |
|-------|-----|---------|
| `id` | yes | unique, `domain.kebab-name` |
| `kind` | yes | `notation` / `definition` / `equation` / `constant` / `result` / `approximation` / `datum` / `pipeline` (a computation of ours, anchored in `PRD.md` with `source = "PRD"`) |
| `domain` | yes | `gr`, `rt`, `opacity`, `eos`, `surface`, `colour`, `const`, `obs`, `notation` |
| `name` | yes | short human name |
| `statement` | yes | what it says, in words, self-contained |
| `latex` | eq/const/def | the relation, in the paper's own notation |
| `symbols` | notation | list of `{sym, meaning, units}`; `units` is `"1"` when dimensionless |
| `value`, `units`, `uncertainty` | constant/datum | numeric value as in the source |
| `source` | yes | key of `refs/<key>.pdf` (must exist) |
| `locator` | yes | `"eq. 3"`, `"p. L86"`, `"Table 1"` |
| `anchor` | yes | verbatim quote, ≥ 24 bytes, byte-identical substring of `refs/<source>.tex` (arXiv source, preferred) or `refs/<source>.txt` (pdftotext) |
| `anchor_tex` | no | a second verbatim anchor, byte-identical in `refs/<source>.tex`, ideally the equation itself; added after review, never replacing `anchor` |
| `ocr` | yes | `"good"` or `"poor"`; poor means the extracted equation text is unreliable |
| `image_verified` | if ocr poor | `true` only after a reviewer read the PDF page image |
| `provenance` | all but notation | `measured` (observation or CODATA), `convention` (a standard's definition), `derivable` (follows from first principles; we recompute it), `fit` (a published fit to someone else's calculation), `table` (a numerical table), `code` (output of external code). `fit`/`table`/`code` entries also need `reimplement` = the physics route or `drop_reason` |
| `ownership` | all but notation | `owned` (recomputed by code in this repo, with a verification test), `stub` (temporarily using the published fit/table/code), `leaf` (measured or convention; owned by definition) |
| `impl` | if owned | `src/<path>` implementing it; `verified_by` = `test/<path>` |
| `status` | yes | `proposed` / `blessed` / `rejected` |
| `reviews` | blessed | list of `{reviewer, model, verdict, note}`; blessed needs ≥ 1 `accept` from `model` starting `gpt` and no open `reject`. For `kind = pipeline` (our own computations) an applied `fix` verdict from a gpt model also suffices |
| `depends_on` | no | ids |
| `v1_ref` | no | `archive/v1/path:line` where this was previously used |
| `validity` | approximation | where it holds |
| `notes` | no | free text |

Symbol conflicts: two notation entries giving one `sym` different meanings
fail the check unless both list each other in `conflicts_with` with a note
explaining which notation this project adopts (`adopted = true` on one).

Anchor rule: the anchor is the greppable proof that the source says this.
For scanned papers with poor OCR, anchor a nearby prose sentence and set
`ocr = "poor"`; the `latex` must then be read off the PDF page by a reviewer.

## The recomputation tree

`depends_on` edges form a DAG rooted at the rendered image. A node is **green**
when its `ownership` is `owned` or `leaf` and every node reachable below it is
green; otherwise **yellow** (a stub somewhere beneath) or **red** (ownership
unknown). `check.py --render` writes `truth/TREE.md` with one line per node and
the count of yellow leaves blocking each root. Breadth-first: stub first, then
turn leaves green from the top down.
