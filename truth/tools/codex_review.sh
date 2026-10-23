#!/usr/bin/env bash
# Adversarial review of one truth/entries/*.toml file by codex (GPT-6.1 sol, xhigh).
# Usage: truth/tools/codex_review.sh truth/entries/<domain>.toml truth/reviews/<domain>-r<N>.md [extra instructions]
# Read-only sandbox: codex may open refs/*.txt, refs/*.pdf (pdftotext -f/-l), archive/v1/**, truth/**.
set -euo pipefail
ENTRIES="$1"; OUT="$2"; EXTRA="${3:-}"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
PROMPT=$(cat <<PEOF
You are an adversarial scientific reviewer for a provenanced ground-truth database. Your job is to REFUSE anything you cannot confirm from the local sources. Default to reject when uncertain.

File under review: $ENTRIES  (TOML, one [[entry]] per claim). Schema: truth/SCHEMA.md.
Sources: refs/<source>.txt is pdftotext -layout output of refs/<source>.pdf. Many keys also have refs/<source>.tex (the arXiv TeX source, tree under refs/tex/<source>/). Per truth/SCHEMA.md an anchor is valid if it is a byte-identical substring of EITHER refs/<source>.tex OR refs/<source>.txt; a TeX anchor is preferred, not a defect. Verify TeX anchors against the .tex file. For a page image or cleaner text use e.g. 'pdftotext -f N -l N refs/<source>.pdf -' or 'pdftotext -layout -f N -l N ...'.
Mechanical check: run 'python3 -I truth/check.py' first and report its output.

For EVERY entry, verify independently and in this order:
1. ANCHOR: the anchor string is really in refs/<source>.txt at a place that corresponds to the locator (eq. number / page). An anchor from a different section or an unrelated sentence is a reject.
2. TRANSCRIPTION: the latex / statement / value faithfully reproduce what the source says at that locator: every symbol, sign, exponent, factor (2, pi, 4pi, c, h vs hbar), unit system (CGS vs SI), and the conditions under which the source states it. Open the source text around the anchor and compare symbol by symbol. If ocr = "poor", extract the page with pdftotext -f/-l and if still unreadable, say so and demand image_verified by a human/vision pass; do not accept on faith.
3. SEMANTICS: the statement does not over-claim (e.g. calls an approximation exact, drops a validity condition, or misattributes a result to a paper that only cites it).
4. NOTATION: symbols used agree with other entries in truth/entries/*.toml; flag silent collisions (same symbol, different meaning) and silent synonyms (different symbols, same quantity).
5. NEED: is the entry something the reboot in PRD.md actually needs (Milestones 1-3)? Flag clearly unneeded entries as 'drop' (not a correctness verdict).

Output format (markdown, nothing else), one line per entry, then a short list of cross-entry issues:
| id | verdict | reason |
verdict is one of: accept | fix | reject | drop. 'fix' means a precise correction you can state in the reason (quote the corrected latex/value/anchor). 'reject' means not supportable from local sources or wrong beyond a one-line fix; a reject MUST also state the SURVIVING WEAKER STATEMENT (what the local source does support, possibly 'nothing'), so the reject is a work order, not a demolition.
Be specific: cite line numbers in refs/<source>.txt. Do not pad. Do not accept an entry because it looks standard; accept only because the local source says it.
$EXTRA
PEOF
)
mkdir -p "$(dirname "$OUT")"
codex exec -m gpt-6.1-sol -c model_reasoning_effort="xhigh" -s read-only --skip-git-repo-check \
  -o "$OUT" "$PROMPT" >"$OUT.log" 2>&1 || { echo "codex exec failed; see $OUT.log" >&2; exit 1; }
echo "review written to $OUT"
