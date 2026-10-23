#!/usr/bin/env bash
# Analytic-argument run by codex (gpt-6-astra, xhigh) on one research-grade critical path.
# Usage: truth/tools/codex_analysis.sh truth/analysis/briefs/<topic>.md truth/analysis/<topic>.md
# Read-only sandbox: codex may open refs/*.txt, refs/*.tex, refs/*.pdf (pdftotext -f/-l), archive/v1/**, truth/**.
set -euo pipefail
BRIEF="$1"; OUT="$2"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"; cd "$ROOT"
PROMPT="$(cat truth/analysis/briefs/_preamble.md; echo; cat "$BRIEF")"
mkdir -p "$(dirname "$OUT")"
codex exec -m gpt-6-astra -c model_reasoning_effort="xhigh" -s read-only --skip-git-repo-check \
  -o "$OUT" "$PROMPT" </dev/null >"$OUT.log" 2>&1 || { echo "codex exec failed; see $OUT.log" >&2; exit 1; }
echo "analysis written to $OUT"
