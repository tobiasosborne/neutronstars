#!/usr/bin/env python3
"""Validate truth/entries/*.toml and render NOTATION.md / INDEX.md.

Usage: python3 -I truth/check.py [--render] [--all]
  --render  write truth/NOTATION.md and truth/INDEX.md from blessed entries
  --all     render proposed entries too (marked)
Exit 1 on any error. Every error names the entry id and the field.
"""
import sys, tomllib, pathlib, re, collections

ROOT = pathlib.Path(__file__).resolve().parent.parent
ENTRIES = ROOT / "truth" / "entries"
REFS = ROOT / "refs"
KINDS = {"notation", "definition", "equation", "constant", "result", "approximation", "datum", "pipeline"}
DOMAINS = {"gr", "rt", "opacity", "eos", "surface", "colour", "const", "obs", "notation", "pipeline"}
STATUS = {"proposed", "blessed", "rejected"}
MIN_ANCHOR = 24

def load():
    out = []
    for f in sorted(ENTRIES.glob("*.toml")):
        with open(f, "rb") as fh:
            try:
                doc = tomllib.load(fh)
            except tomllib.TOMLDecodeError as e:
                print(f"ERROR {f.name}: TOML parse error: {e}"); sys.exit(1)
        for e in doc.get("entry", []):
            e["_file"] = f.name
            out.append(e)
    return out

_txt_cache = {}
def ref_text(key):
    """Concatenated bytes of refs/<key>.tex (arXiv source, preferred) and refs/<key>.txt (pdftotext)."""
    if key == "PRD":
        return (ROOT / "PRD.md").read_bytes()
    if key not in _txt_cache:
        parts = [ (REFS / f"{key}{ext}").read_bytes() for ext in (".tex", ".txt") if (REFS / f"{key}{ext}").exists() ]
        _txt_cache[key] = b"\n".join(parts) if parts else None
    return _txt_cache[key]

def check(entries):
    errs, warns = [], []
    ids = collections.Counter(e.get("id", "<missing>") for e in entries)
    for i, n in ids.items():
        if n > 1: errs.append(f"{i}: duplicate id ({n}x)")
    byid = {e.get("id"): e for e in entries}
    symtab = collections.defaultdict(list)
    for e in entries:
        i = e.get("id", f"<no id in {e['_file']}>")
        def need(f):
            if f not in e or e[f] in ("", [], None):
                errs.append(f"{i}: missing field '{f}'")
                return False
            return True
        for f in ("id", "kind", "domain", "name", "statement", "source", "locator", "anchor", "ocr", "status"):
            need(f)
        if e.get("kind") not in KINDS: errs.append(f"{i}: bad kind {e.get('kind')!r}")
        if e.get("domain") not in DOMAINS: errs.append(f"{i}: bad domain {e.get('domain')!r}")
        if e.get("status") not in STATUS: errs.append(f"{i}: bad status {e.get('status')!r}")
        if e.get("ocr") not in ("good", "poor"): errs.append(f"{i}: ocr must be good|poor")
        if i and not re.fullmatch(r"[a-z]+\.[a-z0-9-]+", str(i)): errs.append(f"{i}: id must be domain.kebab-name")
        k = e.get("kind")
        if k in ("equation", "constant", "definition") and not e.get("latex"):
            errs.append(f"{i}: {k} needs 'latex'")
        if k in ("constant", "datum") and "value" not in e:
            errs.append(f"{i}: {k} needs 'value'")
        if k == "approximation" and not e.get("validity"):
            errs.append(f"{i}: approximation needs 'validity'")
        if k == "notation":
            syms = e.get("symbols") or []
            if not syms: errs.append(f"{i}: notation needs 'symbols'")
            for s in syms:
                for f in ("sym", "meaning", "units"):
                    if f not in s: errs.append(f"{i}: symbol entry missing '{f}'")
                if "sym" in s: symtab[s["sym"]].append(i)
        # source + anchor
        src = e.get("source")
        if src:
            if src != "PRD" and not any((REFS / f"{src}{ext}").exists() for ext in (".pdf", ".html", ".txt", ".csv", ".tex")):
                errs.append(f"{i}: source refs/{src}.{{pdf,html,txt,csv}} not found")
            txt = ref_text(src)
            anc = e.get("anchor")
            if anc is not None:
                if len(anc.encode()) < MIN_ANCHOR:
                    errs.append(f"{i}: anchor shorter than {MIN_ANCHOR} bytes")
                elif txt is None:
                    errs.append(f"{i}: no refs/{src}.tex or .txt to check anchor against")
                elif anc.encode() not in txt:
                    errs.append(f"{i}: anchor not byte-identical in refs/{src}.tex/.txt: {anc[:60]!r}")
            at = e.get("anchor_tex")
            if at is not None:
                texp = REFS / f"{src}.tex"
                if len(at.encode()) < MIN_ANCHOR: errs.append(f"{i}: anchor_tex shorter than {MIN_ANCHOR} bytes")
                elif not texp.exists(): errs.append(f"{i}: anchor_tex given but refs/{src}.tex does not exist")
                elif at.encode() not in texp.read_bytes(): errs.append(f"{i}: anchor_tex not byte-identical in refs/{src}.tex: {at[:60]!r}")
        if e.get("ocr") == "poor" and e.get("status") == "blessed" and not e.get("image_verified"):
            errs.append(f"{i}: ocr=poor and blessed but image_verified is not true")
        # reviews gate
        if e.get("status") == "blessed":
            rv = e.get("reviews") or []
            ok_verdicts = ("accept", "fix") if k == "pipeline" else ("accept",)
            acc = [r for r in rv if r.get("verdict") in ok_verdicts and str(r.get("model", "")).startswith("gpt")]
            rej = [r for r in rv if r.get("verdict") == "reject"]
            if not acc: errs.append(f"{i}: blessed without an 'accept' review from a gpt model")
            if rej and not e.get("reject_resolved"):
                errs.append(f"{i}: blessed with an unresolved 'reject' review")
        PROV = {"measured", "convention", "derivable", "fit", "table", "code"}
        if k != "notation":
            pv = e.get("provenance")
            if pv is None:
                warns.append(f"{i}: no 'provenance' class yet (required before final blessing)")
            elif pv not in PROV:
                errs.append(f"{i}: bad provenance {pv!r}")
            elif pv in ("fit", "table", "code") and not (e.get("reimplement") or e.get("drop_reason")):
                errs.append(f"{i}: provenance={pv} needs 'reimplement' (physics route) or 'drop_reason'")
        OWN = {"owned", "stub", "leaf"}
        if k != "notation":
            ow = e.get("ownership")
            if ow is None:
                warns.append(f"{i}: no 'ownership' yet (owned|stub|leaf)")
            elif ow not in OWN:
                errs.append(f"{i}: bad ownership {ow!r}")
            elif ow == "owned" and not e.get("impl"):
                errs.append(f"{i}: ownership=owned needs 'impl' (src path)")
            elif ow == "leaf" and e.get("provenance") not in (None, "measured", "convention"):
                errs.append(f"{i}: ownership=leaf only for provenance measured|convention")
        for d in e.get("depends_on", []) or []:
            if d not in byid: errs.append(f"{i}: depends_on unknown id {d!r}")
            elif d == i: errs.append(f"{i}: depends_on itself")
    # symbol conflicts
    for sym, owners in symtab.items():
        if len(owners) > 1:
            es = [byid[o] for o in owners]
            meanings = set()
            for o in es:
                for s in o["symbols"]:
                    if s.get("sym") == sym: meanings.add(s.get("meaning", "").strip().lower())
            if len(meanings) > 1:
                ok = all(set(owners) - {o["id"]} <= set(o.get("conflicts_with", [])) for o in es) \
                     and sum(1 for o in es if o.get("adopted")) == 1
                (warns if ok else errs).append(f"symbol {sym!r} has conflicting meanings in {owners}" + ("" if ok else " (declare conflicts_with + one adopted=true)"))
    return errs, warns

def colour_tree(entries):
    """Return {id: (colour, n_stubs_below)} with colour in green|yellow|red. Detects cycles."""
    byid = {e["id"]: e for e in entries if e.get("id")}
    memo, onstack, cyc = {}, set(), []
    def own_colour(e):
        if e.get("kind") == "notation": return "green", 0
        ow = e.get("ownership")
        if ow in ("owned", "leaf"): return "green", 0
        if ow == "stub": return "yellow", 1
        return "red", 0
    def visit(i):
        if i in memo: return memo[i]
        if i in onstack: cyc.append(i); return ("red", 0)
        onstack.add(i)
        e = byid[i]; col, stubs = own_colour(e)
        for d in e.get("depends_on", []) or []:
            if d not in byid or byid[d].get("kind") == "notation": continue
            c2, s2 = visit(d); stubs += s2
            if c2 == "red" or col == "red": col = "red"
            elif c2 == "yellow": col = "yellow"
        onstack.discard(i); memo[i] = (col, stubs); return memo[i]
    for i in byid: visit(i)
    return memo, cyc

def render_tree(entries):
    memo, cyc = colour_tree(entries)
    byid = {e["id"]: e for e in entries if e.get("id")}
    parents = collections.defaultdict(set)
    for e in entries:
        for d in e.get("depends_on", []) or []: parents[d].add(e["id"])
    roots = [i for i in byid if not parents[i] and byid[i].get("kind") != "notation"]
    dot = {"green": "🟢", "yellow": "🟡", "red": "🔴"}
    lines = ["# Recomputation tree (rendered by truth/check.py; do not edit)", "",
             f"roots: {len(roots)}; nodes: {len(byid)}; " + ", ".join(f"{c}: {sum(1 for v in memo.values() if v[0]==c)}" for c in ("green","yellow","red")), ""]
    if cyc: lines += ["**CYCLES involving:** " + ", ".join(sorted(set(cyc))), ""]
    def walk(i, depth, seen):
        c, n = memo[i]; e = byid[i]
        lines.append("  " * depth + f"- {dot[c]} **{i}** [{e.get('ownership','?')}/{e.get('provenance','?')}] {e['name']}" + (f" — {n} stub(s) below" if n and depth == 0 else ""))
        if i in seen: lines.append("  " * (depth+1) + "- (see above)"); return
        seen.add(i)
        for d in e.get("depends_on", []) or []:
            if d in byid: walk(d, depth + 1, seen)
    for r in sorted(roots, key=lambda i: (memo[i][0] != "red", memo[i][0] != "yellow", i)):
        walk(r, 0, set()); lines.append("")
    (ROOT / "truth" / "TREE.md").write_text("\n".join(lines) + "\n")
    return len(roots)

def render(entries, include_proposed=False):
    sel = [e for e in entries if e["status"] == "blessed" or (include_proposed and e["status"] == "proposed")]
    tag = lambda e: "" if e["status"] == "blessed" else " *(proposed)*"
    lines = ["# Notation dictionary (rendered by truth/check.py; do not edit)", ""]
    lines += ["| symbol | meaning | units | entry | source |", "|---|---|---|---|---|"]
    rows = []
    for e in sel:
        if e["kind"] != "notation": continue
        for s in e["symbols"]:
            rows.append((s["sym"], s["meaning"], s["units"], e["id"] + tag(e), f"{e['source']} ({e['locator']})"))
    for r in sorted(rows, key=lambda r: r[0].lower()):
        lines.append("| " + " | ".join(str(x).replace("|", "\\|") for x in r) + " |")
    (ROOT / "truth" / "NOTATION.md").write_text("\n".join(lines) + "\n")
    idx = ["# Ground-truth index (rendered by truth/check.py; do not edit)", ""]
    for dom in sorted({e["domain"] for e in sel}):
        idx.append(f"## {dom}"); idx.append("")
        for e in sorted((e for e in sel if e["domain"] == dom), key=lambda e: e["id"]):
            idx.append(f"- **{e['id']}** [{e['kind']}]{tag(e)} — {e['name']}. {e['statement']}  ")
            if e.get("latex"): idx.append(f"  `{e['latex']}`  ")
            idx.append(f"  src: `{e['source']}` {e['locator']}; anchor: “{e['anchor'][:80]}”")
        idx.append("")
    (ROOT / "truth" / "INDEX.md").write_text("\n".join(idx) + "\n")
    return len(sel)

if __name__ == "__main__":
    entries = load()
    errs, warns = check(entries)
    for w in warns: print("WARN", w)
    for e in errs: print("ERROR", e)
    c = collections.Counter(e.get("status") for e in entries)
    print(f"{len(entries)} entries: {dict(c)}; {len(errs)} errors, {len(warns)} warnings")
    if "--render" in sys.argv:
        n = render(entries, "--all" in sys.argv)
        nr = render_tree(entries)
        print(f"rendered {n} entries to truth/NOTATION.md and truth/INDEX.md; {nr} tree roots to truth/TREE.md")
    sys.exit(1 if errs else 0)
