#!/usr/bin/env python3
"""Red tests for truth/check.py: prove the checker can FAIL.
Copies truth/ into a temp root (refs/ symlinked), applies one mutation per test,
and asserts a nonzero exit. Usage: python3 -I truth/tools/red_tests.py
"""
import subprocess, sys, tempfile, pathlib, shutil, os, re
ROOT = pathlib.Path(__file__).resolve().parents[2]

def run(root):
    r = subprocess.run([sys.executable, "-I", str(root / "truth" / "check.py")], capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr

def fresh():
    tmp = pathlib.Path(tempfile.mkdtemp(prefix="truth-red-"))
    shutil.copytree(ROOT / "truth", tmp / "truth", ignore=shutil.ignore_patterns("reviews", "harvest", "__pycache__"))
    os.symlink(ROOT / "refs", tmp / "refs")
    shutil.copy(ROOT / "PRD.md", tmp / "PRD.md")
    return tmp

BAD = ("TOML parse error", "Traceback", "FileNotFoundError")

def right_error(out, eid):
    """First ERROR line that is not a parse error/traceback/FileNotFoundError and names eid, else None."""
    if any(b in out for b in ("Traceback", "FileNotFoundError")): return None
    for l in out.splitlines():
        if l.startswith("ERROR") and not any(b in l for b in BAD) and eid in l: return l
    return None

def mutated_id(old, new):
    """Id of the entry block containing the first changed character (the entry actually mutated)."""
    i = next((k for k, (x, y) in enumerate(zip(old, new)) if x != y), min(len(old), len(new)))
    a = new.rfind("[[entry]]") if i >= len(old) else new.rfind("[[entry]]", 0, i + 1)  # pure append: the appended copy
    return re.search(r'id = "([^"]+)"', new[a:]).group(1)

def first_entry_file(tmp):
    files = sorted((tmp / "truth" / "entries").glob("*.toml"))
    if not files: print("no entries to mutate; red tests need at least one entry"); sys.exit(2)
    return files[0]

MUTATIONS = {
  "corrupt anchor":        lambda s: re.sub(r'(anchor = "[^"]{10})', r'\1ZZQX', s, count=1),
  "bless without review":  lambda s: re.sub(r'status = "(proposed|rejected)"', 'status = "blessed"', s, count=1),
  "duplicate id":          lambda s: s + "\n" + s[s.index("[[entry]]"):s.index("[[entry]]", s.index("[[entry]]") + 1)] if s.count("[[entry]]") > 1 else s + "\n" + s[s.index("[[entry]]"):],
  "fit without route":     lambda s: _set_key(s, "provenance", "fit"),
  "owned without impl":    lambda s: _set_key(s, "ownership", "owned"),
  "unknown kind":          lambda s: re.sub(r'kind = "[a-z]+"', 'kind = "bogus"', s, count=1),
  "self-dependency":       lambda s: _self_dep(s),
}

def _set_key(s, key, val):
    """Set `key = "val"` in the first entry: replace an existing top-level line, else insert after [[entry]]."""
    a = s.index("[[entry]]"); b = s.find("[[entry]]", a + 1); b = len(s) if b < 0 else b
    block = s[a:b]
    if re.search(rf'(?m)^{key} = .*$', block):
        block = re.sub(rf'(?m)^{key} = .*$', f'{key} = "{val}"', block, count=1)
    else:
        block = block.replace("[[entry]]\n", f'[[entry]]\n{key} = "{val}"\n', 1)
    return s[:a] + block + s[b:]

def _self_dep(s):
    """Point the first entry's depends_on at itself (replace an existing list, else insert one)."""
    m = re.search(r'id = "([a-z]+\.[a-z0-9-]+)"', s)
    if not m: return s
    i = m.group(1); nxt = s.find("[[entry]]", m.end()); nxt = len(s) if nxt < 0 else nxt
    block = s[m.end():nxt]
    if "depends_on" in block:
        block2 = re.sub(r'depends_on = \[[^\]]*\]', f'depends_on = ["{i}"]', block, count=1)
    else:
        block2 = f'\ndepends_on = ["{i}"]' + block
    return s[:m.end()] + block2 + s[nxt:]

def main():
    base = fresh(); code0, out0 = run(base); shutil.rmtree(base)
    if code0 != 0: print("baseline is not clean; red tests meaningless"); sys.exit(2)
    print(f"baseline exit {code0}: {out0.strip().splitlines()[-1]}")
    failed = []
    for name, mut in MUTATIONS.items():
        tmp = fresh(); f = first_entry_file(tmp); s = f.read_text(); m = mut(s)
        if m == s: print(f"SKIP  {name}: mutation did not apply"); shutil.rmtree(tmp); continue
        eid = mutated_id(s, m)
        f.write_text(m); code, out = run(tmp); shutil.rmtree(tmp)
        err = right_error(out, eid)
        ok = code != 0 and err is not None
        print(f"{'RED ' if ok else 'MISS'}  {name}: exit {code} — {(err or next((l for l in out.splitlines() if l.startswith('ERROR')), out.strip().splitlines()[-1]))[:140]}")
        if not ok: failed.append(name)
    if failed: print("checker FAILED to go red for:", failed); sys.exit(1)
    print("all mutations caught")

if __name__ == "__main__": main()
