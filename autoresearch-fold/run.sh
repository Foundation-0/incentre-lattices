#!/bin/bash
# One experiment: lock check -> search.py (3300 s cap) -> lock check -> verify.py -> lock check -> results.tsv line -> keep best.
# The cap stays above search.py's WALL_LIMIT (3000 s), so a slow run ends INCOMPLETE (exit 2) instead of being killed.
# The lock is a tripwire against accidental edits, not a sandbox: search.py can write anything in this folder.
cd "$(dirname "$0")" || exit 1
PY=.venv/bin/python
CAP=3300
lock_check() {  # verify.py must match verify.sha256 and, in a git checkout, neither file may differ from HEAD
  if ! shasum -a 256 -c verify.sha256 >/dev/null 2>&1 ||
     { git rev-parse -q --verify HEAD >/dev/null 2>&1 && ! git diff --quiet HEAD -- verify.py verify.sha256; }; then
    echo "$1: verify.py or verify.sha256 does not match the pin or git HEAD (they are locked; never edit them)" >&2
    exit 3
  fi
}
lock_check "REFUSING TO RUN"
mkdir -p out best
rm -f out/path.npz out/verify.json
if command -v timeout >/dev/null 2>&1; then TO="timeout $CAP"; elif command -v gtimeout >/dev/null 2>&1; then TO="gtimeout $CAP"; else TO="perl -e alarm($CAP);exec(@ARGV) --"; fi
$TO $PY search.py > out/search.log 2>&1
rc=$?
tail -n 3 out/search.log
lock_check "REFUSING TO JUDGE (changed while search.py ran)"
if [ -f out/path.npz ]; then
  $PY -I verify.py out/path.npz > out/verify.json 2>out/verify.err  # -I: neither this folder nor the cwd is on sys.path
else
  echo '{"valid": false, "reasons": ["no path written (search exit code '$rc')"]}' > out/verify.json
fi
lock_check "REFUSING TO LOG OR KEEP"
$PY -I - "$rc" <<'PYEOF'
import json, re, subprocess, sys, shutil, time, os, hashlib
rc = sys.argv[1]
rep = json.load(open("out/verify.json"))
src = open("search.py").read()
try:
    h = subprocess.check_output(["git", "hash-object", "search.py"], stderr=subprocess.DEVNULL).decode().strip()[:7]
except Exception:
    h = hashlib.sha1(b"blob %d\0" % len(src.encode()) + src.encode()).hexdigest()[:7]
mm = re.search(r'^DESC\s*=\s*(["\'])(.*?)\1\s*$', src, re.M)
desc = (mm.group(2) if mm else "no DESC").replace("\t", " ")
def f(k, fmt="%.6f"):
    v = rep.get(k)
    return fmt % v if isinstance(v, (int, float)) and not isinstance(v, bool) else "NA"
valid = bool(rep.get("valid"))
flex = rep.get("flex_bars_plus_contacts", "NA")
line = "\t".join([time.strftime("%Y-%m-%dT%H:%M:%S"), h, str(valid).lower(), f("final_density"), f("max_density"),
                  str(rep.get("contacts_per_ball_mean", "NA") if "contacts_per_ball_mean" not in rep else
                      "%.2f(min%d,max%d)" % (rep["contacts_per_ball_mean"], rep["contacts_per_ball_min"], rep["contacts_per_ball_max"])),
                  str(flex), "m=%s" % rep.get("m", "NA"), "exit=%s" % rc, desc])
hdr = "timestamp\tsearch_hash\tvalid\tfinal_density\tmax_density\tend_contacts\tflex_count\tm\texit\tdesc"
new = not os.path.exists("results.tsv")
with open("results.tsv", "a") as fh:
    if new:
        fh.write(hdr + "\n")
    fh.write(line + "\n")
print(line)
if not valid:
    print("INVALID:", "; ".join(rep.get("reasons", []))[:600])
else:
    prev = -1.0
    if os.path.exists("best/verify.json"):
        prev = json.load(open("best/verify.json")).get("max_density") or -1.0
    if rep["max_density"] > prev + 1e-9:
        os.makedirs("best", exist_ok=True)
        shutil.copy("out/path.npz", "best/path.npz")
        shutil.copy("search.py", "best/search.py")
        shutil.copy("out/verify.json", "best/verify.json")
        print("NEW BEST: max_density %.6f (was %.6f) -> best/" % (rep["max_density"], prev))
    else:
        print("not better than best (%.6f)" % prev)
PYEOF
