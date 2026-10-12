"""Download all COD entries for MB6 (Pm-3m) and MB12 (Fm-3m) and extract the boron free parameter x."""
import json, os, re, sys, time, ssl, certifi, urllib.request, urllib.parse
CTX = ssl.create_default_context(cafile=certifi.where())
METALS = "Li Na K Rb Cs Mg Ca Sr Ba Sc Y La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Th U Np Pu Zr Hf Ti Am".split()
os.makedirs("cod", exist_ok=True)
def get(url):
    for _ in range(3):
        try:
            return urllib.request.urlopen(url, timeout=40, context=CTX).read().decode("utf-8", "replace")
        except Exception as e:
            time.sleep(2)
    return ""
out = []
for M in METALS:
    for nB, sg in ((6, 221), (12, 225)):
        f = f"B{nB} {M}"
        js = get("https://www.crystallography.net/cod/result?" + urllib.parse.urlencode({"formula": f, "format": "json"}))
        try:
            rows = json.loads(js) if js.strip() else []
        except Exception:
            rows = []
        for r in rows:
            if str(r.get("sgNumber")) != str(sg):
                continue
            cid = r["file"]; path = f"cod/{cid}.cif"
            if not os.path.exists(path):
                txt = get(f"https://www.crystallography.net/cod/{cid}.cif")
                open(path, "w").write(txt)
            out.append(dict(M=M, n=nB, cod=cid, a=float(r["a"]), T=r.get("celltemp") or r.get("diffrtemp"),
                            P=r.get("cellpressure") or r.get("diffrpressure")))
        print(M, nB, len([o for o in out if o["M"] == M and o["n"] == nB]), file=sys.stderr)
json.dump(out, open("cod_index.json", "w"), indent=1)
