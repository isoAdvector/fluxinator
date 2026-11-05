"""Splits data/browser.json (browserData.py) into the four data files of the results browser
(dashboard artifact "Wave Solver Explorer"): browser_runs.json (one row per run, rounded
metrics), browser_series.json (time series), browser_snaps.json (profiles and surfaces at
10, 20, 30 s), browser_exact.json (exact surface and crest profile).
Usage: python3 browserFiles.py"""
import json
import os
from common import DATA

SHORT = {"interFoam": "interFoam", "interIsoFoam": "interIsoFoam", "interIsoGoam": "interIsoGoam",
         "explicit + LW + interfaceUpwind": "final"}
CORES = {760: 4}


def rnd(v, n):
    return None if v is None else round(float(v), n)


B = json.load(open(os.path.join(DATA, "browser.json")))
runs, series, snaps = [], [], []
for c in B["cases"]:
    s = SHORT[c["solver"]]
    key = f"{c['scan']}|{c['nx']}|{c['par']:g}|{s}"
    cores = CORES.get(c["nx"], 1)
    runs.append(dict(
        key=key, scan=c["scan"], mesh=f"{c['nx']}x{c['nz']}", nx=c["nx"], par=c["par"], solver=s,
        dH=rnd(c.get("dHfit"), 2), dc=rnd(c.get("dc"), 3), Umax=rnd(c.get("UmaxMed"), 3),
        Uair=rnd(c.get("UairMed"), 3), emax=rnd(c.get("emax"), 2), r0=rnd(c.get("r0"), 3),
        ulow=rnd(c.get("ulow"), 4), steps=c["steps"], clock=rnd(c.get("clock"), 0), cores=cores,
        msStep=rnd(None if c.get("clock") is None else c["clock"]/c["steps"]*1000, 2)))
    series.append(dict(key=key, tH=c["tH"], dH=c["dH"], tU=c["tU"], Umax=c["Umax"],
                       tAir=c.get("tair", []), Uair=c.get("Uair", [])))
    for t, v in c.get("snaps", {}).items():
        snaps.append(dict(key=key, t=int(t), z=v["z"], u=v["u"], ue=v["ue"], x=v["x"], h=v["h"]))
order = ["interFoam", "interIsoFoam", "interIsoGoam", "final", "centred"]
runs.sort(key=lambda r: (r["scan"], r["nx"], -r["par"], order.index(r["solver"])))
e = B["exact"]
files = dict(runs=runs, series=series, snaps=snaps,
             exact=[dict(x=e["x"], h=e["h"], z=e["z"], u=e["u"], H=e["H"], T=e["T"], d=e["d"],
                         L=round(e["L"], 3))])
for k, v in files.items():
    p = os.path.join(DATA, f"browser_{k}.json")
    json.dump(v, open(p, "w"), separators=(",", ":"))
    print(p, len(v), "rows", os.path.getsize(p)//1024, "kB")
