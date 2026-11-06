#!/usr/bin/env python3
"""Order of accuracy of the six coupling algorithms of this tutorial.

    python3 plots/analyse_order.py

For each algorithm (explicit and implicit coupling, order 2 to 4) and each of
the outlet pressure, the inlet pressure and the inlet/outlet flow rate, the RMSE
of each time step against the dt = 1e-5 run of the same algorithm, over the
output times common to every run of every algorithm in 1.6 <= t <= 2.4 s, the
final cardiac cycle. The order between neighbouring time steps is
p = ln(e_coarse/e_fine)/ln(dt_coarse/dt_fine); the headline order p_fit is the
least-squares gradient of ln(RMSE) on ln(dt) over every time step shown, which
is the line drawn on the figure.

The records are the per-step function objects of each case (system/controlDict):
the patch pressures times rho, and the patch flow, at 17 digits, from
postProcessing/<name>/*/surfaceFieldValue.dat.

Writes orderAccuracy.json and orderAccuracy.csv next to this script.
Standard library only.
"""
import json
import math
from pathlib import Path

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent
RHO, T0, T1, TICK = 1060.0, 1.6, 2.4, 1e-7
REF = "1e-5"

# (scheme, order) of the six algorithms, in the order of the figure
ALGOS = [(s, o) for o in ("second", "third", "fourth") for s in ("explicit", "implicit")]
# quantity -> (function object, scale, label)
QUANT = {"outletPressure": ("outletPressure", RHO, "Outlet Pressure"),
         "inletPressure": ("boundaryPressure", RHO, "Inlet Pressure"),
         "inletFlow": ("inletFlow", 1.0, "Inlet/Outlet Flow Rate")}


def dts_of(scheme, order):
    """The time steps of an algorithm other than the reference, coarse to fine."""
    found = [p.name[3:] for p in (ROOT/scheme/f"{order}Order").glob("dt_*")]
    return sorted((d for d in found if d != REF), key=float, reverse=True)


def load(scheme, order, dt, fo, scale):
    """tick -> value of a record in the window; the run must have reached T1."""
    case = ROOT/scheme/f"{order}Order"/f"dt_{dt}"
    files = sorted((case/"postProcessing"/fo).glob("*/surfaceFieldValue.dat"))
    if not files:
        raise SystemExit(f"No {fo} record in {case}: run the case first")
    lo, hi = round(T0/TICK), round(T1/TICK)
    out, last = {}, 0.0
    for line in files[0].read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        t, v = line.split()[:2]
        t = float(t)
        last = t
        k = round(t/TICK)
        if lo <= k <= hi:
            out[k] = float(v)*scale
    if last < T1 - 1e-6:
        raise SystemExit(f"{case} stopped at t = {last:g}: it did not reach {T1:g} s")
    return out


def fit(h, e):
    x = [math.log(v) for v in h]
    y = [math.log(v) for v in e]
    mx, my = sum(x)/len(x), sum(y)/len(y)
    p = sum((a - mx)*(b - my) for a, b in zip(x, y))/sum((a - mx)**2 for a in x)
    return p, my - p*mx


def main(verbose=True):
    results = {}
    for qname, (fo, scale, label) in QUANT.items():
        data = {(s, o, d): load(s, o, d, fo, scale)
                for s, o in ALGOS for d in dts_of(s, o) + [REF]}
        common = sorted(set.intersection(*(set(v) for v in data.values())))
        res = {"label": label, "commonTimes": len(common), "algorithms": {}}
        for s, o in ALGOS:
            ref = data[(s, o, REF)]
            levels = []
            for d in dts_of(s, o):
                a = data[(s, o, d)]
                levels.append({"dt": float(d), "rmse": math.sqrt(
                    sum((a[k] - ref[k])**2 for k in common)/len(common))})
            levels.sort(key=lambda r: r["dt"])                    # fine -> coarse
            for i in range(1, len(levels)):
                levels[i]["order"] = (math.log(levels[i]["rmse"]/levels[i-1]["rmse"])
                                      / math.log(levels[i]["dt"]/levels[i-1]["dt"]))
            p, lnC = fit([r["dt"] for r in levels], [r["rmse"] for r in levels])
            res["algorithms"][f"{s}{o.capitalize()}Order"] = {
                "scheme": s, "order": o, "expectedOrder": {"second": 2, "third": 3, "fourth": 4}[o],
                "levels": levels, "fittedOrder": p, "lnC": lnC}
        results[qname] = res

    (OUT/"orderAccuracy.json").write_text(json.dumps(results, indent=2))
    with open(OUT/"orderAccuracy.csv", "w") as f:
        f.write("quantity,algorithm,deltaT,rmse,pairwiseOrder,fittedOrder,expectedOrder\n")
        for q, res in results.items():
            for a, d in res["algorithms"].items():
                for r in d["levels"]:
                    f.write(f"{q},{a},{r['dt']:g},{r['rmse']:.12g},"
                            f"{r.get('order', float('nan')):.6g},{d['fittedOrder']:.6g},"
                            f"{d['expectedOrder']}\n")

    if verbose:
        for q, res in results.items():
            print(f"\n=== {res['label']}: {res['commonTimes']} common times in "
                  f"{T0:g} <= t <= {T1:g} ===")
            print(f"{'algorithm':<22}{'dt':>9}{'RMSE vs 1e-5':>14}{'pairwise p':>12}")
            for a, d in res["algorithms"].items():
                for r in d["levels"]:
                    o = r.get("order")
                    print(f"{a:<22}{r['dt']:>9g}{r['rmse']:>14.4e}"
                          f"{('' if o is None else f'{o:.3f}'):>12}")
                print(f"{'':<22}{'p_fit =':>9} {d['fittedOrder']:.3f}"
                      f"   (theoretical {d['expectedOrder']})")
    return results


if __name__ == "__main__":
    main()
