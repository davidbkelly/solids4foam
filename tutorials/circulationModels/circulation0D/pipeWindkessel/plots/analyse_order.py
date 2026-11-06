#!/usr/bin/env python3
"""Order of accuracy of every case of this tutorial, for plotOrderOfAccuracy.py.

    python3 plots/analyse_order.py

For each scheme (explicit, implicit) and coupling order (first to fourth) the
RMSE of the applied outlet pressure Po against the exact solution (see
../analyse.py) is taken for every time step that has run. The order is
p = ln(e_coarse/e_fine)/ln(dt_coarse/dt_fine) between neighbouring time steps,
and the headline order is the least-squares gradient of ln(RMSE) on ln(dt)
over every level, which is the line drawn on the figure.

Writes orderAccuracy.json and orderAccuracy.csv next to this script.
Standard library only.
"""
import json
import math
import sys
from pathlib import Path

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent
sys.path.insert(0, str(ROOT))
import analyse  # noqa: E402

SCHEMES = ("explicit", "implicit")
ORDERS = ("first", "second", "third", "fourth")
EXPECTED = {name: i + 1 for i, name in enumerate(ORDERS)}


def main(verbose=True):
    results = {}
    for scheme in SCHEMES:
        for order in ORDERS:
            key = f"{scheme.capitalize()}{order.capitalize()}Order"
            rows = []
            for case in (ROOT/scheme/f"{order}Order").glob("dt_*"):
                e = analyse.rmse(case)
                if e:
                    rows.append({"dt": float(case.name[3:]), "rmse": e})
            rows.sort(key=lambda r: -r["dt"])
            for i in range(1, len(rows)):
                r = math.log(rows[i-1]["dt"]/rows[i]["dt"])
                rows[i]["order"] = math.log(rows[i-1]["rmse"]/rows[i]["rmse"])/r
            results[key] = {
                "levels": rows,
                "fittedOrder": (analyse.fit([(r["dt"], r["rmse"]) for r in rows])
                                if len(rows) > 1 else float("nan")),
                "expectedOrder": EXPECTED[order],
            }

    (OUT/"orderAccuracy.json").write_text(json.dumps(results, indent=2))
    with open(OUT/"orderAccuracy.csv", "w") as f:
        f.write("algorithm,deltaT,rmse,pairwiseOrder,fittedOrder,expectedOrder\n")
        for a, d in results.items():
            for r in d["levels"]:
                f.write(f"{a},{r['dt']:g},{r['rmse']:.12g},"
                        f"{r.get('order', float('nan')):.6g},"
                        f"{d['fittedOrder']:.6g},{d['expectedOrder']}\n")

    if verbose:
        print(f"{'algorithm':<22}{'dt':>9}{'RMSE (Pa)':>13}{'pairwise p':>12}")
        for a, d in results.items():
            for r in d["levels"]:
                o = r.get("order")
                print(f"{a:<22}{r['dt']:>9g}{r['rmse']:>13.4e}"
                      f"{('' if o is None else f'{o:.3f}'):>12}")
            print(f"{'':<22}{'fitted p =':>9} {d['fittedOrder']:.3f}"
                  f"  ({len(d['levels'])} levels)   theoretical {d['expectedOrder']}\n")
    return results


if __name__ == "__main__":
    main()
