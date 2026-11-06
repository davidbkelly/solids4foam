#!/usr/bin/env python3
"""Temporal order of accuracy of the coupling, against the exact solution.

    python3 analyse.py

For every case that has run, the RMSE of the applied outlet pressure Po over
the steps 0 < t <= 1 against the closed form

    Pc(t) = A R/(1+(wRC)^2) [sin wt - wRC cos wt + wRC exp(-t/RC)]
    Po(t) = Pc(t) + Rch A sin wt     (A = 1, w = 2 pi, R = 10, Rch = 1, C = 0.01)

Po is read to 17 digits from the
outletPressure function object (p/rho at the outlet, times rho). It prints the
RMSE for each time step and, per scheme and order, the fitted order of
convergence (least squares of ln RMSE against ln dt) and the orders between
neighbouring time steps. Standard library only.
"""
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
R, RCH, C, A, W, RHO = 10.0, 1.0, 0.01, 1.0, 2 * math.pi, 1060.0


def exact_po(t):
    k = W * R * C
    pc = A * R / (1 + k * k) * (math.sin(W * t) - k * math.cos(W * t)
                                + k * math.exp(-t / (R * C)))
    return pc + RCH * A * math.sin(W * t)


def rmse(case):
    files = sorted((case / "postProcessing" / "outletPressure").glob("*/surfaceFieldValue.dat"))
    if not files:
        return None
    err = []
    for line in files[0].read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        t, p = (float(x) for x in line.split()[:2])
        if 0 < t <= 1 + 1e-9:
            err.append(p * RHO - exact_po(t))
    return math.sqrt(sum(e * e for e in err) / len(err)) if err else None


def fit(points):
    x = [math.log(dt) for dt, _ in points]
    y = [math.log(e) for _, e in points]
    mx, my = sum(x) / len(x), sum(y) / len(y)
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / sum((a - mx) ** 2 for a in x)


def main():
    for family in sorted(p for p in HERE.glob("*/*") if p.is_dir()):
        pts = []
        for case in family.glob("dt_*"):
            e = rmse(case)
            if e:
                pts.append((float(case.name[3:]), e))
        if not pts:
            continue
        pts.sort(reverse=True)
        print(f"{family.relative_to(HERE)}")
        for i, (dt, e) in enumerate(pts):
            local = (f"  order {math.log(pts[i-1][1] / e) / math.log(pts[i-1][0] / dt):.3f}"
                     if i else "")
            print(f"  dt = {dt:<8g} RMSE Po = {e:.3e} Pa{local}")
        if len(pts) > 1:
            print(f"  fitted order {fit(pts):.3f}")


if __name__ == "__main__":
    main()
