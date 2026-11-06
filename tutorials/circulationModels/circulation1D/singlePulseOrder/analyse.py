#!/usr/bin/env python3
"""Order of accuracy of the inviscid single pulse against the linear solution.

    python3 analyse.py

For every run <set>/n<N> that has finished, the error of P along the vessel at
t = 0.8 s against

    P(x, t) = Z0 s Qin(t - x/c0),   Qin(t) = 1e-6 exp(-1e4 (t - 0.05)^2),
    Z0 = rho c0/A0,

s the amplitude factor of the set (1e-3 linear, 1 full), in the L1, L2 and
max norms over the cell centres, relative to the peak Z0 s 1e-6. Prints, per
set, the errors, the orders between neighbouring meshes and the least-squares
order over the four finest meshes; writes orders.md and plots/orderOfAccuracy.png
(if matplotlib is available). Needs numpy.
"""
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))

import boileauMetrics as bm  # noqa: E402
from singlePulse import C0, Z0  # noqa: E402

T_END = 0.8
L = 10.0
SETS = {"linear": 1e-3, "full": 1.0}


def qin(t):
    return 1e-6*np.exp(-1e4*(t - 0.05)**2)


def errors(path, scale):
    p = np.loadtxt(path, comments="#")[1:-1]     # cell centres only
    x, P = p[:, 0], p[:, 3]
    exact = Z0*scale*qin(T_END - x/C0)
    e = np.abs(P - exact)/(Z0*scale*1e-6)
    return np.mean(e), math.sqrt(np.mean(e**2)), e.max()


def main():
    results = {}
    for s, scale in SETS.items():
        rows = []
        for d in sorted((HERE / s).glob("n*"), key=lambda p: int(p.name[1:])):
            f = d / "postProcessing" / "circulation1D" / "0" / "profiles" / f"{T_END:g}" / "tube.dat"
            if f.is_file():
                n = int(d.name[1:])
                rows.append((n, L/n) + errors(f, scale))
        if rows:
            results[s] = rows

    if not results:
        sys.exit("No finished runs: run ./Allrun first")

    out = [f"# Single pulse, inviscid: order of accuracy at t = {T_END} s\n",
           "Errors of P against the linear solution, relative to the peak; MUSCL "
           "(limiter none), HLL, SSPRK2, CFL 0.5.\n"]
    for s, rows in results.items():
        out.append(f"## {s} (inflow x {SETS[s]:g})\n")
        out.append("| nCells | dx [mm] | L1 | L2 | Linf | order L1 | order L2 | order Linf |")
        out.append("|---|---|---|---|---|---|---|---|")
        for i, r in enumerate(rows):
            if i:
                q = rows[i - 1]
                o = [math.log(q[k]/r[k])/math.log(q[1]/r[1]) for k in (2, 3, 4)]
                ostr = " | ".join(f"{v:.3f}" for v in o)
            else:
                ostr = "- | - | -"
            out.append(f"| {r[0]} | {1e3*r[1]:.4g} | {r[2]:.3e} | {r[3]:.3e} | "
                       f"{r[4]:.3e} | {ostr} |")
        fine = rows[-4:]
        if len(fine) >= 2:
            lx = np.log([r[1] for r in fine])
            fits = [np.polyfit(lx, np.log([r[k] for r in fine]), 1)[0] for k in (2, 3, 4)]
            out.append(f"\nLeast-squares order over the {len(fine)} finest meshes: "
                       f"L1 {fits[0]:.3f}, L2 {fits[1]:.3f}, Linf {fits[2]:.3f}\n")
    text = "\n".join(out) + "\n"
    (HERE / "orders.md").write_text(text)
    print(text)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return

    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    colours = {"linear": bm.COLOURS["DCG"], "full": bm.COLOURS["LCG"]}
    for s, rows in results.items():
        dx = np.array([r[1] for r in rows])
        ax.loglog(1e3*dx, [r[3] for r in rows], "o-", color=colours[s],
                  linewidth=2, markersize=6, label=f"{s} (inflow x {SETS[s]:g}), L2")
    rows = results.get("linear", next(iter(results.values())))
    dx = np.array([r[1] for r in rows])
    ref = rows[-1][3]*(dx/dx[-1])**2
    ax.loglog(1e3*dx, ref, color=bm.INK_2, linewidth=1, linestyle=(0, (4, 2)),
              label="slope 2")
    ax.set_xlabel("dx [mm]")
    ax.set_ylabel("relative L2 error of P")
    ax.set_title("Single pulse: error against the linear solution at t = 0.8 s",
                 fontsize=10, color=bm.INK)
    bm.style(ax)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    (HERE / "plots").mkdir(exist_ok=True)
    fig.savefig(HERE / "plots" / "orderOfAccuracy.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
