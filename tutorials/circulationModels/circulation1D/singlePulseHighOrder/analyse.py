#!/usr/bin/env python3
"""Order of accuracy of the higher-order 1D scheme, by self-convergence.

    python3 analyse.py

For every set (see Allrun) the profile at t = 0.8 s of each run is compared
with the finest run of the set (32000 cells), averaged onto the coarser
cells, which is exact for cell averages. Errors of the flow Q and of the area
perturbation A - A0, relative to their peaks, in the L1, L2 and max norms;
the order between successive meshes and by least squares over the three
finest. Writes orders.md and plots/highOrder.png. Needs numpy (matplotlib for
the figure).
"""
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SETS = {
    "MUSCL_SSPRK2": 2, "kExact1_SSPRK2": 2, "kExact2_SSPRK3": 3,
    "kExact3_RK4": 4, "kExact3_RK4_visc": 4,
}
REF = 32000
T = "0.8"


def profile(d):
    f = d/"postProcessing"/"circulation1D"/"0"/"profiles"/T/"tube.dat"
    if not f.is_file():
        return None
    a = np.loadtxt(f, comments="#")[1:-1]          # cell averages only
    return a[:, 1], a[:, 2]                        # A, Q


def main():
    lines = ["# Higher-order 1D scheme: order of accuracy (self-convergence)",
             "", f"Single pulse, full amplitude, t = {T} s, CFL 0.5; each "
             f"run against the {REF}-cell run of its set.", ""]
    results = {}
    for name, expect in SETS.items():
        ref = profile(HERE/name/f"n{REF}")
        if ref is None:
            continue
        Ar, Qr = ref
        A0 = Ar.min()
        rows = []
        for d in sorted((HERE/name).glob("n*"), key=lambda p: int(p.name[1:])):
            n = int(d.name[1:])
            if n == REF:
                continue
            p = profile(d)
            if p is None:
                continue
            r = REF//n
            Qc = Qr.reshape(n, r).mean(axis=1)
            Ac = Ar.reshape(n, r).mean(axis=1)
            eQ = np.abs(p[1] - Qc)/np.abs(Qr).max()
            eA = np.abs(p[0] - Ac)/(Ar.max() - A0)
            rows.append((n, eQ.mean(), math.sqrt((eQ**2).mean()), eQ.max(),
                         eA.max()))
        results[name] = rows

        lines += [f"## {name} (expected {expect})", "",
                  "| cells | Q L1 | order | Q L2 | order | Q max | order | "
                  "A max | order |", "|---|---|---|---|---|---|---|---|---|"]
        for i, r in enumerate(rows):
            cells = [str(r[0])]
            for k in (1, 2, 3, 4):
                cells.append(f"{r[k]:.3e}")
                cells.append("" if i == 0 else
                             f"{math.log(rows[i-1][k]/r[k])/math.log(2):.2f}")
            lines.append("| " + " | ".join(cells) + " |")
        fine = rows[-3:]
        if len(fine) == 3:
            h = np.log([1/r[0] for r in fine])
            fits = [np.polyfit(h, np.log([r[k] for r in fine]), 1)[0]
                    for k in (1, 2, 3, 4)]
            lines += ["", "Least-squares order over the three finest meshes: "
                      "Q L1 %.2f, L2 %.2f, max %.2f; A max %.2f." % tuple(fits)]
        lines.append("")

    if not results:
        sys.exit("No finished runs: run ./Allrun first")
    (HERE/"orders.md").write_text("\n".join(lines))
    print("\n".join(lines))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.ticker
    except ImportError:
        return
    (HERE/"plots").mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(5.5, 4))
    for name, rows in results.items():
        n = np.array([r[0] for r in rows])
        ax.loglog(10/n*1e3, [r[2] for r in rows], "o-", ms=4, label=name)
    # Slope guides through the finest point of a set of that order
    for q, name, ls in ((2, "MUSCL_SSPRK2", ":"), (3, "kExact2_SSPRK3", "--"),
                        (4, "kExact3_RK4", "-.")):
        if name in results:
            rows = results[name]
            dx = 10/np.array([rows[-3][0], rows[-1][0]])*1e3
            ax.loglog(dx, rows[-1][2]*(dx/dx[-1])**q*0.6, "k", ls=ls, lw=0.8,
                      label=f"slope {q}")
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xticks([1.25, 2.5, 5, 10, 20, 40])
    ax.set_xticklabels(["1.25", "2.5", "5", "10", "20", "40"])
    ax.set_xlabel("dx [mm]")
    ax.set_ylabel("Q error, L2 / peak")
    ax.set_title("Single pulse: self-convergence")
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(HERE/"plots"/"highOrder.png", dpi=150)


if __name__ == "__main__":
    main()
