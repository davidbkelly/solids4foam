#!/usr/bin/env python3
"""Order of accuracy of the inviscid pulse in a tapered vessel.

    python3 analyse.py

The pressure along the vessel at t = 0.4 s of each run n<N> is compared with
that of the next finer run, averaged onto the coarser cells (pairs of cells).
With d_N the relative L2 norm of that difference, the order from three meshes
N, 2N, 4N is log2(d_N/d_2N). Writes orders.md. Needs numpy.
"""
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
T_END = 0.4


def main():
    runs = {}
    for d in HERE.glob("n*"):
        f = d/"postProcessing"/"circulation1D"/"0"/"profiles"/f"{T_END:g}"/"tube.dat"
        if f.is_file():
            runs[int(d.name[1:])] = np.loadtxt(f, comments="#")[1:-1, 3]
    ns = sorted(runs)
    pairs = [(a, b) for a, b in zip(ns[:-1], ns[1:]) if b == 2*a]
    if len(pairs) < 2:
        raise SystemExit("Need at least three successive meshes: run ./Allrun")

    peak = np.abs(runs[ns[-1]]).max()
    d = {}
    lines = ["# Tapered vessel: order of accuracy at t = 0.4 s\n",
             "| meshes | L2 difference / peak | order |", "|---|---|---|"]
    for a, b in pairs:
        fine = runs[b]
        diff = runs[a] - 0.5*(fine[0::2] + fine[1::2])
        d[a] = math.sqrt(np.mean(diff**2))/peak
    for a, b in pairs:
        order = f"{math.log2(d[a]/d[b]):.2f}" if b in d else "-"
        lines.append(f"| {a} - {b} | {d[a]:.3e} | {order} |")
    text = "\n".join(lines) + "\n"
    (HERE/"orders.md").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
