#!/usr/bin/env python3
"""Fluid at rest: how far the network moves from rest.

    python3 analyse.py

For each variant that has run, over every profile written (t = 0, 0.5, ..., 2 s)
and every vessel: the largest |Q|, relative to the flow scale A c of the widest
cell (c = 4 m/s), and the largest |P - P_rest|, relative to P_rest. A
well-balanced scheme keeps both at round-off. Writes rest.md. Needs numpy.
"""
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
VARIANTS = {"atDiastolicPressure": 10000.0, "atLowerPressure": 4000.0}


def main():
    lines = ["# Fluid at rest\n",
             "| Variant | P_rest [Pa] | profiles | max abs Q [m3/s] | max abs Q / (A c) "
             "| max abs (P - P_rest) [Pa] | relative |",
             "|---|---|---|---|---|---|---|"]
    for v, prest in VARIANTS.items():
        prof = HERE/v/"postProcessing"/"circulation1D"/"0"/"profiles"
        if not prof.is_dir():
            continue
        qmax = dpmax = 0.0
        scale = 0.0
        n = 0
        for f in sorted(prof.glob("*/*.dat")):
            p = np.loadtxt(f, comments="#")
            qmax = max(qmax, np.abs(p[:, 2]).max())
            dpmax = max(dpmax, np.abs(p[:, 3] - prest).max())
            # Flow scale A c, with c = 4 m/s, the wave speed of these vessels
            scale = max(scale, p[:, 1].max()*4.0)
            n += 1
        lines.append(f"| {v} | {prest:g} | {n} | {qmax:.2e} | {qmax/scale:.2e} | "
                     f"{dpmax:.2e} | {dpmax/prest:.2e} |")
    text = "\n".join(lines) + "\n"
    (HERE/"rest.md").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
