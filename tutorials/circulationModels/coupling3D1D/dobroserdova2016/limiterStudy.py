#!/usr/bin/env python3
"""The step inflow in the resonator: limiters and the unlimited high-order
schemes against the exact delay-equation solution.

    python3 limiterStudy.py [--cells 100,400]

The 1D-only hard-coupling case (the pipe as an inertance node, c0 = 3.5 m/s,
step inflow 1e-7 m3/s for 0.1 s) on the linearised equations, of which the
delay equations are the exact solution: the exact flow has jumps (the
step passes to the interface through the inlet's +1 reflection every
round trip) and the exact pressure L dQ/dt is impulsive. Run with the
MUSCL reconstruction and each limiter (minmod, vanLeer, MC, none) with
SSPRK2, and with the unlimited k-exact reconstructions of degree 2
(SSPRK3) and 3 (RK4), at CFL 0.5, on 100 and 400 cells. Reported, for the
interface flow over 0 to 0.5 s against the exact solution: the overshoot
(max Q minus the exact max, relative to the inflow amplitude), the total
variation relative to the exact one, the RMS and max error. Writes
limiters.md and plots/limiters.png.
"""
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import makeCases as mc  # noqa: E402
import analyse as an  # noqa: E402

RUNS = HERE/"runs"
SCHEMES = [
    ("minmod", "MUSCL", None, "minmod", "SSPRK2"),
    ("vanLeer", "MUSCL", None, "vanLeer", "SSPRK2"),
    ("MC", "MUSCL", None, "MC", "SSPRK2"),
    ("none", "MUSCL", None, "none", "SSPRK2"),
    ("kExact2", "kExact", 2, None, "SSPRK3"),
    ("kExact3", "kExact", 3, None, "RK4"),
]
END = 0.5


def edit(s, pattern, value):
    s, n = re.subn(pattern, value, s, count=1, flags=re.M)
    assert n == 1, pattern
    return s


def setup(label, rec, k, limiter, ts, n):
    case = RUNS/f"step_limiter_{label}_n{n}"
    if (case/"postProcessing").is_dir():
        return case
    if case.exists():
        shutil.rmtree(case)
    (case/"constant").mkdir(parents=True)
    (case/"system").mkdir()
    p = mc.properties("hard", 3.5, "step")
    p = edit(p, r"^(\s+nCells\s+)\d+;", rf"\g<1>{n};")
    p = edit(p, r"^(\s+initialPressure\s+0;)",
             "\\1\n        linearised          true;")
    p = edit(p, r"^(\s+CFL\s+)\S+;", r"\g<1>0.5;")
    if rec == "kExact":
        p = edit(p, r"^(\s+timeScheme\s+)\S+;.*$",
                 f"    timeScheme      {ts};\n    reconstruction  kExact;\n"
                 f"    polynomialOrder {k};")
    else:
        p = edit(p, r"^(\s+limiter\s+)\S+;", rf"\g<1>{limiter};")
        p = edit(p, r"^(\s+timeScheme\s+)\S+;.*$", rf"\g<1>{ts};")
    (case/"constant"/"circulation1DProperties").write_text(p)
    c = mc.control_dict("hard", "step")
    c = edit(c, r"^endTime\s+\S+;", f"endTime         {END};")
    c = edit(c, r"^deltaT\s+\S+;.*$", "deltaT          1e-4;")
    (case/"system"/"controlDict").write_text(c)
    with open(case/"log", "w") as f:
        r = subprocess.run(["circulation1DFoam"], cwd=case, stdout=f,
                           stderr=subprocess.STDOUT)
    if r.returncode:
        sys.exit(f"{case.name} failed: see {case/'log'}")
    return case


def tv(x):
    return np.abs(np.diff(x)).sum()


def main():
    cells = [int(x) for x in (sys.argv[sys.argv.index("--cells") + 1]
                              if "--cells" in sys.argv else "100,400").split(",")]
    rows = []
    series = {}
    for n in cells:
        for label, rec, k, limiter, ts in SCHEMES:
            case = setup(label, rec, k, limiter, ts, n)
            a = np.loadtxt(case/"postProcessing"/"circulation1D"/"0"/"probes"
                           /"interface.dat", comments="#")
            t, P, Q = a[:, 0], a[:, 1], a[:, 2]
            Qx, Px = an.exact("step", 3.5, "hard", t, n2=64000)
            q0 = 1e-7
            rows.append((n, label, (Q.max() - Qx.max())/q0, tv(Q)/tv(Qx),
                         math.sqrt(np.mean((Q - Qx)**2))/q0,
                         np.abs(Q - Qx).max()/q0))
            series[(n, label)] = (t, Q, Qx)
    lines = ["# The step inflow in the resonator: limiters and the unlimited "
             "schemes", "",
             f"1D only, linearised, hard coupling (inertance to ground), c0 = "
             f"3.5 m/s, step 1e-7 m3/s for 0.1 s, 0-{END} s, CFL 0.5; the "
             "interface flow against the exact delay-equation solution (whose "
             "total variation over the window is "
             f"{tv(series[(cells[0], 'minmod')][2])/1e-7:.2f} x the inflow).", "",
             "| cells | scheme | overshoot / inflow | TV / TV exact | RMS error "
             "/ inflow | max error / inflow |",
             "|---|---|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r[0]} | {r[1]} | {r[2]:+.3f} | {r[3]:.3f} | {r[4]:.2e} "
                     f"| {r[5]:.2e} |")
    (HERE/"limiters.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, axs = plt.subplots(1, len(cells), figsize=(5.5*len(cells), 3.8),
                            squeeze=False)
    for ax, n in zip(axs[0], cells):
        t, Q, Qx = series[(n, "minmod")]
        ax.plot(t, Qx*1e7, "k", lw=1.8, label="exact")
        for label, *_ in SCHEMES:
            t, Q, _ = series[(n, label)]
            ax.plot(t, Q*1e7, lw=0.9, label=label)
        ax.set_xlim(0.09, 0.2)
        ax.set_xlabel("t [s]")
        ax.set_title(f"{n} cells (zoom on the first reflections)", fontsize=9)
        ax.grid(alpha=0.3)
    axs[0, 0].set_ylabel("Q / inflow")
    axs[0, 0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    (HERE/"plots").mkdir(exist_ok=True)
    fig.savefig(HERE/"plots"/"limiters.png", dpi=150)


if __name__ == "__main__":
    main()
