#!/usr/bin/env python3
"""The Womersley segment at physiological amplitude: a numerical reference.

    python3 physiological.py [run] [--jobs 2]

The alpha = 6 case with an incident wave of 5 kPa (a tenth of the pulse
pressure of the aorta) on the nonlinear 1D equations, over the first second
from rest, with the fluid at BDF 4 and the coupling at order 4, dt = 0.25 ms,
the vessels with kExact 3 + RK4 at CFL 0.5 on 800 cells, on the tutorial's
50 x 100 wedge and on 100 x 200. There is no closed form (the 1D equations
are nonlinear at this amplitude: 5 kPa is 0.12 of rho c0^2); the finer run
is the reference, labelled numerical, and the coarser one bounds its
discretisation error. Writes reference/alpha6_physiological_*.csv from the
finer run and physiological.md.
"""
import os
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import makeCases as mc  # noqa: E402
import spaceStudy as ss  # noqa: E402

P_INC = 5000.0
DT = 2.5e-4
END = 1.0
LEVELS = [1, 2]


def name(level):
    return f"physiological_alpha6_x{level:g}"


def setup(level):
    case = ss.RUNS/name(level)
    if ss.finished(case):
        return case
    mc.make(6.0, DT, True, name=case.name, end=END, pInc=P_INC)
    nr, nz = int(round(50*level)), int(round(100*level))
    ss.edit(case/"system"/"blockMeshDict",
            [(r"\(50 1 100\) simpleGrading \(0\.1 1 1\)",
              f"({nr} 1 {nz}) simpleGrading (0.1 1 1)")])
    ss.edit(case/"system"/"fvSchemes",
            [(r"^(\s+default\s+)backward;.*$", r"\g<1>BDF 4;")])
    ss.edit(case/"system"/"controlDict",
            [(r'^libs\s+\("libcirculation1DCoupling.so"\);',
              'libs            ("libcirculation1DCoupling.so" '
              '"libsolids4FoamModels.so");     // BDF ddt scheme'),
             (rf"timeStart\s+{END - 1/mc.F:g};",
              f"timeStart       {END - 0.5/mc.F:g};"),
             (rf"writeInterval\s+{1/(8*mc.F):g};",
              f"writeInterval   {0.5/(8*mc.F):g};")])
    ss.edit(case/"constant"/"circulation1DProperties",
            [(r"^(\s+order\s+)\S+;", r"\g<1>4;"),
             (r"^(\s+CFL\s+)\S+;", r"\g<1>0.5;"),
             (r"^(\s+nCells\s+)\d+;", r"\g<1>800;"),
             (r"^(\s+timeScheme\s+)\S+;.*$",
              "    timeScheme      RK4;\n    reconstruction  kExact;\n"
              "    polynomialOrder 3;")])
    s = (case/"system"/"fvSolution").read_text()
    i = s.index("residualControl")
    s = s[:i] + re.sub(r"tolerance\s+[\d.e+-]+;", "tolerance 1e-10;", s[i:],
                       count=1)
    (case/"system"/"fvSolution").write_text(s)
    return case


def analyse():
    cases = {lv: ss.RUNS/name(lv) for lv in LEVELS}
    if not all(ss.finished(c) for c in cases.values()):
        sys.exit("not all runs finished")
    t1, p1 = ss.ports(cases[1])
    t2, p2 = ss.ports(cases[2])
    lines = ["# Womersley segment at physiological amplitude (numerical reference)",
             "", f"alpha = 6, incident wave {P_INC:g} Pa, 0-{END} s from rest, "
             f"BDF 4 / coupling order 4 at dt = {DT*1e3:g} ms, nonlinear "
             "vessels (kExact 3, RK4, 800 cells). The 100 x 200 run is the "
             "reference; the 50 x 100 run against it (RMS / RMS of the "
             "signal):", ""]
    for k in ("to3D.P", "to3D.Q", "from3D.P"):
        e = np.sqrt(np.mean((p1[k] - p2[k])**2))/np.sqrt(np.mean(p2[k]**2))
        lines.append(f"- {k}: {e:.2e}")
    pr1, pr2 = ss.profiles(cases[1]), ss.profiles(cases[2])
    times = sorted(pr2)
    umax = max(np.abs(pr2[tk][1]).max() for tk in times)
    e = max(np.abs(np.interp(pr2[tk][0], pr1[tk][0], pr1[tk][1]) - pr2[tk][1]).max()
            for tk in times)/umax
    lines.append(f"- mid-length profile, max over {len(times)} instants: {e:.2e}")
    lines += ["", "Against the linear exact solution scaled to 5 kPa (what the "
              "linear model would give), RMS / RMS:", ""]
    import analyse as an
    ex = an.exact_transient(6.0, t2)
    for k in ("to3D.P", "to3D.Q", "from3D.P"):
        x = ex[k]*P_INC/mc.P_INC
        lines.append(f"- {k}: {np.sqrt(np.mean((p2[k] - x)**2))/np.sqrt(np.mean(x**2)):.2e}")
    (HERE/"physiological.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    outdir = HERE/"reference"
    os.makedirs(outdir, exist_ok=True)
    hdr = (f"# Womersley segment at physiological amplitude: NUMERICAL reference (alpha 6, incident "
           f"{P_INC:g} sin(2 pi t) min(1, t) Pa, nonlinear 1D sqrtElastic vessels; 100 x 200 wedge, "
           f"BDF 4 and coupling order 4 at dt = {DT*1e3:g} ms, kExact 3 + RK4 vessels on 800 cells); "
           "the 50 x 100 run differs by the amounts in physiological.md\n")
    with open(outdir/"alpha6_physiological_ports.csv", "w") as f:
        f.write(hdr + "t [s],P_in [Pa],Q [m3/s],P_out [Pa]\n")
        for row in zip(t2, p2["to3D.P"], p2["to3D.Q"], p2["from3D.P"]):
            f.write(",".join(f"{v:.12g}" for v in row) + "\n")
    r = pr2[times[0]][0]
    with open(outdir/"alpha6_physiological_profiles.csv", "w") as f:
        f.write(hdr + "# u(r, t) [m/s] at the mid-length; r [m] then t = "
                + ", ".join(f"{tk:g}" for tk in times) + " s\n")
        for j in range(len(r)):
            f.write(f"{r[j]:.12g}," + ",".join(f"{pr2[tk][1][j]:.12g}" for tk in times) + "\n")


def main():
    if sys.argv[1:2] == ["run"]:
        from concurrent.futures import ThreadPoolExecutor
        jobs = int(sys.argv[sys.argv.index("--jobs") + 1]) if "--jobs" in sys.argv else 2
        cases = [setup(lv) for lv in LEVELS]
        with ThreadPoolExecutor(jobs) as ex:
            list(ex.map(ss.run, cases))
    analyse()


if __name__ == "__main__":
    main()
