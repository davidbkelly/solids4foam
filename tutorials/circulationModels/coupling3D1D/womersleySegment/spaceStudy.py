#!/usr/bin/env python3
"""Order in space of the Womersley segment against the exact solution.

    python3 spaceStudy.py [run] [--alpha 6] [--levels 1,2,4] [--jobs 2]
    python3 spaceStudy.py            # analyse what has run

The wedge mesh (50 radial cells graded 10:1 to the wall, 100 axial) is
coarsened and refined, level k = nr/50 = nz/100, and the coupled case is run
from rest over the first second of the ramped incident wave with the
vessels linearised at fourth order (kExact 3, RK4), the fluid at BDF 4 with
the coupling at order 4 and dt = 0.5 ms (its time error is far below the
spatial errors measured here: the BDF 4 family's self-convergence at this
step is below 1e-7), so that what is left against the exact transient
solution (analyse.exact_transient, the Womersley transfer function
convolved with the incident wave) is the 3D discretisation error.

Measured: the RMS error over the window of the port pressures and the
flow, and the error of the axial velocity profile at the pipe's
mid-length at eight instants over the second half of the window (the sets
function object), max over the instants of max |u - u_exact| / max |u|,
both against the exact transient; and the observed order between
successive levels. Writes space.md and plots/space.png.
"""
import argparse
import math
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import makeCases as mc  # noqa: E402
import analyse as an  # noqa: E402

RUNS = HERE/"runs"
END = 1.0
DT = 5e-4
LEVELS = [0.5, 1, 2]
FAMILY = "graded"      # --family: graded (the tutorial's 10:1 wall grading,
                       # the wedge of flat wall) or uniform (no grading, the
                       # wedge scaled by sqrt(theta/sin theta) so that its
                       # cross-section is pi R^2 and the exact solution is
                       # that of the pipe of radius R, independent of the mesh)


def name(alpha, level, dt=None):
    return (f"space_alpha{alpha:g}" + ("_u" if FAMILY == "uniform" else "")
            + f"_x{level:g}" + (f"_dt{dt:g}" if dt and dt != DT else ""))


def finished(case):
    log = case/"log.pimpleFoam"
    return log.is_file() and "\nEnd\n" in log.read_text()[-2000:]


def edit(path, pairs):
    s = path.read_text()
    for pattern, value in pairs:
        s, n = re.subn(pattern, value, s, count=1, flags=re.M)
        assert n == 1, (path, pattern)
    path.write_text(s)


def setup(alpha, level, force=False, dt=None):
    dt = dt or DT
    case = RUNS/name(alpha, level, dt)
    if finished(case) and not force:
        return case
    mc.make(alpha, dt, True, name=case.name, end=END)
    nr, nz = int(round(50*level)), int(round(100*level))
    if FAMILY == "uniform":
        edit(case/"system"/"blockMeshDict",
             [(r"\(50 1 100\) simpleGrading \(0\.1 1 1\)",
               f"({nr} 1 {nz}) simpleGrading (1 1 1)")])
    else:
        edit(case/"system"/"blockMeshDict",
             [(r"\(50 1 100\) simpleGrading \(0\.1 1 1\)",
               f"({nr} 1 {nz}) simpleGrading (0.1 1 1)")])
    edit(case/"system"/"fvSchemes",
            [(r"^(\s+default\s+)backward;.*$", r"\g<1>BDF 4;")])
    edit(case/"system"/"controlDict",
            [(r'^libs\s+\("libcirculation1DCoupling.so"\);',
              'libs            ("libcirculation1DCoupling.so" '
              '"libsolids4FoamModels.so");     // BDF ddt scheme'),
             (rf"timeStart\s+{END - 1/mc.F:g};",
              f"timeStart       {END - 0.5/mc.F:g};"),
             (rf"writeInterval\s+{1/(8*mc.F):g};",
              f"writeInterval   {0.5/(8*mc.F):g};")])
    edit(case/"constant"/"circulation1DProperties",
            [(r"^(\s+order\s+)\S+;", r"\g<1>4;"),
             (r"^(\s+CFL\s+)\S+;", r"\g<1>0.5;"),
             (r"^(\s+timeScheme\s+)\S+;.*$",
              "    timeScheme      RK4;\n    reconstruction  kExact;\n"
              "    polynomialOrder 3;"),
             (r"^(\s+initialPressure\s+0;)",
              "\\1\n    linearised          true;")])
    s = (case/"system"/"fvSolution").read_text()
    i = s.index("residualControl")
    s = s[:i] + re.sub(r"tolerance\s+[\d.e+-]+;", "tolerance 1e-10;", s[i:],
                       count=1)
    (case/"system"/"fvSolution").write_text(s)
    return case


def run(case):
    if finished(case):
        return
    with open(case/"log.blockMesh", "w") as f:
        subprocess.run(["blockMesh"], cwd=case, stdout=f,
                       stderr=subprocess.STDOUT, check=True)
    with open(case/"log.pimpleFoam", "w") as f:
        subprocess.run(["pimpleFoam"], cwd=case, stdout=f,
                       stderr=subprocess.STDOUT)
    print(f"{case.name}: {'finished' if finished(case) else 'FAILED'}",
          flush=True)


def ports(case):
    f = case/"postProcessing"/"circulation1DCoupling"/"0"/"coupling.dat"
    names = None
    with open(f) as fh:
        for line in fh:
            if line.startswith("# time"):
                names = line[1:].split()
                break
    a = np.loadtxt(f, comments="#", ndmin=2)
    col = {n: i for i, n in enumerate(names)}
    return a[:, 0], {k: a[:, col[k]] for k in ("to3D.P", "to3D.Q", "from3D.P")}


def profiles(case):
    """{time: (r, u_z)} from the sets written over the window"""
    out = {}
    d = case/"postProcessing"/"profile"
    for td in sorted(d.iterdir(), key=lambda p: float(p.name)):
        f = td/"midLength_U.xy"
        if f.is_file():
            a = np.loadtxt(f)
            out[float(td.name)] = (a[:, 0], a[:, 3])
    return out


def wall_cell(level):
    """radial width of the wall cell of a level: nr cells over R, uniform or
    graded 10:1 towards the wall (simpleGrading 0.1)"""
    nr = int(round(50*level))
    if FAMILY == "uniform":
        return an.R/nr
    q = 0.1**(1.0/(nr - 1))
    return an.R*(1 - q)/(1 - q**nr)*q**(nr - 1)


def wall_cell_r(level):
    """the radius of the inner face of the wall cell of a level: the
    profile is compared up to it. Inside the wall cell the sampled values
    converge at first order only (0.8 on the uniform family, where the wall
    cell is ten times wider than the graded family's; beyond the cell's
    centre the set's cellPoint interpolation runs to the wall's point
    values), while everywhere else the profile converges at second order"""
    return an.R - wall_cell(level)


def analyse(alpha, levels):
    rows = []
    for level in levels:
        case = RUNS/name(alpha, level)
        if not finished(case):
            continue
        t, p = ports(case)
        # the problem the wedge discretises (see analyse.py): the pipe of
        # radius R, for both families
        geo = dict(A=an.A, R=an.R)        # the wedge with flowScale pi/tan(theta/2) is the pipe of radius R
        ex = an.exact_transient(alpha, t, **geo)
        eP = max(math.sqrt(np.mean((p[k] - ex[k])**2))
                 / math.sqrt(np.mean(ex[k]**2)) for k in ("to3D.P", "from3D.P"))
        eQ = (math.sqrt(np.mean((p["to3D.Q"] - ex["to3D.Q"])**2))
              / math.sqrt(np.mean(ex["to3D.Q"]**2)))
        prof = profiles(case)
        times = sorted(prof)
        r = prof[times[0]][0]
        # up to the wall cell (see wall_cell_r)
        inside = r <= wall_cell_r(level)
        r = r[inside]
        uex = an.exact_transient_profile(alpha, r, times, **geo)
        umax = np.abs(uex).max()
        eU = max(np.abs(prof[tk][1][inside] - uex[i]).max()
                 for i, tk in enumerate(times))/umax
        log = (case/"log.pimpleFoam").read_text()
        steps = len(re.findall(r"^Time = ", log, re.M))
        iters = len(re.findall(r"^PIMPLE: iteration", log, re.M))
        ex_t = re.findall(r"ExecutionTime = ([\d.]+)", log)
        rows.append((level, int(round(50*level)), int(round(100*level)),
                     eP, eQ, eU, iters/max(steps, 1),
                     float(ex_t[-1]) if ex_t else math.nan))
    geometry = (
        ("The wedge has no wall grading. " if FAMILY == "uniform" else "")
        + "With its flat face at the distance R from the axis and the patch "
        "flux scaled by pi/tan(theta/2), the discrete problem is the circular "
        "pipe of radius R exactly; the exact solution is Womersley's for R, "
        "independent of the mesh.")
    lines = [f"# Order in space of the Womersley segment, alpha = {alpha:g}"
             + (" (uniform mesh family)" if FAMILY == "uniform" else ""),
             "",
             f"Mesh refined as nr x nz; 0-{END} s from rest, BDF 4 and "
             f"coupling order 4 at dt = {DT*1e3:g} ms, linearised vessels at "
             "fourth order. Errors against the exact transient solution: port "
             "pressures (max of the two) and flow, RMS over the window "
             "relative to the RMS of the exact signal; mid-length velocity "
             "profile, max over eight instants of max |u - u_exact| / max "
             "|u_exact|, from the axis to the wall cell (inside the wall cell "
             "the sampled values converge at first order only). "
             + geometry, "",
             "| nr x nz | P ports | order | Q | order | profile | order | "
             "PIMPLE it./step | time [s] |",
             "|---|---|---|---|---|---|---|---|---|"]
    for i, r in enumerate(rows):
        cells = [f"{r[1]} x {r[2]}"]
        for k in (3, 4, 5):
            cells.append(f"{r[k]:.2e}")
            cells.append("" if i == 0 else
                         f"{math.log2(rows[i-1][k]/r[k]):.2f}")
        cells += [f"{r[6]:.1f}", f"{r[7]:.0f}"]
        lines.append("| " + " | ".join(cells) + " |")
    text = "\n".join(lines) + "\n"
    tag = f"alpha{alpha:g}" + ("_u" if FAMILY == "uniform" else "")
    (HERE/f"space_{tag}.md").write_text(text)
    print(text)
    plot(alpha, rows, tag)


def plot(alpha, rows, tag):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    if len(rows) < 2:
        return
    fig, ax = plt.subplots(figsize=(4.5, 3.6))
    h = np.array([1/r[0] for r in rows])
    for k, lab in ((3, "port pressures"), (4, "flow"), (5, "mid-length profile")):
        ax.loglog(h, [r[k] for r in rows], "o-", label=lab)
    ax.loglog(h, rows[-1][5]*(h/h[-1])**2, "k:", lw=0.8, label="slope 2")
    ax.set_xticks(h)
    ax.set_xticklabels([f"{r[1]}x{r[2]}" for r in rows])
    ax.set_xlabel("mesh")
    ax.set_ylabel("error vs the exact transient")
    ax.set_title(f"Womersley segment, alpha = {alpha:g}: order in space"
                 + (" (uniform mesh)" if FAMILY == "uniform" else ""),
                 fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    (HERE/"plots").mkdir(exist_ok=True)
    fig.savefig(HERE/"plots"/f"space_{tag}.png", dpi=150)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("what", nargs="?", choices=["run"])
    ap.add_argument("--alpha", type=float, default=6.0)
    ap.add_argument("--levels", default=",".join(str(x) for x in LEVELS))
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--dt", type=float, help=f"3D time step in place of {DT} (a time-error check; the runs get a _dt suffix and are not analysed here)")
    ap.add_argument("--family", default="graded", choices=["graded", "uniform"],
                    help="graded: the tutorial's wedge; uniform: no wall "
                    "grading, cross-section pi R^2 (runs space_alpha*_u_x*)")
    args = ap.parse_args()
    global FAMILY
    FAMILY = args.family
    levels = [float(x) for x in args.levels.split(",")]
    if args.what == "run":
        RUNS.mkdir(exist_ok=True)
        cases = [setup(args.alpha, lv, args.force, args.dt) for lv in levels]
        cases.sort(key=lambda c: -c.stat().st_mtime)
        with ThreadPoolExecutor(args.jobs) as ex:
            list(ex.map(run, cases))
        if args.dt:
            return
    analyse(args.alpha, levels)


if __name__ == "__main__":
    main()
