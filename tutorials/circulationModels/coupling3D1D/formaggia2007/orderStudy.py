#!/usr/bin/env python3
"""Order of accuracy in time of the coupled FSI + 1D tube, by self-convergence.

    python3 orderStudy.py [--families Euler,backward,BDF3,BDF4] [--dts 1e-4,5e-5,2.5e-5,1.25e-5]
                          [--endTime 0.012] [--ranks 4] [--jobs 2] [--analyse]

For each time-scheme family the coupled run (5 cm FSI tube + 10 cm 1D vessel,
1D wave speed 4.92 m/s, wall clamped at the interface) is repeated at the
listed time steps, with the fluid ddt scheme, the solid d2dt2 scheme and the
3D-1D coupling order all set to the family's order q:

    family    fluid ddt   solid d2dt2   coupling   q
    Euler     Euler       Euler         1          1
    backward  backward    BDF 2         2          2
    BDF3      BDF 3       BDF 3         3          3
    BDF4      BDF 4       BDF 4         4          4

(OpenFOAM's Euler d2dt2 is the central second difference; BDF 2 applied
twice is the backward d2dt2 scheme, which OpenFOAM.com does not have. The BDF
schemes come from libsolids4FoamModels.) The 1D vessel keeps its own SSPRK3
at CFL 0.5, far below the error being measured.

Observed order between three successive steps h, h/2, h/4 is
log2(|u_h - u_h/2| / |u_h/2 - u_h/4|) in the RMS over the run of the mean
pressure at z = 2.5 cm (3D), the pressure and flow at the port, and the
radial wall displacement at z = 2.5 cm. The finest run of every family is
also compared with the finest BDF4 run: they must agree to within the
coarsest of their errors. Writes orders.md and plots/orders.png.
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
import analyse as an                                        # noqa: E402
import makeCases as mc                                      # noqa: E402

RUNS = HERE/"runs"
C0 = 4.92
FAMILIES = {
    "Euler":    dict(q=1, ddt="Euler",    d2dt2="Euler"),
    "backward": dict(q=2, ddt="backward", d2dt2="BDF 2"),
    "BDF3":     dict(q=3, ddt="BDF 3",    d2dt2="BDF 3"),
    "BDF4":     dict(q=4, ddt="BDF 4",    d2dt2="BDF 4"),
    # Diagnostics: which of the three order-3 ingredients matters
    "BDF3c2":   dict(q=3, ddt="BDF 3",    d2dt2="BDF 3", coupling=2),
    "BDF3s2":   dict(q=3, ddt="BDF 3",    d2dt2="BDF 2"),
    "BDF3f2":   dict(q=3, ddt="backward", d2dt2="BDF 3"),
    # Fluid and coupling at order 3 or 4 with the solid at order 2: the BDF 3
    # and 4 d2dt2 schemes amplify the undamped wall modes and the FSI loop
    # diverges with them (see bdfD2dt2Scheme), so this is the usable pairing
    "BDF4s2":   dict(q=4, ddt="BDF 4",    d2dt2="BDF 2"),
}


SMOOTH_PULSE = ("#{ (arg() < 3e-3) ? 1.3332*sqr(sin(1047.1975511965977*arg()))"
                " : 0 #}")                   # p/rho: 1333.2 Pa sin^2(pi t/3 ms)


def run_name(family, dt, smooth=False):
    return f"order_{family}_dt{dt:g}" + ("_smooth" if smooth else "")


def make(family, dt, end_time, force, smooth=False, tight=False):
    name = run_name(family, dt, smooth)
    case = RUNS/name
    if case.exists() and mc.finished(case) and not force:
        return name
    mc.make_3d(name, mc.L3D, True, C0, dt, True)
    if smooth:
        # A smooth pulse of the same amplitude and duration, in place of the
        # tutorial's step with its 0.1 ms ramp, so that the time error of
        # the schemes is what the study measures
        mc.edit(case/"0"/"fluid"/"p",
                [(r"uniformValue\s+tableFile;",
                  "uniformValue\n        {\n            type        expression;\n"
                  f"            expression  {SMOOTH_PULSE};\n        }}")])
    if tight:
        mc.edit(case/"constant"/"fsiProperties",
                [(r"outerCorrTolerance\s+1e-6;", "outerCorrTolerance  1e-9;")])
        mc.edit(case/"system"/"fluid"/"fvSolution",
                [(r'("p\|pFinal\|pcorr\|pcorrFinal"\s*\{\s*solver\s+GAMG;\s*tolerance\s+)1e-06;',
                  r"\g<1>1e-10;"),
                 (r"(solver\s+PBiCG;\s*preconditioner\s+DILU;\s*tolerance\s+)1e-09;",
                  r"\g<1>1e-12;")])
    f = FAMILIES[family]
    mc.edit(case/"system"/"fluid"/"fvSchemes",
            [(r"(ddtSchemes\s*\{\s*default\s+)Euler;", rf"\g<1>{f['ddt']};")])
    mc.edit(case/"system"/"solid"/"fvSchemes",
            [(r"(d2dt2Schemes\s*\{\s*default\s+)Euler;", rf"\g<1>{f['d2dt2']};")])
    mc.edit(case/"constant"/"fluid"/"circulation1DProperties",
            [(r"order\s+2;", f"order       {f.get('coupling', f['q'])};")])
    mc.edit(case/"system"/"controlDict",
            [(r"endTime\s+[\d.e-]+;", f"endTime         {end_time:g};"),
             (r"writeInterval\s+1e-3;", f"writeInterval   {end_time:g};")])
    return name


def run(name, ranks):
    case = RUNS/name
    log = case/"log.solids4Foam"
    if mc.finished(case):
        return
    cmds = ["blockMesh -region solid", "blockMesh -region fluid",
            "topoSet -region fluid"]
    if ranks > 1:
        cmds += ["decomposePar -region fluid", "decomposePar -region solid"]
    for c in cmds:
        with open(case/("log." + c.split()[0]), "w") as fh:
            subprocess.run(c.split(), cwd=case, stdout=fh, stderr=subprocess.STDOUT,
                           check=True)
    solver = (["mpirun", "-np", str(ranks), "solids4Foam", "-parallel"]
              if ranks > 1 else ["solids4Foam"])
    with open(log, "w") as fh:
        subprocess.run(solver, cwd=case, stdout=fh, stderr=subprocess.STDOUT)
    print(f"{name}: {'finished' if mc.finished(case) else 'FAILED'}", flush=True)


def signals(name):
    """{label: (t, x)} of the compared quantities"""
    s = an.sections_3d(name)
    s.update(an.sections_coupler(name))
    t, P, _ = s[2.5]
    tp, Pp, Qp = s[5]
    td, D = an.displacement(name)
    return {"P(z=2.5 cm)": (t, P), "P(port)": (tp, Pp), "Q(port)": (tp, Qp),
            "wall u_r(z=2.5 cm)": (td, D)}


def diff(a, b):
    """RMS of a - b on a's times, both covering the same span"""
    t, x = a
    tb, xb = b
    m = (t >= tb[0]) & (t <= tb[-1])
    return math.sqrt(np.mean((x[m] - np.interp(t[m], tb, xb))**2))


def analyse(families, dts, end_time, smooth=False):
    dts = sorted(dts, reverse=True)
    lines = [f"# Time order of the coupled FSI + 1D tube (self-convergence, "
             f"0-{1e3*end_time:g} ms"
             + (", smooth sin^2 pulse)" if smooth else ", the tutorial's step pulse)"),
             ""]
    finest = {}
    table = {}
    for fam in families:
        names = [run_name(fam, dt, smooth) for dt in dts]
        done = [n for n in names if mc.finished(RUNS/n)]
        failed = [n for n in names if n not in done]
        if failed:
            lines += [f"{fam}: failed at dt = "
                      + ", ".join(f"{dt:g}" for n, dt in zip(names, dts)
                                  if n in failed), ""]
        # Orders need three consecutive halvings, so the finished runs must
        # be consecutive in the list
        done = done[:next((i for i, n in enumerate(names) if n in failed),
                          len(names))]
        if len(done) < 3:
            lines += [f"{fam}: {len(done)} consecutive finished run(s), need 3",
                      ""]
            continue
        sig = {n: signals(n) for n in done}
        finest[fam] = sig[done[-1]]
        lines += [f"## {fam} (q = {FAMILIES[fam]['q']})", "",
                  "| dt pair | " + " | ".join(sig[done[0]]) + " |",
                  "|---|" + "---|"*len(sig[done[0]])]
        for i in range(len(done) - 2):
            row = []
            for key in sig[done[0]]:
                e1 = diff(sig[done[i]][key], sig[done[i + 1]][key])
                e2 = diff(sig[done[i + 1]][key], sig[done[i + 2]][key])
                p = math.log2(e1/e2) if e2 > 0 else math.nan
                row.append(f"{p:.2f}")
                table.setdefault(fam, {}).setdefault(key, []).append((dts[i + 1], e2))
            lines.append(f"| {dts[i]:g} / {dts[i+1]:g} / {dts[i+2]:g} | "
                         + " | ".join(row) + " |")
        lines.append("")
    if "BDF4" in finest and len(finest) > 1:
        lines += ["## Finest runs against the finest BDF4 run (RMS of the "
                  "difference / RMS of the signal)", ""]
        ref = finest["BDF4"]
        for fam, sig in finest.items():
            if fam == "BDF4":
                continue
            cells = []
            for key in sig:
                t, x = ref[key]
                cells.append(f"{key}: {diff(sig[key], ref[key])/math.sqrt(np.mean(x**2)):.2e}")
            lines.append(f"- {fam}: " + ", ".join(cells))
        lines.append("")
    out = HERE/("orders_smooth.md" if smooth else "orders.md")
    out.write_text("\n".join(lines))
    print("\n".join(lines))
    plot(table, "orders_smooth.png" if smooth else "orders.png")


def plot(table, filename):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    if not table:
        return
    keys = list(next(iter(table.values())))
    fig, axs = plt.subplots(1, len(keys), figsize=(3.6*len(keys), 3.4), squeeze=False)
    cols = {"Euler": "#999999", "backward": "#0072B2", "BDF3": "#009E73", "BDF4": "#D55E00"}
    for ax, key in zip(axs[0], keys):
        for fam, d in table.items():
            if key in d:
                dt, e = zip(*d[key])
                ax.loglog(dt, e, "o-", color=cols.get(fam, "k"), label=fam)
        ax.set_title(key, fontsize=9)
        ax.set_xlabel("dt [s]")
        ax.grid(alpha=0.3, which="both")
    axs[0, 0].set_ylabel("|u(dt) - u(dt/2)|, RMS")
    axs[0, 0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    (HERE/"plots").mkdir(exist_ok=True)
    fig.savefig(HERE/"plots"/filename, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--families", default="Euler,backward,BDF3,BDF4")
    ap.add_argument("--dts", default="1e-4,5e-5,2.5e-5,1.25e-5")
    ap.add_argument("--endTime", type=float, default=0.012)
    ap.add_argument("--ranks", type=int, default=4)
    ap.add_argument("--jobs", type=int, default=2, help="runs at a time")
    ap.add_argument("--smooth", action="store_true",
                    help="smooth sin^2 inlet pulse instead of the step")
    ap.add_argument("--tight", action="store_true",
                    help="FSI tolerance 1e-9, p 1e-10, U 1e-12")
    ap.add_argument("--analyse", action="store_true", help="analyse only")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    families = args.families.split(",")
    dts = [float(x) for x in args.dts.split(",")]
    if not args.analyse:
        RUNS.mkdir(exist_ok=True)
        names = [make(f, dt, args.endTime, args.force, args.smooth, args.tight)
                 for f in families for dt in dts]
        with ThreadPoolExecutor(args.jobs) as pool:
            list(pool.map(lambda n: run(n, args.ranks), names))
    analyse(families, dts, args.endTime, args.smooth)


if __name__ == "__main__":
    main()
