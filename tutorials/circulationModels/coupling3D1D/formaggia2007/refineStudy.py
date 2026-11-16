#!/usr/bin/env python3
"""Mesh refinement of the coupled FSI tube for each end condition of the
wall (clamped, free, spring), against the full tube refined the same way.

    python3 refineStudy.py run [--levels 1,2] [--jobs 4] [--ranks 4]
    python3 refineStudy.py            # analyse what has run

Level 2 doubles the cells of both regions in every direction (fluid
20 x 20 x 160 + ..., solid 8 x 40 x 160 over 5 cm) at the same time step
(25 us: at this step the time error is far below the mesh error, see the
order study), for the coupled runs coupled_c<c0>[_free|_spring]_x2 and the
full tube full3D_x2. Level 1 is the tutorial's mesh (the runs of
makeCases.py: coupled_c<c0>, _free, _spring and full3D, made here if
absent). Reported: RMS(P coupled - P full) / pulse at the sections on
both levels for each end condition, the wall displacement, and the mesh
sensitivity of each run (level 2 against level 1) beside the differences
between end conditions, which is what decides whether the comparison of
end conditions is converged. Writes refine.md and plots/refine.png.
"""
import argparse
import re
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
ENDS = {"clamped": "", "free": "_free", "spring": "_spring"}
SECTIONS = [1.25, 2.5, 3.75, 5, 6.25, 7.5]


def c0():
    return mc.c0_plane_strain()


def base_name(end):
    return "coupled_c" + f"{c0():.3g}".replace(".", "p") + ENDS[end]


def name(end, level):
    return base_name(end) + ("" if level == 1 else f"_x{level:g}")


def ref_name(level):
    return "full3D" + ("" if level == 1 else f"_x{level:g}")


def existing(n):
    """The finished run of that name, or its _lib2 twin (the runs made on
    the Mac with the checkout's own build carry _lib2)"""
    for cand in (n, n + "_lib2", n.replace("_free", "_free_lib2")):
        if an.finished(cand, "log.solids4Foam"):
            return cand
    return None


def refine(case, level):
    for region in ("fluid", "solid"):
        bm = case/"system"/region/"blockMeshDict"
        s = bm.read_text()
        def double(m):
            a, b, c = (int(x) for x in m.group(2).split())
            return f"{m.group(1)}({a*level} {b*level} {c*level})"
        s, n = re.subn(r"(hex\s*\([\d\s]+\)\s*)\((\d+ \d+ \d+)\)", double, s)
        assert n >= 1, bm
        bm.write_text(s)


def setup(level, force=False):
    cases = []
    for end in ENDS:
        n = name(end, level)
        if level == 1 and existing(n):
            cases.append(RUNS/existing(n))
            continue
        if an.finished(n, "log.solids4Foam") and not force:
            cases.append(RUNS/n)
            continue
        mc.make_3d(n, mc.L3D, True, c0(), mc.DELTA_T, True,
                   free_end=(end == "free"), spring_end=(end == "spring"))
        if level != 1:
            refine(RUNS/n, level)
        cases.append(RUNS/n)
    r = ref_name(level)
    if level == 1 and existing(r):
        cases.append(RUNS/existing(r))
    elif an.finished(r, "log.solids4Foam") and not force:
        cases.append(RUNS/r)
    else:
        mc.make_3d(r, mc.LFULL, False, c0(), mc.DELTA_T, True)
        if level != 1:
            refine(RUNS/r, level)
        cases.append(RUNS/r)
    return cases


def run(case, ranks):
    if an.finished(case.name, "log.solids4Foam"):
        return
    def sh(cmd, log):
        with open(case/log, "w") as f:
            r = subprocess.run(cmd, cwd=case, stdout=f, stderr=subprocess.STDOUT)
        if r.returncode:
            sys.exit(f"{case.name}: {' '.join(cmd)} failed, see {case/log}")
    sh(["blockMesh", "-region", "solid"], "log.blockMesh.solid")
    sh(["blockMesh", "-region", "fluid"], "log.blockMesh.fluid")
    sh(["topoSet", "-region", "fluid"], "log.topoSet")
    if ranks > 1 and case.name.startswith("full3D"):
        for region in ("fluid", "solid"):
            dp = case/"system"/region/"decomposeParDict"
            dp.write_text(re.sub(r"numberOfSubdomains\s+\d+;",
                                 f"numberOfSubdomains {ranks};", dp.read_text()))
            sh(["decomposePar", "-region", region], f"log.decomposePar.{region}")
        with open(case/"log.solids4Foam", "w") as f:
            subprocess.run(["mpirun", "--oversubscribe", "-np", str(ranks), "solids4Foam", "-parallel"],
                           cwd=case, stdout=f, stderr=subprocess.STDOUT)
    else:
        with open(case/"log.solids4Foam", "w") as f:
            subprocess.run(["solids4Foam"], cwd=case, stdout=f,
                           stderr=subprocess.STDOUT)
    print(f"{case.name}: {'finished' if an.finished(case.name, 'log.solids4Foam') else 'FAILED'}",
          flush=True)


def sections(n):
    s = an.sections_3d(n)
    if n.startswith("coupled"):
        s.update(an.sections_coupler(n))
    return s


def err_table(a, b, zs):
    """RMS(P_a - P_b)/pulse at the sections both have"""
    out = {}
    for z in zs:
        if z in a and z in b:
            t, P, _ = a[z]
            tb, Pb, _ = b[z]
            out[z] = an.rms(t, P, tb, Pb)/mc.PULSE
    return out


def analyse(levels):
    data, refs = {}, {}
    for level in levels:
        r = existing(ref_name(level)) if level == 1 else ref_name(level)
        if r and an.finished(r, "log.solids4Foam"):
            refs[level] = (r, sections(r))
        for end in ENDS:
            n = existing(name(end, level)) if level == 1 else name(end, level)
            if n and an.finished(n, "log.solids4Foam"):
                data[(end, level)] = (n, sections(n))
    lines = ["# Mesh refinement of the coupled tube, per end condition of the wall",
             "",
             "Level 2: twice the cells in every direction of both regions, "
             f"same step ({mc.DELTA_T*1e6:g} us). RMS(P coupled - P full) over "
             "the window / pulse at the sections (z < 5 cm in the 3D part, "
             "5 cm the port, beyond it the 1D vessel), each level against the "
             "full tube of the same level; D the radial wall displacement at "
             "z = 2.5 cm.", "",
             "| end | level | run / reference | "
             + " | ".join(f"z = {z:g}" for z in SECTIONS)
             + " | D peak full / coupled [um] | RMS dD / peak |",
             "|---|---|---|" + "---|"*len(SECTIONS) + "---|---|"]
    for end in ENDS:
        for level in levels:
            if (end, level) not in data or level not in refs:
                continue
            n, s = data[(end, level)]
            r, sr = refs[level]
            e = err_table(s, sr, SECTIONS)
            tf, Df = an.displacement(r)
            tc, Dc = an.displacement(n)
            lines.append(f"| {end} | {level} | {n} / {r} | "
                         + " | ".join(f"{e[z]:.2e}" if z in e else "-" for z in SECTIONS)
                         + f" | {1e6*Df.max():.1f} / {1e6*Dc.max():.1f} | "
                         f"{an.rms(tc, Dc, tf, Df)/Df.max():.2e} |")
    if len(levels) >= 2 and all(l in refs for l in levels[:2]):
        l1, l2 = levels[:2]
        lines += ["", f"Mesh sensitivity, level {l2} against level {l1} of the "
                  "same run (RMS / pulse), and the differences between end "
                  f"conditions at level {l2} (RMS / pulse):", "",
                  "| what | " + " | ".join(f"z = {z:g}" for z in SECTIONS) + " |",
                  "|---|" + "---|"*len(SECTIONS)]
        e = err_table(refs[l2][1], refs[l1][1], SECTIONS)
        lines.append("| full tube, level 2 - level 1 | "
                     + " | ".join(f"{e[z]:.2e}" if z in e else "-" for z in SECTIONS) + " |")
        for end in ENDS:
            if (end, l1) in data and (end, l2) in data:
                e = err_table(data[(end, l2)][1], data[(end, l1)][1], SECTIONS)
                lines.append(f"| {end}, level 2 - level 1 | "
                             + " | ".join(f"{e[z]:.2e}" if z in e else "-" for z in SECTIONS) + " |")
        ends = [e for e in ENDS if (e, l2) in data]
        for i, a in enumerate(ends):
            for b in ends[i + 1:]:
                e = err_table(data[(a, l2)][1], data[(b, l2)][1], SECTIONS)
                lines.append(f"| {a} - {b}, level {l2} | "
                             + " | ".join(f"{e[z]:.2e}" if z in e else "-" for z in SECTIONS) + " |")
    text = "\n".join(lines) + "\n"
    (HERE/"refine.md").write_text(text)
    print(text)
    plot(data, refs, levels)


def plot(data, refs, levels):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    if not refs:
        return
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    cols = {"clamped": "#D55E00", "free": "#009E73", "spring": "#CC79A7"}
    for ax, z in zip(axs, (2.5, 5)):
        for level, ls in zip(levels, ("--", "-")):
            if level in refs and z in refs[level][1]:
                t, P, _ = refs[level][1][z]
                ax.plot(1e3*t, P, "k", ls=ls, lw=1.2, label=f"full tube, level {level}")
            for end in ENDS:
                if (end, level) in data and z in data[(end, level)][1]:
                    t, P, _ = data[(end, level)][1][z]
                    ax.plot(1e3*t, P, color=cols[end], ls=ls, lw=1.0,
                            label=f"{end}, level {level}")
        ax.set_title(f"P at z = {z:g} cm", fontsize=10)
        ax.set_xlabel("t [ms]")
        ax.grid(alpha=0.3)
    axs[0].set_ylabel("mean pressure [Pa]")
    axs[1].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    (HERE/"plots").mkdir(exist_ok=True)
    fig.savefig(HERE/"plots"/"refine.png", dpi=150)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("what", nargs="?", choices=["run", "setup"])
    ap.add_argument("--levels", default="1,2")
    ap.add_argument("--jobs", type=int, default=4, help="runs at a time")
    ap.add_argument("--ranks", type=int, default=4, help="ranks of the full-tube runs")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    levels = [int(x) for x in args.levels.split(",")]
    if args.what in ("run", "setup"):
        RUNS.mkdir(exist_ok=True)
        cases = [c for lv in levels for c in setup(lv, args.force)]
        if args.what == "setup":
            return
        # the full tubes first, then the finest coupled runs
        cases.sort(key=lambda c: (not c.name.startswith("full3D"), "_x" not in c.name))
        with ThreadPoolExecutor(args.jobs) as ex:
            list(ex.map(lambda c: run(c, args.ranks), cases))
    analyse(levels)


if __name__ == "__main__":
    main()
