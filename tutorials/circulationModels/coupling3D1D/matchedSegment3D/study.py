#!/usr/bin/env python3
"""B5, the matched segment: the refinement study (Allrun-study).

    python3 study.py run [NPROC] [--levels N] [--aitken | --explicit P]
                                                set up and run
    python3 study.py                            analyse what has run

The 3D cells along the pipe and the 3D step are refined together at a fixed
acoustic Courant number c0 dt/dz = 0.49:
    dz = 5, 2.5, 1.25, 0.625 mm  with  dt = 0.4, 0.2, 0.1, 0.05 ms
(--levels N runs the first N; default 3). The vessels stay at 1600 cells
(0.625 mm, k-exact 3 + RK4, linearised), whose error (about 1e-5 of the
peak) is far below the 3D one, so the study measures the 3D solver carrying
the wave and the coupling. Implicit coupling of order 2, BDF2 in the 3D;
--aitken uses the Aitken-relaxed interface (runs/aitken_*), --explicit P
explicit coupling of order P (runs/explicitP_*), as side checks.

Nothing reaches the pipe before 0.16 s, even at 1e-12 of the peak, so each
case runs in two phases (as Allrun): to 0.12 s with one outer corrector, then
from the written state to 0.35 s with 20 (system/fvSolution), which iterate
the implicit coupling to about 1e-5 of the step's change of the port flow.

Writes results.md and plots/matchedSegment3D.png. Cases are serial; NPROC
(default 3) run at a time. Already finished cases are skipped.
"""
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import analyse as an  # noqa: E402

RUNS = HERE/"runs"
LEVELS = [(40, 4e-4), (80, 2e-4), (160, 1e-4), (320, 5e-5)]
NOUTER = 20
PHASE1 = 0.12
END = 0.35
KEYS = ["P1", "Q1", "P2", "Q2", "refl", "trans"]


TAGS = ["dirichlet", "aitken", "explicit1", "explicit2"]


def name(nz, tag):
    return f"{tag}_nz{nz}"


def finished(case):
    log = case/"log.rhoPimpleFoam.2"
    return log.is_file() and log.read_text().rstrip().endswith("End")


def foam_set(case, entry, value, file):
    subprocess.run(["foamDictionary", "-entry", entry, "-set", str(value),
                    file], cwd=case, check=True, stdout=subprocess.DEVNULL)


def setup(nz, dt, tag):
    case = RUNS/name(nz, tag)
    if finished(case):
        return None
    if case.exists():
        shutil.rmtree(case)
    case.mkdir(parents=True)
    for d in ("0", "constant", "system"):
        shutil.copytree(HERE/d, case/d,
                        ignore=shutil.ignore_patterns("polyMesh"))
    # A text edit: foamDictionary would write the blocks back with $nz expanded
    bmd = case/"system"/"blockMeshDict"
    bmd.write_text(re.sub(r"^nz\s+\d+;", f"nz      {nz};", bmd.read_text(),
                          flags=re.M))
    foam_set(case, "deltaT", f"{dt:g}", "system/controlDict")
    props = "constant/circulation1DProperties"
    if tag == "aitken":
        foam_set(case, "coupling.interface", "aitken", props)
    if tag.startswith("explicit"):
        foam_set(case, "coupling.scheme", "explicit", props)
        foam_set(case, "coupling.order", tag[len("explicit"):], props)
    return case


def run(case):
    def solver(log):
        with open(case/log, "w") as f:
            return subprocess.run(["rhoPimpleFoam"], cwd=case, stdout=f,
                                  stderr=subprocess.STDOUT).returncode
    with open(case/"log.blockMesh", "w") as f:
        subprocess.run(["blockMesh"], cwd=case, stdout=f,
                       stderr=subprocess.STDOUT, check=True)
    # Phase 1: the pulse is still in the proximal vessel
    foam_set(case, "endTime", PHASE1, "system/controlDict")
    foam_set(case, "writeInterval", PHASE1, "system/controlDict")
    foam_set(case, "PIMPLE.nOuterCorrectors", 1, "system/fvSolution")
    r = solver("log.rhoPimpleFoam.1")
    # Phase 2: the pulse through the pipe, the coupling iterated
    if r == 0:
        foam_set(case, "startFrom", "latestTime", "system/controlDict")
        foam_set(case, "endTime", END, "system/controlDict")
        foam_set(case, "writeInterval", END, "system/controlDict")
        foam_set(case, "PIMPLE.nOuterCorrectors", NOUTER, "system/fvSolution")
        r = solver("log.rhoPimpleFoam.2")
    print(f"{case.name}: {'done' if r == 0 else 'stopped'}", flush=True)


def order(e0, e1):
    return np.log2(e0/e1) if e0 > 0 and e1 > 0 else float("nan")


def analyse():
    lines = ["# B5, the matched segment: refinement against the exact solution",
             "",
             "3D cells along the pipe and the 3D step refined together "
             "(acoustic Courant 0.49); vessels at 1600 cells, k-exact 3 + RK4, "
             "linearised; implicit coupling of order 2 with 20 outer "
             "correctors per step, BDF2; errors relative "
             "to the incident peak (maximum over the run), the delay relative "
             f"to L/c0 = {an.LPIPE/an.C0*1e3:.3f} ms.", ""]
    results = {}
    titles = {"dirichlet": "Implicit coupling, Dirichlet interface",
              "aitken": "Implicit coupling, Aitken-relaxed interface",
              "explicit1": "Explicit coupling of order 1",
              "explicit2": "Explicit coupling of order 2"}
    for tag in TAGS:
        rows = []
        for nz, dt in LEVELS:
            case = RUNS/name(nz, tag)
            if finished(case):
                rows.append((nz, dt, an.errors(case)))
        if not rows:
            continue
        results[tag] = rows
        lines += [f"## {titles[tag]}",
                  "",
                  "| dz [mm] | dt [ms] | " + " | ".join(
                      f"{k} | order" for k in KEYS) + " | delay | exch./step |",
                  "|" + "---|"*(2 + 2*len(KEYS) + 2)]
        for i, (nz, dt, e) in enumerate(rows):
            cells = [f"{an.LPIPE/nz*1e3:g}", f"{dt*1e3:g}"]
            for k in KEYS:
                cells.append(f"{e[k]:.2e}")
                cells.append("" if i == 0 else
                             f"{order(rows[i-1][2][k], e[k]):.2f}")
            cells += [f"{e['delay']:+.1e}", f"{e['updates']:.1f}"]
            lines.append("| " + " | ".join(cells) + " |")
        lines.append("")
    (HERE/"results.md").write_text("\n".join(lines))
    print("\n".join(lines))
    if results:
        figure(results)


def figure(results):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    tag = "dirichlet" if "dirichlet" in results else next(iter(results))
    rows = results[tag]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    dz = np.array([an.LPIPE/nz*1e3 for nz, _, _ in rows])
    labels = {"P1": "P, inlet", "Q1": "Q, inlet", "P2": "P, outlet",
              "Q2": "Q, outlet", "refl": "reflection", "trans": "transmitted"}
    for k in KEYS:
        ax[0].loglog(dz, [e[k] for _, _, e in rows], "o-", label=labels[k])
    ax[0].loglog(dz, rows[0][2]["P2"]*(dz/dz[0])**2, "k:", label="slope 2")
    ax[0].set_xticks(dz, [f"{x:g}" for x in dz])
    ax[0].minorticks_off()
    ax[0].set_xlabel("dz [mm], with c0 dt/dz = 0.49")
    ax[0].set_ylabel("error / incident peak")
    ax[0].set_title("Against the exact solution")
    ax[0].legend(fontsize=8)
    ax[0].grid(True, which="both", alpha=0.3)

    nz, dt, _ = rows[-1]
    d = an.load(RUNS/name(nz, tag))
    t = d["t"]
    ax[1].plot(t, d["P1"]/an.P_INC, label="inlet")
    ax[1].plot(t, d["P2"]/an.P_INC, label="outlet")
    ax[1].plot(t, an.pulse(t, an.T1), "k:", label="exact")
    ax[1].plot(t, an.pulse(t, an.T2), "k:")
    ax[1].set_xlim(0.17, 0.29)
    ax[1].set_xlabel("t [s]")
    ax[1].set_ylabel("P / incident peak")
    ax[1].set_title(f"Pressure at the pipe's ends, dz = {an.LPIPE/nz*1e3:g} mm")
    ax[1].legend(fontsize=8)
    ax[1].grid(True, alpha=0.3)

    for nz, _, _ in rows:
        d = an.load(RUNS/name(nz, tag))
        m = d["t"] > an.T_PM + 8*an.SIGMA
        ax[2].plot(d["t"][m], d["Ppm"][m]/an.P_INC*1e3,
                   label=f"dz = {an.LPIPE/nz*1e3:g} mm")
    ax[2].axhline(0, color="k", ls=":", label="exact")
    ax[2].set_xlabel("t [s]")
    ax[2].set_ylabel("P / incident peak  x 1000")
    ax[2].set_title("Reflection, at the proximal midpoint")
    ax[2].legend(fontsize=8)
    ax[2].grid(True, alpha=0.3)
    (HERE/"plots").mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(HERE/"plots"/"matchedSegment3D.png", dpi=150)


def main():
    args = sys.argv[1:]
    if args and args[0] == "run":
        nproc = int(args[1]) if len(args) > 1 and args[1].isdigit() else 3
        nlev = 3
        if "--levels" in args:
            nlev = int(args[args.index("--levels") + 1])
        tag = "aitken" if "--aitken" in args else "dirichlet"
        if "--explicit" in args:
            tag = "explicit" + args[args.index("--explicit") + 1]
        cases = [c for nz, dt in LEVELS[:nlev]
                 if (c := setup(nz, dt, tag)) is not None]
        with ThreadPoolExecutor(max_workers=nproc) as pool:
            list(pool.map(run, cases))
    analyse()


if __name__ == "__main__":
    main()
