#!/usr/bin/env python3
"""The rigid-segment 3D-1D benchmark at fourth order in space and time.

    python3 highOrder.py [run] [NPROC]      set up and run, then analyse
    python3 highOrder.py                    analyse what has run

The vessels use the k-exact reconstruction of degree 3 with RK4 (fourth order
in space and time, CFL 0.5), the pipe pimpleFoam with BDF 4; the 3D solution,
plug flow, is exact in space. The mesh and the 3D time step are refined
together, dx = 2.5, 1.25, 0.625, 0.3125 mm with dt = 2, 1, 0.5, 0.25 ms, for
the couplings

    I4   implicit, order 4        E4   explicit, order 4
    I2   implicit, order 2        I1   implicit, order 1

and each run is compared with the implicit order 4 run at dx = 0.15625 mm,
dt = 0.125 ms (self-convergence). The inflow is that of Boileau Eq. 14 at
full amplitude: the 3D solution is exact at any amplitude, and at 1e-3 the
round-off of the vessels (5e-7 of the peak) would hide a fourth-order error.

Measured, relative to their peaks: the 3D flux (max over the steps), the
pressures at the two pipe ends (the vessels' boundary states, point values),
and the profiles along both vessels at t = 0.35 s (cell averages; the
reference averaged onto the coarser cells), which hold the reflected and the
transmitted pulse. Writes highOrder.md and plots/highOrder3D1D.png.
"""
import math
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RUNS = HERE/"runs"
END = 0.35
LEVELS = [(400, 2e-3), (800, 1e-3), (1600, 5e-4), (3200, 2.5e-4)]
REF = (6400, 1.25e-4)
COUPLINGS = [("implicit", 4), ("explicit", 4), ("implicit", 2),
             ("implicit", 1)]


def name(scheme, order, n):
    return f"HO_{scheme[0].upper()}{order}_n{n}"


def edit(path, pairs):
    s = path.read_text()
    for pattern, value in pairs:
        s, k = re.subn(pattern, value, s, count=1, flags=re.M)
        if k != 1:
            sys.exit(f"no match for {pattern} in {path}")
    path.write_text(s)


def finished(case):
    log = case/"log.pimpleFoam"
    return log.is_file() and log.read_text().rstrip().endswith("End")


def setup(scheme, order, n, dt):
    case = RUNS/name(scheme, order, n)
    if finished(case):
        return None
    if case.exists():
        shutil.rmtree(case)
    shutil.copytree(HERE/"base", case)
    edit(case/"system"/"controlDict", [
        (r"^deltaT\s+\S+;", f"deltaT          {dt!r};"),
        (r"^endTime\s+\S+;", f"endTime         {END!r};"),
        (r"^writeControl\s+\S+;", "writeControl    runTime;"),
        (r"^writeInterval\s+\S+;.*$", f"writeInterval   {END!r};"),
        (r'^libs\s+\("libcirculation1DCoupling.so"\);',
         'libs            ("libcirculation1DCoupling.so" '
         '"libsolids4FoamModels.so");     // BDF ddt scheme'),
    ])
    edit(case/"system"/"fvSchemes", [
        (r"^(\s+default\s+)backward;.*$", r"\g<1>BDF 4;")])
    edit(case/"constant"/"circulation1DProperties", [
        (r"^(\s+scheme\s+)\S+;", rf"\g<1>{scheme};"),
        (r"^(\s+order\s+)\S+;", rf"\g<1>{order};"),
        (r"^(\s+CFL\s+)\S+;", r"\g<1>0.5;"),
        (r"^(\s+timeScheme\s+)\S+;.*$",
         "    timeScheme      RK4;\n    reconstruction  kExact;\n"
         "    polynomialOrder 3;"),
        (r"^(\s+nCells\s+)\S+;.*$", rf"\g<1>{n};"),
        (r"^(\s+scale\s+)\S+;.*$", r"\g<1>1;"),
    ])
    return case


def run(case):
    with open(case/"log.blockMesh", "w") as f:
        subprocess.run(["blockMesh"], cwd=case, stdout=f,
                       stderr=subprocess.STDOUT, check=True)
    with open(case/"log.pimpleFoam", "w") as f:
        r = subprocess.run(["pimpleFoam"], cwd=case, stdout=f,
                           stderr=subprocess.STDOUT)
    print(f"{case.name}: {'done' if r.returncode == 0 else 'stopped'}",
          flush=True)


def load(case, explicit):
    """Coupler rows (the 3D flux of explicit coupling is the sample of the
    next row) and the profiles at END"""
    out = case/"postProcessing"/"circulation1DCoupling"/"0"
    a = np.loadtxt(out/"coupling.dat", comments="#", ndmin=2)
    t, Q, P1, P2 = a[:, 0], a[:, 10], a[:, 3], a[:, 5]
    if explicit:
        t, P1, P2, Q = t[:-1], P1[:-1], P2[:-1], a[1:, 10]
    prof = {}
    for v in ("proximal", "distal"):
        f = out/"profiles"/f"{END:g}"/f"{v}.dat"
        prof[v] = np.loadtxt(f, comments="#")[1:-1, 2]      # Q, cell averages
    return t, Q, P1, P2, prof


def analyse():
    refcase = RUNS/name("implicit", 4, REF[0])
    if not finished(refcase):
        sys.exit("No reference run: python3 highOrder.py run")
    tr, Qr, P1r, P2r, profr = load(refcase, False)
    key = np.rint(tr/REF[1]).astype(int)

    def at(tc, f):
        idx = {k: i for i, k in enumerate(key)}
        return f[[idx[int(round(x/REF[1]))] for x in tc]]

    peakQ = np.abs(Qr).max()
    peakP = np.abs(P1r).max()
    peakV = max(np.abs(p).max() for p in profr.values())

    results = {}
    for scheme, order in COUPLINGS:
        rows = []
        for n, dt in LEVELS:
            case = RUNS/name(scheme, order, n)
            if not finished(case):
                continue
            t, Q, P1, P2, prof = load(case, scheme == "explicit")
            eQ = np.abs(Q - at(t, Qr)).max()/peakQ
            eP = max(np.abs(P1 - at(t, P1r)).max(),
                     np.abs(P2 - at(t, P2r)).max())/peakP
            r = REF[0]//n
            eV = max(np.abs(prof[v] - profr[v].reshape(n, r).mean(axis=1))
                     .max() for v in prof)/peakV
            rows.append((n, dt, eQ, eP, eV))
        results[(scheme, order)] = rows

    lines = ["# Rigid segment 3D-1D at fourth order: self-convergence", "",
             "Vessels: kExact (degree 3) + RK4, CFL 0.5; pipe: pimpleFoam, "
             "BDF 4; full amplitude. dx and dt refined together; reference: "
             f"implicit 4 at {REF[0]} cells, dt = {REF[1]*1e3:g} ms. Errors "
             "relative to the peaks: 3D flux, pressures at the pipe ends "
             "(max of the two), profiles of Q along both vessels at "
             f"t = {END} s.", ""]
    for (scheme, order), rows in results.items():
        if not rows:
            continue
        lines += [f"## {scheme} order {order}", "",
                  "| cells | dt [ms] | 3D flux | order | pipe-end pressure | "
                  "order | vessel profiles | order |",
                  "|---|---|---|---|---|---|---|---|"]
        for i, r in enumerate(rows):
            cells = [str(r[0]), f"{r[1]*1e3:g}"]
            for k in (2, 3, 4):
                cells.append(f"{r[k]:.2e}")
                cells.append("" if i == 0 else
                             f"{math.log(rows[i-1][k]/r[k])/math.log(2):.2f}")
            lines.append("| " + " | ".join(cells) + " |")
        lines.append("")
    (HERE/"highOrder.md").write_text("\n".join(lines))
    print("\n".join(lines))
    figure(results)


def figure(results):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.ticker
    except ImportError:
        return
    (HERE/"plots").mkdir(exist_ok=True)
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.8), sharey=True)
    titles = ("3D flux", "pressure at the pipe ends",
              "waves in the vessels (profiles)")
    for ax, k, title in zip(axs, (2, 3, 4), titles):
        for (scheme, order), rows in results.items():
            if rows:
                dx = np.array([1/r[0] for r in rows])*1e3
                ax.loglog(dx, [r[k] for r in rows],
                          "s--" if scheme == "implicit" else "o-", ms=4,
                          label=f"{scheme} {order}")
        # Slope guides through the finest point of the coupling of that order
        guides = [(("implicit", 4), 4, "-.")]
        if k == 4:
            guides += [(("implicit", 2), 2, "--"), (("implicit", 1), 1, ":")]
        for c, q, ls in guides:
            rows = results[c]
            dx = np.array([1/r[0] for r in rows])*1e3
            ax.loglog(dx, rows[-1][k]*(dx/dx[-1])**q*0.5, "k", ls=ls, lw=0.8,
                      label=f"slope {q}")
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xticks([0.3125, 0.625, 1.25, 2.5])
        ax.set_xticklabels(["0.31", "0.62", "1.25", "2.5"])
        ax.set_xlabel("dx [mm], with dt = 0.8 ms per mm")
        ax.set_title(title)
    axs[0].set_ylabel("error / peak")
    axs[2].legend(frameon=False, fontsize=7, loc="lower right", ncol=2)
    fig.tight_layout()
    fig.savefig(HERE/"plots"/"highOrder3D1D.png", dpi=150)


def main():
    if sys.argv[1:2] == ["run"]:
        nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 6
        RUNS.mkdir(exist_ok=True)
        cases = [setup("implicit", 4, *REF)]
        for scheme, order in COUPLINGS:
            for n, dt in LEVELS:
                cases.append(setup(scheme, order, n, dt))
        cases = [c for c in cases if c is not None]
        print(f"running {len(cases)} cases, {nproc} at a time", flush=True)
        with ThreadPoolExecutor(nproc) as pool:
            list(pool.map(run, cases))
    analyse()


if __name__ == "__main__":
    main()
