#!/usr/bin/env python3
"""The rigid-segment 3D-1D benchmark at fourth order in space and time.

    python3 highOrder.py [run] [NPROC]      set up and run, then analyse
    python3 highOrder.py                    analyse what has run
    python3 highOrder.py [run] [NPROC] --linear     against the closed form
    python3 highOrder.py --reference     the finest nonlinear run as a numerical
                                         reference (reference/*.csv)

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

With --linear the vessels solve the linearised 1D equations ("linearised
true"), of which the closed form of BENCHMARKS3D1D.md is the exact solution
at any amplitude, so every run is compared with the closed form (the
profiles with its cell averages, by Gauss quadrature) and the refinement
goes one level further (0.15625 mm, 0.125 ms). The closed form includes the
return of the reflected pulse to the flow inlet, where it reflects with +1
in pressure, which has reached the first centimetres of the proximal vessel
at t = 0.35 s.

Measured, relative to their peaks: the 3D flux (max over the steps), the
pressures at the two pipe ends (the vessels' boundary states, point values),
and the profiles along both vessels at t = 0.35 s (cell averages; the
reference averaged onto the coarser cells), which hold the reflected and the
transmitted pulse. Writes highOrder.md and plots/highOrder3D1D.png
(highOrder_linear.md, plots/highOrder3D1D_linear.png with --linear).
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
sys.path.insert(0, str(HERE.parents[3]/"src"/"circulationModels"/"docs"))
sys.path.insert(0, str(HERE))

from benchmarks3D1D import transmitted  # noqa: E402
from sweep import C0, TAU, Z0  # noqa: E402

RUNS = HERE/"runs"
END = 0.35
LEVELS = [(400, 2e-3), (800, 1e-3), (1600, 5e-4), (3200, 2.5e-4)]
REF = (6400, 1.25e-4)
COUPLINGS = [("implicit", 4), ("explicit", 4), ("implicit", 2),
             ("implicit", 1)]
LINEAR = "--linear" in sys.argv
ROBIN = "--robin" in sys.argv          # the coupler's Robin (impedance) interface

# The closed form: the inflow of Boileau Eq. 14 at full amplitude
L_VESSEL = 1.0
SIGMA = 1/math.sqrt(2e4)
T_PEAK = 0.05
T_ARRIVE = T_PEAK + L_VESSEL/C0         # incident peak at the pipe
P_INC = Z0*1e-6                         # incident peak [Pa]


def name(scheme, order, n):
    return (f"{'HOL' if LINEAR else 'HO'}{'R' if ROBIN else ''}"
            f"_{scheme[0].upper()}{order}_n{n}")


def pinc(T):
    return P_INC*np.exp(-T**2/(2*SIGMA**2))


def pt(T):
    return transmitted(T, P_INC, SIGMA, TAU)


def exact_pipe(t):
    """Q through the pipe, P at its inlet and outlet"""
    T = t - T_ARRIVE
    p = pt(T)
    return p/Z0, 2*pinc(T) - p, p


def exact_profile(vessel, x, t):
    """Q along a vessel at time t, x from its in end (point values)"""
    if vessel == "proximal":
        d = L_VESSEL - x                            # distance to the pipe
        Tr = t - T_ARRIVE - d/C0                    # reflected, backward
        Trr = t - T_ARRIVE - (L_VESSEL + x)/C0      # re-reflected at the
        pr = lambda T: pinc(T) - pt(T)              # flow inlet, +1
        return (pinc(t - x/C0 - T_PEAK) - pr(Tr) + pr(Trr))/Z0
    return pt(t - T_ARRIVE - x/C0)/Z0


GAUSS = np.polynomial.legendre.leggauss(5)


def exact_profile_avg(vessel, x, dx, t):
    """Cell averages of the closed form over [x - dx/2, x + dx/2]"""
    return sum(w*exact_profile(vessel, x + 0.5*dx*g, t)
               for g, w in zip(*GAUSS))/2


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
    ] + ([(r"^(\s+initialPressure\s+0;)",
           "\\1\n    linearised          true;")] if LINEAR else [])
      + ([(r"^(\s+order\s+\S+;)", "\\1\n    interface   robin;")] if ROBIN else []))
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
        a = np.loadtxt(f, comments="#")[1:-1]
        prof[v] = a[:, 2]                                   # Q, cell averages
        prof[v + ".x"] = a[:, 0]
    return t, Q, P1, P2, prof


def analyse():
    levels = LEVELS + [REF] if LINEAR else LEVELS

    if LINEAR:
        peakQ = P_INC/Z0
        peakP = np.abs(exact_pipe(np.linspace(0, END, 20001))[1]).max()
        peakV = peakQ
    else:
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
        peakV = max(np.abs(profr[v]).max() for v in ("proximal", "distal"))

    results = {}
    for scheme, order in COUPLINGS:
        rows = []
        for n, dt in levels:
            case = RUNS/name(scheme, order, n)
            if not finished(case):
                continue
            t, Q, P1, P2, prof = load(case, scheme == "explicit")
            if LINEAR:
                Qx, P1x, P2x = exact_pipe(t)
                eQ = np.abs(Q - Qx).max()/peakQ
                eP = max(np.abs(P1 - P1x).max(), np.abs(P2 - P2x).max())/peakP
                eV = max(np.abs(prof[v] - exact_profile_avg(
                    v, prof[v + ".x"], L_VESSEL/n, END)).max()
                    for v in ("proximal", "distal"))/peakV
            else:
                eQ = np.abs(Q - at(t, Qr)).max()/peakQ
                eP = max(np.abs(P1 - at(t, P1r)).max(),
                         np.abs(P2 - at(t, P2r)).max())/peakP
                r = REF[0]//n
                eV = max(np.abs(prof[v] - profr[v].reshape(n, r).mean(axis=1))
                         .max() for v in ("proximal", "distal"))/peakV
            rows.append((n, dt, eQ, eP, eV))
        results[(scheme, order)] = rows

    if LINEAR:
        lines = ["# Rigid segment 3D-1D at fourth order: against the closed "
                 "form", "",
                 "Vessels: linearised 1D equations, kExact (degree 3) + RK4, "
                 "CFL 0.5; pipe: pimpleFoam, BDF 4; full amplitude. dx and "
                 "dt refined together. Errors against the closed form "
                 f"(tau = {TAU*1e3:.4f} ms with the patch area pi R^2), "
                 "relative to the peaks: 3D flux, pressures at the pipe ends "
                 "(max of the two), profiles of Q along both vessels at "
                 f"t = {END} s (cell averages of the closed form).", ""]
    else:
        lines = ["# Rigid segment 3D-1D at fourth order: self-convergence", "",
                 "Vessels: kExact (degree 3) + RK4, CFL 0.5; pipe: pimpleFoam, "
                 "BDF 4; full amplitude. dx and dt refined together; "
                 f"reference: implicit 4 at {REF[0]} cells, dt = "
                 f"{REF[1]*1e3:g} ms. Errors relative to the peaks: 3D flux, "
                 "pressures at the pipe ends (max of the two), profiles of Q "
                 f"along both vessels at t = {END} s.", ""]
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
    suffix = ("_linear" if LINEAR else "") + ("_robin" if ROBIN else "")
    (HERE/f"highOrder{suffix}.md").write_text("\n".join(lines))
    print("\n".join(lines))
    figure(results, suffix)


def figure(results, suffix=""):
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
        ticks = [0.3125, 0.625, 1.25, 2.5]
        labels = ["0.31", "0.62", "1.25", "2.5"]
        if LINEAR:
            ticks, labels = [0.15625] + ticks, ["0.16"] + labels
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
        ax.set_xlabel("dx [mm], with dt = 0.8 ms per mm")
        ax.set_title(title)
    axs[0].set_ylabel("error / peak" + (" (closed form)" if LINEAR else ""))
    axs[2].legend(frameon=False, fontsize=7, loc="lower right", ncol=2)
    fig.tight_layout()
    fig.savefig(HERE/"plots"/f"highOrder3D1D{suffix}.png", dpi=150)


def write_reference(outdir=None):
    """The physiological-amplitude variant's numerical reference: the finest
    nonlinear run (implicit order 4, 6400 cells per vessel, dt = 0.125 ms,
    kExact 3 + RK4 vessels, BDF 4 pipe, inflow x 1) as CSV: the pipe flow
    and pressures at every step, and the profiles at t = 0.35 s (cell
    averages). Its discretisation error is about 4e-7 of the peak (the
    self-convergence of highOrder.md and the linear study's errors at this
    level)."""
    import os
    outdir = outdir or str(HERE/"reference")
    os.makedirs(outdir, exist_ok=True)
    case = RUNS/f"HO_I4_n{REF[0]}"
    if not finished(case):
        sys.exit(f"{case} has not run: python3 highOrder.py run")
    hdr = ("# rigid-segment benchmark at physiological amplitude (inflow Boileau Eq. 14 x 1, nonlinear 1D "
           "sqrtElastic vessels): NUMERICAL reference, the finest run of highOrder.py (implicit coupling "
           f"order 4, {REF[0]} cells per vessel, dt = {REF[1]*1e3:g} ms, kExact 3 + RK4 vessels at CFL 0.5, "
           "BDF 4 pipe of 3072 cells); discretisation error about 4e-7 of the peak\n")
    out = case/"postProcessing"/"circulation1DCoupling"/"0"
    a = np.loadtxt(out/"coupling.dat", comments="#", ndmin=2)
    with open(os.path.join(outdir, "rigidSegment_physiological_pipe.csv"), "w") as f:
        f.write(hdr + "t [s],Q [m3/s],P1 pipe inlet [Pa],P2 pipe outlet [Pa]\n")
        for row in zip(a[:, 0], a[:, 10], a[:, 3], a[:, 5]):
            f.write(",".join(f"{v:.12g}" for v in row) + "\n")
    for v in ("proximal", "distal"):
        b = np.loadtxt(out/"profiles"/f"{END:g}"/f"{v}.dat", comments="#")[1:-1]
        with open(os.path.join(outdir, f"rigidSegment_physiological_{v}_t0.35.csv"), "w") as f:
            f.write(hdr + "x [m] (cell centre),A [m2],Q [m3/s],P [Pa]   (cell averages at t = 0.35 s)\n")
            for row in b[:, :4]:
                f.write(",".join(f"{x:.12g}" for x in row) + "\n")
    print(f"written to {outdir}/")


def main():
    if "--reference" in sys.argv:
        write_reference()
        return
    args = [a for a in sys.argv[1:] if a not in ("--linear", "--robin")]
    if args[:1] == ["run"]:
        nproc = int(args[1]) if len(args) > 1 else 6
        RUNS.mkdir(exist_ok=True)
        cases = [] if LINEAR else [setup("implicit", 4, *REF)]
        for scheme, order in COUPLINGS:
            for n, dt in LEVELS + ([REF] if LINEAR else []):
                cases.append(setup(scheme, order, n, dt))
        cases = [c for c in cases if c is not None]
        print(f"running {len(cases)} cases, {nproc} at a time", flush=True)
        with ThreadPoolExecutor(nproc) as pool:
            list(pool.map(run, cases))
    analyse()


if __name__ == "__main__":
    main()
