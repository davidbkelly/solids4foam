#!/usr/bin/env python3
"""Coupling accuracy and stability of the 3D-1D ports with the rigid-pipe
surrogate (Test-circulation1DPorts), for a pipe integrated by BDF1 to BDF4.

    python3 sweep.py WORKDIR [accuracy] [stability] [--ddt LIST] [--cfl C]
                     [--nproc N]

    --ddt    pipe time schemes, comma-separated: Euler, backward (BDF2), BDF3,
             BDF4 (default backward)
    --cfl    CFL number of the 1D vessels in the accuracy runs (default 0.5;
             0.1 makes their time error, 5e-8 of the peak with SSPRK3, small
             enough for a fourth-order coupling to show)
    --nproc  cases run at a time (default 6)

OpenFOAM and the circulation applications must be on PATH. Runs are made in
WORKDIR; the results go to WORKDIR/results.md and WORKDIR/coupling_<ddt>.png.

With two equal vessels of characteristic impedance Zc = rho c/A on either side
of a rigid pipe of inertance L = rho l/A, each vessel end answers the port
flow as a resistor Zc, so the pipe sees Pin - Pout = 2 Pinc - 2 Zc Q and the
coupling number is z = 2 Zc dt/L = dt/tau, tau = l/(2c). A model problem (the
pipe's BDF formula with the flow extrapolated to t1 by order p) gives the
explicit stability limits of PREDICTED, in dt/tau; the plain fixed-point
iteration of implicit coupling contracts by dt/(a0 tau), a0 the leading BDF
coefficient.

accuracy   explicit and implicit, orders 1-4, for a range of dt, SSPRK3 in
           the vessels: the error of the pipe flow (max over time relative to
           its peak) against the implicit coupling of the pipe's order at
           dt = 1/8 of the finest, and the difference from the monolithic run
           (inertance node, same 1D mesh and CFL)
stability  explicit orders 1-4 at 0.9, 1.1 and 1.25 times the predicted
           limit: the measured amplification per step against the model's
           root; and the contraction rate of the implicit iteration
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
TUT = HERE.parents[2]/"tutorials"/"circulationModels"/"coupling3D1D"/\
    "rigidSegmentInertance"

RHO = 1050.0
A0 = math.pi*0.01**2
BETA = 4.0/3.0*math.sqrt(math.pi)*400e3*1.5e-3
C0 = math.sqrt(BETA/(2*RHO*A0))*A0**0.25
TAU = 0.08/(2*C0)

DTS = [8e-3, 4e-3, 2e-3, 1e-3, 5e-4, 2.5e-4]

# Leading BDF coefficient, order, and the explicit limits of orders 1-4 in
# dt/tau from the model problem (roots of the step polynomial)
BDF = {
    "Euler":    (1.0, 1, (2.0, 0.6667, 0.2857, 0.1333)),
    "backward": (1.5, 2, (4.0, 1.3333, 0.5714, 0.2667)),
    "BDF3":     (11/6, 3, (4.5, 2.2222, 0.9524, 0.4444)),
    "BDF4":     (25/12, 4, (1.6428, 2.25, 1.5238, 0.7111)),
}
NPROC = 6


def predicted(order, ddt):
    return BDF[ddt][2][order - 1]


def pmap(f, items):
    with ThreadPoolExecutor(NPROC) as pool:
        return list(pool.map(f, items))


def setup(src, dst, edits):
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst, ignore=shutil.ignore_patterns(
        "postProcessing", "order", "plots", "log*", "results.md", "0.*",
        "[1-9]*"))
    for rel, pairs in edits.items():
        f = dst/rel
        s = f.read_text()
        for pattern, value in pairs:
            s, n = re.subn(pattern, value, s, count=1, flags=re.M)
            if n != 1:
                sys.exit(f"setup: no match for {pattern} in {f}")
        f.write_text(s)


def run(case, app):
    with open(case/"log", "w") as log:
        r = subprocess.run([app], cwd=case, stdout=log, stderr=subprocess.STDOUT)
    return r.returncode


def ports_case(work, name, scheme, order, dt, ddt, end=0.6, rk="SSPRK2",
               cfl=0.5):
    case = work/name
    setup(HERE/"case", case, {
        "system/controlDict": [
            (r"^deltaT\s+\S+;", f"deltaT          {dt!r};"),
            (r"^endTime\s+\S+;", f"endTime         {end!r};"),
            (r"^(\s+scheme\s+)\S+;", rf"\g<1>{scheme};"),
            (r"^(\s+order\s+)\S+;", rf"\g<1>{order};"),
        ],
        "constant/circulation1DProperties": [
            (r"^(\s+ddtScheme\s+)\S+;", rf"\g<1>{ddt};"),
            (r"^(\s+timeScheme\s+)\S+;", rf"\g<1>{rk};"),
            (r"^(\s+CFL\s+)\S+;", rf"\g<1>{cfl};"),
        ],
    })
    # Stop at 1000 x the expected peak flow
    s = (case/"system/controlDict").read_text()
    s = s.replace("tolerance   1e-12;", "tolerance   1e-12;\n    maxFlow     1e-6;")
    (case/"system/controlDict").write_text(s)
    status = run(case, "Test-circulation1DPorts")
    data = None
    f = case/"postProcessing"/"surrogate3D.dat"
    if f.is_file() and f.stat().st_size > 0:
        data = np.loadtxt(f, comments="#", ndmin=2)
    return status, data


def reference(work, cfl):
    case = work/f"reference_cfl{cfl:g}"
    if not (case/"postProcessing").is_dir():
        setup(TUT, case, {
            "system/controlDict": [
                (r"^deltaT\s+\S+;", "deltaT          2.5e-4;")],
            "constant/circulation1DProperties": [
                (r"^(\s+timeScheme\s+)\S+;", r"\g<1>SSPRK3;"),
                (r"^(\s+CFL\s+)\S+;", rf"\g<1>{cfl};")],
        })
        if run(case, "circulation1DFoam") != 0:
            sys.exit("reference run failed: see " + str(case/"log"))
    a = np.loadtxt(case/"postProcessing"/"circulation1D"/"0"/"nodes"/
                   "rigid.dat", comments="#")
    return a[:, 0], a[:, 1]


def accuracy(work, ddt, cfl):
    """Every scheme and order at every dt, against the implicit coupling of
    the pipe's own order at dt/8"""
    tr, Qr = reference(work, cfl)
    q = BDF[ddt][1]
    tag = f"{ddt}_cfl{cfl:g}"
    status, d = ports_case(work, f"fine_{tag}", "implicit", max(q, 2),
                           DTS[-1]/8, ddt, rk="SSPRK3", cfl=cfl)
    if status != 0 or d is None:
        sys.exit("fine reference failed: see " + str(work/f"fine_{tag}"/"log"))
    tf, Qf = d[:, 0], d[:, 1]
    peak = Qr.max()

    jobs = [(s, o, dt) for s in ("explicit", "implicit") for o in (1, 2, 3, 4)
            for dt in DTS]

    def one(job):
        s, o, dt = job
        name = f"acc_{tag}_{s}{o}_dt{dt:g}"
        status, d = ports_case(work, name, s, o, dt, ddt, rk="SSPRK3",
                               cfl=cfl)
        if status != 0 or d is None:
            return (dt, float("nan"), float("nan"), 0)
        e = np.abs(d[:, 1] - np.interp(d[:, 0], tf, Qf)).max()/peak
        m = np.abs(d[:, 1] - np.interp(d[:, 0], tr, Qr)).max()/peak
        print(f"{name}: error {e:.3e}, from monolithic {m:.3e}, "
              f"iterations {d[:, 4].mean():.1f}", flush=True)
        return (dt, e, m, d[:, 4].mean())

    res = dict(zip(jobs, pmap(one, jobs)))
    return {(s, o): [res[(s, o, dt)] for dt in DTS]
            for s in ("explicit", "implicit") for o in (1, 2, 3, 4)}


COEFFS = {"Euler": [1, -1], "backward": [1.5, -2, 0.5],
          "BDF3": [11/6, -3, 1.5, -1/3], "BDF4": [25/12, -4, 3, -4/3, 0.25]}
EXTRAP = {1: [1], 2: [2, -1], 3: [3, -3, 1], 4: [4, -6, 4, -1]}


def model_factor(z, ddt, order):
    """The model problem's largest root at dt/tau = z: the per-step
    amplification of explicit coupling of this order"""
    a, e = COEFFS[ddt], EXTRAP[order]
    c = np.zeros(max(len(a), len(e) + 1))
    c[:len(a)] += a
    for k, ek in enumerate(e):
        c[k + 1] += z*ek
    return float(np.abs(np.roots(c)).max())


def factor(work, order, ddt, x):
    """Explicit run at dt = x tau to t = 0.9 s: the amplification per step
    of the step-to-step oscillation of the flow after the pulse, a
    least-squares fit of log|second difference of Q| from 0.45 s (the pulse
    has passed the pipe) to 0.9 s (the pulse reflected at the inlet returns
    at about 1 s), with only the values above 1e-13 of the peak flow (below
    that a decaying mode has reached round-off)."""
    status, d = ports_case(work, f"stab_{ddt}{order}_{x:.6f}", "explicit",
                           order, x*TAU, ddt, end=0.9)
    if d is None or len(d) < 4:
        return float("inf")
    t, Q = d[:, 0], d[:, 1]
    if not np.all(np.isfinite(Q)) or status != 0:
        return float("inf")
    d2 = np.abs(Q[2:] - 2*Q[1:-1] + Q[:-2])
    tt = t[1:-1]
    m = (tt > 0.45) & (tt < 0.9)
    k, y = np.arange(len(tt))[m], d2[m]
    keep = y > 1e-13*np.abs(Q).max()
    if keep.sum() < 3:
        return 0.0
    return float(np.exp(np.polyfit(k[keep], np.log(y[keep]), 1)[0]))


def contraction(work, ddt, x):
    """Median contraction rate of the implicit iteration while the pulse
    passes the pipe"""
    case = f"rate_{ddt}_{x:.4f}"
    status, d = ports_case(work, case, "implicit", 2, x*TAU, ddt)
    if d is None:
        return float("nan")
    Q = np.abs(d[:, 1])
    busy = (Q > 0.05*Q.max()) & (d[:, 5] > 0)
    return float(np.median(d[busy, 5])) if busy.any() else float("nan")


def stability(work, ddts):
    """Explicit orders 1-4 at 0.9, 1.1 and 1.25 times the predicted limit:
    the measured amplification per step against the model's root, which
    crosses 1 at the limit"""
    jobs = [(ddt, o, x) for ddt in ddts for o in (1, 2, 3, 4)
            for x in (0.9, 1.1, 1.25)]

    def one(j):
        ddt, o, x = j
        z = x*predicted(o, ddt)
        g = factor(work, o, ddt, z)
        r = model_factor(z, ddt, o)
        print(f"explicit {ddt} order {o} at {x}x the limit: factor {g:.4f}, "
              f"model {r:.4f}", flush=True)
        return (z, r, g)

    res = dict(zip(jobs, pmap(one, jobs)))

    rjobs = [(ddt, x) for ddt in ddts for x in (0.25, 0.5, 0.75, 0.9)]

    def rate(j):
        ddt, x = j
        k = 1/BDF[ddt][0]
        r = contraction(work, ddt, x)
        print(f"implicit {ddt} dt/tau {x}: contraction {r:.4f}, "
              f"predicted {k*x:.4f}", flush=True)
        return (k*x, r)

    rates = dict(zip(rjobs, pmap(rate, rjobs)))
    return res, rates


def option(name, default):
    if name in sys.argv:
        i = sys.argv.index(name)
        value = sys.argv[i + 1]
        del sys.argv[i:i + 2]
        return value
    return default


def main():
    global NPROC
    ddts = option("--ddt", "backward").split(",")
    cfl = float(option("--cfl", "0.5"))
    NPROC = int(option("--nproc", "6"))
    for d in ddts:
        if d not in BDF:
            sys.exit(f"unknown ddt {d}: {', '.join(BDF)}")
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    work = Path(sys.argv[1]).resolve()
    work.mkdir(parents=True, exist_ok=True)
    what = sys.argv[2:] or ["accuracy", "stability"]

    lines = ["# 3D-1D ports with a rigid-pipe surrogate", "",
             f"c0 = {C0:.4f} m/s, pipe 8 cm, tau = l/(2 c0) = "
             f"{TAU*1e3:.4f} ms; 1D dx = 1.25 mm.", ""]

    if "accuracy" in what:
        for ddt in ddts:
            acc = accuracy(work, ddt, cfl)
            q = BDF[ddt][1]
            lines += [f"## Accuracy, pipe {ddt} (order {q}), 1D CFL {cfl:g}",
                      "",
                      "Error: max |Q - Qfine| / peak, Qfine the implicit "
                      f"order {max(q, 2)} coupling at dt = "
                      f"{DTS[-1]/8*1e3:g} ms.", "",
                      "| dt [ms] | dt/tau | " + " | ".join(
                          f"{s[:3]} {o}" for s in ("explicit", "implicit")
                          for o in (1, 2, 3, 4)) + " |",
                      "|---|---|" + "---|"*8]
            for i, dt in enumerate(DTS):
                cells = [f"{dt*1e3:g}", f"{dt/TAU:.3f}"]
                for s in ("explicit", "implicit"):
                    for o in (1, 2, 3, 4):
                        e = acc[(s, o)][i][1]
                        cells.append("unstable" if not e < 1 else f"{e:.2e}")
                lines.append("| " + " | ".join(cells) + " |")
            lines += ["", "Observed order between the two finest dt:", ""]
            for s in ("explicit", "implicit"):
                lines.append(f"- {s}: " + ", ".join(
                    f"order {o} " + (
                        f"{math.log(acc[(s, o)][-2][1]/acc[(s, o)][-1][1])/math.log(DTS[-2]/DTS[-1]):.2f}"
                        if acc[(s, o)][-1][1] > 0 else "-")
                    for o in (1, 2, 3, 4)))
            lines += ["", "Difference from the monolithic run (inertance "
                      f"node, same 1D mesh and CFL) at dt = {DTS[-1]*1e3:g} "
                      "ms:", ""]
            for s in ("explicit", "implicit"):
                lines.append(f"- {s}: " + ", ".join(
                    f"order {o} {acc[(s, o)][-1][2]:.2e}"
                    for o in (1, 2, 3, 4)))
            lines.append("")
            figure(work, ddt, acc)

    if "stability" in what:
        stab, rates = stability(work, ddts)
        lines += ["## Explicit stability", "",
                  "Amplification per step of the step-to-step oscillation "
                  "after the pulse, measured and from the model (its largest "
                  "root, which crosses 1 at the predicted limit); 0 means it "
                  "decayed to round-off.", "",
                  "| pipe ddt | order | predicted limit dt/tau | x limit | "
                  "dt/tau | model | measured | result |",
                  "|---|---|---|---|---|---|---|---|"]
        for (ddt, o, x), (z, r, g) in stab.items():
            res = "stable" if g < 1 else "unstable"
            lines.append(f"| {ddt} | {o} | {predicted(o, ddt):.4f} | {x} | "
                         f"{z:.4f} | {r:.4f} | {g:.4f} | {res} |")
        lines += ["", "## Implicit iteration: contraction rate per iteration",
                  "", "| pipe ddt | dt/tau | predicted | measured |",
                  "|---|---|---|---|"]
        for (ddt, x), (p, r) in rates.items():
            lines.append(f"| {ddt} | {x} | {p:.4f} | {r:.4f} |")
        lines.append("")

    (work/"results.md").write_text("\n".join(lines))
    print("\n".join(lines))


def figure(work, ddt, acc):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    q = BDF[ddt][1]
    fig, axs = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    x = np.array(DTS)/TAU
    for ax, s in zip(axs, ("explicit", "implicit")):
        for o in (1, 2, 3, 4):
            e = np.array([r[1] for r in acc[(s, o)]])
            e = np.where(e < 1, e, np.nan)
            ax.loglog(x, e, "o-", ms=4, label=f"order {o}")
        e = np.array([r[1] for r in acc[(s, q)]])
        ax.loglog(x, e[-1]*(x/x[-1])**q, "k:", lw=1, label=f"slope {q}")
        ax.set_xlabel("dt / tau")
        ax.set_title(f"{s} coupling, {ddt} pipe")
        ax.legend(frameon=False, fontsize=8)
    axs[0].set_ylabel("max |Q - Q_fine| / peak")
    fig.tight_layout()
    fig.savefig(work/f"coupling_{ddt}.png", dpi=150)


if __name__ == "__main__":
    main()
