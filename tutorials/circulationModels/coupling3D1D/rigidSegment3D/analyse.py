#!/usr/bin/env python3
"""The rigid-segment 3D-1D benchmark: the coupled runs against the closed form.

    python3 analyse.py

A rigid 8 cm pipe (pimpleFoam, slip walls) between two 1D vessels, coupled at
both ends by circulation1DCoupling. Its exact solution is plug flow: an
inertance L = rho l/A3D in series, so for an incident pulse
Pinc(T) = P exp(-T^2/(2 sigma^2)) at the pipe the transmitted and reflected
pulses are (src/circulationModels/docs/BENCHMARKS3D1D.md)

    tau dPt/dT + Pt = Pinc,  Pr = Pinc - Pt,  tau = L/(2 Zc) = 6.52 ms.

From whatever has run in runs/:

  1. the 3D flow and the pressures at both ends against the closed form, for
     every accuracy case, and the monolithic 1D run (inertance node) too;
  2. the error in dt against the finest coupled run ("fine"), and its order,
     for explicit orders 1-4 and implicit coupling;
  3. stability: explicit orders at 0.9 and 1.1 times the predicted limits,

for each CFD time scheme that has run: backward (BDF2), BDF 3, BDF 4 (see
sweep.py).

Writes results.md and plots/rigidSegment3D.png. Needs numpy and scipy
(matplotlib for the figure).
"""
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[3]/"src"/"circulationModels"/"docs"))
sys.path.insert(0, str(HERE))

from benchmarks3D1D import transmitted  # noqa: E402
from sweep import (A0, A3D, C0, DTS, FINE_DT, TAU, Z0, STUDY,  # noqa: E402
                   limit)

RUNS = HERE/"runs"
SIGMA = 1/math.sqrt(2e4)
T_ARRIVE = 0.05 + 1.0/C0            # incident peak at the pipe
P_INC = Z0*1e-3*1e-6                # incident peak [Pa]


def exact(t):
    """Q through the pipe, P at its inlet and outlet"""
    T = t - T_ARRIVE
    pt = transmitted(T, P_INC, SIGMA, TAU)
    pi = P_INC*np.exp(-T**2/(2*SIGMA**2))
    return pt/Z0, 2*pi - pt, pt


def finished(name):
    log = RUNS/name/"log.pimpleFoam"
    return log.is_file() and log.read_text().rstrip().endswith("End")


def load3D(name):
    """The coupler's rows: t updates to3D.Q to3D.P from3D.Q from3D.P
    proximalMid.P proximalMid.Q distalMid.P distalMid.Q to3D.Qsampled
    from3D.Qsampled. The 3D flux at t is the last flux sampled in the step
    ending at t (implicit); explicit coupling samples once, at the start of
    the step, so the flux at t is the sample of the next row."""
    f = RUNS/name/"postProcessing"/"circulation1DCoupling"/"0"/"coupling.dat"
    if not f.is_file() or f.stat().st_size == 0:
        return None
    a = np.loadtxt(f, comments="#", ndmin=2)
    d = {"t": a[:, 0], "updates": a[:, 1], "Q": a[:, 10],
         "Qout": a[:, 11], "P1": a[:, 3], "P2": a[:, 5]}
    if "explicit" in name or "stab_" in name:
        d = {k: v[:-1] for k, v in d.items()}
        d["Q"], d["Qout"] = a[1:, 10], a[1:, 11]
    return d


def loadMonolithic(cfl=0.5):
    name = "monolithic" if cfl == 0.5 else f"monolithic_cfl{cfl:g}"
    f = RUNS/name/"postProcessing"/"circulation1D"/"0"/"nodes"/"rigid.dat"
    if not f.is_file():
        return None
    a = np.loadtxt(f, comments="#")
    return {"t": a[:, 0], "Q": a[:, 1], "P1": a[:, 2], "P2": a[:, 3]}


def closed_form_errors(d):
    Q, P1, P2 = exact(d["t"])
    return (np.abs(d["Q"] - Q).max()*Z0/P_INC,
            np.abs(d["P1"] - P1).max()/P_INC,
            np.abs(d["P2"] - P2).max()/P_INC)


def growth(d):
    """Amplification per step of the step-to-step oscillation of the flow
    after the pulse (0.30-0.50 s; the pulse reflected at the inlet returns
    at 0.53 s): a least-squares fit of log|second difference of Q| over the
    steps, using only values above 1e-9 of their largest (below that a
    decaying mode has reached round-off). Returns the factor and the number
    of steps fitted."""
    t, Q = d["t"], d["Q"]
    d2 = np.abs(Q[2:] - 2*Q[1:-1] + Q[:-2])
    tt = t[1:-1]
    m = (tt > 0.30) & (tt < 0.50)
    k = np.arange(len(tt))[m]
    y = d2[m]
    keep = y > 1e-9*y.max()
    if keep.sum() < 3:
        return float("nan"), int(keep.sum())
    return float(np.exp(np.polyfit(k[keep], np.log(y[keep]), 1)[0])), \
        int(keep.sum())


def model_factor(z, study, order):
    """The model problem's largest root at dt/tau = z: the per-step
    amplification of the explicit coupling of this order"""
    a = {"backward": [1.5, -2, 0.5], "BDF3": [11/6, -3, 1.5, -1/3],
         "BDF4": [25/12, -4, 3, -4/3, 0.25]}[study]
    e = {1: [1], 2: [2, -1], 3: [3, -3, 1], 4: [4, -6, 4, -1]}[order]
    c = np.zeros(max(len(a), len(e) + 1))
    c[:len(a)] += a
    for k, ek in enumerate(e):
        c[k + 1] += z*ek
    return float(np.abs(np.roots(c)).max())


def ok(name):
    d = load3D(name)
    return d if d is not None and finished(name) else None


def fmt(e):
    return f"{e:.2e}" if e < 1 else "unstable"


def main():
    lines = ["# Rigid segment, 3D-1D: pimpleFoam pipe coupled to 1D vessels", "",
             f"c0 = {C0:.4f} m/s, Zc = {Z0:.4e} Pa s/m3; pipe 8 cm, "
             f"A3D = {A3D:.6e} m2 ({A3D/A0:.4f} pi R^2: the 32-sided patch), "
             f"tau = {TAU*1e3:.4f} ms; incident peak {P_INC*1e3:.4f} mPa "
             "(inflow x 1e-3); 1D dx = 1.25 mm, SSPRK3 (CFL 0.5 with "
             "backward, 0.1 with BDF 3 and BDF 4).", ""]
    fig = {"acc": {}}

    # ---- 1. Against the closed form ---------------------------------------
    rows = []
    for cfl in (0.5, 0.1):
        mono = loadMonolithic(cfl)
        if mono:
            rows.append((f"monolithic 1D, CFL {cfl:g}",) +
                        closed_form_errors(mono))
    for study, st in STUDY.items():
        q, pre = st["q"], st["prefix"]
        for s, o in st["accuracy"]:
            if o != max(q, 2) and not (s == "explicit" and o == q):
                continue
            for dt in DTS:
                d = ok(f"{pre}acc_{s}{o}_dt{dt:g}")
                if d is not None:
                    rows.append((f"{study}: {s} {o}, dt {dt*1e3:g} ms",) +
                                closed_form_errors(d))
        d = ok(pre + "fine")
        if d is not None:
            rows.append((f"{study}: implicit {max(q, 2)}, dt "
                         f"{FINE_DT*1e3:g} ms",) + closed_form_errors(d))
    if not rows:
        sys.exit("No finished runs: run ./Allrun or ./Allrun-sweep first")
    lines += ["## 1. Against the closed form (max error / incident peak)", "",
              "| run | Q (3D flux) | P1, pipe inlet | P2, pipe outlet |",
              "|---|---|---|---|"]
    lines += ["| " + " | ".join([r[0]] + [fmt(e) for e in r[1:]]) + " |"
              for r in rows]
    lines += ["", "The monolithic rows are the 1D discretisation error alone "
              "(dx = 1.25 mm); a coupled run cannot do better than them. Only "
              "the coupling of the time scheme's own order is listed.", ""]

    show = load3D("acc_implicit2_dt0.0005")
    if show is not None:
        fig["show"] = show
        qd = show["Q"] - show["Qout"]
        lines += [f"Mass through the rigid pipe (implicit 2, 0.5 ms): "
                  f"max |Qin - Qout| = {np.abs(qd).max()*Z0/P_INC:.1e} of the "
                  f"peak; outer iterations {show['updates'].mean():.1f} per "
                  "step on average.", ""]

    # ---- 2. Convergence in dt against the fine run -------------------------
    for study, st in STUDY.items():
        q, pre, couplings = st["q"], st["prefix"], st["accuracy"]
        ref = ok(pre + "fine")
        if ref is None:
            continue
        acc = {}
        for s, o in couplings:
            errs = []
            for dt in DTS:
                d = ok(f"{pre}acc_{s}{o}_dt{dt:g}")
                if d is None:
                    errs.append(float("nan"))
                    continue
                Qi = np.interp(d["t"], ref["t"], ref["Q"])
                errs.append(np.abs(d["Q"] - Qi).max()/ref["Q"].max())
            acc[(s, o)] = errs
        fig["acc"][study] = acc
        lines += [f"## 2. Convergence in dt, CFD {st['ddt']} "
                  "(max |Q - Q_fine| / peak)", "",
                  f"Q_fine: implicit order {max(q, 2)} at dt = "
                  f"{FINE_DT*1e3:g} ms.", "",
                  "| dt [ms] | dt/tau | " +
                  " | ".join(f"{s} {o}" for s, o in couplings) + " |",
                  "|---|---|" + "---|"*len(couplings)]
        for i, dt in enumerate(DTS):
            lines.append(f"| {dt*1e3:g} | {dt/TAU:.3f} | " + " | ".join(
                fmt(acc[k][i]) if np.isfinite(acc[k][i]) else "-"
                for k in couplings) + " |")
        lines += ["", "Order between the two finest dt, and least squares over "
                  "the three finest:", ""]
        for k in couplings:
            e = np.array(acc[k])
            if np.all(np.isfinite(e[-3:])) and np.all(e[-3:] < 1):
                o2 = math.log(e[-2]/e[-1])/math.log(DTS[-2]/DTS[-1])
                ls = np.polyfit(np.log(DTS[-3:]), np.log(e[-3:]), 1)[0]
                lines.append(f"- {k[0]} {k[1]}: {o2:.2f} ({ls:.2f})")
        mono = loadMonolithic(st["cfl"])
        if mono:
            Qm = np.interp(ref["t"], mono["t"], mono["Q"])
            lines += ["", "The fine coupled run against the monolithic 1D run "
                      f"(CFL {st['cfl']:g}): "
                      f"{np.abs(ref['Q'] - Qm).max()/ref['Q'].max():.2e} of "
                      "the peak."]
        lines.append("")

    # ---- 3. Stability -------------------------------------------------------
    st_rows = []
    for study, st in STUDY.items():
        for o in st["stab"]:
            for x in (0.9, 1.1):
                name = f"{st['prefix']}stab_{o}_{x:g}"
                d = load3D(name)
                if d is not None:
                    z = x*limit(o, study)
                    st_rows.append((st["ddt"], o, limit(o, study), z,
                                    model_factor(z, study, o), growth(d),
                                    finished(name)))
    if st_rows:
        lines += ["## 3. Explicit stability", "",
                  "Amplification per step of the step-to-step oscillation of "
                  "the flow after the pulse: predicted (the model's largest "
                  "root) and measured (a fit over 0.30-0.50 s, the steps "
                  "fitted in brackets).", "",
                  "| CFD ddt | order | predicted limit dt/tau | run at dt/tau "
                  "| predicted factor | measured factor | result |",
                  "|---|---|---|---|---|---|---|"]
        for ddt, o, lim, z, r, (g, n), done in st_rows:
            res = "stable" if g < 1 and done else "unstable"
            lines.append(f"| {ddt} | {o} | {lim:.4f} | {z:.4f} | {r:.4f} | "
                         f"{g:.4f} ({n}) | {res} |")
        lines.append("")

    (HERE/"results.md").write_text("\n".join(lines))
    print("\n".join(lines))
    figure(fig)


def figure(fig):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.ticker
    except ImportError:
        return
    if "show" not in fig:
        return
    (HERE/"plots").mkdir(exist_ok=True)
    n = 1 + len(fig["acc"])
    f, axs = plt.subplots(1, n, figsize=(5.5*n, 4))
    axs = np.atleast_1d(axs)

    d = fig["show"]
    ms = (d["t"] - T_ARRIVE)*1e3
    Q, P1, P2 = exact(d["t"])
    ax = axs[0]
    ax.plot(ms, d["P1"]/P_INC, color="#2a6fb0", label="P1, 3D inlet")
    ax.plot(ms, d["P2"]/P_INC, color="#c0392b", label="P2, 3D outlet")
    ax.plot(ms, d["Q"]*Z0/P_INC, ".", ms=3, color="#555555",
            label="3D flux x Zc")
    ax.plot(ms, P1/P_INC, "k--", lw=1, label="closed form")
    ax.plot(ms, P2/P_INC, "k--", lw=1)
    ax.set_xlim(-40, 60)
    ax.set_xlabel("time from the incident peak at the pipe [ms]")
    ax.set_ylabel("P / incident peak")
    ax.set_title("pimpleFoam pipe, implicit 2, dt = 0.5 ms")
    ax.legend(frameon=False, fontsize=8)

    x = np.array(DTS)/TAU
    for ax, (study, acc) in zip(axs[1:], fig["acc"].items()):
        q = STUDY[study]["q"]
        for (s, o), e in acc.items():
            e = np.where(np.array(e) < 1, e, np.nan)     # unstable: not shown
            ax.loglog(x, e, "o-" if s == "explicit" else "s--", ms=4,
                      label=f"{s} {o}")
        e = np.array(acc[("implicit", max(q, 2))])
        ax.loglog(x, e[-1]*(x/x[-1])**q, "k:", lw=1, label=f"slope {q}")
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlabel("dt / tau")
        ax.set_ylabel("max |Q - Q_fine| / peak")
        ax.set_title(f"CFD {STUDY[study]['ddt']}")
        ax.legend(frameon=False, fontsize=7)

    f.tight_layout()
    f.savefig(HERE/"plots"/"rigidSegment3D.png", dpi=150)


if __name__ == "__main__":
    main()
