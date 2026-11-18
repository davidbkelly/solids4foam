#!/usr/bin/env python3
"""The thoracic aorta II as a 3D segment: the last cycle at the measurement
sites, for the original 1D network, the split network (stubs1D), the
inertance limit (stubsInertance), the coupled pimpleFoam run (rigid3D) and
the in vitro measurements; the cost of the coupling. Writes results.md and
plots/*.png.
"""
import math
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from makeCases import SRC, T_CYCLE, T_STATE, N_CYCLES  # noqa: E402
sys.path.insert(0, str(SRC.parents[1]/"common"))
import boileauMetrics as bm  # noqa: E402

RUNS = HERE/"runs"
T0 = T_STATE + (N_CYCLES - 1)*T_CYCLE                  # the last cycle
T1 = T_STATE + N_CYCLES*T_CYCLE
SITES = ["AorticArchII", "ThoracicAortaII", "RIliacFemoralII"]


def probe(case, name):
    """A probe of a 1D run (its output directory is its start time)"""
    d = RUNS/case/"postProcessing"/"circulation1D"
    if not d.is_dir():
        return None
    for td in d.iterdir():
        f = td/"probes"/f"{name}.dat"
        if f.is_file():
            a = np.loadtxt(f, comments="#")
            return a[:, 0], a[:, 1], a[:, 2]
    return None


def coupled_probe(case, name):
    """A probe of the coupled run, from the coupler's coupling.dat columns"""
    d = RUNS/case/"postProcessing"/"circulation1DCoupling"
    if not d.is_dir():
        return None
    f = next(iter(d.iterdir()))/"coupling.dat"
    names = None
    with open(f) as fh:
        for line in fh:
            if line.startswith("# time"):
                names = line[1:].split()
                break
    if f"{name}.P" not in names:
        return None
    a = np.loadtxt(f, comments="#", ndmin=2)
    col = {n: i for i, n in enumerate(names)}
    return a[:, 0], a[:, col[f"{name}.P"]], a[:, col[f"{name}.Q"]]


def original(name):
    f = SRC/"postProcessing"/"circulation1D"/"0"/"probes"/f"{name}.dat"
    a = np.loadtxt(f, comments="#")
    return a[:, 0], a[:, 1], a[:, 2]


def measured(name):
    p = np.loadtxt(SRC/"referenceData"/f"Exp_{name}_P.txt")
    q = np.loadtxt(SRC/"referenceData"/f"Exp_{name}_Q.txt")
    return p[:, 0], 1e3*p[:, 1], 1e-6*q[:, 1]        # kPa, ml/s


def zone(case, name):
    d = RUNS/case/"postProcessing"
    try:
        a = np.loadtxt(d/f"{name}_P"/f"{T_STATE:g}"/"surfaceFieldValue.dat",
                       comments="#", ndmin=2)
        b = np.loadtxt(d/f"{name}_Q"/f"{T_STATE:g}"/"surfaceFieldValue.dat",
                       comments="#", ndmin=2)
    except OSError:
        return None
    return a[:, 0], 1050.0*a[:, 1], b[:, 1]


def cycle(t, *xs):
    """A run's last complete cycle, shifted to start at 0. Every run's
    cycles start at multiples of the period (the restart time is 18 T), so
    the phases coincide."""
    k = math.floor((t[-1] + 1e-9)/T_CYCLE) - 1
    t0 = k*T_CYCLE
    m = (t >= t0 - 1e-9) & (t <= t0 + T_CYCLE + 1e-9)
    return (t[m] - t0,) + tuple(x[m] for x in xs)


def rms(x):
    return math.sqrt(np.mean(x**2))


BOILEAU_KEYS = [("RMS", "P"), ("MAX", "P"), ("RMS", "Q"), ("MAX", "Q")]


def boileau_errors(site, t, P, Q):
    """Errors [%] against the in vitro data as Boileau et al. (2015) Eqs.
    10-11 measure them: the pressure error relative to the measured pressure
    itself (about 10-15 kPa, not the pulse), the flow error relative to the
    peak measured flow, on the in vitro samples placed with their published
    time shift of the record (as boileauMetrics.analyse_network does)"""
    refdir = SRC/"referenceData"
    exp = {q: bm.load_reference(refdir/f"Exp_{site}_{q}.txt") for q in "PQ"}
    t0_pub = bm.time_shifts(refdir/"timeShifts.txt")[site] - 0.001
    k = math.floor((t[-1] - exp["P"][-1, 0] - t0_pub)/T_CYCLE + 1e-9)
    t0 = t0_pub + k*T_CYCLE
    assert t0 >= t[0] - 1e-9, (site, t0, t[0])
    m = {}
    for q, y in (("P", P/1e3), ("Q", Q*1e6)):           # kPa, ml/s
        m.update(bm.errors(np.interp(exp[q][:, 0] + t0, t, y), exp[q][:, 1], q))
    return m


def main():
    lines = ["# The thoracic aorta II (s17) of the 37-artery network as a 3D segment",
             "", f"Last cycle ({T0:.3f}-{T1:.3f} s), period {T_CYCLE} s.", ""]
    data = {}
    boileau = {}
    for site in SITES:
        d = {}
        raw = {"original": original(site)}
        for case in ("stubs1D", "stubsInertance", "stubsInertance_total"):
            name = site
            if site == "ThoracicAortaII" and case.startswith("stubsInertance"):
                name = "portIn"
            r = probe(case, name)
            if r is not None:
                raw[case] = r
        for case in ("rigid3D", "rigid3D_total_mean", "rigid3D_total_area",
                     "rigid3D_total_flux"):
            r = zone(case, "mid") if site == "ThoracicAortaII" \
                else coupled_probe(case, site)
            if r is not None:
                raw[case] = r
        for case, r in raw.items():
            d[case] = cycle(*r)
        boileau[site] = {case: boileau_errors(site, *raw[case])
                         for case in ("original", "stubsInertance", "rigid3D") if case in raw}
        tm, Pm, Qm = measured(site)
        d["measured"] = (tm, Pm, Qm)
        data[site] = d
        lines += [f"## {site}", "",
                  "| run | peak P [kPa] | min P [kPa] | peak Q [ml/s] | RMS dP vs original / pulse | RMS dQ vs original / peak Q | RMS dP vs measured | RMS dQ vs measured |",
                  "|---|---|---|---|---|---|---|---|"]
        t_o, P_o, Q_o = d["original"]
        pulse = P_o.max() - P_o.min()
        for run, (t, P, Q) in d.items():
            Po = np.interp(t, t_o, P_o)
            Qo = np.interp(t, t_o, Q_o)
            Pm_ = np.interp(t, tm, Pm)
            Qm_ = np.interp(t, tm, Qm)
            lines.append(
                f"| {run} | {P.max()/1e3:.2f} | {P.min()/1e3:.2f} | {Q.max()*1e6:.1f} "
                f"| {rms(P - Po)/pulse:.2e} | {rms(Q - Qo)/Q_o.max():.2e} "
                f"| {rms(P - Pm_)/(Pm.max() - Pm.min()):.2e} | {rms(Q - Qm_)/Qm.max():.2e} |")
        for sfx, isfx in (("", ""), ("_total_mean", "_total"), ("_total_area", "_total"),
                          ("_total_flux", "_total")):
            if "rigid3D" + sfx in d and "stubsInertance" + isfx in d:
                t, P, Q = d["rigid3D" + sfx]
                ti, Pi, Qi = d["stubsInertance" + isfx]
                lines.append(f"\nrigid3D{sfx} against the inertance limit{isfx.replace('_', ' (') + ')' if isfx else ''}: RMS dP / pulse "
                             f"{rms(P - np.interp(t, ti, Pi))/pulse:.2e}, RMS dQ / peak Q "
                             f"{rms(Q - np.interp(t, ti, Qi))/Q_o.max():.2e}")
        for a, b in (("rigid3D_total_mean", "rigid3D"), ("rigid3D_total_area", "rigid3D"),
                     ("rigid3D_total_flux", "rigid3D"), ("stubsInertance_total", "stubsInertance")):
            if a in d and b in d:
                t, P, Q = d[a]
                tb, Pb, Qb = d[b]
                lines.append(f"\n{a} against {b}: RMS dP / pulse "
                             f"{rms(P - np.interp(t, tb, Pb))/pulse:.2e}, RMS dQ / peak Q "
                             f"{rms(Q - np.interp(t, tb, Qb))/Q_o.max():.2e}, "
                             f"max dP / pulse {np.abs(P - np.interp(t, tb, Pb)).max()/pulse:.2e}")
        lines.append("")
    # The same runs in Boileau's measure, beside the range of his six schemes
    lines += ["## Against the in vitro data in Boileau's measure", "",
              "Boileau et al. (2015) Eqs. 10-11 [%]: the pressure error relative to the "
              "measured pressure itself (not the pulse, as in the tables above), the flow "
              "error relative to the peak measured flow, with their time alignment of the "
              "in vitro record. The last row of each site is the range of the six schemes "
              "of their Table X.", "",
              "| site | run | " + " | ".join(f"E^{m}_{q}" for m, q in BOILEAU_KEYS) + " |",
              "|---|---|" + "---|"*len(BOILEAU_KEYS)]
    for site, runs in boileau.items():
        for case, m in runs.items():
            lines.append(f"| {site} | {case} | "
                         + " | ".join(f"{m[k]:.2f}" for k in BOILEAU_KEYS) + " |")
        cols = [bm.TABLE_X_ORDER.index(k) for k in BOILEAU_KEYS]
        lines.append(f"| {site} | six schemes (Boileau Table X) | "
                     + " | ".join(f"{min(bm.TABLE_X[site][s][c] for s in bm.SCHEMES):.2f}-"
                                  f"{max(bm.TABLE_X[site][s][c] for s in bm.SCHEMES):.2f}"
                                  for c in cols) + " |")
    lines.append("")
    # cost
    log = RUNS/"rigid3D"/"log.pimpleFoam"
    if log.is_file():
        txt = log.read_text()
        steps = len(re.findall(r"^Time = ", txt, re.M))
        iters = len(re.findall(r"^PIMPLE: iteration", txt, re.M))
        ex = re.findall(r"ExecutionTime = ([\d.]+)", txt)
        cd = RUNS/"rigid3D"/"postProcessing"/"circulation1DCoupling"/f"{T_STATE:g}"/"coupling.dat"
        upd = np.loadtxt(cd, comments="#", ndmin=2)[:, 1] if cd.is_file() else []
        lines += ["## Cost of the coupled run", "",
                  f"{steps} steps, {iters/max(steps,1):.1f} outer iterations per step, "
                  f"{np.mean(upd) if len(upd) else float('nan'):.1f} network re-integrations per step, "
                  f"{float(ex[-1]) if ex else float('nan'):.0f} s wall time", ""]
    (HERE/"results.md").write_text("\n".join(lines))
    print("\n".join(lines))
    figure(data)


def figure(data):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    (HERE/"plots").mkdir(exist_ok=True)
    style = {"original": ("k", "-", "1D network (original)"),
             "stubs1D": ("#56B4E9", "-", "1D, segment split in three"),
             "stubsInertance": ("#0072B2", "--", "1D + inertance (rigid, inviscid)"),
             "rigid3D": ("#D55E00", "-", "1D + 3D rigid pipe (pimpleFoam)"),
             "measured": ("gray", ":", "in vitro (Matthys et al.)")}
    fig, axs = plt.subplots(2, len(data), figsize=(4.4*len(data), 6.2), squeeze=False)
    for j, (site, d) in enumerate(data.items()):
        for run, (t, P, Q) in d.items():
            if run not in style:
                continue                    # the total-pressure variants
            c, ls, label = style[run]
            axs[0, j].plot(t, P/1e3, color=c, ls=ls, lw=1.6 if run == "rigid3D" else 1.1, label=label)
            axs[1, j].plot(t, Q*1e6, color=c, ls=ls, lw=1.6 if run == "rigid3D" else 1.1)
        axs[0, j].set_title(site + (" (3D mid-plane)" if site == "ThoracicAortaII" else ""), fontsize=10)
        axs[1, j].set_xlabel("t in the cycle [s]")
        for ax in axs[:, j]:
            ax.grid(alpha=0.3)
    axs[0, 0].set_ylabel("P [kPa]")
    axs[1, 0].set_ylabel("Q [ml/s]")
    axs[0, 0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(HERE/"plots"/"sites.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
