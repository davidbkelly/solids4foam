#!/usr/bin/env python3
"""Compare the coupled 3D-1D tube with the full FSI tube (Formaggia et al.
2007, section 5) and with the 1D vessel alone.

    python3 analyse.py

At each section z the mean pressure and the flow are taken from the faceZone
function objects of the 3D runs (quarter tube: the flow is multiplied by 4,
the kinematic pressure by rho), from the coupler's log of the coupled runs
at the interface z = 5 cm (the port) and beyond it (the vessel probes), and
from the probes of the 1D vessel alone. Every runs/coupled* directory is a
coupled run (makeCases.py --c0 names them by their wave speed). The error of
each coupled run is measured against full3D section by section, the wave
speed is fitted to the arrival of the half-peak of the pressure at the
sections, and the wall displacement at z = 2.5 cm is compared. Writes
results.md and plots/*.png.
"""
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from makeCases import L3D, PLANES_CM, PULSE, RHO, zone_name  # noqa: E402

RUNS = HERE/"runs"
FLOW_SCALE = 4                        # quarter tube


def finished(case, log):
    f = RUNS/case/log
    return f.is_file() and "\nEnd\n" in f.read_text()[-2000:]


def read_fo(case, name):
    f = RUNS/case/"postProcessing"/"fluid"/name/"0"/"surfaceFieldValue.dat"
    a = np.loadtxt(f, comments="#", ndmin=2)
    return a[:, 0], a[:, 1]


def sections_3d(case):
    """{z_cm: (t, P [Pa], Q [m3/s])} from the faceZones of a 3D run"""
    out = {}
    for zcm in PLANES_CM:
        n = zone_name(zcm)
        try:
            t, p = read_fo(case, n + "_P")
            _, q = read_fo(case, n + "_Q")
        except OSError:
            continue
        out[zcm] = (t, RHO*p, FLOW_SCALE*q)
    return out


def sections_1d(case):
    """{z_cm: (t, P, Q)} from the probes of a 1D-only run"""
    out = {}
    d = RUNS/case/"postProcessing"/"circulation1D"/"0"/"probes"
    for zcm in PLANES_CM:
        f = d/(zone_name(zcm) + ".dat")
        if f.is_file():
            a = np.loadtxt(f, comments="#")
            out[zcm] = (a[:, 0], a[:, 1], a[:, 2])
    return out


def sections_coupler(case):
    """{z_cm: (t, P, Q)} at the port (z = L3D) and the vessel probes, from
    the coupler's coupling.dat"""
    f = RUNS/case/"postProcessing"/"circulation1DCoupling"/"0"/"coupling.dat"
    names = None
    with open(f) as fh:
        for line in fh:
            if line.startswith("# time"):
                names = line[1:].split()
                break
    a = np.loadtxt(f, comments="#", ndmin=2)
    col = {n: i for i, n in enumerate(names)}
    t = a[:, 0]
    out = {}
    for zcm in PLANES_CM:
        if zcm*1e-2 == L3D:
            out[zcm] = (t, a[:, col["from3D.P"]], a[:, col["from3D.Qsampled"]])
        elif zone_name(zcm) + ".P" in col:
            n = zone_name(zcm)
            out[zcm] = (t, a[:, col[n + ".P"]], a[:, col[n + ".Q"]])
    return out


def displacement(case):
    f = RUNS/case/"postProcessing"/"0"/"solidPointDisplacement_displacement.dat"
    a = np.loadtxt(f, comments="#")
    return a[:, 0], a[:, 2]                       # t, Dy (radial at (0, R, z))


def arrival(t, P, level=0.5*PULSE):
    """First time the pressure reaches level, by linear interpolation"""
    i = np.argmax(P >= level)
    if i == 0:
        return math.nan
    return t[i - 1] + (level - P[i - 1])/(P[i] - P[i - 1])*(t[i] - t[i - 1])


def wave_speed(secs):
    zs, ts = [], []
    for zcm, (t, P, _) in sorted(secs.items()):
        ta = arrival(t, P)
        if not math.isnan(ta):
            zs.append(zcm*1e-2)
            ts.append(ta)
    if len(zs) < 2:
        return math.nan, zs, ts
    c, _ = np.polyfit(ts, zs, 1)
    return c, zs, ts


def rms(t, x, tref, xref):
    """RMS of x - xref over the times both cover"""
    m = t <= tref[-1]
    xi = np.interp(t[m], tref, xref)
    return math.sqrt(np.mean((x[m] - xi)**2))


def c0_of(case):
    """The 1D wave speed a coupled run was made with, from its dictionary"""
    f = RUNS/case/"constant"/"fluid"/"circulation1DProperties"
    for line in f.read_text().splitlines():
        if "c0 =" in line:
            return float(line.split("c0 =")[1].split()[0])
    return math.nan


def main():
    coupled = sorted(p.name for p in RUNS.glob("coupled*") if p.is_dir())
    lines = ["# Formaggia et al. (2007), section 5: the pulse tube continued "
             "by a 1D vessel", ""]
    data = {}
    if (RUNS/"full3D"/"postProcessing"/"fluid").is_dir():
        if not finished("full3D", "log.solids4Foam"):
            lines += ["**full3D has not finished: partial data.**", ""]
        data["full3D"] = sections_3d("full3D")
        c3, zs, ts = wave_speed(data["full3D"])
        lines += [f"3D wave speed (half-peak arrival at z = "
                  f"{', '.join(f'{z*100:g}' for z in zs)} cm): {c3:.3f} m/s", ""]
    for c in coupled:
        if not (RUNS/c/"postProcessing"/"fluid").is_dir():
            continue
        if not finished(c, "log.solids4Foam"):
            lines += [f"**{c} has not finished: partial data.**", ""]
        s = sections_3d(c)
        s.update(sections_coupler(c))
        data[c] = s
    if finished("only1D", "log.circulation1DFoam"):
        data["only1D"] = sections_1d("only1D")
        c1, _, _ = wave_speed(data["only1D"])
        lines += [f"1D wave speed, same measure (1D alone): {c1:.3f} m/s", ""]

    if "full3D" in data:
        for c in coupled:
            if c not in data:
                continue
            lines += [f"## {c}: 1D c0 = {c0_of(c):.3g} m/s", "",
                      "| z [cm] | where | peak P full3D / coupled / 1D alone [Pa] "
                      "| arrival full3D / coupled / 1D alone [ms] "
                      "| RMS(P coupled - full3D) / pulse "
                      "| RMS(Q coupled - full3D) / peak Q |",
                      "|---|---|---|---|---|---|"]
            for zcm in PLANES_CM:
                if zcm not in data["full3D"]:
                    continue
                tf, Pf, Qf = data["full3D"][zcm]
                where = ("3D" if zcm*1e-2 < L3D
                         else "port" if zcm*1e-2 == L3D else "1D")
                cells = [f"{Pf.max():.0f}"]
                arr = [f"{1e3*arrival(tf, Pf):.2f}"]
                err = ["-", "-"]
                for cc in (c, "only1D"):
                    if cc in data and zcm in data[cc]:
                        t, P, Q = data[cc][zcm]
                        cells.append(f"{P.max():.0f}")
                        arr.append(f"{1e3*arrival(t, P):.2f}")
                        if cc == c:
                            err = [f"{rms(t, P, tf, Pf)/PULSE:.2e}",
                                   f"{rms(t, Q, tf, Qf)/max(Qf.max(), 1e-30):.2e}"]
                    else:
                        cells.append("-")
                        arr.append("-")
                lines.append(f"| {zcm:g} | {where} | {' / '.join(cells)} | "
                             f"{' / '.join(arr)} | {err[0]} | {err[1]} |")
            tf, Df = displacement("full3D")
            tc, Dc = displacement(c)
            lines += ["", f"Radial wall displacement at z = 2.5 cm: peak "
                      f"{1e6*Df.max():.2f} um (full3D), {1e6*Dc.max():.2f} um "
                      f"({c}); RMS difference {1e6*rms(tc, Dc, tf, Df):.3f} um",
                      ""]
    (HERE/"results.md").write_text("\n".join(lines))
    print("\n".join(lines))
    figure(data, coupled)


def figure(data, coupled):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    (HERE/"plots").mkdir(exist_ok=True)
    style = {"full3D": ("k", "-", "full 3D FSI (15 cm)"),
             "only1D": ("#0072B2", ":", "1D alone (15 cm)")}
    cols = ["#D55E00", "#009E73", "#CC79A7", "#E69F00"]
    for i, c in enumerate(coupled):
        style[c] = (cols[i % len(cols)], "--",
                    f"3D (5 cm) + 1D (10 cm), c0 = {c0_of(c):.3g} m/s"
                    + (", free wall end" if c.endswith("_free") else ""))
    planes = [z for z in PLANES_CM if any(z in d for d in data.values())]
    for col, ylabel, scale in (("P", "mean pressure [Pa]", 1),
                               ("Q", "flow [ml/s]", 1e6)):
        n = len(planes)
        ncol = 4
        nrow = int(math.ceil(n/ncol))
        fig, axs = plt.subplots(nrow, ncol, figsize=(3.4*ncol, 2.9*nrow),
                                sharex=True, sharey=True, squeeze=False)
        for ax, zcm in zip(axs.flat, planes):
            for case, secs in data.items():
                if zcm in secs:
                    t, P, Q = secs[zcm]
                    c, ls, label = style[case]
                    ax.plot(1e3*t, scale*(P if col == "P" else Q), color=c,
                            ls=ls, lw=1.2, label=label)
            where = ("3D" if zcm*1e-2 < L3D else "interface" if zcm*1e-2 == L3D
                     else "1D") + " in the coupled runs"
            ax.set_title(f"z = {zcm:g} cm ({where})", fontsize=9)
            ax.grid(alpha=0.3)
        for ax in axs.flat[len(planes):]:
            ax.axis("off")
        for ax in axs[-1]:
            ax.set_xlabel("t [ms]")
        for ax in axs[:, 0]:
            ax.set_ylabel(ylabel)
        axs.flat[0].legend(frameon=False, fontsize=6)
        fig.tight_layout()
        fig.savefig(HERE/"plots"/f"sections_{col}.png", dpi=150)
        plt.close(fig)
    if "full3D" in data and coupled:
        fig, ax = plt.subplots(figsize=(6, 3.2))
        for case in ["full3D"] + [c for c in coupled if c in data]:
            t, D = displacement(case)
            c, ls, label = style[case]
            ax.plot(1e3*t, 1e6*D, color=c, ls=ls, lw=1.2, label=label)
        ax.set_xlabel("t [ms]")
        ax.set_ylabel("radial displacement at z = 2.5 cm [um]")
        ax.grid(alpha=0.3)
        ax.legend(frameon=False, fontsize=7)
        fig.tight_layout()
        fig.savefig(HERE/"plots"/"displacement.png", dpi=150)
        plt.close(fig)


if __name__ == "__main__":
    main()
