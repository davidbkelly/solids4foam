#!/usr/bin/env python3
"""Compare the coupled runs of the Businaro (2026) pressure Gaussian pulse
with the full FSI tube, with the 1D vessel alone and with the interface
mismatch reported by Businaro.

    python3 analyse.py [--suffix _smoke]

At each section z (every centimetre) the mean pressure and the flow come
from the faceZone function objects of the 3D runs (quarter tube: the flow
times 4, the kinematic pressure times rho; the flow is phi + meshPhi, the
flux through the moving faces, since the solver leaves phi relative to the
mesh and the sections move axially with the wall; the pressure is pSection,
the mean of the two cells of every face, processor faces included, where p
itself holds the neighbour cell's value; at the outlet, whose points are
fixed, phi and p themselves), from the coupler's log beyond
the interface (the vessel probes) and from the probes of the 1D vessel
alone. Every run is compared with full3D (45 cm, nothing reflected from its
far end reaches z <= 30 cm by t = 0.1 s): RMS and maximum pressure
difference over the window relative to the 2 kPa pulse, RMS flow difference
relative to the reference's peak flow at that section. The interface
mismatch is Businaro's Eq. 7.12, e_P = P_3D - P_1D and e_Q = Q_3D - Q_1D at
the interface: the 3D side from the outlet patch at the end of each step,
the 1D side from the coupler (the pressure the vessel returned and the flow
it was forced with). Writes results.md and plots/*.png.
"""
import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from makeCases import (END_TIME, L3D, LFULL, LTUBE, P_MAX, RHO,  # noqa: E402
                       c0_1d, planes, zone_name)

RUNS = HERE/"runs"
FLOW_SCALE = 4                        # quarter tube
Z_TUBE = list(range(1, int(round(LTUBE*100)) + 1))
Z_PORT = int(round(L3D*100))
SNAPSHOTS = (0.06, 0.08, 0.10)        # Businaro's Figure 7.12

# Businaro (2026), Section 7.4: the largest interface mismatch of his
# explicit coupling and what it is relative to
BUSINARO = {"eP": 4.452e-1, "eP_rel": 2.226e-4, "eQ": 5.757e-6, "eQ_rel": 5.573e-2}
BUSINARO["Qmax"] = BUSINARO["eQ"]/BUSINARO["eQ_rel"]


def finished(case):
    for log in ("log.solids4Foam", "log.circulation1DFoam"):
        f = RUNS/case/log
        if f.is_file():
            return "\nEnd\n" in f.read_text()[-2000:]
    return False


def concat(files, load):
    """The rows of the logs of successive restarts, later ones winning"""
    t, v = np.empty(0), None
    for f in sorted(files, key=lambda f: float(f.parent.name)):
        a = load(f)
        keep = t < a[0, 0] - 1e-12
        t = np.concatenate([t[keep], a[:, 0]])
        v = a[:, 1:] if v is None else np.vstack([v[keep], a[:, 1:]])
    return t, v


def read_fo(case, name):
    files = list((RUNS/case/"postProcessing"/"fluid"/name).glob("*/surfaceFieldValue.dat"))
    t, v = concat(files, lambda f: np.loadtxt(f, comments="#", ndmin=2))
    return t, v[:, 0]


def read_coupling(case):
    """{column: values} of the coupler's coupling.dat"""
    files = list((RUNS/case/"postProcessing"/"circulation1DCoupling").glob("*/coupling.dat"))
    with open(files[0]) as fh:
        names = next(l for l in fh if l.startswith("# time"))[1:].split()
    t, v = concat(files, lambda f: np.loadtxt(f, comments="#", ndmin=2))
    out = {n: v[:, i - 1] for i, n in enumerate(names) if i > 0}
    out["time"] = t
    return out


def sections_3d(case, zs):
    """{z_cm: (t, P [Pa], Q [m3/s])} from the faceZones and the outlet"""
    out = {}
    for z in zs:
        n = zone_name(z)
        try:
            t, p = read_fo(case, n + "_P")
            _, q = read_fo(case, n + "_Q")
        except (OSError, IndexError):
            continue
        out[z] = (t, RHO*p, FLOW_SCALE*q)
    return out


def outlet_3d(case):
    t, p = read_fo(case, "outlet_P")
    _, q = read_fo(case, "outlet_Q")
    return t, RHO*p, FLOW_SCALE*q


def sections_coupled(case):
    """{z_cm: (t, P, Q)}: the 3D part, its outlet at the interface, and the
    vessel probes beyond it"""
    out = sections_3d(case, planes(L3D))
    out[Z_PORT] = outlet_3d(case)
    c = read_coupling(case)
    for z in Z_TUBE:
        n = zone_name(z)
        if n + ".P" in c:
            out[z] = (c["time"], c[n + ".P"], c[n + ".Q"])
    return out


def sections_1d(case):
    out = {}
    d = RUNS/case/"postProcessing"/"circulation1D"/"0"/"probes"
    for z in Z_TUBE:
        f = d/(zone_name(z) + ".dat")
        if f.is_file():
            a = np.loadtxt(f, comments="#")
            out[z] = (a[:, 0], a[:, 1], a[:, 2])
    return out


def interface_mismatch(case):
    """Businaro's Eq. 7.12 at the interface: (t, e_P, e_Q, Q_3D)"""
    t3, P3, Q3 = outlet_3d(case)
    c = read_coupling(case)
    P1 = np.interp(t3, c["time"], c["from3D.P"])
    Q1 = np.interp(t3, c["time"], c["from3D.Q"])
    return t3, P3 - P1, Q3 - Q1, Q3


def on(t, tref, x):
    return np.interp(tref, t, x)


def compare(secs, ref):
    """{z: (RMS dP/P_MAX, max |dP|/P_MAX, RMS dQ/max|Qref|, peak P ratio)}"""
    out = {}
    for z in Z_TUBE:
        if z not in secs or z not in ref:
            continue
        tr, Pr, Qr = ref[z]
        m = tr <= END_TIME + 1e-9
        tr, Pr, Qr = tr[m], Pr[m], Qr[m]
        t, P, Q = secs[z]
        dP = on(t, tr, P) - Pr
        dQ = on(t, tr, Q) - Qr
        out[z] = (np.sqrt(np.mean(dP**2))/P_MAX, np.abs(dP).max()/P_MAX,
                  np.sqrt(np.mean(dQ**2))/np.abs(Qr).max(), on(t, tr, P).max()/Pr.max())
    return out


def arrival(t, P, level=0.5*P_MAX):
    """first time P reaches the level (linear interpolation)"""
    i = np.argmax(P >= level*P.max()/P_MAX) if P.max() > 0 else 0
    if i == 0:
        return np.nan
    return t[i - 1] + (level*P.max()/P_MAX - P[i - 1])*(t[i] - t[i - 1])/(P[i] - P[i - 1])


def wave_speed(secs, zs):
    """least-squares speed of the half-peak arrival over the sections zs"""
    pts = [(z*1e-2, arrival(*secs[z][:2], level=0.5*secs[z][1].max()))
           for z in zs if z in secs]
    pts = [(z, ta) for z, ta in pts if np.isfinite(ta)]
    if len(pts) < 3:
        return np.nan
    z, ta = np.array(pts).T
    return 1/np.polyfit(z, ta, 1)[0]


def snapshot(secs, t0):
    zs = sorted(z for z in secs if z <= LTUBE*100)
    return (np.array(zs)*1e-2, np.array([on(secs[z][0], t0, secs[z][1]) for z in zs]),
            np.array([on(secs[z][0], t0, secs[z][2]) for z in zs]))


LABELS = {"full3D": "3D FSI, 45 cm (reference)",
          "full3D_Rt": "3D FSI, 30 cm, resistance Rt at the outlet",
          "coupled": "3D 10 cm + 1D 20 cm, implicit order 2",
          "coupled_explicit": "3D 10 cm + 1D 20 cm, explicit order 1",
          "only1D": "1D alone, 30 cm"}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--suffix", default="")
    args = ap.parse_args()
    s = args.suffix
    name = {k: k + s for k in LABELS}
    if not finished(name["full3D"]):
        sys.exit(f"{name['full3D']} has not finished")
    ref = sections_3d(name["full3D"], planes(LFULL))
    data = {}
    for k in ("full3D_Rt", "coupled", "coupled_explicit", "only1D"):
        n = name[k]
        if not finished(n):
            continue
        if k == "only1D":
            data[k] = sections_1d(n)
        elif k == "full3D_Rt":
            data[k] = sections_3d(n, planes(LTUBE))
            data[k][int(round(LTUBE*100))] = outlet_3d(n)
        else:
            data[k] = sections_coupled(n)

    zs_show = [5, 9, Z_PORT, 11, 15, 20, 25, 30]
    lines = ["# Businaro (2026) pressure Gaussian pulse: 3D-1D coupling against the full FSI tube",
             "",
             f"Inlet pressure {P_MAX:g} Pa exp(-1e4 (t - 0.05)^2), R = 1 cm, wall 1 mm, "
             "E = 700 kPa, Saint Venant-Kirchhoff, rho = 1060, mu = 4 mPa s; BDF2, "
             f"window 0-{END_TIME:g} s. Reference: full3D, 45 cm. RMS and max pressure "
             "difference over the window / 2 kPa; RMS flow difference / the reference's "
             "peak flow at that section; peak pressure ratio run / reference. z <= "
             f"{Z_PORT} cm is 3D in the coupled runs, beyond it the 1D vessel.", ""]
    lines += ["| run | " + " | ".join(f"z = {z}" for z in zs_show) + " |",
              "|---|" + "---|"*len(zs_show)]
    cmp = {k: compare(v, ref) for k, v in data.items()}
    for k, e in cmp.items():
        lines.append(f"| {LABELS[k]}: RMS dP | "
                     + " | ".join(f"{e[z][0]:.2e}" if z in e else "-" for z in zs_show) + " |")
        lines.append("| max dP | "
                     + " | ".join(f"{e[z][1]:.2e}" if z in e else "-" for z in zs_show) + " |")
        lines.append("| RMS dQ | "
                     + " | ".join(f"{e[z][2]:.2e}" if z in e else "-" for z in zs_show) + " |")
        lines.append("| peak P ratio | "
                     + " | ".join(f"{e[z][3]:.3f}" if z in e else "-" for z in zs_show) + " |")

    lines += ["", "Wave speed from the half-peak arrivals (least squares):", ""]
    lines.append(f"- full3D, z = 2-28 cm: {wave_speed(ref, range(2, 29)):.3f} m/s "
                 f"(the 1D vessel's c0 = {c0_1d():.3f} m/s)")
    for k in data:
        zs = range(2, Z_PORT) if k == "full3D_Rt" else range(Z_PORT + 2, 29)
        if k == "only1D":
            zs = range(2, 29)
        lines.append(f"- {k}, z = {zs.start}-{zs.stop - 1} cm: {wave_speed(data[k], zs):.3f} m/s")

    lines += ["", "Interface mismatch (Businaro Eq. 7.12): largest |e_P| and |e_Q| over the "
              "window, and relative to the pulse and to the largest interface flow:", "",
              "| run | max abs e_P [Pa] | / 2 kPa | max abs e_Q [m3/s] | / max Q | max Q at the interface [m3/s] |",
              "|---|---|---|---|---|---|"]
    mism = {}
    for k in ("coupled", "coupled_explicit"):
        if k not in data:
            continue
        t, eP, eQ, Q3 = interface_mismatch(name[k])
        mism[k] = (t, eP, eQ)
        qm = np.abs(Q3).max()
        lines.append(f"| {LABELS[k]} | {np.abs(eP).max():.3e} | {np.abs(eP).max()/P_MAX:.2e} | "
                     f"{np.abs(eQ).max():.3e} | {np.abs(eQ).max()/qm:.2e} | {qm:.4e} |")
    lines.append(f"| Businaro (2026), explicit (lifex) | {BUSINARO['eP']:.3e} | "
                 f"{BUSINARO['eP_rel']:.2e} | {BUSINARO['eQ']:.3e} | {BUSINARO['eQ_rel']:.2e} | "
                 f"{BUSINARO['Qmax']:.4e} |")
    tr, Pr, Qr = ref[Z_PORT]
    lines.append(f"\nThe full tube's largest flow at z = {Z_PORT} cm: {np.abs(Qr).max():.4e} m3/s, "
                 f"peak pressure {Pr.max():.1f} Pa.")

    text = "\n".join(lines) + "\n"
    (HERE/("results" + s + ".md")).write_text(text)
    print(text)
    figures(ref, data, mism, s)


def figures(ref, data, mism, s):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    (HERE/"plots").mkdir(exist_ok=True)
    style = {"full3D": ("k", "-"), "full3D_Rt": ("tab:gray", "--"),
             "coupled": ("tab:blue", "-"), "coupled_explicit": ("tab:orange", ":"),
             "only1D": ("tab:green", "-.")}
    runs = {"full3D": ref, **data}

    fig, ax = plt.subplots(2, len(SNAPSHOTS), figsize=(13, 6), sharex=True)
    for j, t0 in enumerate(SNAPSHOTS):
        for k, secs in runs.items():
            z, P, Q = snapshot(secs, t0)
            c, ls = style[k]
            ax[0, j].plot(100*z, P, color=c, ls=ls, label=LABELS[k])
            ax[1, j].plot(100*z, 1e6*Q, color=c, ls=ls)
        ax[0, j].set_title(f"t = {t0:g} s")
        ax[1, j].set_xlabel("z [cm]")
        for a in ax[:, j]:
            a.axvline(Z_PORT, color="0.8", lw=0.8)
            a.grid(alpha=0.3)
    ax[0, 0].set_ylabel("P [Pa]")
    ax[1, 0].set_ylabel("Q [ml/s]")
    ax[0, 0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(HERE/"plots"/("snapshots" + s + ".png"), dpi=150)
    plt.close(fig)

    zs = [5, Z_PORT, 20, 30]
    fig, ax = plt.subplots(2, len(zs), figsize=(14, 6), sharex=True)
    for j, z in enumerate(zs):
        for k, secs in runs.items():
            if z not in secs:
                continue
            t, P, Q = secs[z]
            c, ls = style[k]
            ax[0, j].plot(t, P, color=c, ls=ls, label=LABELS[k])
            ax[1, j].plot(t, 1e6*Q, color=c, ls=ls)
        ax[0, j].set_title(f"z = {z} cm")
        ax[1, j].set_xlabel("t [s]")
        for a in ax[:, j]:
            a.grid(alpha=0.3)
    ax[0, 0].set_ylabel("P [Pa]")
    ax[1, 0].set_ylabel("Q [ml/s]")
    ax[0, 0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(HERE/"plots"/("sections" + s + ".png"), dpi=150)
    plt.close(fig)

    if mism:
        fig, ax = plt.subplots(1, 2, figsize=(11, 3.5))
        for k, (t, eP, eQ) in mism.items():
            c, ls = style[k]
            ax[0].plot(t, eP, color=c, ls=ls, label=LABELS[k])
            ax[1].plot(t, eQ, color=c, ls=ls)
        ax[0].set_ylabel("e_P = P_3D - P_1D [Pa]")
        ax[1].set_ylabel("e_Q = Q_3D - Q_1D [m3/s]")
        for a in ax:
            a.set_xlabel("t [s]")
            a.grid(alpha=0.3)
        ax[0].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(HERE/"plots"/("interface" + s + ".png"), dpi=150)
        plt.close(fig)


if __name__ == "__main__":
    main()
