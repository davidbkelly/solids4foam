#!/usr/bin/env python3
"""Single-pulse benchmark of Boileau et al. (2015), section 4.1: analysis and figure.

Used by the plot.py of singlePulseInviscid and singlePulseViscous. For the case:

  * the peak of P(t) at the probes x = 0..9 m (a parabola through the three
    samples around the maximum) and its decay relative to the inlet, and for
    the viscous case the damping of Eq. 15,
        Ppeak,visc(x) = Ppeak,inv exp(-(zeta + 2) pi mu x/(rho c0 A0)),
    taken two ways: against the inlet peak, and against the inviscid run's own
    peak at the same x (if ../singlePulseInviscid has run), which removes the
    numerical dissipation common to both;
  * Fig. 2: P/Ppeak along the vessel at t = 0.1, 0.3, ..., 1.5 s, with the six
    schemes of the paper.

Writes peaks.md and plots/<case>_fig2.png. Table II parameters are used for c0,
A0 and Z0. Needs numpy and matplotlib.
"""
import math
from pathlib import Path

import numpy as np

import boileauMetrics as bm

# Table II
RHO = 1050.0
MU_VISCOUS = 4e-3
ZETA = 9
A0 = math.pi*0.01**2
BETA = 4.0/3.0*math.sqrt(math.pi)*400e3*1.5e-3
C0 = math.sqrt(BETA/(2*RHO*A0))*A0**0.25
Z0 = RHO*C0/A0
QPEAK = 1e-6
PPEAK = Z0*QPEAK                        # [Pa]
GAMMA = (ZETA + 2)*math.pi*MU_VISCOUS/(RHO*C0*A0)   # [1/m]

TIMES = [0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5]
PROBES = list(range(10))


def probe_peak(path):
    """Peak value and time of P at a probe, refined by a parabola."""
    a = np.loadtxt(path, comments="#")
    t, P = a[:, 0], a[:, 1]
    i = int(np.clip(P.argmax(), 1, len(P) - 2))
    y0, y1, y2 = P[i - 1], P[i], P[i + 1]
    den = y0 - 2*y1 + y2
    s = 0.5*(y0 - y2)/den if den != 0 else 0.0
    return y1 - 0.25*(y0 - y2)*s, t[i] + s*(t[1] - t[0])


def peaks(case_dir):
    out = Path(case_dir)/"postProcessing"/"circulation1D"/"0"/"probes"
    return {x: probe_peak(out/f"x{x}.dat") for x in PROBES}


def analyse(case_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    case_dir = Path(case_dir)
    viscous = "Viscous" in case_dir.name
    kind = "viscous" if viscous else "inviscid"

    pk = peaks(case_dir)
    p_in, t_in = pk[0]

    lines = [
        f"# Single pulse, {kind}: peak of P(t) along the vessel\n",
        f"Linear theory: c0 = {C0:.4f} m/s, Z0 Qpeak = {PPEAK:.4f} Pa"
        + (f", Eq. 15 damping {GAMMA:.5f} 1/m" if viscous else "") + "\n",
    ]

    inviscid = None
    inv_dir = case_dir.parent/"singlePulseInviscid"
    if viscous and (inv_dir/"postProcessing").is_dir():
        inviscid = peaks(inv_dir)

    if not viscous:
        lines.append("| x [m] | peak [Pa] | decay from inlet [%] | speed x/(t - t_in) [m/s] |")
        lines.append("|---|---|---|---|")
        for x in PROBES:
            p, t = pk[x]
            speed = f"{x/(t - t_in):.4f}" if x else "-"
            lines.append(f"| {x} | {p:.4f} | {100*(1 - p/p_in):.3f} | {speed} |")
    else:
        head = "| x [m] | peak [Pa] | Eq. 15 x inlet peak [Pa] | error [%] |"
        sep = "|---|---|---|---|"
        if inviscid:
            head += " peak/inviscid peak | Eq. 15 exp(-g x) | error [%] |"
            sep += "---|---|---|"
        lines += [head, sep]
        for x in PROBES:
            p = pk[x][0]
            th = p_in*math.exp(-GAMMA*x)
            row = f"| {x} | {p:.4f} | {th:.4f} | {100*(p/th - 1):+.3f} |"
            if inviscid:
                r = p/inviscid[x][0]
                e = math.exp(-GAMMA*x)
                row += f" {r:.5f} | {e:.5f} | {100*(r/e - 1):+.3f} |"
            lines.append(row)
    text = "\n".join(lines) + "\n"
    (case_dir/"peaks.md").write_text(text)
    print(text)

    # ---- Fig. 2
    prof = case_dir/"postProcessing"/"circulation1D"/"0"/"profiles"
    ref = {s: np.loadtxt(case_dir/"referenceData"/f"{s}.txt", comments="#")
           for s in bm.SCHEMES}
    pk_kpa = PPEAK*1e-3

    fig, (ax, az) = plt.subplots(
        2, 1, figsize=(9, 6.4), gridspec_kw={"height_ratios": [1.3, 1]}
    )
    # solids4foam underneath, so that the thin lines of the schemes stay visible
    for j, t in enumerate(TIMES):
        f = prof/f"{t:g}"/"tube.dat"
        if f.is_file():
            p = np.loadtxt(f, comments="#")
            for a in (ax, az):
                a.plot(p[:, 0], p[:, 3]/PPEAK, color=bm.COLOURS["solids4foam"],
                       linewidth=3.0, alpha=0.85,
                       label="solids4foam" if j == 0 else None)
    for j, t in enumerate(TIMES):
        for s in bm.SCHEMES:
            d = ref[s]
            for a in (ax, az):
                a.plot(d[:, 0], d[:, j + 1]/pk_kpa, color=bm.COLOURS[s],
                       linewidth=0.9, label=s if j == 0 else None)
    if viscous:
        x = np.linspace(0, 10, 200)
        for a in (ax, az):
            a.plot(x, np.exp(-GAMMA*x), color=bm.COLOURS["theory"], linewidth=1.4,
                   linestyle=(0, (4, 2)), label="Eq. 15")
    else:
        for a in (ax, az):
            a.axhline(1.0, color=bm.INK_2, linewidth=0.8, linestyle=(0, (4, 2)),
                      label="linear theory")
    ax.set_xlim(0, 10)
    xc = C0*(1.5 - 0.05)
    az.set_xlim(xc - 0.25, xc + 0.25)
    lo = math.exp(-GAMMA*xc) if viscous else 1.0
    az.set_ylim(lo - 0.05, lo + 0.02)
    ax.set_ylabel(r"$P/P_{peak}$")
    az.set_ylabel(r"$P/P_{peak}$")
    az.set_xlabel("x [m]")
    az.set_title("t = 1.5 s, at the peak", fontsize=9, color=bm.INK_2)
    for a in (ax, az):
        bm.style(a)
    handles, names = ax.get_legend_handles_labels()
    fig.suptitle(f"Single pulse, {kind}: t = 0.1, 0.3, ..., 1.5 s "
                 "(Boileau et al. 2015, Fig. 2)", color=bm.INK, fontsize=11, y=0.99)
    fig.legend(handles, names, loc="upper center", bbox_to_anchor=(0.5, 0.955),
               ncol=8, frameon=False, fontsize=8, labelcolor=bm.INK)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    (case_dir/"plots").mkdir(exist_ok=True)
    fig.savefig(case_dir/"plots"/f"{case_dir.name}_fig2.png", dpi=150)
    plt.close(fig)
