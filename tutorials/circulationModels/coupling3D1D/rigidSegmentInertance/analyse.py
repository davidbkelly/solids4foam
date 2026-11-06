#!/usr/bin/env python3
"""Rigid segment as an inertance: the 1D runs against the linear closed form.

    python3 analyse.py

A rigid pipe of length l between two equal vessels is the inertance
L = rho l/A in series (src/circulationModels/docs/BENCHMARKS3D1D.md). For an
incident pressure pulse Pinc(T) = P exp(-T^2/(2 sigma^2)) at the pipe, T from
its peak there, the transmitted and reflected pulses are

    tau dPt/dT + Pt = Pinc,     Pr = Pinc - Pt,     tau = l/(2 c0)

with Pt in closed form (erfc). This compares

  * the base case (dx = 1.25 mm): the probes at the two midpoints and the
    pipe's own output (Q, P1, P2) against the closed form; figure
    plots/rigidSegment.png;
  * the order/n<N> runs: the error of P along both vessels at t = 0.6 s,
    relative to the peak, and its order.

Both use the inflow of Boileau Eq. 14 x 1e-3, linear to about 1e-6.

Writes results.md. Needs numpy and scipy (matplotlib for the figure).
"""
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[3] / "src" / "circulationModels" / "docs"))
sys.path.insert(0, str(HERE.parents[1] / "circulation1D" / "common"))

from benchmarks3D1D import transmitted  # noqa: E402
from singlePulse import C0, Z0  # noqa: E402

L_VESSEL = 2.0
L_PIPE = 0.08
TAU = L_PIPE/(2*C0)
SIGMA = 1/math.sqrt(2e4)            # exp(-1e4 (t - 0.05)^2)
T_PEAK = 0.05
T_ARRIVE = T_PEAK + L_VESSEL/C0     # incident peak at the pipe
T_END = 0.6
QPEAK = 1e-6
SCALE = 1e-3                        # amplitude factor of the inflow


def pinc(T, P):
    return P*np.exp(-T**2/(2*SIGMA**2))


def pt(T, P):
    return transmitted(T, P, SIGMA, TAU)


def exact_profile(vessel, x, t, P):
    """P along a vessel at time t, x from its in end"""
    if vessel == "proximal":
        d = L_VESSEL - x                          # distance to the pipe
        T = t - T_ARRIVE - d/C0                   # reflected
        return pinc(t - x/C0 - T_PEAK, P) + pinc(T, P) - pt(T, P)
    return pt(t - T_ARRIVE - x/C0, P)


def out_dir(case):
    return case / "postProcessing" / "circulation1D" / "0"


def profile_errors(case, scale):
    P = Z0*scale*QPEAK
    e = []
    for v in ("proximal", "distal"):
        a = np.loadtxt(out_dir(case)/"profiles"/f"{T_END:g}"/f"{v}.dat",
                       comments="#")[1:-1]
        e.append(np.abs(a[:, 3] - exact_profile(v, a[:, 0], T_END, P))/P)
    e = np.concatenate(e)
    return np.mean(e), math.sqrt(np.mean(e**2)), e.max()


def base():
    case = HERE
    P = Z0*SCALE*QPEAK
    out = out_dir(case)
    lines = []
    if not (out/"nodes"/"rigid.dat").is_file():
        return lines, None

    rig = np.loadtxt(out/"nodes"/"rigid.dat", comments="#")
    t, Q, P1, P2 = rig.T
    T = t - T_ARRIVE
    exQ = pt(T, P)/Z0
    exP1 = 2*pinc(T, P) - pt(T, P)
    exP2 = pt(T, P)

    prox = np.loadtxt(out/"probes"/"proximalMid.dat", comments="#")
    dist = np.loadtxt(out/"probes"/"distalMid.dat", comments="#")
    exProx = exact_profile("proximal", 1.0, prox[:, 0], P)
    exDist = exact_profile("distal", 1.0, dist[:, 0], P)

    def rel(a, b, s):
        return np.abs(a - b).max()/s

    lines += [
        "## Base case (dx = 1.25 mm)",
        "",
        f"c0 = {C0:.4f} m/s, Z0 = {Z0:.4e} Pa s/m3, pipe 8 cm: "
        f"tau = {TAU*1e3:.3f} ms, sigma = {SIGMA*1e3:.3f} ms, "
        f"incident peak {P*1e3:.3f} mPa (inflow x {SCALE:g}).",
        "",
        "| quantity | closed form peak | 1D peak | max error / incident peak |",
        "|---|---|---|---|",
        f"| pipe flow Q [ul/s] | {exQ.max()*1e9:.4f} | {Q.max()*1e9:.4f} | "
        f"{rel(Q*Z0, exQ*Z0, P):.2e} |",
        f"| P1, pipe inlet [mPa] | {exP1.max()*1e3:.4f} | {P1.max()*1e3:.4f} | "
        f"{rel(P1, exP1, P):.2e} |",
        f"| P2 = Pt, pipe outlet [mPa] | {exP2.max()*1e3:.4f} | "
        f"{P2.max()*1e3:.4f} | {rel(P2, exP2, P):.2e} |",
        f"| proximal midpoint [mPa] | {exProx.max()*1e3:.4f} | "
        f"{prox[:, 1].max()*1e3:.4f} | {rel(prox[:, 1], exProx, P):.2e} |",
        f"| distal midpoint [mPa] | {exDist.max()*1e3:.4f} | "
        f"{dist[:, 1].max()*1e3:.4f} | {rel(dist[:, 1], exDist, P):.2e} |",
        "",
    ]

    # Reflected peak at the proximal midpoint (after the incident has passed)
    late = prox[:, 0] > T_ARRIVE
    lines.append(
        f"Reflected peak {prox[late, 1].max()*1e3:.4f} mPa against "
        f"{exProx[late].max()*1e3:.4f} mPa; energy (reflected + transmitted)/"
        "incident from the closed form is 1 to round-off "
        "(benchmarks3D1D.py rigid)."
    )
    lines.append("")
    return lines, (t, Q, P1, P2, exQ, exP1, exP2, prox, exProx, dist, exDist, P)


def orders():
    rows = []
    for d in sorted((HERE/"order").glob("n*"), key=lambda p: int(p.name[1:])):
        f = out_dir(d)/"profiles"/f"{T_END:g}"/"distal.dat"
        if f.is_file():
            n = int(d.name[1:])
            rows.append((n, L_VESSEL/n) + profile_errors(d, SCALE))
    lines = []
    if not rows:
        return lines, rows

    lines += [
        "## Order of accuracy (both vessels at t = 0.6 s)",
        "",
        "| cells per vessel | dx [mm] | L1 | order | L2 | order | max | order |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for i, r in enumerate(rows):
        cells = [f"{r[0]}", f"{r[1]*1e3:.4g}"]
        for k in (2, 3, 4):
            cells.append(f"{r[k]:.3e}")
            if i == 0:
                cells.append("")
            else:
                cells.append(
                    f"{math.log(rows[i-1][k]/r[k])/math.log(rows[i-1][1]/r[1]):.2f}"
                )
        lines.append("| " + " | ".join(cells) + " |")
    lines.append("")

    fine = rows[-4:]
    if len(fine) >= 2:
        h = np.log([r[1] for r in fine])
        fits = []
        for k in (2, 3, 4):
            fits.append(np.polyfit(h, np.log([r[k] for r in fine]), 1)[0])
        lines.append(
            "Least-squares order over the %d finest meshes: L1 %.2f, L2 %.2f, "
            "max %.2f." % (len(fine), *fits)
        )
        lines.append("")
    return lines, rows


def figure(data, rows):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    (HERE/"plots").mkdir(exist_ok=True)

    fig, axs = plt.subplots(1, 3 if rows else 2, figsize=(13, 3.8))
    if data:
        t, Q, P1, P2, exQ, exP1, exP2, prox, exProx, dist, exDist, P = data
        ax = axs[0]
        ms = 1e3
        ax.plot((t - T_ARRIVE)*ms, P1/P, color="#2a6fb0", label="P1 (pipe in)")
        ax.plot((t - T_ARRIVE)*ms, P2/P, color="#c0392b", label="P2 (pipe out)")
        ax.plot((t - T_ARRIVE)*ms, exP1/P, "k--", lw=1, label="closed form")
        ax.plot((t - T_ARRIVE)*ms, exP2/P, "k--", lw=1)
        ax.set_xlim(-40, 60)
        ax.set_xlabel("time from incident peak at the pipe [ms]")
        ax.set_ylabel("P / incident peak")
        ax.set_title("At the rigid pipe")
        ax.legend(frameon=False, fontsize=8)

        ax = axs[1]
        ax.plot(prox[:, 0], prox[:, 1]/P, color="#2a6fb0",
                label="proximal midpoint")
        ax.plot(dist[:, 0], dist[:, 1]/P, color="#c0392b",
                label="distal midpoint")
        ax.plot(prox[:, 0], exProx/P, "k--", lw=1, label="closed form")
        ax.plot(dist[:, 0], exDist/P, "k--", lw=1)
        ax.set_xlim(0.15, 0.6)
        ax.set_xlabel("t [s]")
        ax.set_title("Incident, reflected and transmitted")
        ax.legend(frameon=False, fontsize=8)

    if rows:
        ax = axs[-1]
        dx = np.array([r[1] for r in rows])*1e3
        for k, lab in ((2, "L1"), (3, "L2"), (4, "max")):
            ax.loglog(dx, [r[k] for r in rows], "o-", ms=4, label=lab)
        ax.loglog(dx, rows[-1][3]*(dx/dx[-1])**2, "k:", lw=1, label="slope 2")
        ax.set_xlabel("dx [mm]")
        ax.set_ylabel("error / peak")
        ax.set_title("Convergence (linear amplitude)")
        ax.legend(frameon=False, fontsize=8)

    fig.tight_layout()
    fig.savefig(HERE/"plots"/"rigidSegment.png", dpi=150)


def main():
    lines = ["# Rigid segment (inertance) against the closed form", ""]
    b, data = base()
    o, rows = orders()
    if not b and not o:
        sys.exit("No finished runs: run ./Allrun first")
    lines += b + o
    (HERE/"results.md").write_text("\n".join(lines))
    print("\n".join(lines))
    figure(data, rows)


if __name__ == "__main__":
    main()
