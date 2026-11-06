#!/usr/bin/env python3
"""Hard and soft coupling of a rigid 3D pipe to a 1D vessel (Dobroserdova et
al. 2016, sections 4.2-4.3) against the exact linear solution.

    python3 analyse.py

The vessel is a lossless transmission line of impedance Zc = rho c0/S and
transit time td = L1D/c0 (14 ms at c0 = 3.5 m/s), much shorter than the
pulses, so waves reflected at the interface return from the flow-source inlet
(which reflects them whole) many times per period: the vessel is a resonator
between the inlet and the rigid pipe. In the waves f (forward) and g
(backward) at the interface, P = f + g, Q = (f - g)/Zc, the inlet condition
gives the delay relation

    f(t) = Zc Qin(t - td) + g(t - 2 td),

and the rigid pipe with its open end is an inductor L = rho l/S to ground,
so the hard coupling obeys L dQ/dt = P, and the soft coupling feeds the
source through R1 - C - R2 - L - ground (R1 + R2 = Zc, C = l0 S/(rho c0^2)).
These delay equations are integrated here exactly (to 1e-9). The pipe
reflects everything, Gamma = (i w L - Zc)/(i w L + Zc) of modulus 1, with a
phase lag; their reference (a 10 cm vessel with an absorbing end) has no
reflection at all, so the "spurious" reflection of the hard coupling, and the
resonance it drives, are the exact solution of the problem posed; the soft
coupling's R = Zc damps that resonance.

Compared, at the 1D-3D interface: the 1D runs with the lumped pipe (exact
up to the 1D discretisation), the 3D hard-coupling runs, their reference,
and the delay-equation solutions. Writes results.md and plots/*.png.
"""
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from makeCases import A, C0S, INFLOWS, L0, L1D, LIND, RHO  # noqa: E402

RUNS = HERE/"runs"


def qin(t, inflow):
    """The inflow; zero before t = 0 (the runs start from rest)"""
    t = np.asarray(t, dtype=float)
    if inflow == "sine":
        return np.where(t >= 0, 1e-7*np.sin(2*np.pi*t/0.3), 0.0)
    return np.where((t >= 0) & (t <= 0.1), 1e-7, 0.0)


def exact(inflow, c0, config, t):
    """Q and P at the interface from the delay equations, by Heun's method on
    a grid that divides the round-trip time 2 td exactly"""
    Zc = RHO*c0/A
    Cb = L0*A/(RHO*c0**2)
    R1 = R2 = 0.5*Zc
    td = L1D/c0
    n2 = 4000                              # steps per round trip
    dt = 2*td/n2
    N = int(math.ceil(t[-1]/dt)) + 2
    tt = np.arange(N + 1)*dt
    g = np.zeros(N + 1)                    # backward wave at the interface
    Q = np.zeros(N + 1)
    Pc = np.zeros(N + 1)
    Q2 = np.zeros(N + 1)

    def f_at(k):                           # forward wave at grid point k
        gd = g[k - n2] if k >= n2 else 0.0
        return Zc*float(qin(np.array([tt[k] - td]), inflow)[0]) + gd

    def rates(fk, Qk, Pck, Q2k):
        if config == "hard":
            return (2*fk - Zc*Qk)/LIND, 0.0, 0.0
        Q1 = (2*fk - Pck)/(Zc + R1)
        return 0.0, (Q1 - Q2k)/Cb, (Pck - R2*Q2k)/LIND

    for k in range(N):
        fk, fk1 = f_at(k), f_at(k + 1)
        if config == "hard":
            a1 = rates(fk, Q[k], 0, 0)[0]
            qp = Q[k] + dt*a1
            a2 = rates(fk1, qp, 0, 0)[0]
            Q[k + 1] = Q[k] + 0.5*dt*(a1 + a2)
        else:
            _, b1, c1 = rates(fk, 0, Pc[k], Q2[k])
            pp, q2p = Pc[k] + dt*b1, Q2[k] + dt*c1
            _, b2, c2 = rates(fk1, 0, pp, q2p)
            Pc[k + 1] = Pc[k] + 0.5*dt*(b1 + b2)
            Q2[k + 1] = Q2[k] + 0.5*dt*(c1 + c2)
            Q[k + 1] = (2*fk1 - Pc[k + 1])/(Zc + R1)
        g[k + 1] = fk1 - Zc*Q[k + 1]
    P = np.array([f_at(k) for k in range(N + 1)]) + g
    return np.interp(t, tt, Q), np.interp(t, tt, P)


def load_1d(name):
    f = RUNS/name/"postProcessing"/"circulation1D"/"0"/"probes"/"interface.dat"
    if not f.is_file():
        return None
    a = np.loadtxt(f, comments="#")
    return a[:, 0], a[:, 2], a[:, 1]                 # t, Q, P


def load_3d(name):
    f = RUNS/name/"postProcessing"/"circulation1DCoupling"/"0"/"coupling.dat"
    log = RUNS/name/"log.pimpleFoam"
    if not f.is_file() or not log.is_file() or \
            not log.read_text().rstrip().endswith("End"):
        return None
    a = np.loadtxt(f, comments="#", ndmin=2)
    # t updates to3D.Q to3D.P interface.P interface.Q to3D.Qsampled
    return a[:, 0], a[:, 6], a[:, 3]


def main():
    lines = ["# Dobroserdova et al. (2016), sections 4.2-4.3: hard and soft "
             "coupling against the exact linear solution", ""]
    figs = {}
    for inflow in INFLOWS:
        lines += [f"## {inflow} inflow", "",
                  "| c0 [m/s] | config | peak Q / peak Qin | peak Q, exact "
                  "| Q error, RMS / max | P error, RMS / max |",
                  "|---|---|---|---|---|---|"]
        for c0name, c0 in C0S.items():
            Zc = RHO*c0/A
            ref = load_1d(f"{inflow}_{c0name}_reference")
            rows = []
            if ref:
                t, Q, P = ref
                rows.append(("reference (10 cm, absorbing)", Q.max()/1e-7,
                             float("nan"), "", "", t, Q, P, None, None))
            def err(x, xe, scale):
                e = (x - xe)/scale
                return f"{math.sqrt(np.mean(e**2)):.1e} / {np.abs(e).max():.1e}"

            for config in ("hard", "soft"):
                d = load_1d(f"{inflow}_{c0name}_{config}")
                if d is None:
                    continue
                t, Q, P = d
                Qe, Pe = exact(inflow, c0, config, t)
                rows.append((config + " (1D + lumped pipe)", Q.max()/1e-7,
                             Qe.max()/1e-7, err(Q, Qe, 1e-7),
                             err(P, Pe, Zc*1e-7), t, Q, P, Qe, Pe))
            d = load_3d(f"{inflow}_{c0name}_hard3D")
            if d is not None:
                t, Q, P = d
                Qe, Pe = exact(inflow, c0, "hard", t)
                rows.append(("hard (1D + pimpleFoam pipe)", Q.max()/1e-7,
                             Qe.max()/1e-7, err(Q, Qe, 1e-7),
                             err(P, Pe, Zc*1e-7), t, Q, P, Qe, Pe))
            for r in rows:
                lines.append(f"| {c0} | {r[0]} | {r[1]:.3f} | " +
                             (f"{r[2]:.3f} | {r[3]} | {r[4]} |"
                              if not (isinstance(r[2], float) and math.isnan(r[2]))
                              else "- | - | - |"))
            figs[(inflow, c0name)] = rows
        w = 2*np.pi/0.3
        lines += ["", "Hard coupling, sine, periodic state: |Gamma| = 1 at "
                  "every c0; with the inlet's own reflection the interface "
                  "flow is Qin/|cos(w td) + i (w L/Zc) sin(w td)| in amplitude: " +
                  ", ".join(f"{1/abs(math.cos(w*L1D/c0) + 1j*(w*LIND/(RHO*c0/A))*math.sin(w*L1D/c0)):.3f} at {c0} m/s"
                            for c0 in C0S.values()) +
                  "; the runs start from rest and the undamped resonance of "
                  "the vessel (about c0/(4 L1D) Hz) persists on top of it.", ""]
    (HERE/"results.md").write_text("\n".join(lines))
    print("\n".join(lines))
    figure(figs)


def figure(figs):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    (HERE/"plots").mkdir(exist_ok=True)
    for inflow in INFLOWS:
        fig, axs = plt.subplots(1, 3, figsize=(13, 3.6), sharey=True)
        for ax, (c0name, c0) in zip(axs, C0S.items()):
            rows = figs.get((inflow, c0name), [])
            for r in rows:
                label, t, Q, Qe = r[0], r[5], r[6], r[8]
                style = {"reference": ("k", "-"), "hard (1D": ("#D55E00", "-"),
                         "soft": ("#0072B2", "-"),
                         "hard (1D + pimpleFoam": ("#009E73", "--")}
                col, ls = next(v for k, v in style.items() if label.startswith(k)) \
                    if any(label.startswith(k) for k in style) else ("gray", "-")
                if label.startswith("hard (1D + pimpleFoam"):
                    col, ls = "#009E73", "--"
                ax.plot(t, Q*1e6, color=col, ls=ls, lw=1.3, label=label)
                if Qe is not None and not label.startswith("hard (1D + pimple"):
                    ax.plot(t, Qe*1e6, ":", color="k", lw=0.9,
                            label=None if "soft" in label else "exact (delay equation)")
            ax.set_title(f"{inflow}, c0 = {c0} m/s")
            ax.set_xlabel("t [s]")
            if inflow == "step":
                ax.set_xlim(0, 0.3)
        axs[0].set_ylabel("Q at the interface [cm3/s]")
        axs[-1].legend(frameon=False, fontsize=7)
        fig.tight_layout()
        fig.savefig(HERE/"plots"/f"dobroserdova_{inflow}.png", dpi=150)
        plt.close(fig)


if __name__ == "__main__":
    main()
