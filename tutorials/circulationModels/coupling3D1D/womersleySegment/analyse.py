#!/usr/bin/env python3
"""The Womersley segment against the exact periodic solution.

    python3 analyse.py

For each run the periodic state of the ports (pressure and flow at the
pipe's inlet and outlet, from the coupler's coupling.dat) is fitted with a
sinusoid at the forcing frequency over the last two periods and compared,
in amplitude and phase, with the transmission-line solution in which the
pipe is the series impedance of Womersley (1955):

    Z = i w rho l / (A (1 - F10)),   F10 = 2 J1(L)/(L J0(L)),   L = i^(3/2) alpha
    Gamma = Z/(Z + 2 Zc),   T = 2 Zc/(Z + 2 Zc),   Zc = rho c0/A

so that, for an incident wave of complex amplitude P at the inlet port,
P_in = (1 + Gamma) P, Q = T P/Zc and P_out = T P. The incident wave is the
inlet's p(t) delayed by the proximal vessel, L1/c0. The inviscid inertance
(Z = i w rho l/A) is tabulated as well, to show what viscosity changes. The
axial velocity profile at the pipe's mid-length over the last period is
compared with Womersley's profile for the fitted flow. Writes results.md
and plots/*.png.
"""
import cmath
import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import jv

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from makeCases import A, C0, F, L1D, L3D, OMEGA, P_INC, R, RHO, nu  # noqa: E402

RUNS = HERE/"runs"
ZC = RHO*C0/A


def womersley(alpha):
    """Series impedance of the pipe, and Womersley's F10 and Lambda"""
    lam = 1j**1.5*alpha
    f10 = 2*jv(1, lam)/(lam*jv(0, lam))
    return 1j*OMEGA*RHO*L3D/(A*(1 - f10)), f10, lam


def exact(alpha):
    """Complex amplitudes (e^{i w t} convention) of P_in, Q, P_out at the
    ports for the incident wave, viscous and inviscid"""
    out = {}
    for name, Z in (("Womersley", womersley(alpha)[0]),
                    ("inertance", 1j*OMEGA*RHO*L3D/A)):
        gamma = Z/(Z + 2*ZC)
        T = 2*ZC/(Z + 2*ZC)
        # p(t) = P sin(w t) = Re{-i P e^{i w t}}, delayed by L1/c0
        p = -1j*P_INC*cmath.exp(-1j*OMEGA*L1D/C0)
        out[name] = dict(Gamma=gamma, T=T, P_in=(1 + gamma)*p, Q=T*p/ZC,
                         P_out=T*p)
    return out


def fit(t, x):
    """Complex amplitude X of x(t) = Re{X e^{i w t}} + c over the samples"""
    M = np.column_stack([np.cos(OMEGA*t), np.sin(OMEGA*t), np.ones_like(t)])
    a, b, c = np.linalg.lstsq(M, x, rcond=None)[0]
    return a - 1j*b, c


def load_ports(case):
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
    return t, {"P_in": a[:, col["to3D.P"]], "Q": a[:, col["to3D.Qsampled"]],
               "P_out": a[:, col["from3D.P"]]}


def load_profiles(case):
    """{time: (r, u_z)} from the sets written over the last period"""
    out = {}
    d = RUNS/case/"postProcessing"/"profile"
    if not d.is_dir():
        return out
    for td in sorted(d.iterdir(), key=lambda p: float(p.name)):
        f = td/"midLength_U.xy"
        if f.is_file():
            a = np.loadtxt(f)
            out[float(td.name)] = (a[:, 0], a[:, 3])     # x, Ux, Uy, Uz
    return out


def womersley_profile(alpha, Qhat, r, t):
    """u(r, t) of Womersley's solution with flow amplitude Qhat"""
    _, f10, lam = womersley(alpha)
    uhat = (Qhat/A)*(1 - jv(0, lam*r/R)/jv(0, lam))/(1 - f10)
    return np.real(uhat*cmath.exp(1j*OMEGA*t))


def main():
    lines = ["# The Womersley segment: ports and profiles against the exact "
             "periodic solution", "",
             f"Incident wave {P_INC} Pa at {F} Hz; Zc = {ZC:.4g} Pa s/m3; "
             f"c0 = {C0:.4g} m/s.", ""]
    figs = {}
    for case in sorted(p.name for p in RUNS.glob("alpha*") if p.is_dir()):
        alpha = float(case[5:])
        log = RUNS/case/"log.pimpleFoam"
        if not (log.is_file() and "\nEnd\n" in log.read_text()[-2000:]):
            lines += [f"## {case}: not finished", ""]
            continue
        t, ports = load_ports(case)
        T = 1/F
        m = t > t[-1] - 2*T
        ex = exact(alpha)
        Zw = womersley(alpha)[0]
        lines += [f"## {case}: alpha = {alpha:g}, nu = {nu(alpha):.3g} m2/s, "
                  f"Z = {Zw.real:.3g} + {Zw.imag:.3g}i Pa s/m3 "
                  f"(inertance {OMEGA*RHO*L3D/A:.3g}i)", "",
                  f"|Gamma| = {abs(ex['Womersley']['Gamma']):.4f} "
                  f"(inviscid {abs(ex['inertance']['Gamma']):.4f}), "
                  f"|T| = {abs(ex['Womersley']['T']):.4f} "
                  f"(inviscid {abs(ex['inertance']['T']):.4f})", "",
                  "| quantity | amplitude: run / Womersley / inertance "
                  "| phase [deg]: run / Womersley / inertance "
                  "| amplitude error | phase error [deg] |",
                  "|---|---|---|---|---|"]
        fitted = {}
        for key, unit in (("P_in", "Pa"), ("Q", "m3/s"), ("P_out", "Pa")):
            X, c = fit(t[m], ports[key][m])
            fitted[key] = X
            W = ex["Womersley"][key]
            I = ex["inertance"][key]
            ph = lambda z: math.degrees(cmath.phase(z))
            dph = (ph(X) - ph(W) + 180) % 360 - 180
            lines.append(
                f"| {key} [{unit}] | {abs(X):.4g} / {abs(W):.4g} / {abs(I):.4g} "
                f"| {ph(X):.2f} / {ph(W):.2f} / {ph(I):.2f} "
                f"| {abs(X - W)/abs(W):.2e} | {dph:.3f} |")
        lines.append("")
        profiles = load_profiles(case)
        if profiles:
            errs = []
            for tp, (r, u) in profiles.items():
                ue = womersley_profile(alpha, fitted["Q"], r, tp)
                errs.append(np.abs(u - ue).max())
            umax = max(np.abs(u).max() for _, (r, u) in profiles.items())
            lines += [f"Mid-length velocity profile over the last period, "
                      f"{len(profiles)} phases: max |u - u_Womersley| / max |u| "
                      f"= {max(errs)/umax:.2e}", ""]
        figs[case] = (alpha, t, ports, m, fitted, ex, profiles)
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
    n = len(figs)
    if n == 0:
        return
    # Profiles
    fig, axs = plt.subplots(1, n, figsize=(4.2*n, 3.6), sharey=False, squeeze=False)
    for ax, (case, (alpha, t, ports, m, fitted, ex, profiles)) in zip(axs[0], figs.items()):
        for k, (tp, (r, u)) in enumerate(sorted(profiles.items())):
            ue = womersley_profile(alpha, fitted["Q"], r, tp)
            col = plt.cm.viridis(k/max(len(profiles) - 1, 1))
            ax.plot(1e3*r, ue, "-", color=col, lw=1)
            ax.plot(1e3*r[::8], u[::8], "o", color=col, ms=3, mfc="none")
        ax.set_title(f"alpha = {alpha:g}", fontsize=10)
        ax.set_xlabel("r [mm]")
        ax.grid(alpha=0.3)
    axs[0, 0].set_ylabel("u_z at mid-length [m/s]")
    axs[0, 0].plot([], [], "k-", label="Womersley (fitted Q)")
    axs[0, 0].plot([], [], "ko", mfc="none", label="pimpleFoam")
    axs[0, 0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(HERE/"plots"/"profiles.png", dpi=150)
    plt.close(fig)
    # Ports over the last two periods
    fig, axs = plt.subplots(n, 3, figsize=(12, 2.8*n), squeeze=False)
    for row, (case, (alpha, t, ports, m, fitted, ex, profiles)) in zip(axs, figs.items()):
        for ax, (key, scale, unit) in zip(row, (("P_in", 1, "Pa"), ("Q", 1e6, "ml/s"),
                                                 ("P_out", 1, "Pa"))):
            tt = t[m]
            ax.plot(tt, scale*ports[key][m], "k-", lw=1.2, label="run")
            W = ex["Womersley"][key]
            ax.plot(tt, scale*np.real(W*np.exp(1j*OMEGA*tt)), "--", color="#D55E00",
                    lw=1.2, label="Womersley line")
            I = ex["inertance"][key]
            ax.plot(tt, scale*np.real(I*np.exp(1j*OMEGA*tt)), ":", color="#0072B2",
                    lw=1, label="inviscid inertance")
            ax.set_title(f"alpha = {alpha:g}: {key} [{unit}]", fontsize=9)
            ax.grid(alpha=0.3)
        row[0].legend(frameon=False, fontsize=7)
    for ax in axs[-1]:
        ax.set_xlabel("t [s]")
    fig.tight_layout()
    fig.savefig(HERE/"plots"/"ports.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
