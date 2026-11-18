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


# The wedge mesh has its flat outer face at the distance R from the axis
# and its patch flux scaled by pi/tan(theta/2) (makeCases.py), so the
# finite-volume problem it discretises is the circular pipe of radius R
# exactly: the exact solutions below are those of the pipe of area pi R^2.
# (Before 17 November 2026 the edge vertices were at R, which made it the
# pipe of radius R cos(theta/2) with every flux times 2 tan(theta/2)/theta.)
# exact_transient() and exact_transient_profile() take the pipe's area and
# radius as arguments (defaults: pi R^2 and R).


def impedance(w, alpha, A=A, R=R):
    """Womersley's series impedance of the pipe at any angular frequency
    w >= 0, for the viscosity of the run at Womersley number alpha (at the
    forcing frequency); the Poiseuille resistance at w = 0. Exponentially
    scaled Bessel functions keep the ratio finite at high frequency."""
    from scipy.special import jve
    w = np.asarray(w, dtype=float)
    visc = nu(alpha)
    Z = np.empty(w.shape, dtype=complex)
    dc = w == 0
    Z[dc] = 8*RHO*visc*L3D/(math.pi*R**4)
    lam = 1j**1.5*R*np.sqrt(w[~dc]/visc)
    f10 = 2*jve(1, lam)/(lam*jve(0, lam))
    Z[~dc] = 1j*w[~dc]*RHO*L3D/(A*(1 - f10))
    return Z


def exact_transient(alpha, t, dtg=2.5e-5, Tpad=256.0, Tin=2.0, A=A, R=R):
    """The exact response of the ports to the inlet's ramped incident wave
    from rest, {to3D.P, to3D.Q, from3D.P, from3D.Q} at the times t: the
    linear problem is a single reflection and transmission at the pipe,
    with Gamma(w) = Z/(Z + 2 Zc) and T(w) = 2 Zc/(Z + 2 Zc) at every
    frequency, so the ports are the convolution of the incident wave
    (delayed by L1/c0) with the impulse responses, done by FFT on a grid
    dtg over Tpad seconds. The incident wave is switched off after Tin
    seconds, which by causality does not change the response before Tin
    and lets it decay before the periodic window wraps. The grid error is
    below 1e-10 of the signal (dtg = 2.5e-5 against 5e-5 and 1e-4)."""
    n = int(round(Tpad/dtg))
    tg = np.arange(n)*dtg
    tau = tg - L1D/C0
    pinc = np.where((tau >= 0) & (tau <= Tin),
                    P_INC*np.sin(OMEGA*tau)*np.minimum(1.0, F*tau), 0.0)
    w = 2*np.pi*np.fft.rfftfreq(n, dtg)
    Z = impedance(w, alpha, A, R)
    gamma = Z/(Z + 2*ZC)
    T = 2*ZC/(Z + 2*ZC)
    P = np.fft.rfft(pinc)
    out = {}
    for key, H in (("to3D.P", 1 + gamma), ("to3D.Q", T/ZC),
                   ("from3D.P", T), ("from3D.Q", T/ZC)):
        out[key] = np.interp(t, tg, np.fft.irfft(P*H, n))
    return out


def exact_transient_profile(alpha, r, times, dtg=1e-4, Tpad=128.0, Tin=2.0,
                            A=A, R=R):
    """The exact axial velocity u(r, t) at the pipe's mid-length at the
    given times, for the same transient: the flow's spectrum times
    Womersley's profile at each frequency, (1 - J0(L r/R)/J0(L))/(1 - F10)
    per unit mean velocity (Poiseuille's 2 (1 - (r/R)^2) at w = 0).
    Returns an array [len(times), len(r)]."""
    from scipy.special import jve
    n = int(round(Tpad/dtg))
    tg = np.arange(n)*dtg
    tau = tg - L1D/C0
    pinc = np.where((tau >= 0) & (tau <= Tin),
                    P_INC*np.sin(OMEGA*tau)*np.minimum(1.0, F*tau), 0.0)
    w = 2*np.pi*np.fft.rfftfreq(n, dtg)
    Z = impedance(w, alpha, A, R)
    Qhat = np.fft.rfft(pinc)*(2*ZC/(Z + 2*ZC))/ZC
    visc = nu(alpha)
    lam = 1j**1.5*R*np.sqrt(w[1:]/visc)
    f10 = 2*jve(1, lam)/(lam*jve(0, lam))
    out = np.empty((len(times), len(r)))
    for j, rj in enumerate(r):
        phi = np.empty(w.shape, dtype=complex)
        phi[0] = 2*(1 - (rj/R)**2)
        ratio = jve(0, lam*rj/R)/jve(0, lam)*np.exp(
            np.abs((lam*rj/R).imag) - np.abs(lam.imag))
        phi[1:] = (1 - ratio)/(1 - f10)
        u = np.fft.irfft(Qhat*phi/A, n)
        out[:, j] = np.interp(times, tg, u)
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


def wall_cell(level):
    """radial width of the wall cell of the wedge mesh of level k (nr = 50 k
    cells over R graded 10:1 towards the wall, as spaceStudy.py)"""
    nr = int(round(50*level))
    q = 0.1**(1.0/(nr - 1))
    return R*(1 - q)/(1 - q**nr)*q**(nr - 1)


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
        profiles = load_profiles(case)
        if profiles:
            errs = []
            for tp, (r, u) in profiles.items():
                # up to the wall cell (the tutorial's mesh is level 1 of
                # spaceStudy.py): inside it the sampled values converge at
                # first order only (spaceStudy.wall_cell_r)
                inside = r <= R - wall_cell(1.0)
                ue = womersley_profile(alpha, fitted["Q"], r[inside], tp)
                errs.append(np.abs(u[inside] - ue).max())
            umax = max(np.abs(u).max() for _, (r, u) in profiles.items())
            lines += [f"Mid-length velocity profile over the last period, "
                      f"{len(profiles)} phases, from the axis to the wall cell: "
                      f"max |u - u_Womersley| / max |u| "
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


def write_csv(outdir=None):
    """The exact references as CSV for alpha = 3, 6 and 12: the ports'
    transient from rest over the first 6 s (every 0.1 ms), the periodic
    complex amplitudes, and the mid-length velocity profile at eight
    instants of the first second"""
    import os
    outdir = outdir or str(HERE/"reference")
    os.makedirs(outdir, exist_ok=True)
    for alpha in (3, 6, 12):
        hdr = (f"# Womersley segment, exact linear solution: incident wave {P_INC} sin(2 pi {F} t) min(1, {F} t) Pa "
               f"at the proximal inlet (non-reflecting), vessels L1 = L2 = {L1D} m, R {R} m, c0 {C0:.10g} m/s, "
               f"Zc {ZC:.10g} Pa s/m3, rho {RHO}; pipe l {L3D} m, R {R} m, no slip, nu {nu(alpha):.10g} m2/s "
               f"(alpha {alpha} at {F} Hz); from rest")
        t = np.round(np.arange(0, 6 + 1e-9, 1e-4), 6)
        ex = exact_transient(alpha, t, Tin=20.0)
        with open(os.path.join(outdir, f"alpha{alpha}_ports.csv"), "w") as fh:
            fh.write(hdr + "\nt [s],P_in [Pa],Q [m3/s],P_out [Pa]\n")
            for row in zip(t, ex["to3D.P"], ex["to3D.Q"], ex["from3D.P"]):
                fh.write(",".join(f"{v:.12g}" for v in row) + "\n")
        w = exact(alpha)["Womersley"]
        with open(os.path.join(outdir, f"alpha{alpha}_periodic.csv"), "w") as fh:
            fh.write(hdr + "\n# periodic state, x(t) = Re{X exp(i w t)}: quantity,|X|,arg(X) [deg]\n")
            for k in ("P_in", "Q", "P_out"):
                fh.write(f"{k},{abs(w[k]):.12g},{math.degrees(cmath.phase(w[k])):.12g}\n")
            fh.write(f"Gamma,{abs(w['Gamma']):.12g},{math.degrees(cmath.phase(w['Gamma'])):.12g}\n")
            fh.write(f"T,{abs(w['T']):.12g},{math.degrees(cmath.phase(w['T'])):.12g}\n")
        r = np.linspace(0, R, 101)
        times = [0.5 + k/16 for k in range(8)]
        u = exact_transient_profile(alpha, r, times, Tin=20.0)
        with open(os.path.join(outdir, f"alpha{alpha}_profiles.csv"), "w") as fh:
            fh.write(hdr + "\n# axial velocity u(r, t) [m/s] at the pipe's mid-length; r [m] then one column per t = "
                     + ", ".join(f"{tk:g}" for tk in times) + " s\n")
            for j, rj in enumerate(r):
                fh.write(f"{rj:.12g}," + ",".join(f"{u[i, j]:.12g}" for i in range(len(times))) + "\n")
    print(f"written to {outdir}/")


if __name__ == "__main__":
    if "--csv" in sys.argv:
        write_csv()
    else:
        main()
