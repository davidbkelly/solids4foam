#!/usr/bin/env python3
"""The FSI tube with pulses of increasing length: coupled against the full
tube as the pulse grows long against the radius.

    python3 pulseSweep.py run [--durations 3,6,12,24] [--jobs 2] [--ranks 4]
    python3 pulseSweep.py            # analyse what has run

The tutorial's pulse (1333 Pa for 3 ms) is about three radii long
(c0 T = 16 mm at c0 = 5.25 m/s, R = 5 mm), so the 1D long-wave assumption
does not hold in it and the coupled-against-full comparison measures the
model mismatch as much as the coupling. Here the pulse is smooth,
p(t) = PULSE sin^2(pi t/T) for 0 <= t <= T, of duration T = 3, 6, 12 and
24 ms, on the inlet of (i) the full FSI tube, long enough that nothing
reflected from its far end (p = 0) reaches z = 12.5 cm within the window
(15, 20, 20 and 25 cm), and (ii) the coupled 5 cm tube continued by the 1D
vessel (the wall's cut end clamped, c0 from the static plane-strain match,
as the tutorial). The window is T + 30 ms (the tail of the pulse past
z = 12.5 cm at 5.25 m/s plus a margin), dt = 25 us as the tutorial.

Reported per duration: RMS(P coupled - P full) over the window relative
to the pulse at every section, the peak pressures, the wall displacement
at z = 2.5 cm, against the ratio R/(c0 T) of radius to pulse length.
Writes pulse.md and plots/pulse.png. Runs: full3D_T<ms> and
coupled_c<c0>_T<ms> (analyse.py compares them too).
"""
import argparse
import math
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import makeCases as mc  # noqa: E402
import analyse as an  # noqa: E402

RUNS = HERE/"runs"
DURATIONS_MS = [3, 6, 12, 24]
LENGTH_M = {3: 0.15, 6: 0.20, 12: 0.20, 24: 0.25}   # the full tube
MARGIN = 0.030                                       # window = T + margin


def c0():
    return mc.c0_plane_strain()


def names(T_ms):
    c = "coupled_c" + f"{c0():.3g}".replace(".", "p")
    return f"full3D_T{T_ms:g}", f"{c}_T{T_ms:g}"


def pulse_table(T):
    """The inlet's time series: a sin^2 pulse of duration T, kinematic
    pressure (p/rho) as the tutorial's file"""
    n = 200
    t = np.linspace(0, T, n + 1)
    v = mc.PULSE/mc.RHO*np.sin(math.pi*t/T)**2
    lines = ["("]
    lines += [f"    ({ti:.8g} {vi:.8g})" for ti, vi in zip(t, v)]
    lines += [f"    ({T + 1e-6:.8g} 0)", "    (1 0)", ")"]
    return "\n".join(lines) + "\n"


def setup(T_ms, force=False):
    T = T_ms*1e-3
    end = T + MARGIN
    full, coupled = names(T_ms)
    out = []
    for name, length, is_coupled in ((full, LENGTH_M[T_ms], False),
                                     (coupled, mc.L3D, True)):
        case = RUNS/name
        if an.finished(name, "log.solids4Foam") and not force:
            out.append(case)
            continue
        mc.make_3d(name, length, is_coupled, c0(), mc.DELTA_T, True, end=end)
        (case/"0"/"fluid"/"time-series").write_text(pulse_table(T))
        out.append(case)
    return out


def run(case, ranks):
    if an.finished(case.name, "log.solids4Foam"):
        return
    def sh(cmd, log):
        with open(case/log, "w") as f:
            r = subprocess.run(cmd, cwd=case, stdout=f, stderr=subprocess.STDOUT)
        if r.returncode:
            sys.exit(f"{case.name}: {' '.join(cmd)} failed, see {case/log}")
    sh(["blockMesh", "-region", "solid"], "log.blockMesh.solid")
    sh(["blockMesh", "-region", "fluid"], "log.blockMesh.fluid")
    sh(["topoSet", "-region", "fluid"], "log.topoSet")
    if ranks > 1 and not case.name.startswith("coupled"):
        sh(["decomposePar", "-region", "fluid"], "log.decomposePar.fluid")
        sh(["decomposePar", "-region", "solid"], "log.decomposePar.solid")
        with open(case/"log.solids4Foam", "w") as f:
            subprocess.run(["mpirun", "--oversubscribe", "-np", str(ranks), "solids4Foam", "-parallel"],
                           cwd=case, stdout=f, stderr=subprocess.STDOUT)
    else:
        with open(case/"log.solids4Foam", "w") as f:
            subprocess.run(["solids4Foam"], cwd=case, stdout=f,
                           stderr=subprocess.STDOUT)
    print(f"{case.name}: {'finished' if an.finished(case.name, 'log.solids4Foam') else 'FAILED'}",
          flush=True)


def analyse(durations):
    rows = []
    series = {}
    for T_ms in durations:
        full, coupled = names(T_ms)
        if not (an.finished(full, "log.solids4Foam")
                and an.finished(coupled, "log.solids4Foam")):
            continue
        ref = an.sections_3d(full)
        s = an.sections_3d(coupled)
        s.update(an.sections_coupler(coupled))
        T = T_ms*1e-3
        lam = c0()*T
        errs, peaks = {}, {}
        for zcm in mc.PLANES_CM:
            if zcm in ref and zcm in s:
                tf, Pf, _ = ref[zcm]
                t, P, _ = s[zcm]
                errs[zcm] = an.rms(t, P, tf, Pf)/mc.PULSE
                peaks[zcm] = (Pf.max(), P.max())
        tf, Df = an.displacement(full)
        tc, Dc = an.displacement(coupled)
        rows.append((T_ms, LENGTH_M[T_ms], lam, mc.R/lam, errs, peaks,
                     Df.max(), Dc.max(), an.rms(tc, Dc, tf, Df)))
        series[T_ms] = (ref, s)
    zs = [z for z in mc.PLANES_CM if any(z in r[4] for r in rows)]
    lines = ["# The FSI tube with pulses of increasing length", "",
             f"sin^2 pulses of {mc.PULSE:.0f} Pa and duration T on the full "
             "FSI tube (length L) and on the coupled 5 cm tube + 1D vessel "
             f"(clamped wall end, c0 = {c0():.3g} m/s); window T + "
             f"{MARGIN*1e3:g} ms, dt = {mc.DELTA_T*1e6:g} us. RMS(P coupled - "
             "P full) over the window / pulse at each section (z <= 5 cm in "
             "the 3D part, 5 cm the port, beyond it the 1D vessel against the "
             "full tube), and the radial wall displacement at z = 2.5 cm.", "",
             "| T [ms] | L [cm] | c0 T [mm] | R/(c0 T) | "
             + " | ".join(f"z = {z:g}" for z in zs)
             + " | max, 3D part | D peak full / coupled [um] | RMS dD / peak |",
             "|---|---|---|---|" + "---|"*len(zs) + "---|---|---|"]
    for r in rows:
        T_ms, L, lam, ratio, errs, peaks, Df, Dc, dD = r
        e3 = max(v for z, v in errs.items() if z*1e-2 <= mc.L3D)
        lines.append(f"| {T_ms:g} | {L*100:g} | {lam*1e3:.1f} | {ratio:.3f} | "
                     + " | ".join(f"{errs[z]:.2e}" if z in errs else "-" for z in zs)
                     + f" | {e3:.2e} | {1e6*Df:.1f} / {1e6*Dc:.1f} | {dD/Df:.2e} |")
    lines += ["", "Peak pressure, full / coupled [Pa]:", "",
              "| T [ms] | " + " | ".join(f"z = {z:g}" for z in zs) + " |",
              "|---|" + "---|"*len(zs)]
    for r in rows:
        lines.append(f"| {r[0]:g} | " + " | ".join(
            f"{r[5][z][0]:.0f} / {r[5][z][1]:.0f}" if z in r[5] else "-" for z in zs) + " |")
    # The difference split: the full tube's wave speed (half-peak arrival
    # over z = 1.25-5 cm) against the 1D c0; at each section the peak ratio,
    # the time shift that best aligns the coupled signal with the full one,
    # and the RMS difference left after that shift (shape and amplitude)
    lines += ["", "The difference split. c_3D: the full tube's wave speed from "
              "the half-peak arrivals at z = 1.25-5 cm (1D: "
              f"{c0():.3g} m/s). At each section: peak ratio coupled / full; "
              "best-aligning time shift [ms]; RMS / pulse after the shift.", "",
              "| T [ms] | c_3D [m/s] | " + " | ".join(f"z = {z:g}" for z in zs) + " |",
              "|---|---|" + "---|"*len(zs)]
    for T_ms in sorted(series):
        ref, s = series[T_ms]
        c3, _, _ = an.wave_speed({z: ref[z] for z in ref if z*1e-2 <= mc.L3D})
        cells = []
        for z in zs:
            if z not in ref or z not in s:
                cells.append("-")
                continue
            tf, Pf, _ = ref[z]
            t, P, _ = s[z]
            best = min((math.sqrt(np.mean((np.interp(tf, t + sh, P) - Pf)**2)), sh)
                       for sh in np.arange(-5e-3, 5e-3, 2.5e-5))
            cells.append(f"{P.max()/Pf.max():.3f}; {1e3*best[1]:+.2f}; "
                         f"{best[0]/mc.PULSE:.1e}")
        lines.append(f"| {T_ms:g} | {c3:.2f} | " + " | ".join(cells) + " |")
    text = "\n".join(lines) + "\n"
    (HERE/"pulse.md").write_text(text)
    print(text)
    plot(rows, series)


def plot(rows, series):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    if not rows:
        return
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
    ax = axs[0]
    ratios = [r[3] for r in rows]
    for z in (2.5, 5, 7.5):
        e = [r[4].get(z, math.nan) for r in rows]
        ax.loglog(ratios, e, "o-", label=f"P at z = {z:g} cm")
    ax.loglog(ratios, [r[8]/r[6] for r in rows], "s--", label="wall displacement, z = 2.5 cm")
    ax.set_xlabel("R / (c0 T)")
    ax.set_ylabel("RMS difference / pulse (or / peak)")
    ax.set_title("coupled against the full tube", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.3, which="both")
    ax = axs[1]
    for T_ms, (ref, s) in series.items():
        if 5 in ref and 5 in s:
            tf, Pf, _ = ref[5]
            t, P, _ = s[5]
            l, = ax.plot(tf*1e3, Pf, lw=1.2, label=f"T = {T_ms:g} ms, full")
            ax.plot(t*1e3, P, "--", lw=1.0, color=l.get_color())
    ax.set_xlabel("t [ms]")
    ax.set_ylabel("P at z = 5 cm (the port) [Pa]")
    ax.set_title("full (solid) and coupled (dashed)", fontsize=10)
    ax.legend(frameon=False, fontsize=7)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    (HERE/"plots").mkdir(exist_ok=True)
    fig.savefig(HERE/"plots"/"pulse.png", dpi=150)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("what", nargs="?", choices=["run", "setup"])
    ap.add_argument("--durations", default=",".join(str(x) for x in DURATIONS_MS))
    ap.add_argument("--jobs", type=int, default=2, help="runs at a time")
    ap.add_argument("--ranks", type=int, default=4, help="ranks of the full-tube runs")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    durations = [float(x) for x in args.durations.split(",")]
    durations = [int(x) if x == int(x) else x for x in durations]
    if args.what in ("run", "setup"):
        RUNS.mkdir(exist_ok=True)
        cases = [c for T in durations for c in setup(T, args.force)]
        if args.what == "setup":
            return
        cases.sort(key=lambda c: -LENGTH_M.get(int(c.name.split("_T")[1]), 0))
        with ThreadPoolExecutor(args.jobs) as ex:
            list(ex.map(lambda c: run(c, args.ranks), cases))
    analyse(durations)


if __name__ == "__main__":
    main()
