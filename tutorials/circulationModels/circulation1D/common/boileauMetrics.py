#!/usr/bin/env python3
"""Error metrics and figures of Boileau et al. (2015) for the circulation1D tutorials.

    python3 boileauMetrics.py selftest      # recompute Tables IV, VI, VIII and X
                                            # from the published data of the schemes

Imported by the plot.py of each case. The metrics are Eqs. 10-13 of

    E. Boileau et al., A benchmark study of numerical schemes for one-dimensional
    arterial blood flow modelling, Int. J. Numer. Meth. Biomed. Engng. 31 (2015)
    e02732,

on the time points of the reference (1 kHz), over one cardiac cycle:

    RMS   sqrt(mean(((y - r)/N)^2))      MAX   max|y - r|/N
    SYS   |max y - max r|/max r          DIAS  |min y - min r|/N'

with N = r pointwise for the pressure P (N' = min r) and N = N' = max r for the
flow Q, the change in radius dr and the inlet-outlet pressure difference dP.
They are given in per cent, as magnitudes, as in the tables.

Time alignment. The data supplement does not state how the 3D cycle lines up with
the 1D runs, and the "last cycle" windows of its READMEs do not reproduce the
tables. The windows that do, to the two decimals published for the FVM and McC
schemes (and LCG, FEM where their files allow), are

    carotid      3D(tau) <-> 1D(k T + 0.056 + tau),  T = 1.1 s    (k = 8)
    aorta        3D(tau) <-> 1D(k T + 0.002 + tau),  T = 0.955 s  (k = 19)
    bifurcation  3D(tau) <-> 1D(k T + 0.001 + tau),  T = 1.1 s    (k = 26)
    37 arteries  in vitro(tau) <-> 1D(tini - 0.001 + tau), T = 0.827 s,
                 tini of each site in timeShifts.txt

with the inflow formulas of the supplement used as given from t = 0. The phase
is kept and k is the last whole cycle of the run ("selftest" shows the match).
The right anterior tibial window of Table X runs past the end of the published
1D data (15 s), so that row does not reproduce.

Needs numpy and matplotlib.
"""
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

SCHEMES = ["DCG", "LCG", "FEM", "FVM", "McC", "STM"]
QUANTITIES = ["P", "Q", "dr", "dP"]
METRICS = ["RMS", "MAX", "SYS", "DIAS"]

# Units of the reference columns, and the factor from the SI solver output
REF_UNITS = {"P": "kPa", "Q": "ml/s", "dr": "mm", "dP": "kPa"}
FROM_SI = {"P": 1e-3, "Q": 1e6, "dr": 1e3, "dP": 1e-3}

# Tables IV and VI of the paper [%], in the order of SCHEMES
TABLES = {
    "commonCarotid": {
        ("RMS", "P"): [0.20, 0.24, 0.23, 0.28, 0.24, 0.24],
        ("RMS", "Q"): [0.39, 0.26, 0.28, 0.29, 0.30, 0.29],
        ("RMS", "dr"): [0.96, 1.00, 1.00, 1.07, 1.01, 0.98],
        ("RMS", "dP"): [4.42, 4.35, 4.27, 4.21, 4.47, 4.47],
        ("MAX", "P"): [0.34, 0.52, 0.50, 0.66, 0.52, 0.53],
        ("MAX", "Q"): [1.23, 1.02, 1.12, 1.07, 1.20, 1.13],
        ("MAX", "dr"): [1.75, 1.93, 1.88, 2.28, 1.93, 1.96],
        ("MAX", "dP"): [16.69, 16.17, 15.88, 15.60, 16.58, 16.90],
        ("SYS", "P"): [0.27, 0.27, 0.26, 0.29, 0.27, 0.26],
        ("SYS", "Q"): [0.55, 0.46, 0.54, 0.58, 0.62, 0.55],
        ("SYS", "dr"): [1.65, 1.63, 1.61, 1.68, 1.63, 1.63],
        ("SYS", "dP"): [15.05, 15.37, 15.16, 14.83, 15.73, 15.76],
        ("DIAS", "P"): [0.28, 0.28, 0.27, 0.29, 0.28, 0.27],
        ("DIAS", "Q"): [0.25, 0.22, 0.26, 0.24, 0.23, 0.27],
        ("DIAS", "dr"): [0.08, 0.11, 0.08, 0.12, 0.09, 0.08],
        ("DIAS", "dP"): [4.99, 4.76, 5.03, 4.91, 4.45, 4.85],
    },
    "upperThoracicAorta": {
        ("RMS", "P"): [1.09, 1.21, 1.09, 1.14, 1.08, 1.12],
        ("RMS", "Q"): [2.55, 2.19, 2.23, 2.17, 2.22, 2.28],
        ("RMS", "dr"): [2.32, 2.41, 2.33, 2.44, 2.33, 2.90],
        ("RMS", "dP"): [7.70, 7.23, 7.23, 7.14, 7.23, 7.52],
        ("MAX", "P"): [3.54, 3.55, 3.26, 3.18, 3.27, 3.30],
        ("MAX", "Q"): [9.13, 7.20, 7.04, 7.07, 7.04, 7.40],
        ("MAX", "dr"): [7.80, 7.67, 7.32, 7.20, 7.09, 8.04],
        ("MAX", "dP"): [31.83, 29.31, 29.15, 29.13, 29.31, 31.19],
        ("SYS", "P"): [0.57, 0.32, 0.56, 0.71, 0.55, 0.59],
        ("SYS", "Q"): [5.56, 5.36, 5.63, 5.29, 5.58, 5.56],
        ("SYS", "dr"): [2.54, 1.93, 2.50, 2.86, 2.49, 4.48],
        ("SYS", "dP"): [8.92, 8.89, 9.29, 8.53, 9.04, 9.70],
        ("DIAS", "P"): [0.85, 1.12, 0.88, 0.99, 0.87, 0.85],
        ("DIAS", "Q"): [2.75, 2.67, 3.09, 3.48, 2.67, 12.95],
        ("DIAS", "dr"): [1.85, 2.24, 1.89, 2.05, 1.89, 0.00],
        ("DIAS", "dP"): [7.32, 7.24, 7.60, 6.37, 7.32, 7.47],
    },
}

# Table VIII [%]: the aortic bifurcation at three sites, order of SCHEMES
TABLE_VIII = {
    "aortaMid": {
        ("RMS", "P"): [0.37, 0.39, 0.37, 0.40, 0.38, 0.38],
        ("RMS", "Q"): [0.89, 0.71, 0.82, 0.81, 0.83, 0.82],
        ("RMS", "dr"): [2.48, 2.43, 2.55, 2.49, 2.48, 2.45],
        ("MAX", "P"): [0.61, 0.67, 0.65, 0.68, 0.66, 0.66],
        ("MAX", "Q"): [2.51, 2.24, 2.53, 2.44, 2.56, 2.52],
        ("MAX", "dr"): [4.05, 3.90, 4.03, 4.01, 3.98, 3.94],
        ("SYS", "P"): [0.53, 0.54, 0.51, 0.55, 0.53, 0.51],
        ("SYS", "Q"): [2.47, 2.23, 2.52, 2.44, 2.55, 2.52],
        ("SYS", "dr"): [3.96, 3.86, 3.99, 4.00, 3.94, 3.93],
        ("DIAS", "P"): [0.50, 0.51, 0.47, 0.53, 0.49, 0.48],
        ("DIAS", "Q"): [1.16, 1.03, 1.16, 1.16, 1.17, 1.17],
        ("DIAS", "dr"): [1.39, 1.37, 1.49, 1.35, 1.40, 1.41],
    },
    "aortaEnd": {
        ("RMS", "P"): [0.42, 0.43, 0.41, 0.44, 0.42, 0.42],
        ("RMS", "Q"): [1.20, 1.13, 1.12, 1.08, 1.12, 0.46],
        ("RMS", "dr"): [4.07, 4.02, 4.14, 4.08, 4.08, 4.05],
        ("MAX", "P"): [0.73, 0.75, 0.72, 0.80, 0.74, 0.73],
        ("MAX", "Q"): [3.55, 3.54, 3.53, 3.28, 3.53, 1.48],
        ("MAX", "dr"): [6.87, 6.76, 6.89, 6.95, 6.85, 6.98],
        ("SYS", "P"): [0.71, 0.70, 0.69, 0.72, 0.70, 0.69],
        ("SYS", "Q"): [3.47, 3.53, 3.52, 3.26, 3.52, 1.08],
        ("SYS", "dr"): [6.84, 6.73, 6.87, 6.88, 6.83, 6.90],
        ("DIAS", "P"): [0.54, 0.55, 0.53, 0.57, 0.54, 0.53],
        ("DIAS", "Q"): [1.76, 1.78, 1.76, 1.74, 1.77, 0.69],
        ("DIAS", "dr"): [1.93, 1.92, 2.02, 1.89, 1.94, 1.88],
    },
    "iliacMid": {
        ("RMS", "P"): [0.45, 0.45, 0.44, 0.47, 0.45, 0.44],
        ("RMS", "Q"): [0.92, 0.66, 0.65, 0.68, 0.67, 0.66],
        ("RMS", "dr"): [4.29, 4.30, 4.35, 4.31, 4.31, 4.26],
        ("MAX", "P"): [0.84, 0.86, 0.84, 0.92, 0.86, 0.85],
        ("MAX", "Q"): [2.33, 2.00, 2.01, 1.80, 2.01, 2.02],
        ("MAX", "dr"): [7.29, 7.36, 7.37, 7.48, 7.35, 7.34],
        ("SYS", "P"): [0.83, 0.83, 0.81, 0.85, 0.82, 0.81],
        ("SYS", "Q"): [1.68, 1.84, 1.81, 1.58, 1.89, 1.79],
        ("SYS", "dr"): [7.25, 7.26, 7.27, 7.31, 7.25, 7.22],
        ("DIAS", "P"): [0.55, 0.55, 0.54, 0.57, 0.54, 0.54],
        ("DIAS", "Q"): [1.17, 1.19, 1.17, 1.18, 1.19, 1.18],
        ("DIAS", "dr"): [2.36, 2.37, 2.45, 2.35, 2.38, 2.41],
    },
}

# Table X [%]: the 37-artery network against in vitro data. Per site and
# scheme: RMS, MAX, SYS, DIAS of P, then of Q
TABLE_X_ORDER = [(m, q) for q in ("P", "Q") for m in METRICS]
TABLE_X = {
    "AorticArchII": {
        "DCG": [1.78, 3.59, 1.40, 0.42, 12.32, 29.58, 8.78, 17.95],
        "LCG": [1.68, 3.21, 1.02, 0.93, 12.34, 31.83, 10.02, 16.93],
        "FEM": [1.89, 3.78, 1.40, 0.17, 12.02, 29.00, 8.81, 17.08],
        "FVM": [1.87, 3.72, 1.46, 0.30, 12.11, 29.13, 8.75, 17.93],
        "McC": [1.94, 3.97, 1.60, 0.12, 12.11, 29.40, 8.81, 17.74],
        "STM": [1.84, 3.76, 1.39, 0.36, 12.19, 29.45, 8.75, 17.99],
    },
    "ThoracicAortaII": {
        "DCG": [2.36, 5.29, 0.96, 1.66, 25.59, 67.52, 60.70, 39.65],
        "LCG": [2.17, 5.03, 0.79, 1.13, 25.36, 70.75, 61.59, 35.47],
        "FEM": [2.49, 5.66, 0.97, 1.97, 25.26, 64.93, 60.39, 38.29],
        "FVM": [2.44, 5.57, 1.03, 1.81, 25.43, 65.69, 61.59, 39.26],
        "McC": [2.53, 5.70, 1.13, 1.98, 25.62, 66.22, 62.24, 38.74],
        "STM": [2.42, 5.58, 0.98, 1.69, 25.37, 65.32, 61.16, 38.77],
    },
    "LSubclavianI": {
        "DCG": [3.09, 6.13, 3.92, 3.97, 14.31, 38.39, 2.11, 11.57],
        "LCG": [3.12, 6.06, 3.38, 4.63, 13.87, 34.76, 2.92, 11.33],
        "FEM": [3.05, 6.10, 3.91, 3.72, 14.17, 37.89, 2.53, 10.86],
        "FVM": [3.12, 6.29, 3.96, 3.93, 14.24, 38.12, 2.25, 11.05],
        "McC": [3.11, 6.35, 4.09, 3.69, 14.31, 38.25, 2.08, 10.92],
        "STM": [3.11, 6.13, 3.80, 3.97, 14.45, 38.96, 1.18, 10.83],
    },
    "RIliacFemoralII": {
        "DCG": [3.82, 9.00, 2.20, 4.13, 24.49, 59.63, 50.20, 39.72],
        "LCG": [3.97, 9.69, 2.19, 5.57, 24.17, 61.14, 49.23, 36.92],
        "FEM": [3.69, 8.66, 2.32, 3.67, 23.90, 59.28, 49.99, 37.18],
        "FVM": [3.75, 9.09, 2.26, 4.04, 24.19, 60.61, 51.42, 39.06],
        "McC": [3.65, 8.95, 2.33, 3.73, 24.80, 61.61, 52.41, 40.21],
        "STM": [3.75, 9.32, 2.37, 4.17, 24.12, 60.41, 51.29, 38.26],
    },
    "LUlnar": {
        "DCG": [2.65, 7.18, 0.99, 2.68, 12.74, 30.14, 5.37, 17.48],
        "LCG": [2.57, 6.51, 1.45, 3.54, 12.42, 25.91, 4.30, 17.93],
        "FEM": [2.70, 7.29, 1.13, 2.34, 12.43, 27.75, 3.63, 16.14],
        "FVM": [2.74, 7.51, 1.09, 2.54, 12.70, 29.30, 4.60, 17.15],
        "McC": [2.75, 7.42, 1.07, 2.27, 12.78, 29.19, 4.43, 16.80],
        "STM": [2.74, 7.46, 0.81, 2.63, 12.91, 29.55, 5.68, 17.02],
    },
    "RAnteriorTibial": {
        "DCG": [3.25, 9.87, 0.57, 0.13, 10.49, 35.07, 8.37, 15.09],
        "LCG": [3.43, 12.24, 0.57, 0.80, 11.05, 35.80, 12.37, 14.59],
        "FEM": [3.21, 9.12, 0.55, 0.04, 9.88, 31.16, 8.40, 13.25],
        "FVM": [3.30, 9.42, 0.66, 0.14, 10.22, 33.57, 9.79, 14.75],
        "McC": [3.41, 10.33, 0.89, 0.03, 10.24, 34.38, 9.76, 14.43],
        "STM": [3.27, 9.22, 0.72, 0.15, 10.24, 34.54, 10.32, 14.07],
    },
    "RUlnar": {
        "DCG": [2.54, 6.32, 2.18, 3.91, 11.67, 31.30, 16.89, 3.47],
        "LCG": [2.58, 6.88, 1.92, 4.33, 11.47, 28.19, 17.56, 7.48],
        "FEM": [2.66, 6.36, 2.55, 3.96, 11.22, 29.28, 15.77, 4.08],
        "FVM": [2.49, 6.32, 2.20, 3.92, 11.62, 30.93, 16.50, 4.60],
        "McC": [2.42, 6.06, 2.15, 3.66, 11.63, 31.09, 16.58, 4.39],
        "STM": [2.50, 6.50, 2.06, 4.09, 11.73, 31.58, 17.40, 5.34],
    },
    "Splenic": {
        "DCG": [2.35, 6.25, 0.79, 0.14, 9.44, 23.74, 0.09, 5.97],
        "LCG": [2.34, 6.23, 0.09, 0.77, 9.79, 22.82, 3.59, 6.75],
        "FEM": [2.28, 5.95, 0.90, 0.04, 9.02, 22.83, 1.32, 5.46],
        "FVM": [2.33, 6.12, 0.87, 0.04, 9.22, 23.20, 0.58, 5.95],
        "McC": [2.22, 5.57, 0.97, 0.19, 9.04, 23.30, 1.47, 5.21],
        "STM": [2.36, 5.96, 0.38, 0.03, 9.56, 24.47, 1.05, 6.40],
    },
}
NETWORK_SITES = list(TABLE_X)
NETWORK_T = 0.827

# Bifurcation sites: probe name, reference file suffix, title
BIFURCATION_T = 1.1
BIFURCATION_PHASE = 0.001
BIFURCATION_SITES = [
    ("aortaMid", "AortaMid", "midpoint of the aorta"),
    ("aortaEnd", "Junction", "end of the aorta"),
    ("iliacMid", "IliacMid", "midpoint of an iliac artery"),
]

# Period, phase of the 3D cycle in the 1D time, names of the case's probe and
# vessel, and the table number in the paper
CASES = {
    "commonCarotid": dict(T=1.1, phase=0.056, probe="mid", vessel="carotid",
                          table="IV", figure="3", title="Common carotid artery"),
    "upperThoracicAorta": dict(T=0.955, phase=0.002, probe="mid", vessel="aorta",
                               table="VI", figure="4",
                               title="Upper thoracic aorta"),
}

# Figure colours: the six schemes, the 3D reference and solids4foam, in a fixed
# order (a colour-blind-checked categorical palette)
COLOURS = {
    "DCG": "#2a78d6", "LCG": "#eb6834", "FEM": "#1baf7a",
    "FVM": "#eda100", "McC": "#e87ba4", "STM": "#008300",
    "solids4foam": "#6250d6", "3D": "#0b0b0b", "theory": "#0b0b0b",
}
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e4e3df"


# --------------------------------------------------------------------------- data

def load_reference(path):
    """Columns t [s], P [kPa], Q [ml/s], dr [mm], dP [kPa] of a reference file."""
    return np.loadtxt(path, comments="#")


def load_solution(case_dir, probe, vessel):
    """The solver's probe and vessel output in the units of the reference.

    Returns t and a dict of the four quantities."""
    out = Path(case_dir) / "postProcessing" / "circulation1D" / "0"
    p = np.loadtxt(out / "probes" / f"{probe}.dat", comments="#")
    v = np.loadtxt(out / "vessels" / f"{vessel}.dat", comments="#")
    if len(p) != len(v):
        raise SystemExit(f"{out}: probe and vessel files differ in length")
    y = {
        "P": p[:, 1]*FROM_SI["P"],
        "Q": p[:, 2]*FROM_SI["Q"],
        "dr": p[:, 5]*FROM_SI["dr"],
        "dP": v[:, 3]*FROM_SI["dP"],
    }
    return p[:, 0], y


def cycle_start(t_end, T, phase):
    """Start of the last whole cycle k T + phase that ends by t_end."""
    k = math.floor((t_end - phase - T)/T + 1e-9)
    if k < 0:
        raise SystemExit("The run is shorter than one cycle")
    return k*T + phase


# ------------------------------------------------------------------------ metrics

def errors(yi, r, q):
    """Eqs. 10-13 [%] of yi against r on the same points, for quantity q."""
    if q == "P":
        n, n_dias = r, r.min()
    else:
        n = n_dias = r.max()
    return {
        ("RMS", q): 100*math.sqrt(np.mean(((yi - r)/n)**2)),
        ("MAX", q): 100*np.max(np.abs(yi - r)/np.abs(n)),
        ("SYS", q): 100*abs(yi.max() - r.max())/abs(r.max()),
        ("DIAS", q): 100*abs(yi.min() - r.min())/abs(n_dias),
    }


def metrics(t, y, ref, t0, quantities=QUANTITIES):
    """Eqs. 10-13 [%] of the solution (t, y: dict of quantities) against the
    reference array (columns t and the quantities in order), with the
    reference time tau at t = t0 + tau."""
    tr = ref[:, 0]
    out = {}
    for j, q in enumerate(quantities):
        out.update(errors(np.interp(tr + t0, t, y[q]), ref[:, j + 1], q))
    return out


def periodicity(t, P, T, t0, n=1000):
    """Relative RMS change of P from the cycle before [t0 - T, t0] to [t0, t0 + T]."""
    tau = np.linspace(0, T, n, endpoint=False)
    p1 = np.interp(t0 + tau, t, P)
    p0 = np.interp(t0 - T + tau, t, P)
    return math.sqrt(np.mean((p1 - p0)**2))/math.sqrt(np.mean(p1**2))


def verdict(v, pub):
    lo, hi = min(pub), max(pub)
    if v < lo - 0.005:
        return "below range"
    if v <= hi + 0.005:
        return "in range"
    return f"above (+{v - hi:.2f})"


def table_rows(table, s4f, quantities=QUANTITIES):
    """Rows (metric, quantity, six schemes, solids4foam, verdict) of a table
    {(metric, quantity): six values}."""
    rows = []
    for m in METRICS:
        for q in quantities:
            pub = table[(m, q)]
            v = s4f[(m, q)]
            rows.append((m, q, pub, v, verdict(v, pub)))
    return rows


def markdown_table(table, s4f, quantities=QUANTITIES):
    head = "| Error | " + " | ".join(SCHEMES) + " | solids4foam | vs range |"
    sep = "|---" * (len(SCHEMES) + 3) + "|"
    lines = [head, sep]
    for m, q, pub, v, verdict in table_rows(table, s4f, quantities):
        lines.append(
            f"| E^{m}_{q} | " + " | ".join(f"{x:.2f}" for x in pub)
            + f" | **{v:.2f}** | {verdict} |"
        )
    return "\n".join(lines)


# ------------------------------------------------------------------------ figures

def style(ax):
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK_2)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=INK_2, labelsize=8)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)


def waveform_figure(case, case_dir, t, y, t0, ref3d, schemes, path):
    """Boileau Figs. 3-4: P, dP, Q and dr at the midpoint over the cycle, the six
    schemes, the 3D reference and solids4foam. schemes: name -> reference array
    in the 1D time of that scheme."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    c = CASES[case]
    tau = ref3d[:, 0]
    labels = {"P": "P [kPa]", "dP": r"$\Delta P$ [kPa]", "Q": "Q [ml/s]",
              "dr": r"$\Delta r$ [mm]"}
    order = ["P", "dP", "Q", "dr"]
    col = {q: QUANTITIES.index(q) + 1 for q in QUANTITIES}

    fig, axes = plt.subplots(2, 2, figsize=(9, 6.4), sharex=True)
    for ax, q in zip(axes.flat, order):
        # solids4foam underneath, so that the thin lines of the schemes stay
        # visible where they coincide with it
        ax.plot(tau, np.interp(tau + t0, t, y[q]), color=COLOURS["solids4foam"],
                linewidth=3.0, alpha=0.85, label="solids4foam")
        for s, data in schemes.items():
            ax.plot(tau, np.interp(tau + t0, data[:, 0], data[:, col[q]]),
                    color=COLOURS[s], linewidth=0.9, label=s)
        ax.plot(tau, ref3d[:, col[q]], color=COLOURS["3D"], linewidth=1.4,
                linestyle=(0, (4, 2)), label="3D")
        ax.set_ylabel(labels[q])
        style(ax)
    for ax in axes[1]:
        ax.set_xlabel("t [s] (3D cycle time)")
    handles, names = axes[0, 0].get_legend_handles_labels()
    fig.suptitle(f"{c['title']}: midpoint (Boileau et al. 2015, Fig. {c['figure']})",
                 color=INK, fontsize=11, y=0.99)
    fig.legend(handles, names, loc="upper center", bbox_to_anchor=(0.5, 0.955),
               ncol=8, frameon=False, fontsize=8, labelcolor=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------- the cases' entry point

def analyse_case(case, case_dir):
    """Metrics table and figure of a physiological case. Writes plots/ and
    metrics.md in the case and returns the metrics."""
    c = CASES[case]
    case_dir = Path(case_dir)
    refdir = case_dir / "referenceData"

    t, y = load_solution(case_dir, c["probe"], c["vessel"])
    t0 = cycle_start(t[-1], c["T"], c["phase"])
    ref3d = load_reference(refdir / "3D.txt")
    m = metrics(t, y, ref3d, t0)
    per = periodicity(t, y["P"], c["T"], t0)

    schemes = {s: load_reference(refdir / f"{s}.txt") for s in SCHEMES}
    # The schemes' own cycle: the window of Boileau's tables
    k_pub = {"commonCarotid": 8, "upperThoracicAorta": 19}[case]
    t0_pub = k_pub*c["T"] + c["phase"]

    plots = case_dir / "plots"
    plots.mkdir(exist_ok=True)

    # The six schemes are plotted in their own window, solids4foam in its last cycle
    shifted = {}
    for s, data in schemes.items():
        d = data.copy()
        d[:, 0] += t0 - t0_pub
        shifted[s] = d
    waveform_figure(case, case_dir, t, y, t0, ref3d, shifted,
                    plots / f"{case}_fig{c['figure']}.png")

    text = (
        f"# {c['title']}: errors against the 3D solution [%]\n\n"
        f"Boileau et al. (2015) Table {c['table']} with a solids4foam column. "
        f"solids4foam cycle: t = {t0:.4f} to {t0 + c['T']:.4f} s; cycle-to-cycle "
        f"relative RMS change of P at the probe {per:.1e} "
        f"({'periodic' if per < 1e-4 else 'NOT periodic'} at 1e-4).\n\n"
        + markdown_table(TABLES[case], m) + "\n"
    )
    (case_dir / "metrics.md").write_text(text)
    print(text)
    return m, per


# --------------------------------------------------------------- network cases

def load_probe(case_dir, probe):
    """t and P [kPa], Q [ml/s], dr [mm] of a probe."""
    p = np.loadtxt(Path(case_dir) / "postProcessing" / "circulation1D" / "0"
                   / "probes" / f"{probe}.dat", comments="#")
    return p[:, 0], {"P": p[:, 1]*FROM_SI["P"], "Q": p[:, 2]*FROM_SI["Q"],
                     "dr": p[:, 5]*FROM_SI["dr"]}


def site_figure(path, title, rows, columns, xlabel, ref_label):
    """A grid of waveforms: rows of (site title, {column: (tau, ref, s4f,
    {scheme: values})}), columns of (key, axis label)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(len(rows), len(columns),
                             figsize=(3.3*len(columns) + 1, 2.5*len(rows) + 1.2),
                             sharex=True, squeeze=False)
    for i, (site_title, data) in enumerate(rows):
        for j, (key, label) in enumerate(columns):
            ax = axes[i, j]
            tau, ref, s4f, schemes = data[key]
            ax.plot(tau, s4f, color=COLOURS["solids4foam"], linewidth=3.0,
                    alpha=0.85, label="solids4foam")
            for s, v in schemes.items():
                ax.plot(tau, v, color=COLOURS[s], linewidth=0.9, label=s)
            if ref is not None:
                ax.plot(tau, ref, color=COLOURS["3D"], linewidth=1.4,
                        linestyle=(0, (4, 2)), label=ref_label)
            ax.set_ylabel(label)
            if j == 0:
                ax.set_title(site_title, loc="left", fontsize=9, color=INK)
            style(ax)
    for ax in axes[-1]:
        ax.set_xlabel(xlabel)
    handles, names = axes[0, 0].get_legend_handles_labels()
    fig.suptitle(title, color=INK, fontsize=11, y=0.99)
    fig.legend(handles, names, loc="upper center", bbox_to_anchor=(0.5, 0.96),
               ncol=8, frameon=False, fontsize=8, labelcolor=INK)
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.75/(2.5*len(rows) + 1.2)))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def analyse_bifurcation(case_dir):
    """Table VIII and Fig. 5 of the aortic bifurcation. Writes metrics.md and
    plots/."""
    case_dir = Path(case_dir)
    refdir = case_dir / "referenceData"
    T, phase = BIFURCATION_T, BIFURCATION_PHASE
    t0_pub = 26*T + phase
    quantities = ["P", "Q", "dr"]
    labels = {"P": "P [kPa]", "Q": "Q [ml/s]", "dr": r"$\Delta r$ [mm]"}

    text = ["# Aortic bifurcation: errors against the 3D solution [%]\n",
            "Boileau et al. (2015) Table VIII with a solids4foam column."]
    rows = []
    for probe, suffix, title in BIFURCATION_SITES:
        t, y = load_probe(case_dir, probe)
        t0 = cycle_start(t[-1], T, phase)
        ref = load_reference(refdir / f"3D_{suffix}.txt")
        m = metrics(t, y, ref, t0, quantities)
        if probe == "aortaMid":
            per = periodicity(t, y["P"], T, t0)
            text[1] += (f" solids4foam cycle: t = {t0:.4f} to {t0 + T:.4f} s;"
                        f" cycle-to-cycle relative RMS change of P at the aorta"
                        f" midpoint {per:.1e}.\n")
        text.append(f"\n## {title[0].upper() + title[1:]}\n")
        text.append(markdown_table(TABLE_VIII[probe], m, quantities))

        tau = ref[:, 0]
        data = {}
        for j, q in enumerate(quantities):
            schemes = {}
            for s in SCHEMES:
                d = load_reference(refdir / f"{s}_{suffix}.txt")
                schemes[s] = np.interp(tau + t0_pub, d[:, 0], d[:, j + 1])
            data[q] = (tau, ref[:, j + 1], np.interp(tau + t0, t, y[q]), schemes)
        rows.append((title, data))

    (case_dir / "plots").mkdir(exist_ok=True)
    site_figure(case_dir / "plots" / "aorticBifurcation_fig5.png",
                "Aortic bifurcation (Boileau et al. 2015, Fig. 5)", rows,
                [(q, labels[q]) for q in quantities],
                "t [s] (3D cycle time)", "3D")
    out = "\n".join(text) + "\n"
    (case_dir / "metrics.md").write_text(out)
    print(out)


def time_shifts(path):
    """Site -> tini of timeShifts.txt."""
    out = {}
    for line in Path(path).read_text().splitlines():
        f = line.split()
        if f and not line.startswith("#"):
            out[f[0]] = float(f[1])
    return out


def analyse_network(case_dir):
    """Table X and Figs. 6-9 of the 37-artery network against the in vitro
    data. Writes metrics.md and plots/."""
    case_dir = Path(case_dir)
    refdir = case_dir / "referenceData"
    T = NETWORK_T
    shifts = time_shifts(refdir / "timeShifts.txt")
    pub = {s: {q: load_reference(refdir / f"{s}_{q}.txt") for q in "PQ"}
           for s in SCHEMES}

    head = ("| Site | Scheme | " + " | ".join(f"E^{m}_{q}" for m, q in TABLE_X_ORDER)
            + " |")
    lines = [
        "# 37-artery network: errors against the in vitro data [%]\n",
        "Boileau et al. (2015) Table X with solids4foam rows (the range of the six "
        "schemes in brackets where solids4foam is outside it).\n",
        head, "|---" * (len(TABLE_X_ORDER) + 2) + "|",
    ]
    data_rows = {}
    counts = {"in range": 0, "below range": 0, "above": 0}
    for site in NETWORK_SITES:
        t, y = load_probe(case_dir, site)
        exp = {q: load_reference(refdir / f"Exp_{site}_{q}.txt") for q in "PQ"}
        span = exp["P"][-1, 0]
        t0_pub = shifts[site] - 0.001
        k = math.floor((t[-1] - span - t0_pub)/T + 1e-9)
        t0 = t0_pub + k*T
        m = {}
        for q in "PQ":
            m.update(errors(np.interp(exp[q][:, 0] + t0, t, y[q]), exp[q][:, 1], q))

        cells = []
        for i, key in enumerate(TABLE_X_ORDER):
            vals = [TABLE_X[site][s][i] for s in SCHEMES]
            v = verdict(m[key], vals)
            counts["above" if v.startswith("above") else v] += 1
            cells.append(f"**{m[key]:.2f}**" if v == "in range"
                         else f"**{m[key]:.2f}** ({min(vals):.2f}-{max(vals):.2f})")
        for s in SCHEMES:
            lines.append(f"| {site} | {s} | "
                         + " | ".join(f"{x:.2f}" for x in TABLE_X[site][s]) + " |")
        lines.append(f"| {site} | **solids4foam** | " + " | ".join(cells) + " |")

        data = {}
        for j, q in enumerate("PQ"):
            col = NETWORK_SITES.index(site) + 1
            tau = exp[q][:, 0]
            schemes = {s: np.interp(tau + t0_pub, pub[s][q][:, 0], pub[s][q][:, col])
                       for s in SCHEMES}
            data[q] = (tau, exp[q][:, 1], np.interp(tau + t0, t, y[q]), schemes)
        data_rows[site] = data

    lines.append(f"\nsolids4foam: {counts['in range']} in the range of the six "
                 f"schemes, {counts['below range']} below it, {counts['above']} "
                 f"above it (of {sum(counts.values())}).")

    (case_dir / "plots").mkdir(exist_ok=True)
    names = {
        "AorticArchII": "aortic arch II", "ThoracicAortaII": "thoracic aorta II",
        "LSubclavianI": "left subclavian I", "RIliacFemoralII": "right iliac-femoral II",
        "LUlnar": "left ulnar", "RAnteriorTibial": "right anterior tibial",
        "RUlnar": "right ulnar", "Splenic": "splenic",
    }
    for fig, sites in (("6", NETWORK_SITES[0:2]), ("7", NETWORK_SITES[2:4]),
                       ("8", NETWORK_SITES[4:6]), ("9", NETWORK_SITES[6:8])):
        site_figure(case_dir / "plots" / f"thirtySevenArtery_fig{fig}.png",
                    f"37-artery network (Boileau et al. 2015, Fig. {fig})",
                    [(names[s], data_rows[s]) for s in sites],
                    [("P", "P [kPa]"), ("Q", "Q [ml/s]")],
                    "t [s] (in vitro time)", "in vitro")
    out = "\n".join(lines) + "\n"
    (case_dir / "metrics.md").write_text(out)
    print(out)


# ------------------------------------------------------------------------ ADAN56

ADAN56_SITES = [
    ("AorticArchI", "aortic arch I"), ("ThoAortaIII", "thoracic aorta III"),
    ("AbdomAortaV", "abdominal aorta V"), ("RightCCA", "right common carotid"),
    ("RightRenal", "right renal"), ("RightCommonIliac", "right common iliac"),
    ("RightICA", "right internal carotid"), ("RightRadial", "right radial"),
    ("RightInternalIliac", "right internal iliac"),
    ("RightPostInteross", "right posterior interosseous"),
    ("RightFemoralII", "right femoral II"), ("RightAntTibial", "right anterior tibial"),
]
# The five schemes the paper finds consistent with one another
ADAN56_CONSISTENT = ["DCG", "LCG", "FEM", "FVM", "McC"]


def analyse_adan56(case_dir):
    """Figs. 11-16 of the ADAN56 network, and the differences between
    solids4foam and the five consistent schemes, against the differences
    among those schemes. Writes metrics.md and plots/."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    case_dir = Path(case_dir)
    refdir = case_dir / "referenceData"
    T, t0 = 1.0, 9.0
    tau = np.arange(0, T, 0.001)
    pub = {s: {q: load_reference(refdir / f"{s}_{q}.txt") for q in "PQ"}
           for s in SCHEMES}

    def rms(y, r, q):
        return errors(y, r, q)[("RMS", q)]

    lines = [
        "# ADAN56 network: solids4foam against the five consistent schemes\n",
        "Relative RMS differences [%] over the last cycle (t = 9 to 10 s), as in "
        "Boileau Eqs. 10 (P pointwise, Q by the maximum flow): solids4foam "
        "against each of DCG, LCG, FEM, FVM and McC (largest), and the largest "
        "difference between two of those five schemes.\n",
        "| Site | P: solids4foam vs the five | P: among the five | Q: solids4foam "
        "vs the five | Q: among the five |",
        "|---|---|---|---|---|",
    ]
    rows_fig = []
    for j, (site, title) in enumerate(ADAN56_SITES):
        t, y = load_probe(case_dir, site)
        data, cells = {}, []
        for q in "PQ":
            s4f = np.interp(t0 + tau, t, y[q])
            sch = {s: np.interp(t0 + tau, pub[s][q][:, 0], pub[s][q][:, j + 1])
                   for s in SCHEMES}
            vs = max(rms(s4f, sch[s], q) for s in ADAN56_CONSISTENT)
            among = max(rms(sch[a], sch[b], q) for a in ADAN56_CONSISTENT
                        for b in ADAN56_CONSISTENT if a != b)
            cells += [f"{vs:.2f}", f"{among:.2f}"]
            data[q] = (tau, None, s4f, sch)
        lines.append(f"| {title} | " + " | ".join(cells) + " |")
        rows_fig.append((title, data))

    plots = case_dir / "plots"
    plots.mkdir(exist_ok=True)
    for k, fig in enumerate(("11", "12", "13", "14")):
        site_figure(plots / f"ADAN56_fig{fig}.png",
                    f"ADAN56 network (Boileau et al. 2015, Fig. {fig})",
                    rows_fig[3*k:3*k + 3], [("P", "P [kPa]"), ("Q", "Q [ml/s]")],
                    "t [s] (last cycle)", None)

    # ---- along the aorta: Figs. 15 and 16
    probes = sorted((case_dir / "postProcessing" / "circulation1D" / "0" / "probes")
                    .glob("aorta[0-9][0-9][0-9].dat"))
    x = np.array([int(p.stem[5:])*1e-3 for p in probes])
    series = [np.loadtxt(p, comments="#") for p in probes]
    tp = series[0][:, 0]
    cyc = (tp >= t0 - 1e-9) & (tp <= t0 + T + 1e-9)
    P = np.array([s[:, 1]*1e-3 for s in series])          # [kPa], x by time
    aorta = {s: load_reference(refdir / f"{s}_aorta.txt") for s in SCHEMES}

    times = [0.1, 0.2, 0.3, 0.5, 0.7, 0.0]
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.2), sharex=True)
    worst15 = 0.0
    for c, (ax, tm) in enumerate(zip(axes.flat, times)):
        i = np.argmin(np.abs(tp - (t0 + (tm if tm > 0 else T))))
        ax.plot(x, P[:, i], color=COLOURS["solids4foam"], linewidth=3.0, alpha=0.85,
                label="solids4foam")
        for s in SCHEMES:
            ax.plot(aorta[s][:, 0], aorta[s][:, c + 1], color=COLOURS[s],
                    linewidth=0.9, label=s)
        ref = np.interp(x, aorta["FVM"][:, 0], aorta["FVM"][:, c + 1])
        worst15 = max(worst15, np.max(np.abs(P[:, i] - ref)/ref))
        ax.set_title(f"t = {tm:g} T", loc="left", fontsize=9, color=INK)
        ax.set_ylabel("P [kPa]")
        style(ax)
    for ax in axes[1]:
        ax.set_xlabel("distance from the aortic root [m]")
    handles, names = axes[0, 0].get_legend_handles_labels()
    fig.suptitle("ADAN56: pressure along the aorta (Boileau et al. 2015, Fig. 15)",
                 color=INK, fontsize=11, y=0.99)
    fig.legend(handles, names, loc="upper center", bbox_to_anchor=(0.5, 0.955),
               ncol=7, frameon=False, fontsize=8, labelcolor=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(plots / "ADAN56_fig15.png", dpi=150)
    plt.close(fig)

    stats = {"systolic": (P[:, cyc].max(axis=1), 7),
             "mean": (P[:, cyc].mean(axis=1), 9),
             "diastolic": (P[:, cyc].min(axis=1), 8)}
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    worst16 = {}
    for name, (v, col) in stats.items():
        ax.plot(x, v, color=COLOURS["solids4foam"], linewidth=3.0, alpha=0.85,
                label="solids4foam" if name == "systolic" else None)
        for s in SCHEMES:
            ax.plot(aorta[s][:, 0], aorta[s][:, col], color=COLOURS[s],
                    linewidth=0.9, label=s if name == "systolic" else None)
        ref = np.interp(x, aorta["FVM"][:, 0], aorta["FVM"][:, col])
        worst16[name] = np.max(np.abs(v - ref)/ref)
        ax.annotate(name, (x[-1], v[-1]), xytext=(4, 0), textcoords="offset points",
                    fontsize=8, color=INK_2, va="center")
    ax.set_xlabel("distance from the aortic root [m]")
    ax.set_ylabel("P [kPa]")
    ax.set_title("ADAN56: systolic, mean and diastolic pressure along the aorta "
                 "(Boileau et al. 2015, Fig. 16)", fontsize=10, color=INK)
    style(ax)
    ax.legend(frameon=False, fontsize=8, ncol=4)
    fig.tight_layout()
    fig.savefig(plots / "ADAN56_fig16.png", dpi=150)
    plt.close(fig)

    # The same for the other consistent schemes, on the FVM points
    others15, others16 = 0.0, {k: 0.0 for k in stats}
    xf = aorta["FVM"][:, 0]
    for s in ADAN56_CONSISTENT:
        if s == "FVM":
            continue
        for c in range(6):
            v = np.interp(xf, aorta[s][:, 0], aorta[s][:, c + 1])
            others15 = max(others15, np.max(np.abs(v/aorta["FVM"][:, c + 1] - 1)))
        for name, (_, col) in stats.items():
            v = np.interp(xf, aorta[s][:, 0], aorta[s][:, col])
            others16[name] = max(others16[name],
                                 np.max(np.abs(v/aorta["FVM"][:, col] - 1)))

    lines.append(
        f"\nAlong the aorta, the largest difference from FVM, relative: "
        f"{100*worst15:.2f}% over the six instants of Fig. 15; "
        + ", ".join(f"{k} {100*v:.2f}%" for k, v in worst16.items())
        + " (Fig. 16). For DCG, LCG, FEM and McC against FVM: "
        f"{100*others15:.2f}%; "
        + ", ".join(f"{k} {100*v:.2f}%" for k, v in others16.items()) + "."
    )
    out = "\n".join(lines) + "\n"
    (case_dir / "metrics.md").write_text(out)
    print(out)


# ----------------------------------------------------------------------- selftest

def selftest():
    """Recompute Tables IV, VI, VIII and X from the schemes' data in the case
    folders."""
    tutorials = HERE.parent
    selftest_single(tutorials)
    selftest_networks(tutorials)


def selftest_networks(tutorials):
    print("\nAortic bifurcation (Table VIII): recomputed - published, max |d| [%]")
    refdir = tutorials / "aorticBifurcation" / "referenceData"
    t0 = 26*BIFURCATION_T + BIFURCATION_PHASE
    print("                " + " ".join(f"{s:>6}" for s in SCHEMES))
    for probe, suffix, _ in BIFURCATION_SITES:
        ref = load_reference(refdir / f"3D_{suffix}.txt")
        dev = []
        for s in SCHEMES:
            d = load_reference(refdir / f"{s}_{suffix}.txt")
            y = {q: d[:, j + 1] for j, q in enumerate(["P", "Q", "dr"])}
            m = metrics(d[:, 0], y, ref, t0, ["P", "Q", "dr"])
            dev.append(max(abs(m[k] - v[SCHEMES.index(s)])
                           for k, v in TABLE_VIII[probe].items()))
        print(f"{probe:<16}" + " ".join(f"{d:6.2f}" for d in dev))

    print("\n37-artery network (Table X): recomputed - published, max |d| [%]")
    refdir = tutorials / "thirtySevenArtery" / "referenceData"
    shifts = time_shifts(refdir / "timeShifts.txt")
    print("                " + " ".join(f"{s:>6}" for s in SCHEMES))
    for site in NETWORK_SITES:
        col = NETWORK_SITES.index(site) + 1
        dev = []
        for s in SCHEMES:
            m = {}
            for q in "PQ":
                e = load_reference(refdir / f"Exp_{site}_{q}.txt")
                d = load_reference(refdir / f"{s}_{q}.txt")
                m.update(errors(np.interp(e[:, 0] + shifts[site] - 0.001,
                                          d[:, 0], d[:, col]), e[:, 1], q))
            dev.append(max(abs(m[k] - TABLE_X[site][s][i])
                           for i, k in enumerate(TABLE_X_ORDER)))
        print(f"{site:<16}" + " ".join(f"{d:6.2f}" for d in dev))


def selftest_single(tutorials):
    worst = 0.0
    for case, c in CASES.items():
        refdir = tutorials / case / "referenceData"
        ref3d = load_reference(refdir / "3D.txt")
        k = {"commonCarotid": 8, "upperThoracicAorta": 19}[case]
        t0 = k*c["T"] + c["phase"]
        print(f"\n{c['title']} (Table {c['table']}), window t0 = {t0:.4f} s:"
              " recomputed - published [%]")
        print("        " + " ".join(f"{s:>6}" for s in SCHEMES))
        diffs = {s: [] for s in SCHEMES}
        for s in SCHEMES:
            d = load_reference(refdir / f"{s}.txt")
            y = {q: d[:, j + 1] for j, q in enumerate(QUANTITIES)}
            m = metrics(d[:, 0], y, ref3d, t0)
            for key in TABLES[case]:
                diffs[s].append(m[key] - TABLES[case][key][SCHEMES.index(s)])
        print("max|d| " + " ".join(f"{max(map(abs, diffs[s])):6.2f}" for s in SCHEMES))
        for s in ("LCG", "FEM", "FVM", "McC"):
            worst = max(worst, max(map(abs, diffs[s])))
    print(f"\nWorst LCG/FEM/FVM/McC deviation from the published tables: "
          f"{worst:.3f} % (published to two decimals)")
    print("DCG (both cases) and STM (aorta) do not reproduce: their files use"
          " another time origin (DCG starts at t = 0.001, the aorta STM run is"
          " 20.055 s long), which the window cannot recover.")


if __name__ == "__main__":
    if sys.argv[1:] == ["selftest"]:
        selftest()
    else:
        print(__doc__)
