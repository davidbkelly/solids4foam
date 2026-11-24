#!/usr/bin/env python3
"""B5, the matched segment: coupled runs against the exact solution.

    python3 analyse.py [case ...]       errors of the given cases (default .)

The 3D pipe (rhoPimpleFoam, slip walls, inviscid) holds a compressible fluid
with the vessels' wave speed c0 and impedance Zc, so the whole system is one
uniform lossless line: the incident pulse Pinc(t) = Zc Qin(t - x/c0) passes
through the pipe with no reflection, delayed by L/c0. Measured against that,
relative to the incident peak:
    P1, Q1    the pressure and the 3D flux at the pipe's inlet (to3D)
    P2, Q2    the same at its outlet (from3D)
    refl      the largest wave in the proximal vessel's midpoint once the
              incident pulse has passed it: the spurious reflection
    trans     the distal vessel's midpoint against the delayed pulse
    delay     the delay between the pressure peaks at the two ports, against
              L/c0
The coupler's rows (coupling.dat) of every start time are joined, so a run
continued from a later time (Allrun's two phases) reads as one.
"""
import math
import re
import sys
from pathlib import Path

import numpy as np

RHO = 1050.0
A0 = math.pi*0.01**2
BETA = 4.0/3.0*math.sqrt(math.pi)*400e3*1.5e-3
C0 = math.sqrt(BETA/(2*RHO*A0))*A0**0.25         # 6.1721 m/s
Z0 = RHO*C0/A0
SCALE = 1e-3                                     # inflow amplitude factor
Q_INC = SCALE*1e-6                               # incident peak flow [m3/s]
P_INC = Z0*Q_INC                                 # incident peak pressure [Pa]
SIGMA = 1/math.sqrt(2e4)                         # Qin = exp(-1e4 (t - 0.05)^2)
L1D = 1.0                                        # each vessel [m]
LPIPE = 0.2                                      # the 3D pipe [m]
T1 = 0.05 + L1D/C0                               # incident peak at the inlet
T2 = T1 + LPIPE/C0                               # and at the outlet
T_PM = 0.05 + 0.5*L1D/C0                         # at the proximal midpoint
T_DM = T2 + 0.5*L1D/C0                           # at the distal midpoint


def pulse(t, t0):
    """The incident pulse, relative to its peak, centred at t0"""
    return np.exp(-(t - t0)**2/(2*SIGMA**2))


def load(case):
    """The coupler's rows of every start time, joined: t updates to3D.Q
    to3D.P from3D.Q from3D.P proximalMid.P proximalMid.Q distalMid.P
    distalMid.Q to3D.Qsampled from3D.Qsampled"""
    root = Path(case)/"postProcessing"/"circulation1DCoupling"
    parts = sorted((float(d.name), d) for d in root.iterdir()
                   if (d/"coupling.dat").is_file())
    rows = [np.loadtxt(d/"coupling.dat", comments="#", ndmin=2)
            for _, d in parts]
    a = np.vstack(rows)
    a = a[np.argsort(a[:, 0], kind="stable")]
    keep = np.r_[np.diff(a[:, 0]) > 0, True]     # a repeated time: the later
    a = a[keep]
    d = {"t": a[:, 0], "updates": a[:, 1], "P1": a[:, 3], "P2": a[:, 5],
         "Q1": a[:, 10], "Q2": a[:, 11], "Ppm": a[:, 6], "Pdm": a[:, 8]}
    # Explicit coupling samples the 3D flux once, at the start of a step: the
    # flux at t is the sample of the next row (as rigidSegment3D/analyse.py)
    props = Path(case)/"constant"/"circulation1DProperties"
    if re.search(r"^\s*scheme\s+explicit\s*;", props.read_text(), re.M):
        d = {k: v[:-1] for k, v in d.items()}
        d["Q1"], d["Q2"] = a[1:, 10], a[1:, 11]
    return d


def peak_time(t, y):
    """Time of the maximum, by a parabola through the three largest samples"""
    i = int(np.argmax(y))
    i = min(max(i, 1), len(y) - 2)
    y0, y1, y2 = y[i - 1], y[i], y[i + 1]
    h = t[i + 1] - t[i]
    d = y0 - 2*y1 + y2
    return t[i] + (0.5*h*(y0 - y2)/d if d != 0 else 0.0)


def errors(case):
    d = load(case)
    t = d["t"]
    e = {
        "P1": np.max(np.abs(d["P1"] - P_INC*pulse(t, T1)))/P_INC,
        "Q1": np.max(np.abs(d["Q1"] - Q_INC*pulse(t, T1)))/Q_INC,
        "P2": np.max(np.abs(d["P2"] - P_INC*pulse(t, T2)))/P_INC,
        "Q2": np.max(np.abs(d["Q2"] - Q_INC*pulse(t, T2)))/Q_INC,
        "trans": np.max(np.abs(d["Pdm"] - P_INC*pulse(t, T_DM)))/P_INC,
    }
    after = t > T_PM + 8*SIGMA
    e["refl"] = np.max(np.abs(d["Ppm"][after]))/P_INC
    delay = peak_time(t, d["P2"]) - peak_time(t, d["P1"])
    e["delay"] = (delay - LPIPE/C0)/(LPIPE/C0)
    e["updates"] = float(np.mean(d["updates"][t > T1 - 6*SIGMA]))
    return e


def main():
    cases = sys.argv[1:] or ["."]
    print("errors / incident peak; delay relative to L/c0 = "
          f"{LPIPE/C0*1e3:.3f} ms")
    for c in cases:
        e = errors(c)
        print(f"{c}: P1 {e['P1']:.2e}  Q1 {e['Q1']:.2e}  P2 {e['P2']:.2e}  "
              f"Q2 {e['Q2']:.2e}  reflection {e['refl']:.2e}  "
              f"transmitted {e['trans']:.2e}  delay {e['delay']:+.2e}  "
              f"exchanges/step {e['updates']:.1f}")


if __name__ == "__main__":
    main()
