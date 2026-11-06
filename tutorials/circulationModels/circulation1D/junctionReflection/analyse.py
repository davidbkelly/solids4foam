#!/usr/bin/env python3
"""Reflection and transmission of a linear pulse at a 1 -> 2 junction.

    python3 analyse.py

For each variant that has run, the pressure peaks at the probes:

    incident     the parent probe (x = 2 m), on the way to the junction
    reflected    the parent probe, on the way back
    transmitted  each daughter probe (x = 1 m)

against linear theory for a junction of vessels with admittances
Y = A/(rho c) (from the solver's log, Y = 1/Z):

    R = (Y0 - sum Yi)/(Y0 + sum Yi),     T = 1 + R

The pulse is Boileau's Eq. 14 scaled by 1e-3, so it is linear to about 1e-6.
Writes reflection.md. Needs numpy.
"""
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
VARIANTS = ["stiffDaughters", "matchedDaughters", "compliantDaughters"]
X_PARENT, L_PARENT, X_DAUGHTER, T_PULSE = 2.0, 4.0, 1.0, 0.05


def vessels_from_log(path):
    """Vessel -> (c, Z) of the initial state, from the solver's summary."""
    out, name = {}, None
    for line in Path(path).read_text().splitlines():
        m = re.match(r"\s*Vessel (\w+):", line)
        if m:
            name = m.group(1)
        m = re.search(r"initial state: .* c ([\d.eE+-]+) m/s, Z = rho c/A ([\d.eE+-]+)", line)
        if m and name:
            out[name] = (float(m.group(1)), float(m.group(2)))
    return out


def peak(t, p, centre, half):
    """Signed extreme of p within centre +- half, refined by a parabola."""
    sel = np.where(np.abs(t - centre) < half)[0]
    i = sel[np.argmax(np.abs(p[sel]))]
    y0, y1, y2 = p[i - 1], p[i], p[i + 1]
    den = y0 - 2*y1 + y2
    return y1 - 0.25*(y0 - y2)*((0.5*(y0 - y2)/den) if den else 0.0)


def main():
    lines = ["# Junction reflection: linear pulse at a 1 -> 2 junction\n",
             "| Variant | sum Y/Y0 | R theory | R measured | T1 theory | T1 measured "
             "| T2 theory | T2 measured |",
             "|---|---|---|---|---|---|---|---|"]
    for v in VARIANTS:
        d = HERE/v
        log = d/"log.circulation1DFoam"
        if not log.is_file():
            continue
        ves = vessels_from_log(log)
        c0, Z0 = ves["parent"]
        Y0 = 1/Z0
        Ys = [1/ves[n][1] for n in ("daughter1", "daughter2")]
        R = (Y0 - sum(Ys))/(Y0 + sum(Ys))
        out = d/"postProcessing"/"circulation1D"/"0"/"probes"
        pp = np.loadtxt(out/"parent.dat", comments="#")
        inc = peak(pp[:, 0], pp[:, 1], T_PULSE + X_PARENT/c0, 0.1)
        ref = peak(pp[:, 0], pp[:, 1], T_PULSE + (2*L_PARENT - X_PARENT)/c0, 0.1)
        trans = []
        for n in ("daughter1", "daughter2"):
            pd = np.loadtxt(out/f"{n}.dat", comments="#")
            ci = ves[n][0]
            trans.append(peak(pd[:, 0], pd[:, 1],
                              T_PULSE + L_PARENT/c0 + X_DAUGHTER/ci, 0.1)/inc)
        lines.append(
            f"| {v} | {sum(Ys)/Y0:.4f} | {R:+.4f} | {ref/inc:+.4f} | {1 + R:.4f} | "
            f"{trans[0]:.4f} | {1 + R:.4f} | {trans[1]:.4f} |"
        )
    lines.append(
        "\nThe measured ratios include the numerical decay of the peak over the "
        "extra distance it travels (about 0.1% per metre at dx = 1.25 mm)."
    )
    text = "\n".join(lines) + "\n"
    (HERE/"reflection.md").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
