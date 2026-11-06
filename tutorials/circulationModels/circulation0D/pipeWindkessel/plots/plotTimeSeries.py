#!/usr/bin/env python3
"""The outlet pressure of every algorithm against the exact solution, in time.

    python3 plots/plotTimeSeries.py [dt ...]        (default: 0.1 0.01 0.001)

For each time step two figures are written, from the cases in
../explicit and ../implicit:

    outletPressure_dt<dt>.png   Po(t) of the eight algorithms (coupling scheme
                                and order 1 to 4) on one axis, with the exact
                                solution as the dashed black line
    pointwiseError_dt<dt>.png   the pointwise error Po(t) - Po_exact(t), one
                                panel per algorithm; each panel has its own
                                vertical scale, as the errors of the first and
                                fourth order differ by many orders of magnitude

Colour and marker shape give the order of the coupling; filled markers and a
solid line are explicit, hollow markers and a dashed line implicit. A marker is
drawn every few steps only, so that the fine time steps stay readable; the
lines pass through every step. The figures are drawn as SVG and converted with
rsvg-convert; without it the SVG files are kept instead.
"""

import math
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT))
import figure_style as fs                                    # noqa: E402
import analyse                                               # noqa: E402
from plotOrderOfAccuracy import ORDER_STYLE, ORDER_WORD, marker, sup  # noqa: E402

ORDERS = ("first", "second", "third", "fourth")
GRID, GUIDE, EXACT = "#e3e3e3", "#8a8a8a", "#000000"
MARKERS = 20          # markers per curve at most


def read_series(scheme, order, dt):
    """(t, Po) of a case, from t = 0 (Po = 0) to the end of the run."""
    case = ROOT/scheme/f"{order}Order"/f"dt_{dt}"
    files = sorted((case/"postProcessing"/"outletPressure").glob("*/surfaceFieldValue.dat"))
    if not files:
        sys.exit(f"No results in {case}: run the case first")
    points = [(0.0, 0.0)]
    for line in files[0].read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            t, p = (float(x) for x in line.split()[:2])
            if t <= 1 + 1e-9:
                points.append((t, p*analyse.RHO))
    return points


def nice_ticks(lo, hi, count=5):
    raw = (hi - lo)/count
    magnitude = 10**math.floor(math.log10(raw))
    step = min((m*magnitude for m in (1, 2, 2.5, 5, 10)), key=lambda c: abs(c - raw))
    value = math.ceil(lo/step)*step
    ticks = []
    while value <= hi + step*1e-9:
        v = round(value, 12)
        ticks.append(0.0 if v == 0 else v)
        value += step
    return ticks


def style(scheme, order):
    colour, shape = ORDER_STYLE[ORDERS.index(order) + 1]
    return colour, shape, scheme == "implicit"


def curve(points, X, Y, colour, shape, implicit, r, stride):
    dash = ' stroke-dasharray="9 6"' if implicit else ""
    coords = " ".join(f"{X(t):.2f},{Y(v):.2f}" for t, v in points)
    out = [f'<polyline points="{coords}" fill="none" stroke="{colour}" '
           f'stroke-width="2.4" stroke-linejoin="round"{dash}/>']
    # The shared initial condition is left to the exact or zero line.
    for t, v in points[1::stride]:
        out.append(marker(shape, X(t), Y(v), colour, implicit, r))
    return out


def write(svg_text, name):
    svg = HERE/f"{name}.svg"
    svg.write_text(svg_text)
    conv = shutil.which("rsvg-convert")
    if conv:
        subprocess.run([conv, "-f", "png", "-z", "1.6",
                        "-o", str(svg.with_suffix(".png")), str(svg)], check=True)
        svg.unlink()
        print(f"wrote {name}.png in {HERE}")
    else:
        print(f"rsvg-convert not found - kept {svg.name} in {HERE}")


def pressure_figure(series, dt, n):
    """All eight Po(t) and the exact solution on one axis."""
    W, H = 1240, 1060
    L, R, T, PB = 175, 45, 45, 745
    s = fs.sizes(W)
    exact = [(i/400, analyse.exact_po(i/400)) for i in range(401)]
    ys = [v for pts in series.values() for _, v in pts] + [v for _, v in exact]
    ymin, ymax = min(ys), max(ys)
    pad = 0.10*(ymax - ymin)
    ymin, ymax = ymin - pad, ymax + pad

    def X(t): return L + t*(W - L - R)
    def Y(v): return T + (ymax - v)/(ymax - ymin)*(PB - T)

    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
         f'viewBox="0 0 {W} {H}">', f'<style>{fs.css(W)}</style>',
         f'<rect width="{W}" height="{H}" fill="#ffffff"/>']
    for v in nice_ticks(ymin, ymax, 6):
        y = Y(v)
        o.append(f'<line x1="{L}" y1="{y:.2f}" x2="{W-R}" y2="{y:.2f}" '
                 f'stroke="{GRID}" stroke-width="1"/>')
        o.append(f'<text class="tick" x="{L-16}" y="{y+s["tick"]*0.36:.2f}" '
                 f'text-anchor="end">{v:g}</text>')
    for i in range(6):
        x = X(i/5)
        o.append(f'<line x1="{x:.2f}" y1="{T}" x2="{x:.2f}" y2="{PB}" '
                 f'stroke="{GRID}" stroke-width="1"/>')
        o.append(f'<text class="tick" x="{x:.2f}" y="{PB+s["tick"]*1.5:.0f}" '
                 f'text-anchor="middle">{i/5:g}</text>')
    o.append(f'<polyline points="' + " ".join(f"{X(t):.2f},{Y(v):.2f}" for t, v in exact)
             + f'" fill="none" stroke="{EXACT}" stroke-width="2.6" stroke-dasharray="9 6"/>')
    stride = max(1, n//MARKERS)
    for (scheme, order), pts in series.items():
        colour, shape, implicit = style(scheme, order)
        o += curve(pts, X, Y, colour, shape, implicit, 7.0, stride)
    o.append(f'<rect x="{L}" y="{T}" width="{W-L-R}" height="{PB-T}" fill="none" '
             f'stroke="#000000" stroke-width="{s["axis_line"]}"/>')
    o.append(f'<text class="axlab" x="{(L+W-R)/2:.0f}" y="820" '
             f'text-anchor="middle">Time [s]</text>')
    o.append(f'<text class="axlab" transform="translate(42,{(T+PB)/2:.0f}) '
             f'rotate(-90)" text-anchor="middle">Outlet Pressure [Pa]</text>')
    o.append('<line x1="215" y1="866" x2="257" y2="866" stroke="#000000" '
             'stroke-width="2.6" stroke-dasharray="9 6"/>')
    o.append(f'<text class="leg" x="273" y="{866+s["legend"]*0.36:.0f}">Exact</text>')
    o.append(f'<text class="leg" x="{W-R}" y="{866+s["legend"]*0.36:.0f}" '
             f'text-anchor="end">&#916;t = {dt}</text>')
    for k, order in enumerate(ORDERS):
        colour, shape = ORDER_STYLE[k + 1]
        y = 902 + k*34
        for scheme, x in (("explicit", 215), ("implicit", 720)):
            implicit = scheme == "implicit"
            dash = ' stroke-dasharray="9 6"' if implicit else ""
            o.append(f'<line x1="{x}" y1="{y}" x2="{x+42}" y2="{y}" '
                     f'stroke="{colour}" stroke-width="2.4"{dash}/>')
            o.append(marker(shape, x + 21, y, colour, implicit, 7.0))
            o.append(f'<text class="leg" x="{x+58}" y="{y+s["legend"]*0.36:.0f}">'
                     f'{scheme.capitalize()}, {order} order</text>')
    o.append("</svg>")
    return "\n".join(o)


def error_figure(series, dt, n):
    """The pointwise error Po - Po_exact, one panel per algorithm."""
    W, H = 1400, 1560
    OUTER_L, OUTER_R, TOP, BOTTOM = 110, 30, 30, 90
    AX_L, AX_R, TITLE, GAP_Y = 150, 30, 55, 85
    s = fs.sizes(W)
    cw = (W - OUTER_L - OUTER_R)/2
    ph = (H - TOP - BOTTOM - 3*GAP_Y)/4
    stride = max(1, n//MARKERS)
    small = round(s["legend"]*0.9)

    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
         f'viewBox="0 0 {W} {H}">', f'<style>{fs.css(W)}</style>',
         f'<rect width="{W}" height="{H}" fill="#ffffff"/>']
    for r, order in enumerate(ORDERS):
        for c, scheme in enumerate(("explicit", "implicit")):
            pts = [(t, v - analyse.exact_po(t)) for t, v in series[scheme, order]]
            colour, shape, implicit = style(scheme, order)
            emax = max(abs(v) for _, v in pts)
            lo, hi = min(min(v for _, v in pts), 0), max(max(v for _, v in pts), 0)
            pad = 0.10*(hi - lo)
            lo, hi = lo - pad, hi + pad
            px = OUTER_L + c*cw + AX_L
            pw = cw - AX_L - AX_R
            py = TOP + r*(ph + GAP_Y) + TITLE
            pah = ph - TITLE

            def X(t): return px + t*pw
            def Y(v): return py + (hi - v)/(hi - lo)*pah

            for v in nice_ticks(lo, hi, 4):
                y = Y(v)
                o.append(f'<line x1="{px:.2f}" y1="{y:.2f}" x2="{px+pw:.2f}" '
                         f'y2="{y:.2f}" stroke="{GRID}" stroke-width="1"/>')
                o.append(f'<text class="tick" x="{px-12:.2f}" '
                         f'y="{y+s["tick"]*0.36:.2f}" text-anchor="end">{v:.2g}</text>')
            for t in (0, 0.5, 1):
                o.append(f'<line x1="{X(t):.2f}" y1="{py:.2f}" x2="{X(t):.2f}" '
                         f'y2="{py+pah:.2f}" stroke="{GRID}" stroke-width="1"/>')
                o.append(f'<text class="tick" x="{X(t):.2f}" '
                         f'y="{py+pah+s["tick"]*1.3:.0f}" text-anchor="middle">{t:g}</text>')
            o.append(f'<line x1="{px:.2f}" y1="{Y(0):.2f}" x2="{px+pw:.2f}" '
                     f'y2="{Y(0):.2f}" stroke="{GUIDE}" stroke-width="1.6"/>')
            o += curve(pts, X, Y, colour, shape, implicit, 5.5, stride)
            o.append(f'<rect x="{px:.2f}" y="{py:.2f}" width="{pw:.2f}" '
                     f'height="{pah:.2f}" fill="none" stroke="#000000" '
                     f'stroke-width="{s["axis_line"]}"/>')
            o.append(f'<text class="leg" x="{px:.2f}" y="{py-16:.2f}">'
                     f'{scheme.capitalize()}, {ORDER_WORD[ORDERS.index(order)+1]}'
                     f'</text>')
            o.append(f'<text class="leg" x="{px+pw:.2f}" y="{py-16:.2f}" '
                     f'text-anchor="end" fill="#555555" font-size="{small}">'
                     f'max |e| = {emax:.2e} Pa</text>')
    o.append(f'<text class="axlab" x="{OUTER_L + AX_L + (W-OUTER_L-OUTER_R-AX_L-AX_R)/2:.0f}" '
             f'y="{H-24}" text-anchor="middle">Time [s]</text>')
    o.append(f'<text class="axlab" transform="translate(40,{H/2:.0f}) rotate(-90)" '
             f'text-anchor="middle">Outlet Pressure Error [Pa]  (&#916;t = {dt})</text>')
    o.append("</svg>")
    return "\n".join(o)


def main():
    dts = sys.argv[1:] or ["0.1", "0.01", "0.001"]
    for dt in dts:
        series = {(scheme, order): read_series(scheme, order, dt)
                  for scheme in ("explicit", "implicit") for order in ORDERS}
        n = len(next(iter(series.values()))) - 1
        tag = "dt" + dt.replace(".", "p")
        write(pressure_figure(series, dt, n), f"outletPressure_{tag}")
        write(error_figure(series, dt, n), f"pointwiseError_{tag}")


if __name__ == "__main__":
    main()
