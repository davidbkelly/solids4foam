#!/usr/bin/env python3
"""Convergence figure of all eight coupling algorithms of this tutorial.

    python3 plots/plotOrderOfAccuracy.py

Log-log RMSE of the outlet pressure against the exact solution versus time step,
for explicit and implicit coupling of orders one to four (the cases in
../explicit and ../implicit, analysed by analyse_order.py). Eight series exceed
what hue alone can separate, so identity is carried by two channels:

    hue and marker shape  the order of the coupling (1, 2, 3, 4)
    filled / hollow       explicit / implicit coupling (solid / dashed fit)

Two figures are written, the second without the reference slopes:

    orderOfAccuracy.png           with a reference slope of each theoretical
                                  order, dt^1 to dt^4, drawn between the
                                  explicit and implicit series of that order
    orderOfAccuracyNoGuides.png   the data and fits only

The figures are drawn as SVG and converted with rsvg-convert; without it the
SVG files are kept instead. Text styling comes from figure_style.py.
"""

import math
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import figure_style as fs                                    # noqa: E402
from analyse_order import main as analyse                    # noqa: E402

W, H = 1240, 1000
L, R, T, B = 165, 45, 45, 115

# One hue per order of the coupling, colour-vision separable.
ORDER_STYLE = {1: ("#0072B2", "circle"), 2: ("#D55E00", "square"),
               3: ("#009E73", "triangle"), 4: ("#E69F00", "diamond")}
ORDER_WORD = {1: "1st", 2: "2nd", 3: "3rd", 4: "4th"}
GRID, GUIDE = "#e3e3e3", "#8a8a8a"


def marker(shape, x, y, c, hollow, r=7.5):
    """Filled with a white ring for explicit; white with a coloured ring for
    implicit, so the coupling reads without relying on the line style."""
    if hollow:
        common = f'fill="#ffffff" stroke="{c}" stroke-width="2.6"'
    else:
        common = f'fill="{c}" stroke="#ffffff" stroke-width="2"'
    if shape == "circle":
        return f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r:.2f}" {common}/>'
    if shape == "square":
        s = r*1.75
        return (f'<rect x="{x-s/2:.2f}" y="{y-s/2:.2f}" width="{s:.2f}" '
                f'height="{s:.2f}" {common}/>')
    if shape == "triangle":
        p = (f"{x:.2f},{y-r*1.15:.2f} {x-r:.2f},{y+r*0.8:.2f} "
             f"{x+r:.2f},{y+r*0.8:.2f}")
        return f'<polygon points="{p}" {common}/>'
    p = (f"{x:.2f},{y-r*1.2:.2f} {x+r*1.2:.2f},{y:.2f} "
         f"{x:.2f},{y+r*1.2:.2f} {x-r*1.2:.2f},{y:.2f}")
    return f'<polygon points="{p}" {common}/>'


def sup(text, size):
    return f'<tspan dy="-{size*0.45:.0f}" font-size="{size*0.72:.0f}">{text}</tspan>'


def draw(series, guides):
    """The SVG text of the figure for the (key, order, implicit, data) series."""
    s = fs.sizes(W)
    dts = sorted({r["dt"] for *_, d in series for r in d["levels"]})
    errs = [r["rmse"] for *_, d in series for r in d["levels"]]

    x0, x1 = math.log10(min(dts)), math.log10(max(dts))
    # One decade of headroom below the data holds the legend.
    y0 = math.floor(math.log10(min(errs))) - 1
    y1 = math.ceil(math.log10(max(errs)))
    pad = 0.06*(x1 - x0)
    # The guide labels sit left of the finest level, so give that side room.
    x0, x1 = x0 - (2.2*pad if guides else pad), x1 + pad

    def X(dt): return L + (math.log10(dt) - x0)/(x1 - x0)*(W - L - R)
    def Y(e):  return T + (y1 - math.log10(e))/(y1 - y0)*(H - T - B)

    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
         f'viewBox="0 0 {W} {H}">', f'<style>{fs.css(W)}</style>',
         f'<rect width="{W}" height="{H}" fill="#ffffff"/>']

    for d in range(y0, y1 + 1):
        yy = Y(10.0**d)
        o.append(f'<line x1="{L}" y1="{yy:.2f}" x2="{W-R}" y2="{yy:.2f}" '
                 f'stroke="{GRID}" stroke-width="1"/>')
        if (d - y0) % 2 == 0 or y1 - y0 <= 8:
            o.append(f'<text class="tick" x="{L-16}" y="{yy+s["tick"]*0.36:.2f}" '
                     f'text-anchor="end">10{sup(d, s["tick"])}</text>')
    # Decade ticks only: the ladder has too many levels to label individually.
    for d in range(math.floor(x0), math.floor(x1) + 1):
        for m in (1, 2, 5):
            v = m*10.0**d
            if not (10**x0 <= v <= 10**x1):
                continue
            x = X(v)
            o.append(f'<line x1="{x:.2f}" y1="{T}" x2="{x:.2f}" y2="{H-B}" '
                     f'stroke="{GRID}" stroke-width="1"/>')
            o.append(f'<text class="tick" x="{x:.2f}" y="{H-B+s["tick"]*1.5:.0f}" '
                     f'text-anchor="middle">{v:g}</text>')

    if guides:
        # Scale of each guide: the geometric mean of the explicit and implicit
        # constants C in e = C dt^p at the theoretical p, so the guide runs
        # between the two series it describes. Drawn first, under the data.
        for order in sorted({o_ for _, o_, _, _ in series}):
            pair = [d for _, o_, _, d in series if o_ == order]
            lnCs = [sum(math.log(r["rmse"]) - order*math.log(r["dt"])
                        for r in d["levels"])/len(d["levels"]) for d in pair]
            scale = math.exp(sum(lnCs)/len(lnCs))
            lo = min(r["dt"] for d in pair for r in d["levels"])
            hi = max(r["dt"] for d in pair for r in d["levels"])
            ya, yb = scale*lo**order, scale*hi**order
            o.append(f'<line x1="{X(lo):.2f}" y1="{Y(ya):.2f}" '
                     f'x2="{X(hi):.2f}" y2="{Y(yb):.2f}" stroke="{GUIDE}" '
                     f'stroke-width="1.6" stroke-dasharray="3 5"/>')
            o.append(f'<text class="leg" x="{X(lo)-14:.2f}" '
                     f'y="{Y(ya)+s["legend"]*0.36:.2f}" text-anchor="end" '
                     f'fill="{GUIDE}">&#916;t{sup(order, s["legend"])}</text>')

    for key, order, implicit, d in series:
        lv = d["levels"]
        colour, shape = ORDER_STYLE[order]
        p_ = d["fittedOrder"]
        lnC = sum(math.log(r["rmse"]) - p_*math.log(r["dt"]) for r in lv)/len(lv)
        lo, hi = min(r["dt"] for r in lv), max(r["dt"] for r in lv)
        line = lambda t: math.exp(lnC)*t**p_
        dash = ' stroke-dasharray="9 6"' if implicit else ""
        o.append(f'<line x1="{X(lo):.2f}" y1="{Y(line(lo)):.2f}" '
                 f'x2="{X(hi):.2f}" y2="{Y(line(hi)):.2f}" '
                 f'stroke="{colour}" stroke-width="2.2"{dash}/>')
        for r in lv:
            o.append(marker(shape, X(r["dt"]), Y(r["rmse"]), colour, implicit))

    o.append(f'<rect x="{L}" y="{T}" width="{W-L-R}" height="{H-T-B}" '
             f'fill="none" stroke="#000000" stroke-width="{s["axis_line"]}"/>')
    o.append(f'<text class="axlab" x="{(L+W-R)/2:.0f}" y="{H-26}" '
             f'text-anchor="middle">Timestep &#916;t [s]</text>')
    o.append(f'<text class="axlab" transform="translate(40,{(T+H-B)/2:.0f}) '
             f'rotate(-90)" text-anchor="middle">'
             f'Outlet Pressure RMSE [Pa]</text>')

    # Legend in the lower right: at the coarse time steps every error is still
    # large, so the bottom-right corner is the one region no series reaches.
    lx = W - R - 330
    ly = H - B - 7*34 - 28
    for i, (key, order, implicit, d) in enumerate(series):
        colour, shape = ORDER_STYLE[order]
        yy = ly + i*34
        dash = ' stroke-dasharray="9 6"' if implicit else ""
        o.append(f'<line x1="{lx}" y1="{yy}" x2="{lx+42}" y2="{yy}" '
                 f'stroke="{colour}" stroke-width="2.4"{dash}/>')
        o.append(marker(shape, lx + 21, yy, colour, implicit))
        label = (f'{"Implicit" if implicit else "Explicit"} '
                 f'{ORDER_WORD[order]},  p = {d["fittedOrder"]:.2f}')
        o.append(f'<text class="leg" x="{lx+58}" y="{yy+s["legend"]*0.36:.0f}">'
                 f'{label}</text>')

    o.append("</svg>")
    return "\n".join(o)


def write(svg_text, name):
    svg = HERE/f"{name}.svg"
    svg.write_text(svg_text)
    conv = shutil.which("rsvg-convert")
    if conv:
        subprocess.run([conv, "-f", "png", "-z", "1.6",
                        "-o", str(svg.with_suffix(".png")), str(svg)],
                       check=True)
        svg.unlink()
        print(f"wrote {name}.png in {HERE}")
    else:
        print(f"rsvg-convert not found - kept {svg.name} in {HERE}")


def main():
    results = analyse(verbose=False)
    series = [(key, d["expectedOrder"], key.startswith("Implicit"), d)
              for key, d in results.items() if d["levels"]]
    series.sort(key=lambda s: (s[1], s[2]))
    if len(series) != 8:
        sys.exit("Expected results for all eight algorithms: run every case first")
    levels = [{round(r["dt"], 12) for r in d["levels"]} for *_, d in series]
    if any(lv != levels[0] for lv in levels):
        sys.exit("Every algorithm needs the same time steps: run every case first")

    write(draw(series, guides=True), "orderOfAccuracy")
    write(draw(series, guides=False), "orderOfAccuracyNoGuides")


if __name__ == "__main__":
    main()
