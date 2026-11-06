#!/usr/bin/env python3
"""Order-of-accuracy figure of the six coupling algorithms of this tutorial.

    python3 plots/plotOrderOfAccuracy.py

All six algorithms (explicit and implicit coupling, order 2 to 4) on one set of
axes, in three stacked panels sharing the time step axis: (a) the outlet
pressure, (b) the inlet pressure and (c) the inlet/outlet flow rate, as analysed
by analyse_order.py. Markers are the RMSE against each algorithm's own
dt = 1e-5 run over the times common to all runs in 1.6 <= t <= 2.4 s; the line
is the least-squares fit over every time step shown, whose slope is the quoted
p_fit. Hue and marker shape give the order; a filled marker and a solid line are
explicit, a hollow marker and a dashed line are implicit. The grey dotted lines
are the dt^2, dt^3 and dt^4 reference slopes, each running between its explicit
and implicit series.

Writes orderOfAccuracy.png, drawn as SVG and converted with rsvg-convert;
without it the SVG is kept instead. Text styling comes from figure_style.py.
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

W, H = 1240, 2050
L, R = 200, 40
HEAD = 175                    # legend
PH, GAP = 560, 34             # panel height, gap between panels
FOOT = 130
GRID, GUIDE = "#e3e3e3", "#8a8a8a"
# One hue per order of the coupling; filled/solid explicit, hollow/dashed implicit
ORDER_STYLE = {2: ("#D55E00", "square"), 3: ("#009E73", "triangle"),
               4: ("#E69F00", "diamond")}
IMPLICIT_DASH = ' stroke-dasharray="9 6"'
ORDERNAME = {2: "second", 3: "third", 4: "fourth"}
PANELS = [("outletPressure", "(a) Outlet Pressure", "Outlet Pressure RMSE [Pa]"),
          ("inletPressure", "(b) Inlet Pressure", "Inlet Pressure RMSE [Pa]"),
          ("inletFlow", "(c) Inlet/Outlet Flow Rate", "Inlet/Outlet Flow Rate RMSE [m{}")]


def marker(shape, x, y, c, hollow, r=7.5):
    """Filled with a white ring for explicit; white with a coloured ring for
    implicit, so the coupling reads without relying on the line style."""
    if hollow:
        common = f'fill="#ffffff" stroke="{c}" stroke-width="2.6"'
    else:
        common = f'fill="{c}" stroke="#ffffff" stroke-width="2"'
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


def main():
    results = analyse(verbose=False)
    s = fs.sizes(W)
    pw = W - L - R

    # (scheme, order) in the legend order: rows second/third/fourth, columns
    # explicit/implicit
    algos = [(f"{sch}{ORDERNAME[o].capitalize()}Order", sch == "implicit", o)
             for o in (2, 3, 4) for sch in ("explicit", "implicit")]
    dts = sorted({r["dt"] for d in results["outletPressure"]["algorithms"].values()
                  for r in d["levels"]})
    x0, x1 = math.log10(min(dts)), math.log10(max(dts))
    pad = 0.07*(x1 - x0)
    x0, x1 = x0 - 2.2*pad, x1 + pad              # guide labels sit left of the finest level

    def X(dt): return L + (math.log10(dt) - x0)/(x1 - x0)*pw

    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
         f'viewBox="0 0 {W} {H}">', f'<style>{fs.css(W)}</style>',
         f'<rect width="{W}" height="{H}" fill="#ffffff"/>']

    # Legend: columns explicit | implicit, rows second / third / fourth
    for i, (key, implicit, order) in enumerate(algos):
        colour, shape = ORDER_STYLE[order]
        dash = IMPLICIT_DASH if implicit else ""
        lx = L + 20 + (1 if implicit else 0)*pw/2
        ly = 38 + (i//2)*38
        o.append(f'<line x1="{lx}" y1="{ly}" x2="{lx+56}" y2="{ly}" '
                 f'stroke="{colour}" stroke-width="2.4"{dash}/>')
        o.append(marker(shape, lx + 28, ly, colour, implicit))
        o.append(f'<text class="leg" x="{lx+74}" y="{ly+s["legend"]*0.34:.1f}">'
                 f'{"Implicit" if implicit else "Explicit"}, {ORDERNAME[order]} order</text>')

    for n, (qkey, plabel, ylab) in enumerate(PANELS):
        res = results[qkey]["algorithms"]
        py = HEAD + n*(PH + GAP)
        errs = [r["rmse"] for d in res.values() for r in d["levels"]]
        y0 = math.floor(math.log10(min(errs)) - 1.0)
        y1 = math.ceil(math.log10(max(errs)) + 0.2)
        bottom = py + PH

        def Y(v): return py + (y1 - math.log10(v))/(y1 - y0)*PH

        for dec in range(y0, y1 + 1):
            yy = Y(10.0**dec)
            o.append(f'<line x1="{L}" y1="{yy:.2f}" x2="{L+pw}" y2="{yy:.2f}" '
                     f'stroke="{GRID}" stroke-width="1"/>')
            o.append(f'<text class="tick" x="{L-12}" y="{yy+s["tick"]*0.36:.2f}" '
                     f'text-anchor="end">10{sup(dec, s["tick"])}</text>')
        for v in dts:
            xx = X(v)
            o.append(f'<line x1="{xx:.2f}" y1="{py}" x2="{xx:.2f}" y2="{bottom}" '
                     f'stroke="{GRID}" stroke-width="1"/>')
            if n == len(PANELS) - 1:
                o.append(f'<text class="tick" x="{xx:.2f}" y="{bottom+s["tick"]*1.35:.0f}" '
                         f'text-anchor="middle">{f"{v:.5f}".rstrip("0")}</text>')

        # Reference slopes of the theoretical order, drawn under the data: scaled
        # to the geometric mean of the explicit and implicit constants C in
        # e = C dt^p, so each runs between the two series it describes; labelled
        # in the margin left of the finest level
        for order in (2, 3, 4):
            pair = [d for d in res.values() if d["expectedOrder"] == order]
            lnCs = [sum(math.log(r["rmse"]) - order*math.log(r["dt"]) for r in d["levels"])
                    / len(d["levels"]) for d in pair]
            scale = math.exp(sum(lnCs)/len(lnCs))
            lo = min(r["dt"] for d in pair for r in d["levels"])
            hi = max(r["dt"] for d in pair for r in d["levels"])
            ya, yb = scale*lo**order, scale*hi**order
            o.append(f'<line x1="{X(lo):.2f}" y1="{Y(ya):.2f}" x2="{X(hi):.2f}" y2="{Y(yb):.2f}" '
                     f'stroke="{GUIDE}" stroke-width="1.6" stroke-dasharray="3 5"/>')
            o.append(f'<text class="leg" x="{X(lo)-14:.2f}" y="{Y(ya)+s["legend"]*0.36:.2f}" '
                     f'text-anchor="end" fill="{GUIDE}">&#916;t{sup(order, s["legend"])}</text>')

        for key, implicit, order in algos:
            d = res[key]
            colour, _ = ORDER_STYLE[order]
            dash = IMPLICIT_DASH if implicit else ""
            lo = min(r["dt"] for r in d["levels"])
            hi = max(r["dt"] for r in d["levels"])
            fitline = lambda t: math.exp(d["lnC"])*t**d["fittedOrder"]
            o.append(f'<line x1="{X(lo):.2f}" y1="{Y(fitline(lo)):.2f}" '
                     f'x2="{X(hi):.2f}" y2="{Y(fitline(hi)):.2f}" '
                     f'stroke="{colour}" stroke-width="2.2"{dash}/>')
        for key, implicit, order in algos:
            colour, shape = ORDER_STYLE[order]
            for r in res[key]["levels"]:
                o.append(marker(shape, X(r["dt"]), Y(r["rmse"]), colour, implicit))

        # Fitted orders, bottom right
        bx = L + pw - 255
        for i, (key, implicit, order) in enumerate(algos):
            colour, shape = ORDER_STYLE[order]
            by = bottom - 24 - (len(algos) - 1 - i)*36
            o.append(marker(shape, bx, by, colour, implicit))
            o.append(f'<text class="leg" x="{bx+22}" y="{by+s["legend"]*0.34:.1f}">p<tspan dy="6" '
                     f'font-size="{s["legend"]*0.72:.0f}">fit</tspan><tspan dy="-6" dx="6">'
                     f'=&#160;{res[key]["fittedOrder"]:.2f}</tspan></text>')

        o.append(f'<rect x="{L}" y="{py}" width="{pw}" height="{PH}" fill="none" '
                 f'stroke="#000000" stroke-width="{s["axis_line"]}"/>')
        o.append(f'<text class="tick" x="{L+16}" y="{py+s["tick"]*1.5:.0f}">{plabel}</text>')
        if "{}" in ylab:
            ylab = ylab.replace("{}", sup("3", s["axis_label"])
                                + f'<tspan dy="{s["axis_label"]*0.45:.0f}">/s]</tspan>')
        o.append(f'<text class="axlab" transform="translate(52,{py+PH/2:.0f}) rotate(-90)" '
                 f'text-anchor="middle">{ylab}</text>')

    o.append(f'<text class="axlab" x="{L+pw/2:.0f}" y="{H-30}" text-anchor="middle">'
             f'Timestep &#916;t [s]</text>')
    o.append("</svg>")

    svg = HERE/"orderOfAccuracy.svg"
    svg.write_text("\n".join(o))
    conv = shutil.which("rsvg-convert")
    if conv:
        subprocess.run([conv, "-f", "png", "-z", "1.6", "-o",
                        str(svg.with_suffix(".png")), str(svg)], check=True)
        svg.unlink()
        print(f"wrote {svg.stem}.png in {HERE}")
    else:
        print(f"rsvg-convert not found - kept {svg.name} in {HERE}")


if __name__ == "__main__":
    main()
