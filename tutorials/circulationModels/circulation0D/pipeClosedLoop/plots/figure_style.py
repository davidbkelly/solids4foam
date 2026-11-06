"""Shared text and line styling for the order-of-accuracy figures.

The sizes are fractions of the canvas width, and each plotting script resolves
them against its own canvas via `sizes()`, so figures drawn on different canvas
sizes keep the same visual text size once scaled to the same printed width.

Keeping the font stack, weights, and line widths in one place is what stops the
figures drifting apart.
"""

from __future__ import annotations

# Serif, to match journal body text and because the axis labels carry maths
# (Delta t, RMSE, superscripts).
FONT_STACK = '"Times New Roman", Times, "DejaVu Serif", serif'
FONT_FILL = "#000000"

# Axis labels and tick labels are bold; legend text stays regular so it does
# not compete with them.
WEIGHT_AXIS_LABEL = "700"
WEIGHT_TICK = "700"
WEIGHT_LEGEND = "400"

# Fractions of canvas width. Calibrated at 26/30/24 px on a 1240 px canvas.
FRAC_TICK = 26/1240
FRAC_AXIS_LABEL = 30/1240
FRAC_LEGEND = 24/1240

# Times is narrower than DejaVu Sans, so legend layout needs its own metric.
CHAR_W = 0.48

# Line widths, also as fractions of canvas width so they scale together.
FRAC_AXIS_LINE = 3.2/1240      # the frame; "bold" axes
FRAC_TICK_MARK = 2.6/1240


def sizes(canvas_width: int) -> dict[str, float]:
    """Resolve the fractional sizes against a canvas width, in px."""
    return {
        "tick": round(FRAC_TICK*canvas_width),
        "axis_label": round(FRAC_AXIS_LABEL*canvas_width),
        "legend": round(FRAC_LEGEND*canvas_width),
        "axis_line": round(FRAC_AXIS_LINE*canvas_width, 2),
        "tick_mark": round(FRAC_TICK_MARK*canvas_width, 2),
    }


def css(canvas_width: int) -> str:
    """The shared `<style>` body: font stack plus the three text classes."""
    s = sizes(canvas_width)
    return (
        f"text{{font-family:{FONT_STACK};fill:{FONT_FILL}}}"
        f".tick{{font-size:{s['tick']}px;font-weight:{WEIGHT_TICK}}}"
        f".axlab{{font-size:{s['axis_label']}px;font-weight:{WEIGHT_AXIS_LABEL}}}"
        f".leg{{font-size:{s['legend']}px;font-weight:{WEIGHT_LEGEND}}}"
    )
