#!/usr/bin/env python3
"""Write constant/circulation1DProperties of the 37-artery network from its mesh table.

    python3 makeNetwork.py                  # tapered vessels, dx 2 mm
    python3 makeNetwork.py --taper pieces --step 0.01 --dx 1e-3 \
        --output other/constant/circulation1DProperties

The table is referenceData/mesh.txt of the data supplement of Boileau et al.
(2015), section 4.5 (Matthys et al. 2007): one row per segment,

    segment  first-node  last-node  length [m]  initial radius [m]  final radius [m]
    wall thickness [m]  PWV start [m/s]  PWV end [m/s]  terminal resistance

The model of the paper (Table IX): rho 1050, mu 2.5 mPa s, zeta 9, E 1.2 MPa,
Pd = Pext = 0 so that Ad = A0 is the area of the table, Pout 432.6 Pa, starting
from rest at P = 0. E and the wall thickness reproduce the PWV columns to 2%.

Two decisions this script makes, which the data supplement does not settle:

  * The terminal resistance column is labelled Pa s/m3 but holds values of
    2.6 to 5.7: they are in 1e9 Pa s/m3 (with those, the mean pressure is
    12 kPa, as measured; as labelled it would be 0.4 kPa).
  * 22 segments taper. The radius is taken linear along the segment. By
    default (--taper native) each is one tapered vessel (radiusProximal,
    radiusDistal). With --taper pieces a tapered segment is instead split
    into uniform pieces joined by two-vessel junctions, so that no radius
    jump exceeds --step: the approach of phase 2, kept for comparison; the
    pieces converge to the tapered vessel as --step goes to zero.

Probes are written at the midpoint of every segment: the eight measured sites
under their names, the others as seg<NN>.
"""
import argparse
import math
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent

RHO = 1050
MU = 2.5e-3
ZETA = 9
E = 1.2e6
POUT = 432.6
R_UNIT = 1e9

# Segment of each measurement site of the paper (Figs. 6-9, Table X). The
# table has no names: each site is the segment whose midpoint P and Q match the
# published FVM waveforms of that site (relative RMS 0.2-0.4%, against 4% or
# more for any other segment), consistent with the topology, lengths and radii.
SITES = {
    10: "AorticArchII",
    17: "ThoracicAortaII",
    11: "LSubclavianI",
    29: "RIliacFemoralII",
    14: "LUlnar",
    34: "RAnteriorTibial",
    7: "RUlnar",
    20: "Splenic",
}

BANNER = r"""/*--------------------------------*- C++ -*----------------------------------*\
| solids4foam: solid mechanics and fluid-solid interaction simulations        |
| Version:     v2.4                                                           |
| Web:         https://solids4foam.github.io                                  |
| Disclaimer:  This offering is not approved or endorsed by OpenCFD Limited,  |
|              producer and distributor of the OpenFOAM software via          |
|              www.openfoam.com, and owner of the OPENFOAM® and OpenCFD®      |
|              trade marks.                                                   |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      circulation1DProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""


def read_mesh(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        f = line.split()
        if len(f) < 10 or line.lstrip().startswith("#"):
            continue
        rows.append(dict(
            seg=int(f[0]), a=int(f[1]), b=int(f[2]), L=float(f[3]),
            r0=float(f[4]), r1=float(f[5]), h=float(f[6]), R=float(f[9]),
        ))
    return rows


def inflow_expression(path):
    """The Fourier series of inflow.txt [ml/s] as an expression in m3/s."""
    s = Path(path).read_text().replace("\r", "")
    period = re.search(r"T=([\d.]+)", s).group(1)
    e = re.search(r"inflow=(.*)", s).group(1).strip().replace("PI", "pi()")
    e = re.sub(r"\bt\b", "arg()", e)
    terms = re.split(r"(?<![eE])(?=[+-][\d.]+(?:[eE][+-]?\d+)?\*(?:sin|cos))", e)
    lines, cur = [], ""
    for t in terms:
        if len(cur) + len(t) > 68:
            lines.append(cur)
            cur = ""
        cur += t
    lines.append(cur)
    return period, "1e-6*" + "\n                ".join(lines)


def pieces(row, step):
    """Uniform pieces of a segment: (length, radius at the piece midpoint)."""
    r0, r1, L = row["r0"], row["r1"], row["L"]
    rel = abs(r1 - r0)/(0.5*(r0 + r1))
    n = max(1, math.ceil(rel/step - 1e-9)) if rel > 0 else 1
    return [(L/n, r0 + (r1 - r0)*(k + 0.5)/n) for k in range(n)]


def vessel_name(seg, k, n):
    return f"s{seg:02d}" if n == 1 else f"s{seg:02d}_{k + 1}"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--taper", choices=["native", "pieces"], default="native",
                    help="tapered vessels, or uniform pieces joined by junctions")
    ap.add_argument("--step", type=float, default=0.02,
                    help="largest relative radius step between pieces")
    ap.add_argument("--dx", type=float, default=2e-3, help="target cell size [m]")
    ap.add_argument("--output", default=str(HERE/"constant"/"circulation1DProperties"))
    args = ap.parse_args()

    rows = read_mesh(HERE/"referenceData"/"mesh.txt")
    period, expr = inflow_expression(HERE/"referenceData"/"inflow.txt")

    vessels, nodes, probes = [], [], []
    ends_at = {}            # mesh node -> list of vessel ends
    n_pieces = 0
    n_internal = 0

    for row in rows:
        native = args.taper == "native"
        ps = [(row["L"], None)] if native else pieces(row, args.step)
        n = len(ps)
        n_pieces += n
        for k, (L, r) in enumerate(ps):
            name = vessel_name(row["seg"], k, n)
            nc = max(2, round(L/args.dx))
            radius = (
                f"radiusProximal {row['r0']:.6g}; radiusDistal {row['r1']:.6g};"
                if native and row["r0"] != row["r1"]
                else f"radius {row['r0'] if native else r:.6g};"
            )
            vessels.append(
                f"    {name:<8} {{ length {L:.6g}; nCells {nc}; {radius} "
                f"wallThickness {row['h']:.6g}; }}"
            )
            if k > 0:
                prev = vessel_name(row["seg"], k - 1, n)
                nodes.append(
                    f"    {name}_j {{ type junction; connections ({prev}.out {name}.in); }}"
                )
                n_internal += 1
        ends_at.setdefault(row["a"], []).append(vessel_name(row["seg"], 0, n) + ".in")
        ends_at.setdefault(row["b"], []).append(vessel_name(row["seg"], n - 1, n) + ".out")

        # Probe at the segment midpoint
        x, k = 0.5*row["L"], 0
        while x > ps[k][0] + 1e-12:
            x -= ps[k][0]
            k += 1
        label = SITES.get(row["seg"], f"seg{row['seg']:02d}")
        probes.append(
            f"    {label:<16} {{ vessel {vessel_name(row['seg'], k, n)}; x {x:.6g}; }}"
        )

    by_seg_end = {r["b"]: r for r in rows}
    for node in sorted(ends_at):
        ends = ends_at[node]
        if node == 0:
            nodes.insert(0, (
                "    inlet\n    {\n        type        flowInlet;\n"
                f"        connection  {ends[0]};\n"
                "        flowRate\n        {\n            type        expression;\n"
                f"            variables   ( \"T = {period}\" );\n"
                f"            expression\n            #{{\n                {expr}\n            #}};\n"
                "        }\n    }"
            ))
        elif len(ends) == 1:
            r = by_seg_end[node]
            if r["R"] <= 0:
                raise SystemExit(f"Node {node} ends segment {r['seg']} without a resistance")
            nodes.append(
                f"    out{r['seg']:02d} {{ type resistance; connection {ends[0]}; "
                f"R {r['R']*R_UNIT:.6g}; Pout {POUT}; }}"
            )
        else:
            nodes.append(
                f"    n{node:02d} {{ type junction; connections ({' '.join(ends)}); }}"
            )

    how = (
        f"{len(rows)} segments, tapered ones as tapered vessels"
        if args.taper == "native" else
        f"{len(rows)} segments as {n_pieces} uniform pieces with radius\n"
        f"// steps of at most {args.step:g}, joined by {n_internal} two-vessel junctions"
    )
    text = BANNER + f"""
// The 37-artery network of Boileau et al. (2015) section 4.5, written by
// makeNetwork.py from referenceData/mesh.txt, dx {args.dx:g} m:
// {how}.
// Do not edit by hand: change makeNetwork.py and run it again.

fluid
{{
    rho             {RHO};
    mu              {MU};
    profileOrder    {ZETA};
    alpha           1;
}}

numerics
{{
    CFL             0.5;
    flux            HLL;
    limiter         MC;
    timeScheme      SSPRK2;
}}

// Every vessel: E {E:g} Pa, Pd = Pext = 0 (A = Ad at P = 0), from rest
vesselDefaults
{{
    youngsModulus       {E:g};
    diastolicPressure   0;
    externalPressure    0;
    initialPressure     0;
}}

vessels
{{
{chr(10).join(vessels)}
}}

nodes
{{
{chr(10).join(nodes)}
}}

probes
{{
{chr(10).join(probes)}
}}
"""
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(text + "\n// ************************************************************************* //\n")
    print(f"{args.output}: {n_pieces} vessels, {n_internal} taper junctions, "
          f"{len(nodes)} nodes")


if __name__ == "__main__":
    main()
