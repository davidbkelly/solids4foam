#!/usr/bin/env python3
"""Write constant/circulation1DProperties of the ADAN56 network from its mesh table.

    python3 makeNetwork.py                  # dx 2 mm
    python3 makeNetwork.py --dx 1e-3 --output other/constant/circulation1DProperties

The table is referenceData/mesh.txt of the data supplement of Boileau et al.
(2015), section 4.6 (Blanco et al. 2015), in cgs units: one row per segment,

    segment  first-node  last-node  length [cm]  initial radius [cm]
    final radius [cm]  R1 [dyn s/cm5]  R2 [dyn s/cm5]  C [cm5/dyn]  name

with R1, R2, C the three-element Windkessel of a terminal segment (-1
otherwise). The model of the paper (Table XI): rho 1040, mu 4 mPa s, zeta 2,
E 225 kPa everywhere, Pd = P0 = 10 kPa with the radii of the table at P0
(Ad = A0), Pext = Pout = 0, the wall thickness of Eq. 17, and the initial state
(A0, 0, P0). Everything is converted to SI: 1 dyn s/cm5 = 1e5 Pa s/m3,
1 cm5/dyn = 1e-5 m3/Pa. The radius is linear along a tapered segment.

Probes are written at the midpoint of the twelve segments of Figs. 11-14,
under the names of the data supplement, and every 5 mm along the aorta (for
Figs. 15-16), as aortaNNN with NNN the distance from the aortic root in mm.
"""
import argparse
from pathlib import Path

HERE = Path(__file__).resolve().parent

RHO = 1040
MU = 4e-3
ZETA = 2
E = 225e3
P0 = 10e3
CM = 1e-2
R_SI = 1e5          # dyn s/cm5 -> Pa s/m3
C_SI = 1e-5         # cm5/dyn -> m3/Pa

# The twelve sites of the data supplement and their segments
SITES = {
    "AorticArchI": "aortic_arch_I",
    "ThoAortaIII": "thoracic_aorta_III",
    "AbdomAortaV": "abdominal_aorta_V",
    "RightCCA": "common_carotid_R",
    "RightRenal": "renal_R",
    "RightCommonIliac": "common_iliac_R",
    "RightICA": "internal_carotid_R",
    "RightRadial": "radial_R",
    "RightInternalIliac": "internal_iliac_R",
    "RightPostInteross": "posterior_interosseous_R",
    "RightFemoralII": "femoral_R_II",
    "RightAntTibial": "anterior_tibial_R",
}

# The aorta, from the root to the iliac bifurcation (Figs. 15-16)
AORTA = [
    "aortic_arch_I", "aortic_arch_II", "aortic_arch_III", "aortic_arch_IV",
    "thoracic_aorta_I", "thoracic_aorta_II", "thoracic_aorta_III",
    "thoracic_aorta_IV", "thoracic_aorta_V", "thoracic_aorta_VI",
    "abdominal_aorta_I", "abdominal_aorta_II", "abdominal_aorta_III",
    "abdominal_aorta_IV", "abdominal_aorta_V",
]
AORTA_PROBE_SPACING = 5e-3

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
    for line in Path(path).read_text().replace("\r", "").splitlines():
        f = line.split()
        if len(f) < 10:
            continue
        rows.append(dict(
            seg=f[0], a=int(f[1]), b=int(f[2]), L=float(f[3])*CM,
            r0=float(f[4])*CM, r1=float(f[5])*CM,
            R1=float(f[6]), R2=float(f[7]), C=float(f[8]), name=f[9],
        ))
    return rows


def read_inflow(path):
    """(t [s], Q [m3/s]) pairs of the inflow table [ml/s]."""
    pts = []
    for line in Path(path).read_text().replace("\r", "").splitlines():
        f = line.split()
        if len(f) == 2:
            pts.append((float(f[0]), float(f[1])*1e-6))
    return pts


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dx", type=float, default=2e-3, help="target cell size [m]")
    ap.add_argument("--output", default=str(HERE/"constant"/"circulation1DProperties"))
    args = ap.parse_args()

    rows = read_mesh(HERE/"referenceData"/"mesh.txt")
    by_name = {r["name"]: r for r in rows}
    inflow = read_inflow(HERE/"referenceData"/"inflow.txt")

    vessels, nodes, probes = [], [], []
    ends_at = {}
    for r in rows:
        nc = max(2, round(r["L"]/args.dx))
        radius = (
            f"radius {r['r0']:.8g};" if r["r0"] == r["r1"] else
            f"radiusProximal {r['r0']:.8g}; radiusDistal {r['r1']:.8g};"
        )
        vessels.append(f"    {r['name']:<28} {{ length {r['L']:.8g}; nCells {nc}; {radius} }}")
        ends_at.setdefault(r["a"], []).append(r["name"] + ".in")
        ends_at.setdefault(r["b"], []).append(r["name"] + ".out")

    table = "\n".join(f"                ({t:g} {q:.8g})" for t, q in inflow)
    by_end = {r["b"]: r for r in rows}
    for node in sorted(ends_at):
        ends = ends_at[node]
        if node == 0:
            nodes.append(
                "    inlet\n    {\n        type        flowInlet;\n"
                f"        connection  {ends[0]};\n"
                "        flowRate                        // referenceData/inflow.txt, T = 1 s\n"
                "        {\n            type        table;\n"
                "            outOfBounds repeat;\n            values\n            (\n"
                f"{table}\n            );\n        }}\n    }}"
            )
        elif len(ends) == 1:
            r = by_end[node]
            if r["R1"] <= 0:
                raise SystemExit(f"Node {node} ends {r['name']} without a Windkessel")
            nodes.append(
                f"    {r['name']}_wk\n    {{\n        type        lumped0D;\n"
                f"        connection  {ends[0]};\n"
                "        circuit\n        {\n"
                "            components\n            {\n"
                "                port { type externalPort; direction outlet; }\n"
                f"                wk   {{ type windkessel3; Rch {r['R1']*R_SI:.8g}; "
                f"C {r['C']*C_SI:.8g}; R {r['R2']*R_SI:.8g}; initialPressure {P0:g}; }}\n"
                "            }\n            connections ( (port wk) );\n        }\n    }"
            )
        else:
            nodes.append(f"    n{node:02d} {{ type junction; connections ({' '.join(ends)}); }}")

    for site, seg in SITES.items():
        probes.append(f"    {site:<20} {{ vessel {seg}; x {0.5*by_name[seg]['L']:.8g}; }}")

    # Along the aorta, checking that the segments follow one another
    for prev, nxt in zip(AORTA[:-1], AORTA[1:]):
        if by_name[prev]["b"] != by_name[nxt]["a"]:
            raise SystemExit(f"{prev} does not lead to {nxt}")
    total = sum(by_name[s]["L"] for s in AORTA)
    k = 0
    while k*AORTA_PROBE_SPACING <= total + 1e-12:
        dist = k*AORTA_PROBE_SPACING
        d = dist
        for s in AORTA:
            L = by_name[s]["L"]
            if d <= L + 1e-12:
                probes.append(f"    aorta{round(dist*1000):03d}           "
                              f"{{ vessel {s}; x {min(d, L):.8g}; }}")
                break
            d -= L
        k += 1

    text = BANNER + f"""
// The ADAN56 network of Boileau et al. (2015) section 4.6, written by
// makeNetwork.py from referenceData/mesh.txt, dx {args.dx:g} m: {len(rows)} segments
// (tapered ones as tapered vessels), cgs converted to SI.
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

// Every vessel: E {E:g} Pa, Pd = P0 = 10 kPa (Ad = A0, the area of the table),
// the wall thickness of Boileau Eq. 17, from (A0, 0, P0)
vesselDefaults
{{
    youngsModulus       {E:g};
    diastolicPressure   {P0:g};
    externalPressure    0;
    initialPressure     {P0:g};
    wallThicknessModel  ADAN56;
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
    print(f"{args.output}: {len(rows)} vessels, {len(nodes)} nodes, "
          f"{len(probes)} probes, aorta {total:.4f} m")


if __name__ == "__main__":
    main()
