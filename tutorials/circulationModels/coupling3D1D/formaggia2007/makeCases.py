#!/usr/bin/env python3
"""Write the runs of the Formaggia et al. (2007) 3D-1D test, section 5: the
pressure-pulse tube of the solids4foam 3dTube tutorial (Formaggia et al. 2001,
Greenshields & Weller 2005) continued by a 1D vessel.

    python3 makeCases.py [--c0 C0] [--dt DT] [--force]

Runs (in runs/):
    full3D    the "exact solution": the FSI tube three times as long (15 cm),
              so that nothing reflected from its far end returns within 20 ms
    coupled   the 5 cm FSI tube with a 10 cm 1D vessel on its outlet
              (circulation1DPressure), absorbing at the far end
    only1D    the same pulse in a 15 cm 1D vessel (pressureInlet)
    coupled_c<c0>[_free]  with --c0 / --freeEnd: another coupled run with
              that 1D wave speed, and/or the wall's end at the interface free
              (zero traction) instead of clamped

The 1D wall is matched to the 3D one through its static compliance: a thick
cylinder in plane strain (fixed ends) with inner radius a = 5 mm, outer
b = 6 mm, E = 300 kPa, nu = 0.3 has u(a)/p = a (1 + nu) [(1 - 2 nu) a^2 +
b^2] / (E (b^2 - a^2)), so dA/dp = 2 pi a u/p and c0 = sqrt(A/(rho dA/dp)) =
5.25 m/s; the sqrtElastic law gives c0^2 = 2 E1D h/(3 rho R), so E1D = 3 rho R
c0^2/(2 h) (the law's 4/3 factor is the nu = 1/2 of Formaggia's 1D model).
--c0 overrides this, e.g. with a wave speed measured from full3D.

The 3D cases are quarter tubes (two symmetry planes): the coupled port has
flowScale 4. Sections of the fluid at z = 1.25, 2.5, ... cm are faceZones
(topoSet) on which surfaceFieldValue writes the mean pressure and the flow.
"""
import argparse
import math
import re
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = HERE/"base"
RUNS = HERE/"runs"

RHO = 1000.0
MU = 3e-3
R = 5e-3
H = 1e-3
E = 3e5
NU = 0.3
A = math.pi*R**2
L3D = 0.05                      # the tutorial tube
LFULL = 0.15                    # the reference tube
L1D = LFULL - L3D               # the vessel continuing the coupled tube
CELLS_PER_M = 1600              # 80 axial cells per 5 cm, as the tutorial
PULSE = 1333.2                  # Pa, for 3 ms
END_TIME = 0.02
DELTA_T = 2.5e-5                # the tutorial's finer step (3 FSI iterations per
                                # step; at 1e-4 it needs 40 and is no cheaper)
PLANES_CM = [1.25, 2.5, 3.75, 5, 6.25, 7.5, 10, 12.5]   # z of the sections


def c0_plane_strain():
    a, b = R, R + H
    u_over_p = a*(1 + NU)*((1 - 2*NU)*a*a + b*b)/(E*(b*b - a*a))
    dAdp = 2*math.pi*a*u_over_p
    return math.sqrt(A/(RHO*dAdp))


def youngs1D(c0):
    """E of the sqrtElastic law for wave speed c0 at the reference area"""
    return 1.5*RHO*R*c0**2/H


HEADER = """/*--------------------------------*- C++ -*----------------------------------*\\
| solids4foam: solid mechanics and fluid-solid interaction simulations        |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "%s";
    object      %s;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""
TAIL = "\n// ************************************************************************* //\n"


def zone_name(zcm):
    return "z" + f"{zcm:g}".replace(".", "p")


def control_dict(coupled, planes, dt):
    libs ='libs            ("libcirculation1DCoupling.so");\n\n' if coupled else ""
    fos = ""
    for zcm in planes:
        n = zone_name(zcm)
        for suffix, op, field in (("P", "areaAverage", "p"), ("Q", "sum", "phi")):
            fos += f"""    {n}_{suffix}
    {{
        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");
        region          fluid;
        regionType      faceZone;
        name            {n};
        operation       {op};
        fields          ({field});
        writeControl    timeStep;
        writeInterval   1;
        writeFields     false;
        log             false;
    }}
"""
    body = f"""application     solids4Foam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {END_TIME};
deltaT          {dt:g};
writeControl    runTime;
writeInterval   1e-3;
purgeWrite      0;
writeFormat     binary;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable yes;
adjustTimeStep  no;

{libs}functions
{{
    // Radial displacement of the wall at z = 2.5 cm, as the tutorial
    displacement
    {{
        type            solidPointDisplacement;
        region          solid;
        point           (0 0.005 0.025);
    }}
{fos}}}
"""
    return HEADER % ("system", "controlDict") + body + TAIL


def topo_set_dict(planes):
    body = "actions\n(\n"
    for zcm in planes:
        n = zone_name(zcm)
        z = zcm*1e-2
        body += f"""    {{
        name    {n}Faces;
        type    faceSet;
        action  new;
        source  boxToFace;
        box     (-1 -1 {z - 1e-4:.6g}) (1 1 {z + 1e-4:.6g});
    }}
    {{
        name    {n}Cells;
        type    cellSet;
        action  new;
        source  boxToCell;
        box     (-1 -1 -1) (1 1 {z:.6g});
    }}
    {{
        name    {n};
        type    faceZoneSet;
        action  new;
        source  setsToFaceZone;
        faceSet {n}Faces;
        cellSet {n}Cells;
    }}
"""
    body += ");\n"
    return HEADER % ("system", "topoSetDict") + body + TAIL


def properties_coupled(c0):
    probes = "".join(f"    {zone_name(z)} {{ vessel distal; x {z*1e-2 - L3D:.6g}; }}\n"
                     for z in PLANES_CM if z*1e-2 > L3D)
    body = f"""// The 1D continuation of the 5 cm FSI tube (Formaggia et al. 2007, section
// 5): a 10 cm vessel whose wall has the static compliance of the 3D one,
// c0 = {c0:.4g} m/s, Poiseuille friction as their 1D model
fluid
{{
    rho             {RHO};
    mu              {MU};
    profileOrder    2;
    alpha           1;
}}
numerics
{{
    CFL             0.5;
    flux            HLL;
    limiter         MC;
    timeScheme      SSPRK3;
}}
vessels
{{
    distal
    {{
        length              {L1D};
        nCells              {int(round(L1D*CELLS_PER_M))};
        radius              {R};
        wallThickness       {H};
        youngsModulus       {youngs1D(c0):.6e};
        diastolicPressure   0;
        externalPressure    0;
        initialPressure     0;
    }}
}}
nodes
{{
    from3D                              // fed by the 3D outlet; quarter tube
    {{
        type        coupledFlow;
        connection  distal.in;
        patch       outlet;
        flowScale   4;
    }}
    outlet
    {{
        type        absorbing;
        connection  distal.out;
    }}
}}
probes                                  // x along the vessel = z - 5 cm
{{
{probes}}}
coupling
{{
    scheme      implicit;
    order       2;
    seed        zero;
}}
"""
    return HEADER % ("constant/fluid", "circulation1DProperties") + body + TAIL


def properties_only1D(c0):
    probes = "".join(f"    {zone_name(z)} {{ vessel tube; x {z*1e-2:.6g}; }}\n"
                     for z in PLANES_CM)
    body = f"""// The pulse tube as a 1D vessel alone, c0 = {c0:.4g} m/s
fluid
{{
    rho             {RHO};
    mu              {MU};
    profileOrder    2;
    alpha           1;
}}
numerics
{{
    CFL             0.5;
    flux            HLL;
    limiter         MC;
    timeScheme      SSPRK3;
}}
vessels
{{
    tube
    {{
        length              {LFULL};
        nCells              {int(round(LFULL*CELLS_PER_M))};
        radius              {R};
        wallThickness       {H};
        youngsModulus       {youngs1D(c0):.6e};
        diastolicPressure   0;
        externalPressure    0;
        initialPressure     0;
    }}
}}
nodes
{{
    inlet
    {{
        type        pressureInlet;
        connection  tube.in;
        pressure
        {{
            type        table;
            values      ((0 {PULSE}) (3e-3 {PULSE}) (3.1e-3 0) (1 0));
        }}
    }}
    outlet
    {{
        type        absorbing;
        connection  tube.out;
    }}
}}
probes
{{
{probes}}}
"""
    return HEADER % ("constant", "circulation1DProperties") + body + TAIL


def control_dict_1d():
    body = f"""application     circulation1DFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {END_TIME};
deltaT          1e-5;
writeControl    runTime;
writeInterval   1e-3;
writeFrequency  1;
purgeWrite      0;
writeFormat     ascii;
writePrecision  8;
timeFormat      general;
timePrecision   6;
runTimeModifiable yes;
"""
    return HEADER % ("system", "controlDict") + body + TAIL


def edit(path, pairs):
    text = path.read_text()
    for old, new in pairs:
        text, n = re.subn(old, new, text)
        if n == 0:
            raise SystemExit(f"{path}: nothing matched {old!r}")
    path.write_text(text)


def finished(case):
    log = case/"log.solids4Foam"
    if not log.is_file():
        log = case/"log.circulation1DFoam"
    return log.is_file() and "\nEnd\n" in log.read_text()[-2000:]


def make_3d(name, length, coupled, c0, dt, force, free_end=False):
    case = RUNS/name
    if case.exists() and finished(case) and not force:
        print(f"{name}: finished, kept")
        return
    if case.exists():
        shutil.rmtree(case)
    shutil.copytree(BASE, case)
    n = int(round(length*CELLS_PER_M))
    for region in ("fluid", "solid"):
        edit(case/"system"/region/"blockMeshDict",
             [(r"(\(\s*[-\d.e]+\s+[-\d.e]+\s+)0\.05\)", rf"\g<1>{length:g})"),
              (r" 80\)", f" {n})")])
    planes = [z for z in PLANES_CM if z*1e-2 < length]
    (case/"system"/"controlDict").write_text(control_dict(coupled, planes, dt))
    (case/"system"/"fluid"/"topoSetDict").write_text(topo_set_dict(planes))
    if coupled:
        edit(case/"0"/"fluid"/"p",
             [(r"outlet\s*\{\s*type\s+fixedValue;\s*value\s+uniform 0;\s*\}",
               "outlet\n    {\n        type        circulation1DPressure;\n"
               "        value       uniform 0;\n    }")])
        (case/"constant"/"fluid"/"circulation1DProperties").write_text(
            properties_coupled(c0))
    if free_end:
        # The wall's end at the interface free to breathe (zero traction)
        # instead of clamped
        edit(case/"0"/"solid"/"D",
             [(r"outlet\s*\{\s*type\s+fixedDisplacement;\s*value\s+uniform \(0 0 0\);\s*\}",
               "outlet\n    {\n        type            solidTraction;\n"
               "        traction        uniform ( 0 0 0 );\n"
               "        pressure        uniform 0;\n"
               "        value           uniform (0 0 0);\n    }")])
        edit(case/"0"/"solid"/"pointD",
             [(r"outlet\s*\{\s*type\s+fixedValue;",
               "outlet\n    {\n        type            calculated;")])
    print(f"{name}: {length*100:g} cm, {n} axial cells, dt = {dt:g} s"
          + (f", 1D c0 = {c0:.4g} m/s" if coupled else "")
          + (", free wall end at the interface" if free_end else ""))


def make_1d(name, c0, force):
    case = RUNS/name
    if case.exists() and finished(case) and not force:
        print(f"{name}: finished, kept")
        return
    if case.exists():
        shutil.rmtree(case)
    (case/"system").mkdir(parents=True)
    (case/"constant").mkdir()
    (case/"system"/"controlDict").write_text(control_dict_1d())
    (case/"constant"/"circulation1DProperties").write_text(properties_only1D(c0))
    print(f"{name}: {LFULL*100:g} cm 1D vessel, c0 = {c0:.4g} m/s")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--c0", type=float, default=None,
                    help="1D wave speed [m/s] (default: static plane-strain match)")
    ap.add_argument("--dt", type=float, default=DELTA_T,
                    help=f"3D time step [s] (default {DELTA_T:g})")
    ap.add_argument("--freeEnd", action="store_true",
                    help="coupled run with the wall's end free at the interface")
    ap.add_argument("--force", action="store_true",
                    help="rewrite finished runs too")
    args = ap.parse_args()
    c0 = args.c0 if args.c0 else c0_plane_strain()
    print(f"plane-strain static match: c0 = {c0_plane_strain():.4g} m/s; "
          f"using {c0:.4g} m/s, E1D = {youngs1D(c0):.4g} Pa")
    RUNS.mkdir(exist_ok=True)
    make_3d("full3D", LFULL, False, c0, args.dt, args.force)
    if args.c0 or args.freeEnd:
        # Another coupled run beside the static-match one, named by its c0
        # and end condition
        name = ("coupled_c" + f"{c0:.3g}".replace(".", "p")
                + ("_free" if args.freeEnd else ""))
        make_3d(name, L3D, True, c0, args.dt, args.force, args.freeEnd)
    else:
        make_3d("coupled", L3D, True, c0, args.dt, args.force)
        make_1d("only1D", c0, args.force)


if __name__ == "__main__":
    main()
