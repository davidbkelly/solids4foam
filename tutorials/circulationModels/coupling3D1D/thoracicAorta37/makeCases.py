#!/usr/bin/env python3
"""A 3D segment inside the 37-artery network of Boileau et al. (2015):
the thoracic aorta II (segment s17, 72 mm, radius 11 to 9.26 mm, an in vitro
measurement site at its midpoint) as a rigid no-slip pipe in pimpleFoam,
between 8 mm 1D stubs, at physiological amplitude.

    python3 makeCases.py            # the 1D runs
    python3 makeCases.py --restart  # after stubs1D: the runs restarting from
                                    # its periodic state (stubsInertance, rigid3D)
    python3 makeCases.py --dt 5e-4 [--force]

Runs (in runs/):
    stubs1D         the network with s17 split into s17a (8 mm), s17m (56 mm)
                    and s17b (8 mm), all 1D, 18 cycles from rest: should equal
                    the original network, and supplies the periodic state
    stubsInertance  s17m replaced by an inertance node, the exact inviscid
                    rigid limit, restarted from that state, two cycles
    rigid3D         s17m replaced by the pimpleFoam cone (O-grid, 28 axial
                    cells, no-slip, mu = 2.5 mPa s), restarted from that state,
                    two cycles; the 3D starts as plug flow at the state's flow
The original network's results are read from the Boileau case itself.
"""
import argparse
import math
import re
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parents[1]/"circulation1D"/"boileau2015"/"thirtySevenArtery"
RUNS = HERE/"runs"

RHO = 1050.0
MU = 2.5e-3
T_CYCLE = 0.827
L_SEG = 0.072
R_IN, R_OUT = 0.011, 0.00926                 # s17's radii (mesh.txt)
H_WALL = 0.00034
L_STUB = 0.008
L_3D = L_SEG - 2*L_STUB
DX = 0.002
T_STATE = 14.886                              # the 18th cycle's start
N_CYCLES = 2
PROBE_X = 0.036                               # ThoracicAortaII along s17


def radius(x):
    return R_IN + (R_OUT - R_IN)*x/L_SEG


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


def vessel_line(name, x0, x1):
    n = max(2, int(round((x1 - x0)/DX)))
    return (f"    {name:8s} {{ length {x1 - x0:g}; nCells {n}; "
            f"radiusProximal {radius(x0):.6g}; radiusDistal {radius(x1):.6g}; "
            f"wallThickness {H_WALL}; }}")


def properties(kind, dt3d=None):
    """The network dictionary for kind in stubs1D, inertance, coupled"""
    t = (SRC/"constant"/"circulation1DProperties").read_text()
    old = re.search(r"^    s17 +\{[^\n]*\}\n", t, re.M)
    if not old:
        raise SystemExit("vessel s17 not found in the Boileau dictionary")
    lines = [vessel_line("s17a", 0, L_STUB)]
    if kind == "stubs1D":
        lines.append(vessel_line("s17m", L_STUB, L_SEG - L_STUB))
    lines.append(vessel_line("s17b", L_SEG - L_STUB, L_SEG))
    t = t[:old.start()] + "\n".join(lines) + "\n" + t[old.end():]
    t = t.replace("s17.in", "s17a.in").replace("s17.out", "s17b.out")
    n17 = re.search(r"^    n17 +\{[^\n]*\}\n", t, re.M)
    if kind == "stubs1D":
        extra = ("    j17a { type junction; connections (s17a.out s17m.in); }\n"
                 "    j17b { type junction; connections (s17m.out s17b.in); }\n")
    elif kind == "inertance":
        # A cone's inertance rho int dz/A is rho l/(pi R0 R1): radius sqrt(R0 R1)
        r = math.sqrt(radius(L_STUB)*radius(L_SEG - L_STUB))
        extra = (f"    seg3D {{ type inertance; connections (s17a.out s17b.in); "
                 f"length {L_3D:g}; radius {r:.6g}; }}\n")
    else:
        extra = ("    to3D   { type coupledFlow; connection s17a.out; patch inlet; }\n"
                 "    from3D { type coupledFlow; connection s17b.in; patch outlet; }\n")
    t = t[:n17.end()] + extra + t[n17.end():]
    probe = re.search(r"^    ThoracicAortaII +\{[^\n]*\}\n", t, re.M)
    rep = ("    ThoracicAortaII { vessel s17m; x %g; }\n" % (PROBE_X - L_STUB)
           if kind == "stubs1D" else "")
    rep += "    portIn  { vessel s17a; x out; }\n    portOut { vessel s17b; x in; }\n"
    t = t[:probe.start()] + rep + t[probe.end():]
    if kind == "coupled":
        t = t.rstrip()
        t = t[:t.rfind("// ****")].rstrip() + f"""

coupling
{{
    scheme      implicit;
    order       2;
    seedDeltaT  {dt3d:g};
    seed        {{ to3D (SEEDQ); from3D (SEEDQ); }}   // set by --restart
}}
""" + TAIL
    return t


def control_1d(start, end):
    return HEADER % ("system", "controlDict") + f"""application     circulation1DFoam;
startFrom       startTime;
startTime       {start:g};
stopAt          endTime;
endTime         {end:g};
deltaT          1e-3;           // output interval; the solver sub-steps at its CFL step
writeControl    runTime;
writeInterval   {T_CYCLE};    // profiles and restart state
writeFrequency  1;
purgeWrite      0;
writeFormat     ascii;
writePrecision  10;
timeFormat      general;
timePrecision   6;
runTimeModifiable yes;
""" + TAIL


def control_3d(dt, start=T_STATE, end=None):
    if end is None:
        end = T_STATE + N_CYCLES*T_CYCLE
    fos = ""
    for name, z in (("mid", PROBE_X - L_STUB), ("in", DX), ("out", L_3D - DX)):
        for suffix, op, field in (("P", "areaAverage", "p"), ("Q", "sum", "phi")):
            fos += f"""    {name}_{suffix}
    {{
        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");
        regionType      faceZone;
        name            {name};
        operation       {op};
        fields          ({field});
        writeControl    timeStep;
        writeInterval   1;
        writeFields     false;
        log             false;
    }}
"""
    return HEADER % ("system", "controlDict") + f"""application     pimpleFoam;
startFrom       startTime;
startTime       {start:g};
stopAt          endTime;
endTime         {end:g};
deltaT          {dt:g};
writeControl    runTime;
writeInterval   {T_CYCLE};
purgeWrite      0;
writeFormat     ascii;
writePrecision  10;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable yes;
adjustTimeStep  no;

libs            ("libcirculation1DCoupling.so");

functions
{{
{fos}}}
""" + TAIL


def topo_set(zs):
    body = "actions\n(\n"
    for name, z in zs:
        body += f"""    {{ name {name}Faces; type faceSet; action new; source boxToFace;
      box (-1 -1 {z - 1e-4:.6g}) (1 1 {z + 1e-4:.6g}); }}
    {{ name {name}Cells; type cellSet; action new; source boxToCell;
      box (-1 -1 -1) (1 1 {z:.6g}); }}
    {{ name {name}; type faceZoneSet; action new; source setsToFaceZone;
      faceSet {name}Faces; cellSet {name}Cells; }}
"""
    return HEADER % ("system", "topoSetDict") + body + ");\n" + TAIL


def block_mesh(nz):
    """O-grid cone: inner square of half-width 0.3 R, four outer blocks graded
    to the wall, radius R(z) from the stub ends"""
    r0, r1 = radius(L_STUB), radius(L_SEG - L_STUB)
    c = math.cos(math.pi/4)
    verts = []
    for z, r in ((0, r0), (L_3D, r1)):
        verts += [(r, 0, z), (0, r, z), (-r, 0, z), (0, -r, z)]
    for z, r in ((0, r0), (L_3D, r1)):
        s = 0.3*r
        verts += [(s, 0, z), (0, s, z), (-s, 0, z), (0, -s, z)]
    # indices: outer ring z0 0-3, outer ring z1 4-7, inner z0 8-11, inner z1 12-15
    v = "\n".join(f"    ({x:.7g} {y:.7g} {z:.7g})" for x, y, z in verts)
    arcs = ""
    for base, (z, r) in ((0, (0, r0)), (4, (L_3D, r1))):
        for k, (sx, sy) in enumerate(((1, 1), (-1, 1), (-1, -1), (1, -1))):
            a, b = base + k, base + (k + 1) % 4
            arcs += f"    arc {a} {b} ({sx*r*c:.7g} {sy*r*c:.7g} {z:.7g})\n"
    return HEADER % ("system", "blockMeshDict") + f"""// Rigid cone of the thoracic aorta II, radius {r0*1e3:.3f} to {r1*1e3:.3f} mm over
// {L_3D*1e3:g} mm, as an O-grid: 8 x 8 core, four outer blocks of 8 cells
// graded 4:1 towards the wall, {nz} axial cells

scale 1;

vertices
(
{v}
);

blocks
(
    hex ( 8  0  1  9  12  4  5 13) (8 8 {nz}) simpleGrading (0.25 1 1)
    hex ( 9  1  2 10  13  5  6 14) (8 8 {nz}) simpleGrading (0.25 1 1)
    hex (10  2  3 11  14  6  7 15) (8 8 {nz}) simpleGrading (0.25 1 1)
    hex (11  3  0  8  15  7  4 12) (8 8 {nz}) simpleGrading (0.25 1 1)
    hex ( 8  9 10 11  12 13 14 15) (8 8 {nz}) simpleGrading (1 1 1)
);

edges
(
{arcs});

boundary
(
    inlet  {{ type patch; faces ((8 0 1 9) (9 1 2 10) (10 2 3 11) (11 3 0 8) (8 9 10 11)); }}
    outlet {{ type patch; faces ((12 4 5 13) (13 5 6 14) (14 6 7 15) (15 7 4 12) (12 13 14 15)); }}
    wall   {{ type wall;  faces ((0 4 5 1) (1 5 6 2) (2 6 7 3) (3 7 4 0)); }}
);
""" + TAIL


def field_p():
    return (HEADER % ("14.886", "p")).replace("dictionary", "volScalarField") + """
dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    inlet   { type circulation1DPressure; value uniform 0; }
    outlet  { type circulation1DPressure; value uniform 0; }
    wall    { type zeroGradient; }
}
""" + TAIL


def field_U(uz):
    return (HEADER % ("14.886", "U")).replace("dictionary", "volVectorField") + f"""
dimensions      [0 1 -1 0 0 0 0];

// Plug flow at the restart state's flow through the segment
internalField   uniform (0 0 {uz:.6g});

boundaryField
{{
    inlet   {{ type pressureInletOutletVelocity; value uniform (0 0 {uz:.6g}); }}
    outlet  {{ type pressureInletOutletVelocity; value uniform (0 0 {uz:.6g}); }}
    wall    {{ type noSlip; }}
}}
""" + TAIL


def transport():
    return HEADER % ("constant", "transportProperties") + f"""
transportModel  Newtonian;
nu              [0 2 -1 0 0 0 0] {MU/RHO:.6g};
""" + TAIL


def turbulence():
    return HEADER % ("constant", "turbulenceProperties") + "\nsimulationType  laminar;\n" + TAIL


def fv_solution():
    return HEADER % ("system", "fvSolution") + """
solvers
{
    "p|pFinal|pcorr|pcorrFinal"
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-12;
        relTol          0;
        maxIter         5000;
    }
    "(U|UFinal)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-12;
        relTol          0;
        maxIter         1000;
    }
}

// As womersleySegment: no momentum predictor on a viscous, wall-graded
// mesh; the implicit coupling re-integrates the network every outer
// iteration and residualControl on p ends the step
PIMPLE
{
    momentumPredictor       no;
    nOuterCorrectors        50;
    nCorrectors             2;
    nNonOrthogonalCorrectors 0;
    residualControl
    {
        p { tolerance 1e-7; relTol 0; }
    }
}
""" + TAIL


def state_flow(state_file, vessel="s17a"):
    """The last-cell flow of a vessel in a circulation1DState file"""
    t = state_file.read_text()
    m = re.search(r"\n\s+%s\s*\{.*?Q\s+\d+\s*\(([^)]*)\)" % vessel, t, re.S)
    return float(m.group(1).split()[-1])


def finished_1d(case):
    log = case/"log.circulation1DFoam"
    return log.is_file() and "\nEnd\n" in log.read_text()[-2000:]


def make_stubs1d(force):
    case = RUNS/"stubs1D"
    if case.exists() and finished_1d(case) and not force:
        print("stubs1D: finished, kept")
        return
    if case.exists():
        shutil.rmtree(case)
    (case/"constant").mkdir(parents=True)
    (case/"system").mkdir()
    (case/"constant"/"circulation1DProperties").write_text(properties("stubs1D"))
    (case/"system"/"controlDict").write_text(control_1d(0, 15.1))
    print("stubs1D: 18 cycles from rest")


def make_restart_runs(dt, force):
    state = RUNS/"stubs1D"/f"{T_STATE:g}"/"circulation1DState"
    if not state.is_file():
        raise SystemExit(f"no restart state {state}: run stubs1D first")
    q0 = state_flow(state)
    # inertance
    case = RUNS/"stubsInertance"
    if not (case.exists() and finished_1d(case) and not force):
        if case.exists():
            shutil.rmtree(case)
        (case/"constant").mkdir(parents=True)
        (case/"system").mkdir()
        (case/f"{T_STATE:g}").mkdir()
        shutil.copy(state, case/f"{T_STATE:g}"/"circulation1DState")
        (case/"constant"/"circulation1DProperties").write_text(properties("inertance"))
        (case/"system"/"controlDict").write_text(
            control_1d(T_STATE, T_STATE + N_CYCLES*T_CYCLE))
        print(f"stubsInertance: restart at {T_STATE} s, {N_CYCLES} cycles")
    # coupled
    case = RUNS/"rigid3D"
    log = case/"log.pimpleFoam"
    if case.exists() and log.is_file() and "\nEnd\n" in log.read_text()[-2000:] \
            and not force:
        print("rigid3D: finished, kept")
        return
    if case.exists():
        shutil.rmtree(case)
    write_coupled(case, dt, T_STATE, T_STATE + N_CYCLES*T_CYCLE, q0, state)
    r_mean = math.sqrt(radius(L_STUB)*radius(L_SEG - L_STUB))
    print(f"rigid3D: restart at {T_STATE} s with Q = {q0*1e6:.2f} ml/s "
          f"(plug {q0/(math.pi*r_mean**2):.3f} m/s), dt = {dt:g} s, "
          f"{N_CYCLES} cycles, {int(round(L_3D/DX))} axial cells")


def write_coupled(case, dt, start, end, q0=0.0, state=None):
    """The coupled pimpleFoam case from start to end: restarted from a
    circulation1DState file with the 3D as plug flow at q0 (the tutorial), or
    from rest (state None, the order study)"""
    for d in ("constant", "system", f"{start:g}"):
        (case/d).mkdir(parents=True)
    if state is not None:
        shutil.copy(state, case/f"{start:g}"/"circulation1DState")
    props = properties("coupled", dt).replace("SEEDQ", f"{q0:.9g}")
    (case/"constant"/"circulation1DProperties").write_text(props)
    (case/"constant"/"transportProperties").write_text(transport())
    (case/"constant"/"turbulenceProperties").write_text(turbulence())
    nz = int(round(L_3D/DX))
    (case/"system"/"blockMeshDict").write_text(block_mesh(nz))
    (case/"system"/"controlDict").write_text(control_3d(dt, start, end))
    (case/"system"/"topoSetDict").write_text(
        topo_set([("mid", PROBE_X - L_STUB), ("in", DX), ("out", L_3D - DX)]))
    (case/"system"/"fvSolution").write_text(fv_solution())
    shutil.copy(HERE.parent/"rigidSegment3D"/"base"/"system"/"fvSchemes",
                case/"system"/"fvSchemes")
    r_mean = math.sqrt(radius(L_STUB)*radius(L_SEG - L_STUB))
    (case/f"{start:g}"/"p").write_text(field_p())
    (case/f"{start:g}"/"U").write_text(field_U(q0/(math.pi*r_mean**2)))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--restart", action="store_true",
                    help="write the runs that restart from stubs1D's state")
    ap.add_argument("--dt", type=float, default=1e-3, help="3D time step")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    RUNS.mkdir(exist_ok=True)
    if args.restart:
        make_restart_runs(args.dt, args.force)
    else:
        make_stubs1d(args.force)


if __name__ == "__main__":
    main()
