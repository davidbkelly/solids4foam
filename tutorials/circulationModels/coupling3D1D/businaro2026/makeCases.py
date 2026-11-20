#!/usr/bin/env python3
"""Write the runs of the "pressure Gaussian pulse" 3D-1D test of Businaro
(2026, MSc thesis, Politecnico di Milano, Section 7.4 and Table F.1), a
compliant FSI tube continued by a 1D vessel, computed there with lifex.

    python3 makeCases.py [--dt DT] [--linearElastic] [--withRt] [--force]

Runs (in runs/):
    full3D            the reference: the FSI tube 45 cm long, p = 0 at its far
                      end, long enough that nothing reflected there reaches
                      z <= 30 cm within the 0.1 s window
    full3D_Rt         with --withRt: Businaro's own 3D reference, the tube
                      30 cm long with the terminal resistance Rt at its outlet
                      (through a 1 cm 1D stub ending in a resistance node)
    coupled           10 cm of FSI tube with a 20 cm 1D vessel on its outlet,
                      absorbing at the far end; implicit coupling of order 2
    coupled_explicit  the same with explicit coupling of order 1: the 1D
                      vessel advanced with the 3D flow of the previous step,
                      the 3D tube with the 1D pressure at the end of the step,
                      Businaro's algorithm (his Eq. 5.21)
    only1D            the 30 cm tube as a 1D vessel alone

The tube (Table F.1): R = 1 cm, wall 1 mm, rho_s = 1000 kg/m3, E = 700 kPa,
nu = 0.3, Saint Venant-Kirchhoff; blood rho = 1060 kg/m3, mu = 4 mPa s; the
inlet pressure 2 kPa exp(-1e4 (t - 0.05)^2) (his Eq. 7.11); the wall's ends
fixed axially and free radially (his Eq. 5.12), at the inlet and at the
interface. BDF2 in the fluid and the wall, dt = 0.1 ms, t <= 0.1 s.

The 1D vessel is his (Eq. 2.3): P = beta (sqrt(A) - sqrt(A0))/A0 with
beta = sqrt(pi) h E/(1 - nu^2), so c0^2 = h E/(2 rho R (1 - nu^2)),
c0 = 6.02 m/s; the sqrtElastic law has c0^2 = 2 E1D h/(3 rho R), so
E1D = 3 E/(4 (1 - nu^2)). Friction KR = 2 pi (gamma + 2) nu with his blunt
profile gamma = 9, but alpha = 1 (his alpha = 1.1 is not available here).

The 3D cases are quarter tubes (two symmetry planes): the coupled port has
flowScale 4. Sections of the fluid every centimetre are faceZones (topoSet)
on which surfaceFieldValue writes the mean pressure and the flow, the flow
through the moving faces (phi + meshPhi) and the pressure as the mean of
the two cells on every face, processor faces included (see
SECTION_FIELDS); the 1D vessels carry probes at the same z.
"""
import argparse
import math
import re
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = HERE/"base"
RUNS = HERE/"runs"

RHO = 1060.0
MU = 4e-3
R = 1e-2
H = 1e-3
E = 7e5
NU = 0.3
A = math.pi*R**2
L3D = 0.10                      # the 3D part of the coupled tube
LTUBE = 0.30                    # Businaro's tube
LFULL = 0.45                    # the reference tube
CELLS_PER_M = 1000              # 1 mm axial cells, as Businaro
P_MAX = 2000.0                  # Pa
A_GAUSS = 1e4                   # 1/s2
T0 = 0.05                       # s
RT = 2.0494e7                   # Pa s/m3, Table F.1
GAMMA = 9                       # blunt profile
END_TIME = 0.1
DELTA_T = 1e-4
DX_1D = 1e-3
Z_DISP = 0.05                   # wall displacement monitored at z = 5 cm


def c0_1d():
    """Businaro's 1D wave speed at the reference area"""
    return math.sqrt(H*E/(2*RHO*R*(1 - NU**2)))


def youngs1D():
    """E of the sqrtElastic law with the same beta as Businaro's law"""
    return 0.75*E/(1 - NU**2)


def planes(length):
    """z [cm] of the sections inside a 3D tube of that length"""
    return list(range(1, int(round(length*100))))


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
    return f"z{zcm:g}"


# The fields the sections are measured on, registered by this object ahead
# of the others. The flow: the solver leaves phi relative to the moving
# mesh, and the faceZones move axially with the wall (its axial waves run
# ahead of the pulse), so their sum(phi) is not the flow past the section;
# phiAbsolute = phi + meshPhi is (meshPhi itself is not registered). The
# pressure: surfaceFieldValue takes the mean of the two cells on an interior
# face, but on a processor face the stored value, which is the neighbour
# cell's, so a section on a processor boundary would be read half a cell
# away; pSection is the mean of the two cells on every face. The outlet's
# points are fixed and it is a boundary, so there phi and p are measured
# directly.
SECTION_FIELDS = """    sectionFields
    {
        type            coded;
        libs            (utilityFunctionObjects);
        name            sectionFields;
        region          fluid;
        codeExecute
        #{
            const fvMesh& m = mesh();
            const surfaceScalarField& phi =
                m.lookupObject<surfaceScalarField>("phi");
            surfaceScalarField* ptr =
                m.getObjectPtr<surfaceScalarField>("phiAbsolute");
            if (!ptr)
            {
                ptr = new surfaceScalarField
                (
                    IOobject
                    (
                        "phiAbsolute",
                        m.time().timeName(),
                        m,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        IOobject::REGISTER
                    ),
                    phi
                );
                ptr->store();
            }
            *ptr = phi;
            if (m.moving())
            {
                *ptr += m.phi();
            }

            const volScalarField& p = m.lookupObject<volScalarField>("p");
            surfaceScalarField* ps =
                m.getObjectPtr<surfaceScalarField>("pSection");
            if (!ps)
            {
                ps = new surfaceScalarField
                (
                    IOobject
                    (
                        "pSection",
                        m.time().timeName(),
                        m,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        IOobject::REGISTER
                    ),
                    m,
                    dimensionedScalar(p.dimensions(), Zero)
                );
                ps->store();
            }
            const labelUList& own = m.owner();
            const labelUList& nei = m.neighbour();
            scalarField& psi = ps->primitiveFieldRef();
            forAll(nei, facei)
            {
                psi[facei] = 0.5*(p[own[facei]] + p[nei[facei]]);
            }
            auto& psb = ps->boundaryFieldRef();
            forAll(psb, patchi)
            {
                const fvPatchScalarField& pp = p.boundaryField()[patchi];
                if (pp.coupled())
                {
                    const scalarField mean
                    (
                        0.5*(pp.patchInternalField() + pp.patchNeighbourField())
                    );
                    psb[patchi] = mean;
                }
                else
                {
                    psb[patchi] = pp;
                }
            }
        #};
    }
"""


def control_dict(coupled, zs, dt):
    libs = 'libs            ("libcirculation1DCoupling.so");\n\n' if coupled else ""
    # The 3D side of the outlet (the interface of the coupled runs)
    fos = SECTION_FIELDS
    for suffix, op, field in (("P", "areaAverage", "p"), ("Q", "sum", "phi")):
        fos += f"""    outlet_{suffix}
    {{
        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");
        region          fluid;
        regionType      patch;
        name            outlet;
        operation       {op};
        fields          ({field});
        writeControl    timeStep;
        writeInterval   1;
        writeFields     false;
        log             false;
    }}
"""
    for zcm in zs:
        n = zone_name(zcm)
        for suffix, op, field in (("P", "areaAverage", "pSection"), ("Q", "sum", "phiAbsolute")):
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
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {END_TIME:g};
deltaT          {dt:g};
writeControl    runTime;
writeInterval   0.01;
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
    // Radial displacement of the wall at z = {Z_DISP*100:g} cm
    displacement
    {{
        type            solidPointDisplacement;
        region          solid;
        point           (0 {R:g} {Z_DISP:g});
    }}
{fos}}}
"""
    return HEADER % ("system", "controlDict") + body + TAIL


def topo_set_dict(zs):
    body = "actions\n(\n"
    for zcm in zs:
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


FLUID_1D = f"""fluid
{{
    rho             {RHO};
    mu              {MU};
    profileOrder    {GAMMA};
    alpha           1;
}}
numerics
{{
    CFL             0.5;
    flux            HLL;
    limiter         MC;
    timeScheme      SSPRK3;
}}
"""


def vessel(name, length):
    return f"""    {name}
    {{
        length              {length:g};
        nCells              {int(round(length/DX_1D))};
        radius              {R};
        wallThickness       {H};
        youngsModulus       {youngs1D():.6e};
        diastolicPressure   0;
        externalPressure    0;
        initialPressure     0;
    }}
"""


def properties_coupled(scheme, order):
    probes = "".join(f"    {zone_name(z)} {{ vessel distal; x {z*1e-2 - L3D:.6g}; }}\n"
                     for z in range(int(round(L3D*100)) + 1, int(round(LTUBE*100)) + 1))
    body = f"""// Businaro (2026): the 20 cm 1D vessel continuing the 10 cm FSI tube,
// his tube law (c0 = {c0_1d():.4g} m/s) and friction (gamma = {GAMMA})
{FLUID_1D}vessels
{{
{vessel("distal", LTUBE - L3D)}}}
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
probes                                  // x along the vessel = z - {L3D*100:g} cm
{{
{probes}}}
coupling
{{
    scheme      {scheme};
    order       {order};
    seed        zero;
}}
"""
    return HEADER % ("constant/fluid", "circulation1DProperties") + body + TAIL


def properties_rt():
    body = f"""// Businaro's 3D reference outlet, P = Rt Q (his Eq. 5.8): a 1 cm 1D stub
// with his tube law ending in the resistance Rt = {RT:g} Pa s/m3
{FLUID_1D}vessels
{{
{vessel("stub", 0.01)}}}
nodes
{{
    from3D
    {{
        type        coupledFlow;
        connection  stub.in;
        patch       outlet;
        flowScale   4;
    }}
    outlet
    {{
        type        resistance;
        connection  stub.out;
        R           {RT:g};
        Pout        0;
    }}
}}
coupling
{{
    scheme      implicit;
    order       2;
    seed        zero;
}}
"""
    return HEADER % ("constant/fluid", "circulation1DProperties") + body + TAIL


def properties_only1D():
    probes = "".join(f"    {zone_name(z)} {{ vessel tube; x {z*1e-2:.6g}; }}\n"
                     for z in range(1, int(round(LTUBE*100)) + 1))
    body = f"""// Businaro (2026): the 30 cm tube as a 1D vessel alone, c0 = {c0_1d():.4g} m/s
{FLUID_1D}vessels
{{
{vessel("tube", LTUBE)}}}
nodes
{{
    inlet
    {{
        type        pressureInlet;
        connection  tube.in;
        pressure
        {{
            type        expression;
            expression  #{{ {P_MAX:g}*exp(-{A_GAUSS:g}*sqr(arg() - {T0:g})) #}};
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
endTime         {END_TIME:g};
deltaT          1e-5;
writeControl    runTime;
writeInterval   0.01;
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


def make_3d(name, length, outlet, dt, force, linear_elastic=False):
    """outlet: "fixed" (p = 0), "rt" (the resistance through a 1D stub) or a
    coupling (scheme, order) for the 20 cm 1D vessel"""
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
             [(r"\nL\s+[\d.e-]+;", f"\nL       {length:g};"),
              (r"\nnz\s+\d+;", f"\nnz      {n};")])
    zs = planes(length)
    coupled = outlet != "fixed"
    (case/"system"/"controlDict").write_text(control_dict(coupled, zs, dt))
    (case/"system"/"fluid"/"topoSetDict").write_text(topo_set_dict(zs))
    if coupled:
        edit(case/"0"/"fluid"/"p",
             [(r"outlet\s*\{\s*type\s+fixedValue;\s*value\s+uniform 0;\s*\}",
               "outlet\n    {\n        type        circulation1DPressure;\n"
               "        value       uniform 0;\n    }")])
        (case/"constant"/"fluid"/"circulation1DProperties").write_text(
            properties_rt() if outlet == "rt" else properties_coupled(*outlet))
    if linear_elastic:
        edit(case/"constant"/"solid"/"mechanicalProperties",
             [(r"type\s+StVenantKirchhoffElastic;", "type            linearElastic;")])
        edit(case/"constant"/"solid"/"solidProperties",
             [(r"\nsolidModel\s+nonLinearGeometryTotalLagrangianTotalDisplacement;",
               "\nsolidModel linearGeometryTotalDisplacement;")])
    print(f"{name}: {length*100:g} cm, {n} axial cells, dt = {dt:g} s, outlet "
          + (outlet if isinstance(outlet, str) else f"{outlet[0]} coupling of order {outlet[1]}")
          + (", linear elastic" if linear_elastic else ""))


def make_1d(name, force):
    case = RUNS/name
    if case.exists() and finished(case) and not force:
        print(f"{name}: finished, kept")
        return
    if case.exists():
        shutil.rmtree(case)
    (case/"system").mkdir(parents=True)
    (case/"constant").mkdir()
    (case/"system"/"controlDict").write_text(control_dict_1d())
    (case/"constant"/"circulation1DProperties").write_text(properties_only1D())
    print(f"{name}: {LTUBE*100:g} cm 1D vessel, c0 = {c0_1d():.4g} m/s")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dt", type=float, default=DELTA_T,
                    help=f"3D time step [s] (default {DELTA_T:g}, Businaro's)")
    ap.add_argument("--linearElastic", action="store_true",
                    help="the wall linear elastic (small strain) instead of Saint Venant-Kirchhoff")
    ap.add_argument("--withRt", action="store_true",
                    help="also full3D_Rt, Businaro's 30 cm reference with the resistance outlet")
    ap.add_argument("--suffix", default="",
                    help="appended to every run name (e.g. _smoke)")
    ap.add_argument("--force", action="store_true", help="rewrite finished runs too")
    args = ap.parse_args()
    s = args.suffix + ("_le" if args.linearElastic else "")
    print(f"1D: c0 = {c0_1d():.4g} m/s, E1D = {youngs1D():.4g} Pa, "
          f"rho c0/A = {RHO*c0_1d()/A:.4g} Pa s/m3 (Rt = {RT:g})")
    RUNS.mkdir(exist_ok=True)
    make_3d("full3D" + s, LFULL, "fixed", args.dt, args.force, args.linearElastic)
    if args.withRt:
        make_3d("full3D_Rt" + s, LTUBE, "rt", args.dt, args.force, args.linearElastic)
    make_3d("coupled" + s, L3D, ("implicit", 2), args.dt, args.force, args.linearElastic)
    make_3d("coupled_explicit" + s, L3D, ("explicit", 1), args.dt, args.force,
            args.linearElastic)
    make_1d("only1D" + args.suffix, args.force)


if __name__ == "__main__":
    main()
