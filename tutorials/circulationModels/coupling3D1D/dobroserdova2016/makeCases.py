#!/usr/bin/env python3
"""Set up the hard- and soft-coupling tests of Dobroserdova, Olshanskii and
Simakov (2016), sections 4.2 and 4.3, in SI.

    python3 makeCases.py            writes the cases into runs/ (a case that
                                    has run is left alone; --force rewrites all)

Their set-up (cgs): a 1D vessel of 5 cm, radius 1 cm, joined to a rigid 3D
cylinder of 5 cm and radius 1 cm whose far end is open (zero traction);
rho = 1 g/cm3, nu = 0.04 cm2/s; the inflow Q = 0.1 sin(2 pi t/0.3) cm3/s
(section 4.2) or 0.1 for t <= 0.1 s, then 0 (section 4.3); c0 = 350, 700 and
1050 cm/s. The vessel here is inviscid (theirs has Poiseuille friction, a
damping time of about 3 s), so that the exact lossless solution applies. Hard coupling: continuity of pressure and flow. Soft coupling: an
elastic sphere (a compliance C) between two resistances at the interface,
R1D0D + R0D3D = rho c0/S (their Eq. 13) and C = l S/(rho c0^2), l = 0.1 cm
(Eq. 14). Their reference is a 10 cm 1D vessel with a non-reflecting outlet.

The rigid 3D pipe with an open end is, exactly, an inductor L = rho l/S to
ground (plug flow; with slip walls there is no viscous term), so the 1D-only
cases below are the exact linear solutions of their multiscale problems:

    reference   a 10 cm vessel, absorbing outlet
    hard        the 5 cm vessel ending in an open inertance node (Pout = 0)
    soft        the same with (port R1 C R2 L ground), R1 = R2 = Zc/2

and hard3D is the 5 cm vessel coupled to a pimpleFoam pipe of 5 cm with
slip walls and p = 0 at its outlet (the 3D run of the hard coupling).
"""
import math
import re
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = HERE/"runs"

RHO = 1000.0
MU = 0.0                        # their nu = 0.04 cm2/s; inviscid here, so that
                                # the exact (lossless-line) solution applies
R = 0.01
A = math.pi*R**2
L1D = 0.05
L3D = 0.05
L0 = 1e-3                       # l in Eq. 14: 0.1 cm
H = 1e-3                        # wall thickness of the sqrtElastic law
C0S = {"c350": 3.5, "c700": 7.0, "c1050": 10.5}
INFLOWS = {
    "sine": ("expression", "#{ 1e-7*sin(20.943951023931955*arg()) #}", 1.2, 1e-3),
    "step": ("expression", "#{ (arg() <= 0.1) ? 1e-7 : 0 #}", 0.5, 5e-4),
}
LIND = RHO*L3D/A                # inductance of the rigid pipe [kg/m4]


def youngs(c0):
    """E of the sqrtElastic law with wave speed c0 at the reference area:
    c0^2 = beta/(2 rho sqrt(Ad)), beta = 4/3 sqrt(pi) E h"""
    beta = 2*RHO*c0**2*math.sqrt(A)
    return beta/(4.0/3.0*math.sqrt(math.pi)*H)


HEADER = """/*--------------------------------*- C++ -*----------------------------------*\\
| solids4foam: solid mechanics and fluid-solid interaction simulations        |
| Version:     v2.4                                                           |
| Web:         https://solids4foam.github.io                                  |
| Disclaimer:  This offering is not approved or endorsed by OpenCFD Limited,  |
|              producer and distributor of the OpenFOAM software via          |
|              www.openfoam.com, and owner of the OPENFOAM® and OpenCFD®      |
|              trade marks.                                                   |
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


def properties(config, c0, inflow):
    kind, expr, _, _ = INFLOWS[inflow]
    Zc = RHO*c0/A
    length = 2*L1D if config == "reference" else L1D
    vessel = f"""vessels
{{
    vessel
    {{
        length              {length};
        nCells              {int(round(length/5e-4))};
        radius              {R};
        wallThickness       {H};
        youngsModulus       {youngs(c0):.6e};
        diastolicPressure   0;
        externalPressure    0;
        initialPressure     0;
    }}
}}
"""
    if config == "reference":
        outlet = """    outlet { type absorbing; connection vessel.out; }
"""
    elif config == "hard":
        outlet = f"""    // The rigid 3D pipe with an open end, exactly: an inertance to Pout = 0
    outlet
    {{
        type        inertance;
        connection  vessel.out;
        length      {L3D};
        radius      {R};
        Pout        0;
    }}
"""
    elif config == "soft":
        Cb = L0*A/(RHO*c0**2)
        outlet = f"""    // Their absorbing 0D model (Eqs. 7-10, I = R0 = 0) between the vessel
    // and the rigid pipe: R1 + R2 = Zc (Eq. 13), C = l S/(rho c0^2) (Eq. 14)
    outlet
    {{
        type        lumped0D;
        connection  vessel.out;
        circuit
        {{
            units       SI;
            components
            {{
                port    {{ type externalPort; direction outlet; }}
                R1      {{ type resistor;   R {0.5*Zc:.6e}; }}
                C       {{ type compliance; C {Cb:.6e}; initialPressure 0; }}
                R2      {{ type resistor;   R {0.5*Zc:.6e}; }}
                L       {{ type inductor;   L {LIND:.6e}; }}
                ground  {{ type ground; }}
            }}
            connections ( (port R1 C R2 L ground) );
        }}
    }}
"""
    else:  # hard3D
        outlet = """    to3D
    {
        type        coupledFlow;
        connection  vessel.out;
        patch       inlet;
    }
"""
    coupling = """
coupling
{
    scheme      implicit;
    order       2;
    seed        zero;
}
""" if config == "hard3D" else ""
    body = f"""// Dobroserdova, Olshanskii & Simakov (2016), section 4.{2 if inflow == 'sine' else 3}:
// the {inflow} inflow, c0 = {c0} m/s, {config} coupling (see makeCases.py).

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
    limiter         {"MC" if inflow == "step" else "none"};{"             // TVD: no overshoot at the jump" if inflow == "step" else ""}
    timeScheme      SSPRK3;
}}

{vessel}
nodes
{{
    inlet
    {{
        type        flowInlet;
        connection  vessel.in;
        flowRate
        {{
            type        {kind};
            expression  {expr};
        }}
    }}

{outlet}}}

probes
{{
    interface { "{" } vessel vessel; x out; { "}" }
}}
{coupling}"""
    return HEADER % ("constant", "circulation1DProperties") + body + TAIL


def control_dict(config, inflow):
    _, _, end, dt = INFLOWS[inflow]
    if config == "hard3D":
        return None
    body = f"""application     circulation1DFoam;

startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {end};
deltaT          {dt};          // output interval; the solver sub-steps at its CFL step

writeControl    runTime;
writeInterval   {end};

writePrecision  12;
timePrecision   6;
"""
    return HEADER % ("system", "controlDict") + body + TAIL


def finished(case):
    log = case/"log.pimpleFoam" if "3D" in case.name else case/"log"
    return log.is_file() and log.read_text().rstrip().endswith("End")


def make_1d(config, c0name, inflow):
    case = RUNS/f"{inflow}_{c0name}_{config}"
    if finished(case) and "--force" not in sys.argv:
        return case
    if case.exists():
        shutil.rmtree(case)
    (case/"constant").mkdir(parents=True)
    (case/"system").mkdir()
    (case/"constant"/"circulation1DProperties").write_text(
        properties(config, C0S[c0name], inflow))
    (case/"system"/"controlDict").write_text(control_dict(config, inflow))
    return case


def edit(path, pairs):
    s = path.read_text()
    for pattern, value in pairs:
        s, n = re.subn(pattern, value, s, count=1, flags=re.M)
        assert n == 1, (path, pattern)
    path.write_text(s)


def make_3d(c0name, inflow):
    """The 3D run of the hard coupling: the rigidSegment3D base case with a
    5 cm pipe, one 1D vessel at its inlet and p = 0 at its outlet"""
    case = RUNS/f"{inflow}_{c0name}_hard3D"
    if finished(case) and "--force" not in sys.argv:
        return case
    if case.exists():
        shutil.rmtree(case)
    base = HERE.parent/"rigidSegment3D"/"base"
    shutil.copytree(base, case)
    (case/"constant"/"circulation1DProperties").write_text(
        properties("hard3D", C0S[c0name], inflow))
    _, _, end, dt = INFLOWS[inflow]
    edit(case/"system"/"controlDict", [
        (r"^endTime\s+\S+;", f"endTime         {end};"),
        (r"^deltaT\s+\S+;.*$", f"deltaT          {dt};"),
    ])
    # pipe length 8 cm -> 5 cm: the z coordinate of the far vertices and arcs,
    # and 10 cells along the pipe instead of 16
    s = (case/"system"/"blockMeshDict").read_text()
    assert s.count(" 0.08)") == 12 and s.count(" 16)") == 5, "blockMeshDict changed"
    s = s.replace(" 0.08)", f" {L3D})").replace(" 16)", " 10)")
    s = s.replace("length 8 cm", "length 5 cm").replace("*16)", "*10)").replace("*16 ", "*10 ")
    (case/"system"/"blockMeshDict").write_text(s)
    # outlet: zero traction, p = 0
    edit(case/"0"/"p", [
        (r"outlet\n    \{\n        type            circulation1DPressure;\n        value           uniform 0;\n    \}",
         "outlet\n    {\n        type            fixedValue;\n        value           uniform 0;\n    }")])
    return case


def main():
    RUNS.mkdir(exist_ok=True)
    names = []
    for inflow in INFLOWS:
        for c0name in C0S:
            for config in ("reference", "hard", "soft"):
                names.append(make_1d(config, c0name, inflow).name)
    for inflow in INFLOWS:
        names.append(make_3d("c350", inflow).name)
    print("\n".join(names))


if __name__ == "__main__":
    main()
