#!/usr/bin/env python3
"""Write the runs of the Womersley segment: a rigid no-slip pipe between two
1D vessels, with a single-frequency incident wave, at three Womersley numbers.

    python3 makeCases.py [--alphas 3,6,12] [--dt 1e-3] [--force]

The pipe (radius R = 1 cm, length 10 cm, an axisymmetric 5-degree wedge in
pimpleFoam) joins two 1 m vessels of the rigid-segment benchmark (c0 = 6.17
m/s, inviscid). The proximal vessel's inlet is a waveInlet: a semi-infinite
line sending a wave p(t) = P sin(w t), ramped over the first period, and
absorbing whatever the segment reflects; the distal end is absorbing. At
f = 1 Hz the viscosity sets the Womersley number alpha = R sqrt(w/nu):
nu = w R^2/alpha^2. The exact periodic solution is the transmission line
with the segment as the series impedance of Womersley (1955),

    Z(w) = i w rho l / (A (1 - F10)),  F10 = 2 J1(L)/(L J0(L)),  L = i^(3/2) alpha

(the inertance i w rho l/A of the inviscid benchmark as alpha -> infinity),
so the reflection and transmission of the incident wave are
Gamma = Z/(Z + 2 Zc) and T = 2 Zc/(Z + 2 Zc); analyse.py compares them with
the fitted periodic state of the ports, and the velocity profile at the
pipe's mid-length with Womersley's.

Each run lasts long enough for the viscous transient to decay: three times
R^2/(j01^2 nu) (the slowest radial mode) plus three periods.
"""
import argparse
import math
import re
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = HERE/"base"
RUNS = HERE/"runs"

RHO = 1050.0
R = 0.01
A = math.pi*R**2
L3D = 0.1
L1D = 1.0
WEDGE_DEG = 5.0
FLOW_SCALE = 360.0/WEDGE_DEG
F = 1.0                                 # Hz
OMEGA = 2*math.pi*F
P_INC = 50.0                            # Pa, incident wave amplitude
H_WALL = 1.5e-3                         # the rigid-segment vessels
E_WALL = 400e3
C0 = math.sqrt(2*E_WALL*H_WALL/(3*RHO*R))
ALPHAS = [3, 6, 12]
J01 = 2.404825557695773                 # first zero of J0


def nu(alpha):
    return OMEGA*R**2/alpha**2


def end_time(alpha):
    transient = 3*R**2/(J01**2*nu(alpha))
    return math.ceil(transient/F + 3)/F


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


def control_dict(alpha, dt, end=None):
    T = end if end is not None else end_time(alpha)
    body = f"""application     pimpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {T:g};
deltaT          {dt:g};
writeControl    runTime;
writeInterval   {T:g};
purgeWrite      0;
writeFormat     ascii;
writePrecision  10;
writeCompression off;
timeFormat      general;
timePrecision   8;
runTimeModifiable yes;
adjustTimeStep  no;

libs            ("libcirculation1DCoupling.so");

functions
{{
    // The axial velocity profile at the pipe's mid-length, eight times
    // over the last period
    profile
    {{
        type            sets;
        libs            ("libsampling.so");
        writeControl    runTime;
        writeInterval   {1/(8*F):g};
        timeStart       {T - 1/F:g};
        interpolationScheme cellPoint;
        setFormat       raw;
        fields          (U);
        sets
        {{
            midLength
            {{
                type    uniform;
                axis    x;
                start   (0 0 {L3D/2:g});
                end     ({R*0.9995:g} 0 {L3D/2:g});
                nPoints 200;
            }}
        }}
    }}
}}
"""
    return HEADER % ("system", "controlDict") + body + TAIL


def properties(pInc=P_INC):
    body = f"""// The 1D side of the Womersley segment: the two vessels of the
// rigid-segment benchmark, joined through the viscous pipe. The proximal
// inlet sends the incident wave in and absorbs what comes back.
fluid
{{
    rho             {RHO};
    mu              0;
    profileOrder    9;
    alpha           1;
}}
numerics
{{
    CFL             0.5;
    flux            HLL;
    limiter         none;
    timeScheme      SSPRK3;
}}
vesselDefaults
{{
    length              {L1D};
    nCells              400;            // dx = 2.5 mm; the wavelength is 6 m
    radius              {R};
    wallThickness       {H_WALL};
    youngsModulus       {E_WALL:g};
    diastolicPressure   0;
    externalPressure    0;
    initialPressure     0;
}}
vessels
{{
    proximal {{}}
    distal   {{}}
}}
nodes
{{
    inlet                               // incident wave, non-reflecting
    {{
        type        waveInlet;
        connection  proximal.in;
        pressure
        {{
            type        expression;
            expression  #{{ {pInc}*sin({OMEGA:.15g}*arg())*min(1, {F}*arg()) #}};
        }}
    }}
    to3D                                // feeds the pipe's inlet (a wedge)
    {{
        type        coupledFlow;
        connection  proximal.out;
        patch       inlet;
        flowScale   {FLOW_SCALE:g};
    }}
    from3D
    {{
        type        coupledFlow;
        connection  distal.in;
        patch       outlet;
        flowScale   {FLOW_SCALE:g};
    }}
    outlet
    {{
        type        absorbing;
        connection  distal.out;
    }}
}}
probes
{{
    proximalMid {{ vessel proximal; x midpoint; }}
    distalMid   {{ vessel distal;   x midpoint; }}
}}
coupling
{{
    scheme      implicit;
    order       2;
    seed        zero;
}}
"""
    return HEADER % ("constant", "circulation1DProperties") + body + TAIL


def edit(path, pairs):
    text = path.read_text()
    for old, new in pairs:
        text, n = re.subn(old, new, text)
        if n == 0:
            raise SystemExit(f"{path}: nothing matched {old!r}")
    path.write_text(text)


def finished(case):
    log = case/"log.pimpleFoam"
    return log.is_file() and "\nEnd\n" in log.read_text()[-2000:]


def run_name(alpha):
    return f"alpha{alpha:g}"


def make(alpha, dt, force, name=None, end=None, pInc=P_INC):
    """name and end override the run's name and end time (the order study);
    pInc the incident amplitude (the physiological-amplitude variant)"""
    case = RUNS/(name or run_name(alpha))
    if case.exists() and finished(case) and not force:
        print(f"{case.name}: finished, kept")
        return case
    if case.exists():
        shutil.rmtree(case)
    shutil.copytree(BASE, case)
    edit(case/"constant"/"transportProperties",
         [(r"nu\s+\[0 2 -1 0 0 0 0\]\s+[\d.e+-]+;",
           f"nu              [0 2 -1 0 0 0 0] {nu(alpha):.6g};")])
    (case/"system"/"controlDict").write_text(control_dict(alpha, dt, end))
    (case/"constant"/"circulation1DProperties").write_text(properties(pInc))
    print(f"{case.name}: nu = {nu(alpha):.4g} m2/s, "
          f"{end if end is not None else end_time(alpha):g} s, dt = {dt:g} s")
    return case


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--alphas", default=",".join(str(a) for a in ALPHAS))
    ap.add_argument("--dt", type=float, default=1e-3)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    print(f"c0 = {C0:.4g} m/s, Zc = {RHO*C0/A:.4g} Pa s/m3, incident {P_INC} Pa "
          f"at {F} Hz")
    RUNS.mkdir(exist_ok=True)
    for a in (float(x) for x in args.alphas.split(",")):
        make(a, args.dt, args.force)


if __name__ == "__main__":
    main()
