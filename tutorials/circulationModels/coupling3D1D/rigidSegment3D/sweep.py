#!/usr/bin/env python3
"""Set up and run the cases of the rigid-segment 3D-1D study (Allrun-sweep).

    python3 sweep.py [NPROC] [--ddt LIST]   the study for the CFD time schemes
                                            in LIST (default backward,BDF3,BDF4)
    python3 sweep.py --one                  the monolithic run and implicit
                                            order 2 at dt = 0.5 ms (Allrun)

Each case is a copy of base/ in runs/<name> that differs only in the coupling
scheme and order (constant/circulation1DProperties), deltaT and endTime
(system/controlDict), and for BDF3 and BDF4 the ddt scheme (system/fvSchemes,
"BDF 3", from libsolids4FoamModels) and the CFL number of the vessels (0.1, so
that their time error, 5e-8 of the peak, stays below a third- or fourth-order
coupling's). Cases are serial; NPROC (default 6) run at a time. For each time
scheme, with <p> its prefix ("" for backward, "BDF3_", "BDF4_"):

    <p>acc_<scheme><order>_dt<dt>  explicit 1-4 and implicit orders (STUDY),
                                   dt = 4 ... 0.25 ms
    <p>fine                        implicit, the scheme's order, at dt = 0.25/8
                                   ms: the reference
    <p>stab_<order>_<x>            explicit orders at dt = x * the predicted
                                   limit, x = 0.9 and 1.1, to t = 0.5 s
    monolithic, monolithic_cfl0.1  the 1D network with the inertance node of
                                   the same pipe (area A3D) in place of the
                                   ports: circulation1DFoam, no 3D

Already finished cases (a log.pimpleFoam ending in "End") are skipped.
"""
import math
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = HERE/"runs"

RHO = 1050.0
A0 = math.pi*0.01**2
BETA = 4.0/3.0*math.sqrt(math.pi)*400e3*1.5e-3
C0 = math.sqrt(BETA/(2*RHO*A0))*A0**0.25
Z0 = RHO*C0/A0
A3D = 0.5*32*0.01**2*math.sin(2*math.pi/32)     # the 32-sided patch polygon
TAU = RHO*0.08/A3D/(2*Z0)

DTS = [4e-3, 2e-3, 1e-3, 5e-4, 2.5e-4]
FINE_DT = DTS[-1]/8

# Per CFD time scheme: its order, the fvSchemes entry, the vessels' CFL, the
# case prefix, the couplings of the accuracy study, the explicit orders whose
# stability is checked, and the predicted explicit limits of orders 1-4 in
# dt/tau (the model problem: the pipe's BDF formula, the flow extrapolated to
# t1 by order p, each vessel end a resistor Zc)
EXPLICIT = [("explicit", o) for o in (1, 2, 3, 4)]
STUDY = {
    "backward": dict(q=2, ddt="backward", cfl=0.5, prefix="",
                     accuracy=EXPLICIT + [("implicit", 2)], stab=(2, 3, 4),
                     limits=(4.0, 4/3, 4/7, 4/15)),
    "BDF3": dict(q=3, ddt="BDF 3", cfl=0.1, prefix="BDF3_",
                 accuracy=EXPLICIT + [("implicit", 1), ("implicit", 3)],
                 stab=(3, 4), limits=(4.5, 20/9, 20/21, 4/9)),
    "BDF4": dict(q=4, ddt="BDF 4", cfl=0.1, prefix="BDF4_",
                 accuracy=EXPLICIT + [("implicit", 1), ("implicit", 4)],
                 stab=(1, 3, 4), limits=(1.6428, 9/4, 32/21, 32/45)),
}
# The study of the backward scheme, used by analyse.py
ACCURACY = STUDY["backward"]["accuracy"]


def limit(order, ddt="backward"):
    """Predicted explicit stability limit in dt/tau"""
    return STUDY[ddt]["limits"][order - 1]


def edit(path, pairs):
    s = path.read_text()
    for pattern, value in pairs:
        s, n = re.subn(pattern, value, s, count=1, flags=re.M)
        if n != 1:
            sys.exit(f"no match for {pattern} in {path}")
    path.write_text(s)


def setup(name, scheme, order, dt, end=0.35, study="backward"):
    st = STUDY[study]
    case = RUNS/(st["prefix"] + name)
    log = case/"log.pimpleFoam"
    if log.is_file() and log.read_text().rstrip().endswith("End"):
        return None
    if case.exists():
        shutil.rmtree(case)
    shutil.copytree(HERE/"base", case)
    edit(case/"system"/"controlDict", [
        (r"^deltaT\s+\S+;", f"deltaT          {dt!r};"),
        (r"^endTime\s+\S+;", f"endTime         {end!r};"),
    ])
    edit(case/"constant"/"circulation1DProperties", [
        (r"^(\s+scheme\s+)\S+;", rf"\g<1>{scheme};"),
        (r"^(\s+order\s+)\S+;", rf"\g<1>{order};"),
        (r"^(\s+CFL\s+)\S+;", rf"\g<1>{st['cfl']};"),
    ])
    if st["ddt"] != "backward":
        edit(case/"system"/"fvSchemes", [
            (r"^(\s+default\s+)backward;.*$", rf"\g<1>{st['ddt']};")])
        edit(case/"system"/"controlDict", [
            (r'^libs\s+\("libcirculation1DCoupling.so"\);',
             'libs            ("libcirculation1DCoupling.so" '
             '"libsolids4FoamModels.so");     // BDF ddt scheme')])
    return case


def run3D(case):
    with open(case/"log.blockMesh", "w") as f:
        subprocess.run(["blockMesh"], cwd=case, stdout=f,
                       stderr=subprocess.STDOUT, check=True)
    with open(case/"log.pimpleFoam", "w") as f:
        r = subprocess.run(["pimpleFoam"], cwd=case, stdout=f,
                           stderr=subprocess.STDOUT)
    print(f"{case.name}: {'done' if r.returncode == 0 else 'stopped'}",
          flush=True)


def monolithic(cfl=0.5):
    case = RUNS/("monolithic" if cfl == 0.5 else f"monolithic_cfl{cfl:g}")
    if (case/"postProcessing").is_dir():
        return
    if case.exists():
        shutil.rmtree(case)
    case.mkdir(parents=True)
    shutil.copytree(HERE/"base"/"constant", case/"constant")
    (case/"system").mkdir()
    s = (HERE/"base"/"system"/"controlDict").read_text()
    s = s.split("libs")[0] + """application     circulation1DFoam;

startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         0.35;
deltaT          3.125e-5;       // every step of the finest 3D run

writeControl    runTime;
writeInterval   1;

writePrecision  17;
timePrecision   12;

// ************************************************************************* //
"""
    (case/"system"/"controlDict").write_text(s)
    p = case/"constant"/"circulation1DProperties"
    edit(p, [(r"^(\s+CFL\s+)\S+;", rf"\g<1>{cfl};")])
    s = p.read_text()
    ports = s[s.index("    to3D"):s.index("    outlet\n")]
    s = s.replace(ports, f"""    rigid                               // the 3D pipe as plug flow
    {{
        type        inertance;
        connections (proximal.out distal.in);
        length      0.08;
        area        {A3D!r};
        pressure    static;
    }}

""")
    p.write_text(s)
    with open(case/"log.circulation1DFoam", "w") as f:
        subprocess.run(["circulation1DFoam"], cwd=case, stdout=f,
                       stderr=subprocess.STDOUT, check=True)
    print(f"{case.name}: done", flush=True)


def main():
    RUNS.mkdir(exist_ok=True)
    monolithic()

    if sys.argv[1:] == ["--one"]:
        case = setup("acc_implicit2_dt0.0005", "implicit", 2, 5e-4)
        if case is not None:
            run3D(case)
        return

    args = sys.argv[1:]
    studies = ["backward", "BDF3", "BDF4"]
    if "--ddt" in args:
        i = args.index("--ddt")
        studies = args[i + 1].split(",")
        del args[i:i + 2]
    nproc = int(args[0]) if args else 6

    cases = []
    for study in studies:
        st = STUDY[study]
        monolithic(st["cfl"])
        cases.append(setup("fine", "implicit", st["q"], FINE_DT, study=study))
        for s, o in st["accuracy"]:
            for dt in DTS:
                cases.append(setup(f"acc_{s}{o}_dt{dt:g}", s, o, dt,
                                   study=study))
        for o in st["stab"]:
            for x in (0.9, 1.1):
                dt = round(x*limit(o, study)*TAU, 10)
                cases.append(setup(f"stab_{o}_{x:g}", "explicit", o, dt,
                                   end=0.5, study=study))

    cases = [c for c in cases if c is not None]
    print(f"tau = {TAU*1e3:.4f} ms; running {len(cases)} cases, {nproc} at a "
          "time", flush=True)
    with ThreadPoolExecutor(nproc) as pool:
        list(pool.map(run3D, cases))


if __name__ == "__main__":
    main()
