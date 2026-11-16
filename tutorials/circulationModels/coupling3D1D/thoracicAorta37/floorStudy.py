#!/usr/bin/env python3
"""Bisecting the dt-independent floor of the BDF 3/4 time-order study.

    python3 floorStudy.py [--variants a,b,...] [--jobs 3] [--analyse]

The BDF4 runs of ../orderStudy.py (0-0.3 s from the smooth start, dt = 1,
0.5 and 0.25 ms) are repeated with one change each, and the self-convergence
of the port flow and pressure is measured again. The floor shows as a
difference between successive steps that stops falling at about 3e-6 of
the signal; a variant that removes it converges at order 4 (the
uniform-vessel cases reach 3.8-3.9 on the same triple).

    ddtCorrOff     pimpleFoam's flux correction fvc::ddtCorr(U, phi) off
    tightSolvers   linear solvers to 1e-15 (maxIter 20000), PIMPLE p residual
                   1e-12 (from 1e-12 and 1e-10)
    cylinder       the 3D cone replaced by a cylinder of its mean radius
                   sqrt(r0 r1) (the mesh is then orthogonal); the tapered
                   stubs unchanged
    uniformStubs   the two tapered 1D stubs made uniform (their proximal
                   radius); the cone unchanged
    nonOrthCorr    one non-orthogonal corrector in the pressure equation

Writes floor.md and plots/floor.png.
"""
import argparse
import math
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import orderStudy as os_  # noqa: E402
import makeCases as mc  # noqa: E402

CASE = "thoracicAorta37"
FAMILY = "BDF4"
DTS = [1e-3, 5e-4, 2.5e-4]
END = 0.3


def ddt_corr_off(case):
    os_.edit(case/"system"/"fvSolution",
             [(r"^(\s+momentumPredictor\s+no;)", "\\1\n    ddtCorr                 no;")])


def tight_solvers(case):
    s = (case/"system"/"fvSolution").read_text()
    s = s.replace("tolerance       1e-12;", "tolerance       1e-15;")
    s = s.replace("maxIter         5000;", "maxIter         20000;")
    s = s.replace("maxIter         1000;", "maxIter         20000;")
    s = re.sub(r"p \{ tolerance 1e-10; relTol 0; \}",
               "p { tolerance 1e-12; relTol 0; }", s)
    assert "1e-15" in s and "1e-12; relTol" in s
    (case/"system"/"fvSolution").write_text(s)


def cylinder(case):
    r = math.sqrt(mc.radius(mc.L_STUB)*mc.radius(mc.L_SEG - mc.L_STUB))
    orig = mc.radius
    mc.radius = lambda x: r
    try:
        nz = int(round(mc.L_3D/mc.DX))
        (case/"system"/"blockMeshDict").write_text(mc.block_mesh(nz))
    finally:
        mc.radius = orig


def uniform_stubs(case):
    p = case/"constant"/"circulation1DProperties"
    s = p.read_text()
    for name in ("s17a", "s17b"):
        m = re.search(rf"^    {name} +\{{.*radiusProximal (\S+); radiusDistal (\S+);.*\}}$",
                      s, re.M)
        assert m, name
        line = m.group(0).replace(f"radiusDistal {m.group(2)};",
                                  f"radiusDistal {m.group(1)};")
        s = s[:m.start()] + line + s[m.end():]
    p.write_text(s)


def non_orth_corr(case):
    os_.edit(case/"system"/"fvSolution",
             [(r"^(\s+nNonOrthogonalCorrectors\s+)0;", r"\g<1>1;")])


VARIANTS = {
    "ddtCorrOff": ddt_corr_off,
    "tightSolvers": tight_solvers,
    "cylinder": cylinder,
    "uniformStubs": uniform_stubs,
    "nonOrthCorr": non_orth_corr,
}


def make_variant(variant):
    os_.SUFFIX = f"_{variant}"
    names = []
    for dt in DTS:
        case = HERE/"runs"/os_.run_name(FAMILY, dt)
        fresh = not os_.finished(case)
        name = os_.make(CASE, FAMILY, dt, END, False)
        if fresh:
            VARIANTS[variant](case)
        names.append(name)
    return names


def errors(names):
    """(e(dt0, dt1), e(dt1, dt2), order) of to3D.Q and to3D.P, relative to
    the RMS of the finest signal, and the finest run's signals"""
    sig = [os_.signals(CASE, n) for n in names]
    out = {}
    for key in ("to3D.Q", "to3D.P"):
        e1 = os_.diff(sig[0][key], sig[1][key])/os_.rms(sig[2][key])
        e2 = os_.diff(sig[1][key], sig[2][key])/os_.rms(sig[2][key])
        out[key] = (e1, e2, math.log2(e1/e2))
    return out, sig[2]


def analyse(variants):
    os_.SUFFIX = ""
    base_names = [os_.run_name(FAMILY, dt) for dt in DTS]
    rows = {"baseline": errors(base_names)}
    for v in variants:
        os_.SUFFIX = f"_{v}"
        names = [os_.run_name(FAMILY, dt) for dt in DTS]
        if all(os_.finished(HERE/"runs"/n) for n in names):
            rows[v] = errors(names)
    os_.SUFFIX = ""
    lines = [f"# The BDF 4 floor of {CASE}: one change at a time", "",
             f"BDF 4 at dt = {', '.join(f'{dt*1e3:g}' for dt in DTS)} ms, "
             f"0-{END} s from the smooth start; differences between "
             "successive steps, RMS over the window relative to the RMS of "
             "the finest run (to3D.Q, to3D.P), and the order they give; the "
             "last column is the finest run against the finest baseline run.",
             "",
             "| variant | Q: e(1, 0.5) | e(0.5, 0.25) | order | P: e(1, 0.5) "
             "| e(0.5, 0.25) | order | vs baseline (Q) |",
             "|---|---|---|---|---|---|---|---|"]
    base_sig = rows["baseline"][1]
    for v, (e, sig) in rows.items():
        q, p = e["to3D.Q"], e["to3D.P"]
        d = os_.diff(sig["to3D.Q"], base_sig["to3D.Q"])/os_.rms(base_sig["to3D.Q"])
        lines.append(f"| {v} | {q[0]:.2e} | {q[1]:.2e} | {q[2]:.2f} | "
                     f"{p[0]:.2e} | {p[1]:.2e} | {p[2]:.2f} | {d:.2e} |")
    (HERE/"floor.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    plot(rows)


def plot(rows):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, ax = plt.subplots(figsize=(4.5, 3.6))
    for v, (e, _) in rows.items():
        q = e["to3D.Q"]
        ax.loglog(DTS[:2], [q[0], q[1]], "o-", label=v)
    ax.loglog(DTS[:2], [rows["baseline"][0]["to3D.Q"][0],
                        rows["baseline"][0]["to3D.Q"][0]/16], "k:", lw=0.8,
              label="slope 4")
    ax.set_xticks(DTS[:2])
    ax.set_xticklabels(["1", "0.5"])
    ax.set_xlabel("dt [ms] (difference with the run at dt/2)")
    ax.set_ylabel("|Q(dt) - Q(dt/2)| / |Q|, RMS")
    ax.set_title("thoracicAorta37, BDF 4: the floor", fontsize=10)
    ax.legend(frameon=False, fontsize=7)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    (HERE/"plots").mkdir(exist_ok=True)
    fig.savefig(HERE/"plots"/"floor.png", dpi=150)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--variants", default=",".join(VARIANTS))
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--analyse", action="store_true")
    args = ap.parse_args()
    variants = args.variants.split(",")
    if not args.analyse:
        jobs = []
        for v in variants:
            jobs += [(v, n) for n in make_variant(v)]
        jobs.sort(key=lambda x: -len(x[1]))

        def run(job):
            os_.SUFFIX = f"_{job[0]}"
            os_.run(CASE, job[1])

        with ThreadPoolExecutor(args.jobs) as ex:
            list(ex.map(run, jobs))
    analyse(variants)


if __name__ == "__main__":
    main()
