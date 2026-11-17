#!/usr/bin/env python3
"""Order in time of the coupled pimpleFoam + 1D cases, by self-convergence.

    python3 orderStudy.py --case dobroserdova2016|womersleySegment|thoracicAorta37
                          [--families Euler,backward,BDF3,BDF4] [--dts ...] [--endTime ...]
                          [--jobs 3] [--analyse] [--smoke]
    python3 orderStudy.py --summary          # one figure over the cases done

Each case is run from rest over a short window at four time steps, halving,
for four time-scheme families, with the fluid ddt scheme and the 3D-1D
coupling order both at the family's order q:

    family    fluid ddt   coupling (implicit)   q
    Euler     Euler       1                     1
    backward  backward    2                     2
    BDF3      BDF 3       3                     3
    BDF4      BDF 4       4                     4

(BDF 3 and BDF 4 are the bdfDdtScheme of libsolids4FoamModels.) The 1D
vessels sub-step with SSPRK3 at CFL 0.1 and the PIMPLE loop iterates the
coupling to a pressure residual of 1e-10, both far below the error being
measured, so that the self-convergence measures the 3D scheme and the
coupling. The window is the start of the transient: the sine inflow's first
period (dobroserdova2016, 0.3 s), the ramped incident wave (womersleySegment,
alpha = 6, 0.6 s) and the first systole of the 37-artery network from rest
(thoracicAorta37, 0.3 s).

Observed order between three successive steps h, h/2, h/4 is
log2(|u_h - u_h/2| / |u_h/2 - u_h/4|) in the RMS over the window of the
port pressures and flows (coupling.dat) and of the 3D mid-plane pressure and
flow where the case samples them. The finest run of every family is also
compared with the finest BDF4 run. Writes <case>/orders.md and
<case>/plots/orders.png; --summary writes plots/timeOrders.png here.

--exact (dobroserdova2016 and womersleySegment, whose linear problems have
exact solutions: the delay equations of the resonator, and the Womersley
transfer function convolved with the incident wave) runs the vessels with
the linearised 1D equations ("linearised true", of which those solutions
are exact at any amplitude) at fourth order in space and time (kExact 3,
RK4, twice the cells), and reports, besides the self-convergence, the RMS
error of every run against the exact solution and the order it gives:
<case>/orders_linear.md, plots/orders_linear.png. The runs carry the suffix
_linear.
"""
import argparse
import math
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

FAMILIES = {
    "Euler":    dict(q=1, ddt="Euler"),
    "backward": dict(q=2, ddt="backward"),
    "BDF3":     dict(q=3, ddt="BDF 3"),
    "BDF4":     dict(q=4, ddt="BDF 4"),
}
COLOURS = {"Euler": "#999999", "backward": "#0072B2", "BDF3": "#009E73",
           "BDF4": "#D55E00"}
CASES = {
    "dobroserdova2016": dict(dts=[2e-3, 1e-3, 5e-4, 2.5e-4, 1.25e-4], end=0.3,
                             title="Dobroserdova: rigid pipe, sine inflow"),
    "womersleySegment": dict(dts=[4e-3, 2e-3, 1e-3, 5e-4], end=0.6,
                             title="Womersley segment, alpha = 6"),
    "thoracicAorta37":  dict(dts=[2e-3, 1e-3, 5e-4, 2.5e-4], end=0.3,
                             title="Thoracic aorta II in the 37-artery network"),
}
CFL = 0.1
P_TOLERANCE = 1e-10
SUFFIX = ""            # run-name suffix of a --limiter/--timeScheme variant
LIMITER = None
TIME_SCHEME = None
EXACT = False          # --exact: linearised vessels, errors against the exact
INFLOW = "sine"        # --inflow: the resonator's inflow (sine, step, smooth)
CYLINDER = False       # --cylinder: the aorta's cone replaced by a cylinder


def run_name(family, dt):
    return f"order_{family}_dt{dt:g}" + SUFFIX


def edit(path, pairs):
    s = path.read_text()
    for pattern, value in pairs:
        s, n = re.subn(pattern, value, s, count=1, flags=re.M)
        assert n == 1, (path, pattern)
    path.write_text(s)


def finished(case):
    log = case/"log.pimpleFoam"
    return log.is_file() and "\nEnd\n" in log.read_text()[-2000:]


def make(case_name, family, dt, end, force, suffix=""):
    """The case's generator writes the run, then the family's scheme and
    order, the tight tolerances and the window are set"""
    case_dir = HERE/case_name
    name = run_name(family, dt) + suffix
    case = case_dir/"runs"/name
    if finished(case) and not force:
        return name
    if case.exists():
        shutil.rmtree(case)
    sys.path.insert(0, str(case_dir))
    for m in [m for m in sys.modules if m in ("makeCases", "analyse")]:
        del sys.modules[m]
    import makeCases as mc                                  # noqa: E402
    sys.path.pop(0)
    (case_dir/"runs").mkdir(exist_ok=True)
    if case_name == "dobroserdova2016":
        mc.make_3d("c350", INFLOW, name=name.replace("order_", "order3D_"),
                   dt=dt, end=end, force=True)
        (case_dir/"runs"/name.replace("order_", "order3D_")).rename(case)
    elif case_name == "womersleySegment":
        mc.make(6.0, dt, True, name=name, end=end)
    elif case_name == "thoracicAorta37":
        mc.write_coupled(case, dt, 0.0, end)
        # A smooth start: the tutorial's network starts at P = 0 against
        # terminal pressures of 432.6 Pa with the inflow's Fourier series at
        # its t = 0 value, and the fronts this launches (one 1D cell wide,
        # 0.3 ms at the port) are features at the scale of dt that hide the
        # schemes' orders. Here the network starts in equilibrium at the
        # terminal pressure and the inflow rises under a sin^2 ramp of 50 ms.
        edit(case/"constant"/"circulation1DProperties",
             [(r"^(\s+initialPressure\s+)0;", r"\g<1>432.6;"),
              (r"1e-6\*1e6\*\(3\.1199",
               "((arg() < 0.05) ? sqr(sin(pi()*arg()/0.1)) : 1)*(3.1199")])
        if CYLINDER:
            # The cone replaced by the cylinder of its geometric-mean radius
            # (floorStudy.py): the non-orthogonal cone mesh caps BDF 3/4 at
            # first order in pimpleFoam, the cylinder does not
            sys.path.insert(0, str(case_dir))
            import floorStudy
            sys.path.pop(0)
            floorStudy.cylinder(case)
    else:
        raise SystemExit(f"unknown case {case_name}")
    f = FAMILIES[family]
    edit(case/"system"/"fvSchemes",
         [(r"^(\s+default\s+)backward;.*$", rf"\g<1>{f['ddt']};")])
    if f["ddt"] not in ("Euler", "backward"):
        edit(case/"system"/"controlDict",
             [(r'^libs\s+\("libcirculation1DCoupling.so"\);',
               'libs            ("libcirculation1DCoupling.so" '
               '"libsolids4FoamModels.so");     // BDF ddt scheme')])
    edit(case/"constant"/"circulation1DProperties",
         [(r"^(\s+order\s+)\S+;", rf"\g<1>{f['q']};"),
          (r"^(\s+CFL\s+)\S+;", rf"\g<1>{CFL};")])
    if LIMITER:
        edit(case/"constant"/"circulation1DProperties",
             [(r"^(\s+limiter\s+)\S+;", rf"\g<1>{LIMITER};")])
    if TIME_SCHEME:
        edit(case/"constant"/"circulation1DProperties",
             [(r"^(\s+timeScheme\s+)\S+;", rf"\g<1>{TIME_SCHEME};")])
    if EXACT:
        # Linearised vessels at fourth order in space and time, on twice the
        # cells, so that the 1D error is far below the time error measured
        p = case/"constant"/"circulation1DProperties"
        n = int(re.search(r"^\s+nCells\s+(\d+);", p.read_text(), re.M)[1])
        edit(p, [(r"^(\s+initialPressure\s+0;)",
                  "\\1\n        linearised          true;"),
                 (r"^(\s+nCells\s+)\d+;", rf"\g<1>{2*n};"),
                 (r"^(\s+timeScheme\s+)\S+;.*$",
                  "    timeScheme      RK4;\n    reconstruction  kExact;\n"
                  "    polynomialOrder 3;")])
    s = (case/"system"/"fvSolution").read_text()
    i = s.index("residualControl")
    s = s[:i] + re.sub(r"tolerance\s+[\d.e+-]+;", f"tolerance {P_TOLERANCE:g};",
                       s[i:], count=1)
    (case/"system"/"fvSolution").write_text(s)
    # Every step to coupling.dat and the function objects; no field output
    # but the last time
    by_step = re.search(r"^writeControl\s+timeStep;",
                        (case/"system"/"controlDict").read_text(), re.M)
    edit(case/"system"/"controlDict",
         [(r"^endTime\s+\S+;", f"endTime         {end:g};"),
          (r"^writeInterval\s+\S+;",
           f"writeInterval   {int(round(end/dt)) if by_step else end:g};")])
    return name


def run(case_name, name):
    case = HERE/case_name/"runs"/name
    if finished(case):
        return
    cmds = ["blockMesh"]
    if (case/"system"/"topoSetDict").is_file():
        cmds.append("topoSet")
    for c in cmds:
        with open(case/f"log.{c}", "w") as fh:
            subprocess.run([c], cwd=case, stdout=fh, stderr=subprocess.STDOUT,
                           check=True)
    with open(case/"log.pimpleFoam", "w") as fh:
        subprocess.run(["pimpleFoam"], cwd=case, stdout=fh,
                       stderr=subprocess.STDOUT)
    print(f"{case_name}/{name}: {'finished' if finished(case) else 'FAILED'}",
          flush=True)


def signals(case_name, name):
    """{label: (t, x)}: the port pressures and flows from coupling.dat, and
    the mid-plane zone values where the case has them"""
    case = HERE/case_name/"runs"/name
    d = case/"postProcessing"/"circulation1DCoupling"
    f = next(iter(d.iterdir()))/"coupling.dat"
    names = None
    with open(f) as fh:
        for line in fh:
            if line.startswith("# time"):
                names = line[1:].split()
                break
    a = np.loadtxt(f, comments="#", ndmin=2)
    out = {}
    for i, n in enumerate(names):
        if re.fullmatch(r"(to3D|from3D)\.(P|Q)", n):
            out[n] = (a[:, 0], a[:, i])
    for suffix in ("P", "Q"):
        dd = case/"postProcessing"/f"mid_{suffix}"
        if dd.is_dir():
            b = np.loadtxt(next(iter(dd.iterdir()))/"surfaceFieldValue.dat",
                           comments="#", ndmin=2)
            out[f"mid.{suffix}"] = (b[:, 0], b[:, 1])
    return out


def diff(a, b):
    """RMS of a - b on a's times over the common span"""
    t, x = a
    tb, xb = b
    m = (t >= tb[0] - 1e-12) & (t <= tb[-1] + 1e-12)
    return math.sqrt(np.mean((x[m] - np.interp(t[m], tb, xb))**2))


def rms(a):
    return math.sqrt(np.mean(a[1]**2))


def exact_signals(case_name, t):
    """{label: (t, x)} of the exact solution of the linear problem at the
    times t, for the port signals the case logs"""
    case_dir = HERE/case_name
    sys.path.insert(0, str(case_dir))
    for m in [m for m in sys.modules if m in ("makeCases", "analyse")]:
        del sys.modules[m]
    import analyse as an                                    # noqa: E402
    sys.path.pop(0)
    if case_name == "dobroserdova2016":
        Q, P = an.exact(INFLOW, 3.5, "hard", t, pipeArea=an.A3D, n2=64000)
        return {"to3D.Q": (t, Q), "to3D.P": (t, P)}
    if case_name == "womersleySegment":
        # the circular pipe of the wedge's cross-section (see its README)
        ex = an.exact_transient(6.0, t, A=an.A_WEDGE, R=an.R_WEDGE)
        return {k: (t, x) for k, x in ex.items()}
    raise SystemExit(f"{case_name} has no exact solution")


def analyse(case_name, families, dts, end):
    dts = sorted(dts, reverse=True)
    case_dir = HERE/case_name
    lines = [f"# Time order of {case_name} (self-convergence, 0-{end:g} s, "
             f"from rest; 1D at CFL {CFL}, p residual {P_TOLERANCE:g}"
             + (f", 1D limiter {LIMITER}" if LIMITER else "")
             + (f", 1D {TIME_SCHEME}" if TIME_SCHEME else "")
             + (", linearised vessels at fourth order (kExact 3, RK4, "
                "twice the cells)" if EXACT else "") + ")", ""]
    if EXACT:
        lines += ["Besides the self-convergence, every run is compared with "
                  "the exact solution of the linear problem (RMS error / RMS "
                  "of the signal), with the order between successive "
                  "steps.", ""]
    finest, table, orders, exact_table = {}, {}, {}, {}
    for fam in families:
        names = [run_name(fam, dt) for dt in dts]
        done = [n for n in names if finished(case_dir/"runs"/n)]
        failed = [n for n in names if n not in done]
        if failed:
            lines += [f"{fam}: failed or not run at dt = "
                      + ", ".join(f"{dt:g}" for n, dt in zip(names, dts)
                                  if n in failed), ""]
        done = done[:next((i for i, n in enumerate(names) if n in failed),
                          len(names))]
        if len(done) < 3:
            lines += [f"{fam}: {len(done)} consecutive finished run(s), need 3",
                      ""]
            continue
        sig = {n: signals(case_name, n) for n in done}
        finest[fam] = sig[done[-1]]
        keys = list(sig[done[0]])
        lines += [f"## {fam} (q = {FAMILIES[fam]['q']})", "",
                  "| dt pair | " + " | ".join(keys) + " |",
                  "|---|" + "---|"*len(keys)]
        for i in range(len(done) - 2):
            row = []
            for key in keys:
                e1 = diff(sig[done[i]][key], sig[done[i + 1]][key])
                e2 = diff(sig[done[i + 1]][key], sig[done[i + 2]][key])
                p = math.log2(e1/e2) if e2 > 0 else math.nan
                row.append(f"{p:.2f}")
                orders.setdefault(fam, {}).setdefault(key, []).append(p)
                tab = table.setdefault(fam, {}).setdefault(key, [])
                if i == 0:
                    tab.append((dts[0], e1/rms(sig[done[-1]][key])))
                tab.append((dts[i + 1], e2/rms(sig[done[-1]][key])))
            lines.append(f"| {dts[i]:g} / {dts[i+1]:g} / {dts[i+2]:g} | "
                         + " | ".join(row) + " |")
        if EXACT:
            lines += ["", f"Against the exact solution ({fam}):", "",
                      "| dt | " + " | ".join(f"{k} | order" for k in keys)
                      + " |", "|---|" + "---|---|"*len(keys)]
            prev = None
            for n, dt in zip(done, dts):
                ex = exact_signals(case_name, sig[n][keys[0]][0])
                errs = [diff(sig[n][k], ex[k])/rms(ex[k]) for k in keys]
                cells = [f"{dt:g}"]
                for j, (k, e) in enumerate(zip(keys, errs)):
                    cells.append(f"{e:.2e}")
                    cells.append("" if prev is None else
                                 f"{math.log2(prev[j]/e):.2f}")
                    exact_table.setdefault(fam, {}).setdefault(k, []).append(
                        (dt, e))
                lines.append("| " + " | ".join(cells) + " |")
                prev = errs
        # iterations per step and cost of the finest run
        log = (case_dir/"runs"/done[-1]/"log.pimpleFoam").read_text()
        steps = len(re.findall(r"^Time = ", log, re.M))
        iters = len(re.findall(r"^PIMPLE: iteration", log, re.M))
        ex = re.findall(r"ExecutionTime = ([\d.]+)", log)
        lines += ["", f"Finest run {done[-1]}: {steps} steps, "
                  f"{iters/max(steps, 1):.1f} PIMPLE iterations per step, "
                  f"{float(ex[-1]) if ex else math.nan:.0f} s", ""]
    if "BDF4" in finest and len(finest) > 1:
        lines += ["## Finest runs against the finest BDF4 run (RMS of the "
                  "difference / RMS of the signal)", ""]
        for fam, sig in finest.items():
            if fam == "BDF4":
                continue
            lines.append(f"- {fam}: " + ", ".join(
                f"{k}: {diff(sig[k], finest['BDF4'][k])/rms(finest['BDF4'][k]):.2e}"
                for k in sig))
        lines.append("")
    (case_dir/f"orders{SUFFIX}.md").write_text("\n".join(lines))
    print("\n".join(lines))
    plot(case_name, table, exact_table)
    return table, orders


def plot(case_name, table, exact_table=None):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    if not table:
        return
    keys = list(next(iter(table.values())))
    fig, axs = plt.subplots(1, len(keys), figsize=(3.3*len(keys), 3.4),
                            squeeze=False)
    for ax, key in zip(axs[0], keys):
        for fam, d in table.items():
            if key in d:
                dt, e = zip(*d[key])
                ax.loglog(dt, e, "o-", color=COLOURS[fam], label=fam)
            if exact_table and fam in exact_table and key in exact_table[fam]:
                dt, e = zip(*exact_table[fam][key])
                ax.loglog(dt, e, "s--", color=COLOURS[fam], ms=4,
                          label=f"{fam} vs exact")
        slopes(ax, table, key)
        ax.set_title(key, fontsize=9)
        ax.set_xlabel("dt [s]")
        ax.grid(alpha=0.3, which="both")
        ticks(ax, table, key)
    axs[0, 0].set_ylabel("|u(dt) - u(dt/2)| / |u|, RMS"
                         + ("; vs exact: |u - exact| / |exact|" if exact_table
                            else ""))
    axs[0, 0].legend(frameon=False, fontsize=7)
    fig.suptitle(CASES[case_name]["title"], fontsize=10)
    fig.tight_layout()
    (HERE/case_name/"plots").mkdir(exist_ok=True)
    fig.savefig(HERE/case_name/"plots"/f"orders{SUFFIX}.png", dpi=150)
    plt.close(fig)


def ticks(ax, table, key):
    """The time steps as the x ticks, no minor labels"""
    from matplotlib.ticker import NullFormatter
    dts = sorted({dt for d in table.values() if key in d for dt, _ in d[key]})
    ax.set_xticks(dts)
    ax.set_xticklabels([f"{dt*1e3:g}" for dt in dts])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("dt [ms]")


def slopes(ax, table, key):
    """Reference slopes dt^q through the coarsest point of each family"""
    for fam, d in table.items():
        if key not in d:
            continue
        q = FAMILIES[fam]["q"]
        dt0, e0 = d[key][0]
        dt1 = d[key][-1][0]
        ax.loglog([dt0, dt1], [e0, e0*(dt1/dt0)**q], ":", color=COLOURS[fam],
                  lw=0.8)


def summary():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    done = [c for c in CASES if (HERE/c/"orders.md").is_file()]
    fig, axs = plt.subplots(1, len(done), figsize=(3.6*len(done), 3.6),
                            squeeze=False)
    rows = ["| case | signal | " + " | ".join(FAMILIES) + " |",
            "|---|---|" + "---|"*len(FAMILIES)]
    for ax, c in zip(axs[0], done):
        table, orders = analyse(c, list(FAMILIES), CASES[c]["dts"], CASES[c]["end"])
        key = "to3D.P"
        for fam, d in table.items():
            if key in d:
                dt, e = zip(*d[key])
                ax.loglog(dt, e, "o-", color=COLOURS[fam], label=fam)
        slopes(ax, table, key)
        ax.set_title(CASES[c]["title"], fontsize=9)
        ax.set_xlabel("dt [s]")
        ax.grid(alpha=0.3, which="both")
        ticks(ax, table, key)
        for k in next(iter(orders.values()), {}):
            rows.append(f"| {c} | {k} | " + " | ".join(
                f"{orders[f][k][-1]:.2f}" if f in orders and k in orders[f]
                else "-" for f in FAMILIES) + " |")
    axs[0, 0].set_ylabel("|P(dt) - P(dt/2)| / |P|, RMS, at the inlet port")
    axs[0, 0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    (HERE/"plots").mkdir(exist_ok=True)
    fig.savefig(HERE/"plots"/"timeOrders.png", dpi=150)
    plt.close(fig)
    text = ("# Observed orders in time (the finest triple of steps) of the "
            "coupled pimpleFoam cases\n\n" + "\n".join(rows) + "\n")
    (HERE/"timeOrders.md").write_text(text)
    print(text)


def main():
    global SUFFIX, LIMITER, TIME_SCHEME, CFL, EXACT, INFLOW, CYLINDER
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--case", choices=list(CASES))
    ap.add_argument("--families", default=",".join(FAMILIES))
    ap.add_argument("--dts")
    ap.add_argument("--endTime", type=float)
    ap.add_argument("--jobs", type=int, default=3, help="runs at a time")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--analyse", action="store_true", help="no runs")
    ap.add_argument("--smoke", action="store_true",
                    help="three steps of the coarsest run of each family")
    ap.add_argument("--summary", action="store_true")
    ap.add_argument("--limiter", help="1D limiter in place of the case's (e.g. none)")
    ap.add_argument("--timeScheme", help="1D time scheme in place of the case's")
    ap.add_argument("--cfl", type=float, help=f"1D CFL in place of {CFL}")
    ap.add_argument("--exact", action="store_true",
                    help="linearised vessels; errors against the exact solution")
    ap.add_argument("--inflow", default="sine", choices=["sine", "step", "smooth"],
                    help="dobroserdova2016: the inflow (smooth: the Gaussian "
                    "pulse, whose exact solution is smooth)")
    ap.add_argument("--cylinder", action="store_true",
                    help="thoracicAorta37: the cone replaced by a cylinder")
    args = ap.parse_args()
    LIMITER, TIME_SCHEME, EXACT = args.limiter, args.timeScheme, args.exact
    SUFFIX = (f"_{LIMITER}" if LIMITER else "") + (f"_{TIME_SCHEME}" if TIME_SCHEME else "")
    if args.cylinder:
        if args.case != "thoracicAorta37":
            ap.error("--cylinder is for thoracicAorta37")
        CYLINDER = True
        SUFFIX += "_cylinder"
        CASES[args.case]["title"] = "Thoracic aorta II (as a cylinder) in the 37-artery network"
    if args.inflow != "sine":
        if args.case != "dobroserdova2016":
            ap.error("--inflow is for dobroserdova2016")
        INFLOW = args.inflow
        SUFFIX += f"_{INFLOW}"
        CASES[args.case]["title"] = f"Dobroserdova: rigid pipe, {INFLOW} inflow"
    if args.cfl:
        CFL = args.cfl
        SUFFIX += f"_cfl{CFL:g}"
    if EXACT:
        if args.case not in ("dobroserdova2016", "womersleySegment"):
            ap.error("--exact is for dobroserdova2016 and womersleySegment")
        SUFFIX += "_linear"
    if args.summary:
        summary()
        return
    if not args.case:
        ap.error("--case or --summary")
    c = CASES[args.case]
    families = args.families.split(",")
    dts = [float(x) for x in args.dts.split(",")] if args.dts else c["dts"]
    end = args.endTime or c["end"]
    if args.smoke:
        for fam in families:
            name = make(args.case, fam, dts[0], 3*dts[0], True, "_smoke")
            run(args.case, name)
            print(signals(args.case, name).keys())
        return
    if not args.analyse:
        names = [(fam, make(args.case, fam, dt, end, args.force))
                 for fam in families for dt in dts]
        # the expensive runs first
        names.sort(key=lambda x: -len(x[1]))
        with ThreadPoolExecutor(args.jobs) as ex:
            list(ex.map(lambda x: run(args.case, x[1]), names))
    analyse(args.case, families, dts, end)


if __name__ == "__main__":
    main()
