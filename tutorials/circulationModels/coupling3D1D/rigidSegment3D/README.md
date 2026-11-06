# Rigid segment, 3D-1D: `rigidSegment3D`

Prepared by David Kelly

---

## Tutorial Aims

- Verify the 3D-1D coupling (`circulation1DCoupling`) against an exact
  solution: a rigid 3D pipe between two 1D vessels.
- Measure the order of the coupling in the time step, for explicit coupling of
  orders 1-4 and implicit coupling, and its stability limits, with the CFD
  time schemes `backward` (BDF2), `BDF 3` and `BDF 4`.

---

## Case Overview

A Gaussian pressure pulse travels down a 1D vessel (1 m, the single-pulse
vessel of Boileau et al. 2015: `c0` = 6.17 m/s), through a rigid 3D pipe
(8 cm long, radius 1 cm, pimpleFoam) and on down a second 1D vessel with an
absorbing end. The 3D pipe takes its pressure at both ends from the vessels
(`circulation1DPressure`), and the vessels take the 3D flow through their
`coupledFlow` ports.

With slip walls the exact 3D solution is plug flow, an inertance
`L = rho l/A3D` in series, and in the linear limit the pipe is a first-order
low-pass filter:

```text
tau dPt/dt + Pt = Pinc,   Pr = Pinc - Pt,   tau = L/(2 Zc) = 6.52 ms
```

with a closed form for a Gaussian `Pinc` (see
`src/circulationModels/docs/BENCHMARKS3D1D.md`). Finite volumes represent plug
flow exactly, so what is left is the coupling's error and the 1D
discretisation's. The inflow of Boileau Eq. 14 is scaled by 1e-3 so that the
pulse is linear. The 3D patches are 32-sided polygons, 0.64% smaller than
`pi R^2`; `A3D` is their area.

The coupling number is `dt/tau`: within one step each vessel end answers the
3D flow as a resistor `Zc`. A model problem (the CFD's BDF formula for the
pipe, the flow extrapolated to the end of the step by the coupling's order
`p`) gives the explicit stability limits in `dt/tau`:

| CFD ddt | p = 1 | 2 | 3 | 4 |
| --- | --- | --- | --- | --- |
| `backward` (BDF2) | 4 | 1.333 | 0.571 | 0.267 |
| `BDF 3` | 4.5 | 2.222 | 0.952 | 0.444 |
| `BDF 4` | 1.643 | 2.25 | 1.524 | 0.711 |

(BDF 4's order 1 loses stability through a complex pair of roots, so before
order 2.) The implicit iteration contracts by `dt/(a0 tau)` per outer
corrector, `a0` the leading BDF coefficient (3/2, 11/6, 25/12).

`BDF 3` and `BDF 4` are the `bdfDdtScheme` of `libsolids4FoamModels`; the
cases that use them load it in `system/controlDict`.

| Path | Purpose |
| --- | --- |
| `base/` | the case: mesh (blockMesh, 3072 cells), BCs, `constant/circulation1DProperties` (the vessels, the ports, `coupling`) |
| `sweep.py` | sets up and runs the cases in `runs/` |
| `analyse.py` | the comparison with the closed form, the orders, the stability checks |
| `Allrun` | one case (implicit order 2, dt = 0.5 ms) and the 1D-only run, about 2 minutes |
| `Allrun-sweep` | the whole study for backward, BDF 3 and BDF 4: about 100 cases, about an hour with 6 at a time (`DDT=backward` for one scheme) |
| `Allrun-highOrder` | fourth order in space and time (`highOrder.py`): 17 cases, about 15 minutes |

---

## Running the Case

Needs OpenFOAM.com and the `circulationModels` module built
(`libcirculation1DCoupling`).

```bash
./Allrun            # results.md, plots/rigidSegment3D.png
./Allrun-sweep      # NPROC=6 by default
```

---

## Set-up notes

- **Convergence of the implicit coupling.** Each outer corrector of PIMPLE is
  one fixed-point iteration on the interface flow; `residualControl` on `p`
  alone ends the step. The initial residual of `U` levels off whatever the
  coupling (the corrected `U` does not satisfy the predictor's equation), so a
  tolerance on it is never met.
- **Solver tolerances** are on the normalised residual: below about 1e-12 the
  solvers run to `maxIter` on round-off while the pipe is still at rest.
- **SSPRK3 in the vessels**: its time error at CFL 0.5 is 6e-6 of the peak;
  SSPRK2's, 1e-3, would hide the coupling's. With BDF 3 and BDF 4 the vessels
  run at CFL 0.1 (time error about 5e-8), or the 1D would cap a third- or
  fourth-order coupling.

---

## Results

From `./Allrun-sweep` (`results.md`, `plots/rigidSegment3D.png`):

- **Against the closed form.** The 3D flux and the pressures at both ends of
  the pipe agree with the closed form to the 1D discretisation error: 1.1e-3
  of the incident peak for the flux and outlet pressure, 3.5e-3 for the inlet
  pressure (implicit order 2, dt = 0.5 ms), against 1.2e-3 and 4.1e-3 for the
  1D-only run with the pipe as an inertance node. Mass through the pipe is
  conserved to 1e-12 of the peak.
- **Consistency.** The finest coupled run (dt = 0.03125 ms) differs from the
  1D-only run by 3.4e-6 of the peak.
- **Order in dt** (against each scheme's finest run, two finest dt):

  | CFD ddt | explicit 1 | explicit 2 | explicit 3 | explicit 4 | implicit |
  | --- | --- | --- | --- | --- | --- |
  | `backward` | 1.00 | 1.99 | 2.14 | 2.01 | 2.00 (order 2) |
  | `BDF 3` | 1.02 | 1.97 | 2.99 | 3.17 | 2.98 (orders 1 and 3) |
  | `BDF 4` | 1.02 | 2.00 | 2.97 | 3.89 | 3.78 (orders 1 and 4) |

  These orders are of the 3D flux. Explicit coupling of order `p` with a CFD
  scheme of order `q` converges at order `min(p, q)`; implicit coupling of
  orders 1 and `q` both at `q`. That is the resistor at work: within a step
  the pipe sees the vessels' pressure at the end of the step only, which
  depends on the flow there and not on the shape of the forcing over the
  step. The waves the ports send into the vessels do carry that shape, and
  there implicit order `p` gives `min(p, q)` (see the fourth-order study
  below). The BDF 4 orders fall a little short of
  4 on the finest step, where the error (5e-7 of the peak) nears the vessels'
  time error. At dt = 0.25 ms the error is 2.2e-4 with backward, 8.8e-6 with
  BDF 3, 5.5e-7 with BDF 4; with BDF 4 the coupled run is at the 1D
  discretisation error against the closed form from dt = 1 ms on, with
  backward only from 0.25 ms. The finest runs of BDF 3 and BDF 4 differ from
  the 1D run with the inertance node by 1.6e-8 and 4.4e-9 of the peak.
- **Stability.** Every explicit run at 0.9 times its predicted limit decays and
  every one at 1.1 times it grows (16 runs: backward orders 2-4, BDF 3 orders
  3-4, BDF 4 orders 1, 3, 4). The amplification per step of the step-to-step
  oscillation after the pulse, fitted from the runs, matches the model's
  largest root: for example 1.1721 against 1.1721 (backward, order 2, 1.1x),
  1.1508 against 1.1508 (BDF 3, order 3), 1.1373 against 1.1371 (BDF 4,
  order 4); BDF 4's order 1, the complex pair, 1.035 against 1.025.
- **Cost.** The implicit iteration takes 11 outer correctors per step on
  average at dt = 0.5 ms (19 while the pulse passes); the case runs in about
  2 minutes in serial.

In explicit coupling the coupler samples the flux once, at the start of the
step, so the flux at the end of a step is the sample logged on the next row;
`analyse.py` reads it so.

## Fourth order in space and time

`./Allrun-highOrder` (`highOrder.py`, 17 cases, about 15 minutes with 6 at a
time): the vessels with the k-exact reconstruction of degree 3 and RK4, the
pipe with `BDF 4`, at full amplitude (the 3D solution is exact at any
amplitude, and at 1e-3 the vessels' round-off, 5e-7 of the peak, would hide a
fourth-order error). The mesh and the 3D step are refined together, from
2.5 mm and 2 ms to 0.3125 mm and 0.25 ms, each run against implicit order 4 at
0.15625 mm and 0.125 ms. Orders between successive levels
(`highOrder.md`, `plots/highOrder3D1D.png`):

| coupling | 3D flux | pipe-end pressures | waves in the vessels |
| --- | --- | --- | --- |
| implicit 4 | 3.90, 3.98, 4.08 | 3.89, 3.99, 4.09 | 3.95, 4.00, 4.09 |
| explicit 4 | 3.90, 3.97, 3.99 | 3.90, 3.98, 3.97 | 3.96, 4.00, 3.98 |
| implicit 2 | 3.90, 3.98, 4.08 | 3.89, 3.99, 4.08 | 1.60, 1.89, 1.97 |
| implicit 1 | 3.89, 3.91, 2.43 | 3.79, 3.34, 1.87 | 0.98, 1.00, 1.00 |

The coupled 3D-1D problem is fourth order in space and time with
fourth-order coupling: at the finest level the error is 4.4e-7 of the peak
with implicit order 4 and 1.4e-6 with explicit order 4. A lower coupling
order spoils the waves in the vessels at once (I2 second order, I1 first) and
the 3D side later, with a much smaller constant: I1's first-order wave error
reaches the 3D flux at the finest level, I2's not yet. In the plot the
implicit 4 and implicit 2 curves coincide on the 3D flux and the pressures.
