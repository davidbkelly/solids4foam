# circulation1D

A one-dimensional (1D) model of blood flow in networks of elastic vessels, for
OpenFOAM v2412, solved by the mesh-free application `circulation1DFoam`. It is
part of the solids4foam `circulationModels` module, see [../README.md](../README.md),
and terminates its vessels with lumped circuits of
[circulation0D](../circulation0D/README.md). The benchmarks are those of
Boileau et al. (2015), in
[tutorials/circulationModels/circulation1D](../../../tutorials/circulationModels/circulation1D).

This is phases 1 to 3 of [circulation1D_spec.md](circulation1D_spec.md):
vessels with inflow, absorbing, resistance and lumped 0D terminals, joined into
networks by junctions, with properties that may vary along a vessel (tapering,
a wall thickness that follows the radius) in a well-balanced scheme. Phase 4,
the coupling to 3D: the `coupledFlow` ports and the `inertance` node (below),
and the coupler that binds the ports to CFD patches,
[../circulation1DCoupling](../circulation1DCoupling) (boundary condition
`circulation1DPressure`), verified by the rigid-segment benchmark
(`tutorials/circulationModels/coupling3D1D/rigidSegment3D`).

## Equations

In each vessel, for the area `A` and the flow `Q`,

```text
dA/dt + dQ/dx = 0
dQ/dt + d/dx(alpha Q^2/A + fluxPressure(A)) = -KR Q/A,   KR = 2 pi (zeta + 2) mu/rho
```

with the tube law of Boileau Eq. 3 (`sqrtElastic`),
`P = Pext + Pd + beta/Ad (sqrt(A) - sqrt(Ad))`, `beta = 4/3 sqrt(pi) E h`.
When `beta` and `Ad` vary along the vessel the momentum equation is
`dQ/dt + d/dx(alpha Q^2/A) + (A/rho) dP/dx = -KR Q/A` (Boileau Eq. 1), whose
state at rest is a uniform pressure with no flow (`A = Ad` only if `Pd + Pext`
is uniform).

## Numerics

Finite volumes on a uniform mesh per vessel: MUSCL reconstruction of `(A, Q)`
with a limiter, an approximate Riemann solver at the interior faces, and an
explicit SSP Runge-Kutta scheme. Every stage, each node sets the boundary state
of its vessel ends from the outgoing characteristic (the reconstructed face
value of the end cell) and its own condition, and the ODE states of the nodes
(for example a Windkessel's pressure) are advanced in the same stage. The time
step is `CFL min(dx/max|lambda|)` over the network; `advance(t0, t1)` lands on
`t1` exactly.

| `numerics` entry | Options | Default |
|---|---|---|
| `flux` | `HLL`, `Rusanov` | `HLL` |
| `reconstruction` | `MUSCL` (with `limiter`), `kExact` (with `polynomialOrder`) | `MUSCL` |
| `limiter` | `minmod`, `vanLeer`, `MC`, `none` (unlimited central slopes) | `minmod` |
| `polynomialOrder` | `1`, `2`, `3`: order 2, 3, 4 in space | |
| `timeScheme` | `SSPRK2`, `SSPRK3`, `RK4` | `SSPRK2` |
| `CFL` | | 0.9 |

**Higher order.** `reconstruction kExact` is the 1D form of the k-exact
least-squares reconstruction of solids4foam's higher-order finite volumes
(`higherOrderHelpers/kExactLeastSquares`): in each cell a polynomial of degree
`k = polynomialOrder` whose average is the cell average, fitted by least squares
to the averages of a stencil of 3 cells (k = 1, 2) or 5 (k = 3), shifted
one-sided at the vessel ends. Its face values are of order k + 1, and the
friction source is integrated over the cell by two-point Gauss quadrature of
the reconstruction. With a time scheme of the same order (`SSPRK3` for k = 2,
`RK4`, the classical scheme, for k = 3) the scheme is of order 3 or 4: see the
`singlePulseHighOrder` tutorial. It is unlimited, for smooth flows, and for
uniform vessels only (a tapered vessel's wall is uniform within each cell, a
second-order representation). A degree-4 reconstruction (fifth order) was
tried and went unstable on fine meshes with its one-sided end stencils; it is
not offered.

A limited reconstruction is TVD only for `CFL <= 0.5` (the solver warns above
it): beyond that the limiter erodes smooth peaks a little more every step, so
the error grows as the mesh is refined. Use `none` for smooth problems and
convergence tests; limiters clip extrema. The scheme is second order: a pulse
needs about 15 cells per standard deviation of its width to keep its peak to
1% over 1500 cell widths (see the single-pulse tutorials).

**Tapered vessels.** A vessel whose wall varies along it has the wall of each
cell's centre, uniform within the cell. A face between two cells is then a
two-vessel junction (mass conserved, total pressure continuous, on the
outgoing characteristics of both sides), solved by Newton: the mass flux is
common and each side has its own momentum flux. The reconstruction is of
`(P, Q)`, which are uniform at rest, so a network at rest (uniform pressure, no
flow) stays at rest to round-off whatever the taper: the scheme is
well-balanced, and needs no geometric source term. It stays second order (see
the `fluidAtRest` and `taperedOrder` tutorials). Uniform vessels keep the
`(A, Q)` reconstruction and the Riemann flux above.

## `constant/circulation1DProperties` (SI)

```c++
fluid     { rho 1060; mu 4e-3; profileOrder 2; alpha 1; }  // alpha must be 1 for now
numerics  { CFL 0.5; flux HLL; limiter MC; timeScheme SSPRK2; }

vesselDefaults  { youngsModulus 1.2e6; }    // optional: entries for every vessel

vessels
{
    carotid
    {
        length              0.126;
        nCells              126;        // or dx 1e-3;
        tubeLaw             sqrtElastic;    // default
        radius              3e-3;       // at the diastolic pressure; or
                                        // radiusProximal, radiusDistal (linear)
        wallThickness       0.3e-3;     // or wallThicknessModel ADAN56 (Eq. 17)
        youngsModulus       700e3;
        diastolicPressure   10933;      // default 0
        externalPressure    0;          // default 0
        initialPressure     0;          // uniform; default 0
        initialFlowRate     0;          // uniform; default 0
    }
}

nodes                   // every vessel end <vessel>.in (x = 0) or .out (x = L) once
{
    inlet   { type flowInlet; connection carotid.in; flowRate <Function1>; }
    outlet  { type lumped0D;  connection carotid.out; circuit { ... } }
    // in a network:
    split   { type junction;  connections (aorta.out iliacL.in iliacR.in); }
}

probes  { mid { vessel carotid; x midpoint; } }     // x: metres, in, midpoint, out
```

## Nodes

| Type | Condition |
|---|---|
| `flowInlet` | Flow `flowRate` (a `Function1` of time, m3/s, positive into the vessel), by Newton on the outgoing characteristic |
| `pressureInlet` | Pressure `pressure` (a `Function1` of time, Pa, the tube law's pressure): the area is `A(P)` and the velocity follows from the outgoing characteristic; reflects an incoming wave with Γ = -1 |
| `absorbing` | Non-reflecting: the incoming Riemann invariant keeps its initial value |
| `resistance` | `P - Pout = R Q` (entries `R`, `Pout`) |
| `lumped0D` | A circulation0D circuit (`circuit`, a circulationProperties of its own) with one `externalPort` for the vessel end; for example a three-element Windkessel with the `windkessel3` preset |
| `junction` | Two or more vessel ends: mass conserved and total pressure `P + rho U^2/2` continuous (Boileau Eqs. 7-8), by Newton with the analytic Jacobian |
| `inertance` | Two vessel ends joined by a rigid straight pipe in plug flow, `L dQ/dt = P1 - P2 - R Q`, `L = rho length/area` (entries `length`, `radius` or `area`, `R`, `pressure static` or `total`): the exact 0D reduction of a rigid 3D segment |
| `coupledFlow` | A port to a 3D domain: the flow is set from outside over each step (`network1D::setPortForcing`) and the end returns its pressure (entries `patch`, and `flowScale` when the patch is a symmetric part of the cross-section, e.g. 4 for a quarter) |

`lumped0D` exchanges flow and pressure through the port as a 3D domain does: the
vessel's flow on its characteristic goes into the port and the circuit returns
the port pressure, solved together by Newton every stage. The circuit's states
are advanced in the vessels' Runge-Kutta stages from
`circulationNetwork::derivatives()`; the circuit's own ODE solver is not used.
A circuit needs at least one state, so a purely resistive terminal is the
`resistance` node.

Periodic tables: a `table` or `tableFile` `Function1` with `outOfBounds repeat`.
A formula: `type expression; expression #{ ... arg() ... #};` with `arg()` the
time.

## Coupling interface

`network1D` has `advance(t0, t1)`, `saveState()` / `restoreState()` (a full
snapshot, so a coupling step can be repeated) and `interfacePressure(node)`.

Its `coupledFlow` nodes are its ports, with the interface of a circulation0D
network's, so that the 3D-0D coupling schemes apply unchanged:

```c++
network1D net(dict, t0);
y0 = net.state();                       // the whole state as one field
// Q = c0 + c1 a + c2 a^2 + c3 a^3, a = (t - t0)/(t1 - t0)
net.setPortForcing("to3D", t0, t1, {c0, c1, c2, c3});
net.advance(t0, t1);                    // evaluated at every Runge-Kutta stage
scalar P = net.portPressure("to3D");    // Pa, at t1
net.setState(t0, y0);                   // to repeat the step (implicit)
```

The port flow is the flow along the vessel at its end (positive from `in` to
`out`), so a port at an `out` end feeds a 3D inlet (`portInlet()` true) and a
port at an `in` end is fed by a 3D outlet: the sign convention of the 0D
coupler's ports. `portNames()`, `portPatch(name)` and `nEqns()` complete it.

Two checks of milestone M1:

- `inertance` against the linear closed form of the rigid-segment benchmark
  ([../docs/BENCHMARKS3D1D.md](../docs/BENCHMARKS3D1D.md)): second order, see
  the `coupling3D1D/rigidSegmentInertance` tutorial;
- the ports, driven as the coupler will by `Test-circulation1DPorts` with a
  rigid-pipe surrogate for the CFD: the coupled flow against the inertance
  node, and the stability limits of the explicit coupling (see
  `applications/test/circulation1DPorts/sweep.py`).

## Classes

| Directory | Contents |
|---|---|
| `tubeLaws` | `tubeLaw1D` (run-time selectable), `sqrtElastic` |
| `numerics` | `slopeLimiter1D`, `riemannFlux1D` (run-time selectable), `rungeKutta1D` (SSPRK2, SSPRK3, RK4), `kExactReconstruction1D` |
| `vessel1D` | `vessel1D`: fields, the wall of each cell, reconstruction, fluxes (Riemann or coupled), sources, sampling |
| `nodes1D` | `node1D` (run-time selectable) and the nodes above |
| `network1D` | `network1D`: vessels, nodes, probes, stepping, state, volume balance |
