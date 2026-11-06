# circulationModels

Models of the circulatory system for solids4foam, from a mesh-free network of
lumped (0D) components, such as chambers, valves, resistances and compliances,
to the coupling of that network with a 3D CFD domain. It is currently only
available for OpenFOAM.com (tested with v2412), where the libraries are built
with solids4foam, and to networks of 1D vessels.

## Contents

- `circulation0D` (library `libcirculation0D`): Lumped-parameter (0D)
  circulation network: components, waveforms, per-cycle metrics, presets, and
  the ODE solution of the network. See
  [circulation0D/README.md](circulation0D/README.md)
- `circulationCoupling` (library `libcirculationCoupling`): Coupling of a
  circulation network to a CFD domain: the `circulationPressure` boundary
  condition and the coupler (explicit/implicit, orders 1 to 4)
- `circulation1D` (library `libcirculation1D`): One-dimensional vessels
  (finite volumes, MUSCL, Riemann fluxes, SSP Runge-Kutta) joined by nodes,
  with circulation0D circuits as terminals. See
  [circulation1D/README.md](circulation1D/README.md)
- `circulation1DCoupling` (library `libcirculation1DCoupling`): Coupling of a
  1D network to a CFD domain through its `coupledFlow` ports: the
  `circulation1DPressure` boundary condition and its coupler, for now a copy
  of `circulationCoupling`'s (the same schemes), until the two share one
- `tools` (no library): `circulationPlot` (plots result columns, needs
  gnuplot), `circulation1DTree` (draws the vessel network of a
  circulation1DFoam case, needs matplotlib) and `addHeader` (adds the
  solids4foam banner to case dictionaries)
- `docs` (no library): The description of the 0D/3D coupling,
  [docs/COUPLING.md](docs/COUPLING.md), and the 3D-1D benchmarks,
  [docs/BENCHMARKS3D1D.md](docs/BENCHMARKS3D1D.md)

The `circulationCoupling` library is not specific to the 0D models: the
coupler talks to the network through ports.

## Time scheme

Coupling orders 3 and 4 want a higher-order time scheme in the CFD, so that the
CFD time error does not dominate. OpenFOAM and solids4foam have no
first-derivative `ddtScheme` above second order (the solids4foam `bdfD2dt2`
scheme is a second-derivative scheme, `d2dt2Schemes`, for solid dynamics), so
the tutorials use `backward`, which limits the coupled solution to second
order. The coupler prints a warning when the CFD time scheme is of lower order
than the coupling. A `BDF` ddt scheme of order 1 to 6 is being proposed
separately.

The applications are `applications/solvers/circulationFoam`, which solves a
circulation network with no mesh, and `applications/solvers/circulation1DFoam`,
which solves a network of 1D vessels.

The tutorials are in `tutorials/circulationModels`: `circulation0D` for the
0D circuits and their coupling to CFD, `circulation1D` for the 1D vessels,
`coupling3D1D` for the 3D-1D coupling benchmarks.

## Building

The libraries are built by `src/Allwmake` and the applications by
`applications/Allwmake`, along with the rest of solids4foam. To build only this
module, after sourcing OpenFOAM:

```bash
(cd src/circulationModels && ./Allwmake)
(cd applications/solvers/circulationFoam && ./Allwmake)
(cd applications/solvers/circulation1DFoam && ./Allwmake)
```

The libraries are written to `$FOAM_MODULE_LIBBIN`, like the other solids4foam
libraries. For a CFD case, load the coupling library in `system/controlDict`:

```c++
libs ("libcirculationCoupling.so");       // a 0D circuit
libs ("libcirculation1DCoupling.so");     // a 1D network
```

## Quick start

```bash
cd tutorials/circulationModels/circulation0D/healthyHuman
./Allrun    # 20 cardiac cycles of a 4-chamber closed loop in under a second
# results: postProcessing/circulation/0/circulation.dat
```

## Status

The step-by-step tutorial and the verification scripts and reference data for
these models are not included yet. The tutorials in
`tutorials/circulationModels` are self-contained.
