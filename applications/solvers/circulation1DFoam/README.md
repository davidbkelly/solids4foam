---
sort: 6
---

# `circulation1DFoam` solver

`circulation1DFoam` solves a network of one-dimensional (1D) arterial vessels:
the area and flow along each vessel, with inflows, absorbing and resistance
terminals and lumped 0D circuits at the vessel ends. It needs no mesh. It is
part of the [`circulationModels`](../../../src/circulationModels/README.md)
module and is currently only built for OpenFOAM.com.

The network is a dictionary, `constant/circulation1DProperties`, described in
[`circulation1D`](../../../src/circulationModels/circulation1D/README.md).

---

## Usage

```bash
circulation1DFoam
```

A case needs only two files, as there is no mesh:

* `constant/circulation1DProperties`: the fluid, the numerics, the vessels, the
  nodes and the probes.
* `system/controlDict`: `startTime`, `endTime`, `deltaT` and the write times.

`deltaT` is the **output interval**: the solver takes its own CFL-limited steps
inside it and lands on every output time. See
[`tutorials/circulationModels/circulation1D`](../../../tutorials/circulationModels/circulation1D).

To see the network of a case (vessels as long and wide as they are, junctions,
terminals and probes), draw it with `circulation1DTree`, in
[`src/circulationModels/tools`](../../../src/circulationModels/tools) (it
needs matplotlib):

```bash
circulation1DTree                   # the case in the current directory, in a window
circulation1DTree -o tree.png boileau2015/ADAN56
circulation1DTree --topology boileau2015/ADAN56 # the branching only, every vessel one length
```

---

## Output

Every `deltaT`, a row is appended to each file in
`postProcessing/circulation1D/<startTime>/` (SI; the first lines name the
columns and the units):

* `probes/<probe>.dat`: `t P Q A U deltaR` at each probe, `deltaR` the
  change in radius from the diastolic radius;
* `vessels/<vessel>.dat`: `t Pin Pout dP Qin Qout volume` of each vessel;
* `nodes/<node>.dat`: the circuit of each `lumped0D` node, as
  `circulationFoam` writes it;
* `balance.dat`: `t volume netInflow drift`, the volume balance, which stays at
  round-off.

At each write time the profiles `x A Q P U` along every vessel are written to
`profiles/<time>/<vessel>.dat`, and the state, at full precision, to
`<time>/circulation1DState`. If this file exists in the start time, the solver
reads it, so a run can be continued.

---

## Building

The solver is built by `applications/Allwmake` along with the rest of
solids4foam. It requires `libcirculation1D` and `libcirculation0D`, built by
`src/Allwmake`. On other OpenFOAM flavours, the build is skipped. To build it
alone, after sourcing OpenFOAM and building `src/circulationModels`:

```bash
cd applications/solvers/circulation1DFoam
./Allwmake
```
