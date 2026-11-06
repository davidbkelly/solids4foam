# circulation1D — implementation spec for solids4foam

Spec for a 1D arterial blood-flow library (`src/circulationModels/circulation1D`), a standalone solver (`applications/solvers/circulation1DFoam`) and tutorials reproducing the six benchmark cases of Boileau et al. (2015), *A benchmark study of numerical schemes for one-dimensional arterial blood flow modelling*, IJNMBE 31(10):e02732. The paper and its supporting data are already in `src/circulationModels/circulation1D/`.

The library will later be coupled to 3D CFD and FSI through `circulationCoupling`. The coupling itself is **out of scope** for now, but the API in §8 must be designed in from the start.

---

## 0. Ground rules

- **Do not modify `circulationFoam`.** It stays a pure 0D ODE/LPM solver. The 1D solver is a new application.
- **Before writing code**, read `src/circulationModels/circulation0D` and `src/circulationModels/circulationCoupling` and follow their conventions:
  - class layout and runtime selection (`declareRunTimeSelectionTable` / `addToRunTimeSelectionTable`);
  - dictionary names and how dictionaries are read;
  - `Make/files` and `Make/options`, and `Allwmake` integration;
  - the OpenFOAM-version `#ifdef` handling used elsewhere in solids4foam;
  - coding style and file headers.
- **Reuse `circulation0D` where it fits.** In particular, its ODE integrators and its 0D network objects for terminal boundary conditions (§5.4). If reuse is awkward, implement a self-contained Windkessel boundary condition first and leave a TODO.
- **Use SI units internally.** Convert input data on read. ADAN56 tables use cgs: 1 dyn·s/cm⁵ = 1e5 Pa·s/m³ and 1 cm⁵/dyn = 1e-5 m³/Pa.
- **No hard-coded case parameters.** Everything comes from dictionaries.
- Inspect the benchmark data files before designing the input format, and check their units, column layout and periodicity.

---

## 1. Governing equations

### 1.1 Per-vessel system

Per vessel, with axial coordinate $x \in [0, L]$ and conservative variables $\mathbf U = (A, Q)$:

$$
\partial_t \begin{pmatrix} A \\ Q \end{pmatrix}
+ \partial_x \begin{pmatrix} Q \\ \alpha \dfrac{Q^2}{A} + \dfrac{\beta}{3\rho A_d} A^{3/2} \end{pmatrix}
= \begin{pmatrix} 0 \\ -K_R \dfrac{Q}{A} + S_{\text{geom}} \end{pmatrix}
$$

| Symbol | Definition |
|---|---|
| $A$ | lumen area |
| $Q = AU$ | volumetric flow rate |
| $\rho$, $\mu$, $\nu = \mu/\rho$ | density, dynamic and kinematic viscosity |
| $\zeta$ | velocity-profile order (Boileau Eq. 2) |
| $K_R = 2\pi(\zeta+2)\nu$ | friction coefficient ($\zeta=2 \Rightarrow 8\pi\nu$, $\zeta=9 \Rightarrow 22\pi\nu$) |
| $\alpha$ | momentum-flux (Coriolis) coefficient |
| $S_{\text{geom}}$ | zero for uniform $\beta$, $A_d$ (see §1.4) |

**About $\alpha$:** the profile gives $\alpha = (\zeta+2)/(\zeta+1)$, but Boileau's Eq. (1) is written with $U\,\partial_x U$, which is equivalent to $\alpha = 1$. Make `alpha` a dictionary entry with **default 1**. All benchmark tutorials use 1.

### 1.2 Tube law (`sqrtElastic`, Boileau Eq. 3)

$$
P = P_{ext} + P_d + \frac{\beta}{A_d}\left(\sqrt A - \sqrt{A_d}\right), \qquad \beta = \tfrac{4}{3}\sqrt{\pi}\, E h
$$

Make the tube law runtime-selectable (abstract `tubeLaw1D` plus `sqrtElastic`) so viscoelastic laws can be added later. Each tube law must provide:

- $P(A)$
- $\partial P/\partial A = \beta/(2A_d\sqrt A)$
- $c(A) = \sqrt{\dfrac{A}{\rho}\dfrac{\partial P}{\partial A}} = \sqrt{\dfrac{\beta}{2\rho A_d}}\,A^{1/4}$
- $A(P)$, the inverse
- the flux pressure term $\dfrac{\beta}{3\rho A_d}A^{3/2}$

Zero-pressure area (Boileau Eq. 16): $A_0 = A_d\left(1 - P_d\sqrt{A_d}/\beta\right)^2$.

### 1.3 Characteristics and invariants

Eigenvalues: $\lambda_\pm = \alpha U \pm \sqrt{c^2 + \alpha(\alpha-1)U^2}$, which reduces to $U \pm c$ for $\alpha = 1$. The flow is subcritical, so one boundary condition is needed at each vessel end.

Riemann invariants for $\alpha = 1$ and locally uniform properties:

$$W_f = U + 4c, \qquad W_b = U - 4c$$

Inverse: $U = (W_f + W_b)/2$, $c = (W_f - W_b)/8$, then $A$ from $c$ via $A = \left(c\,\sqrt{2\rho A_d/\beta}\right)^4$.

Useful identity: $\dfrac{dc}{dA} = \dfrac{c}{4A}$.

### 1.4 Variable properties (phase 3 only)

When $\beta(x)$ and $A_d(x)$ vary along the vessel (tapering in ADAN56), write $k(x) = \beta/A_d$ and keep the flux above with local $k$. The geometric source is then

$$
S_{\text{geom}} = \frac{1}{\rho}\Big[\partial_x k\,\big(A\sqrt{A_d} - \tfrac23 A^{3/2}\big) + A\,k\,\partial_x\sqrt{A_d} - A\,\partial_x(P_{ext}+P_d)\Big].
$$

This exactly balances the flux gradient for the state at rest ($A = A_d$, $Q = 0$) at the continuous level. Verify this derivation before relying on it.

A naive discretisation is **not** well-balanced. Implement a well-balanced variant: hydrostatic-reconstruction type, following Müller, Parés & Toro (JCP 2013), or treat property changes as two-vessel junctions between piecewise-uniform cells. Test it with the fluid-at-rest problem (§7.2).

---

## 2. Spatial and temporal discretisation

- **Mesh.** Each vessel is a uniform 1D cell array (`nCells` or `dx` from the dictionary). Do **not** use an OpenFOAM `fvMesh`: a network of many vessels joined at junctions does not map onto one mesh. Use OpenFOAM containers (`scalarField`, `List`, `dictionary`, `Function1`) for consistency. If, after reading the existing code, there is a strong convention argument for `fvMesh`, flag it rather than silently choosing.
- **Reconstruction.** MUSCL in $(A, Q)$, with a selectable limiter: `minmod`, `vanLeer`, `MC`, and `none` (unlimited central slopes). `none` is needed for clean second-order convergence tests on smooth problems, because limiters clip at extrema.
- **Riemann flux.** Selectable: `HLL` (default) and `Rusanov`. For HLL, $S_L = \min(\lambda_-^L, \lambda_-^R)$ and $S_R = \max(\lambda_+^L, \lambda_+^R)$.
- **Time integration.** SSP-RK2 (Heun) by default, SSP-RK3 optional. Boundary and junction states are recomputed at **every stage**. Friction and $S_{\text{geom}}$ enter as explicit sources in each stage. Warn if $K_R\,\Delta t / A > 0.1$ anywhere.
- **Time step.** $\Delta t = \text{CFL}\,\min_{\text{vessels, cells}} \Delta x / \max|\lambda_\pm|$, global across the network, with CFL from the dictionary (default 0.9).
- **Boundary fluxes.** At a vessel end, the face flux is $\mathbf F(A_b, Q_b)$ evaluated at the boundary state from §5. The outgoing invariant comes from the reconstructed face value of the adjacent cell at the current stage.
- **ODE states inside boundary conditions** (Windkessel $P_c$, any 0D network states) are advanced with the **same RK stages** as the PDE, by appending them to the stage update, so the coupled 1D–0D system keeps its formal order.

---

## 3. Library structure (proposed; adapt to existing conventions)

```
src/circulationModels/circulation1D/
    Make/{files,options}
    tubeLaws/            tubeLaw1D (abstract, runtime-selectable), sqrtElastic
    vessel1D/            geometry, properties, cell fields A, Q, reconstruction, flux, sources
    numerics/            limiters, Riemann fluxes (runtime-selectable)
    nodes1D/             node1D (abstract) and concrete node types, see §5
    network1D/           owns vessels + nodes; stepping; probes; output; coupling API (§8)
    benchmarkData/       existing paper + data (consider moving reference data into tutorials)
```

A **node** connects one or more vessel ends: `(vesselName, proximal|distal)`.

- Nodes with one connection are boundary conditions.
- Nodes with two or more connections are junctions. Two-connection junctions occur in ADAN56, for example at the split segments 7a/7b/7c.

---

## 4. Input format (example; finalise after reading existing dictionaries)

`constant/circulation1DProperties`:

```
fluid      { rho 1060; mu 4e-3; profileOrder 2; alpha 1; }
numerics   { CFL 0.9; flux HLL; limiter minmod; timeScheme SSPRK2; }

vessels
{
    carotid
    {
        length            0.126;
        nCells            126;
        radiusProximal    3e-3;     // reference (diastolic) radius
        radiusDistal      3e-3;     // optional; defaults to radiusProximal
        wallThickness     0.3e-3;   // or wallThicknessModel ADAN56 (Boileau Eq. 17)
        youngsModulus     700e3;
        diastolicPressure 10933;
        externalPressure  0;
        tubeLaw           sqrtElastic;
        initialPressure   0;        // initial A from tube law; Q = 0
    }
}

nodes
{
    inlet  { type flowInlet;     connections ((carotid proximal)); flowRate table/csvFile/coded Function1; period 1.1; }
    outlet { type windkesselRCR; connections ((carotid distal)); R1 2.4875e8; C 1.7529e-10; R2 1.8697e9; Pout 0; }
}

probes  { carotidMid { vessel carotid; x 0.063; } }   // also: "midpoint", "proximal", "distal"
```

`system/controlDict` controls `endTime` and output. `deltaT` in `controlDict` is the **output/coupling interval**. The solver subcycles internally with the CFL step and lands exactly on each `controlDict` time (§8).

---

## 5. Boundary and junction algorithms

Notation:

- At a **proximal** end the outgoing invariant is $W_b$ and $U = W_b + 4c(A)$.
- At a **distal** end the outgoing invariant is $W_f$ and $U = W_f - 4c(A)$.

Each case reduces to a scalar (or small) Newton solve in $A$, with the adjacent cell's $A$ as the initial guess and a relative tolerance of about 1e-12.

### 5.1 `flowInlet` (prescribed $Q_{in}(t)$, Boileau Eq. 6)

Solve $f(A) = A\,(W_b + 4c(A)) - Q_{in} = 0$, with $f'(A) = W_b + 5c$.

Read $Q_{in}$ through `Function1` so tables, CSV and analytic expressions all work. The single-pulse case uses
$Q_{in} = 10^{-6}\exp(-10^4 (t - 0.05)^2)$ m³/s.
The physiological inflows are periodic tables from the supporting data, so provide periodic extension and interpolation (linear, with cubic optional).

### 5.2 `absorbing` (non-reflecting outlet)

Hold the incoming invariant at its initial value, e.g. $W_b = W_b^0 = -4c(A_0)$ at a distal end. Then
$U = (W_f + W_b^0)/2$, $c = (W_f - W_b^0)/8$, and $A$ follows from $c$. This is closed form; no Newton solve is needed.

### 5.3 `resistance` (37-artery terminals)

Solve $g(A) = P(A) - P_{out} - R\,A\,(W_f - 4c(A)) = 0$, with
$g'(A) = \partial P/\partial A - R\,(W_f - 5c)$.

### 5.4 `windkesselRCR` (Boileau Eq. 9, Fig. 1)

State: capacitor pressure $P_c$, with $C\,\dot P_c = Q_b - (P_c - P_{out})/R_2$.

At each stage, using the current $P_c$, solve $P(A) - P_c - R_1 A\,(W_f - 4c(A)) = 0$ (same derivative as §5.3 with $R = R_1$). Then advance $P_c$ in the RK stage with $Q_b$.

Preferred design: a generic `lumped0D` node that wraps any `circulation0D` network, using flow in and pressure out as in the 3D–0D framework, with RCR as one instance.

### 5.5 `junction` (N ≥ 2 vessel ends, Boileau Eqs. 7–8)

For connection $i$:

- $s_i = +1$ if the junction is at the vessel's distal end (outgoing invariant $W_f$);
- $s_i = -1$ if it is at the proximal end (outgoing invariant $W_b$).

Then $U_i = W_i^{out} - 4 s_i\, c_i(A_i)$, which eliminates all $U_i$. Solve for $(A_1, \dots, A_N)$ with Newton (analytic Jacobian):

- mass: $\sum_i s_i A_i U_i = 0$;
- total pressure: $P_1 + \tfrac12\rho U_1^2 - P_i - \tfrac12\rho U_i^2 = 0$ for $i = 2, \dots, N$.

Check: in the linear limit, the reflection coefficient for a parent with daughters is $(Y_0 - \sum Y_i)/(Y_0 + \sum Y_i)$, with $Y = A/(\rho c)$.

### 5.6 `coupledFlow` (stub now, used later for 3D–1D)

$Q_b(t)$ is supplied externally as a callable over $[t^n, t^{n+1}]$, solved as in §5.1. The node returns $P(A_b)$. See §8.

---

## 6. Solver: `applications/solvers/circulation1DFoam`

Minimal structure, mirroring how `circulationFoam` creates `Time` and reads dictionaries (no mesh needed):

1. construct `network1D` from `constant/circulation1DProperties`;
2. loop over `controlDict` times, calling `network.advance(t^n, t^{n+1})` with internal CFL subcycling;
3. write probe time series and optional spatial snapshots;
4. print the global volume-balance diagnostic (§7.4).

### Output

- Probe files, one per probe: columns `t, P, Q, A, U, deltaR`, where $\Delta r = \sqrt{A/\pi} - \sqrt{A_d/\pi}$ at that point.
- Per-vessel inlet–outlet $\Delta P$ (needed for Boileau Figs. 3–4).
- Optional spatial profiles $P(x)$ along a named vessel path at given times, for ADAN56 Figs. 15–16 (systolic, mean and diastolic pressure along the aorta).

---

## 7. Verification and benchmarks

### 7.1 Tutorials (`tutorials/circulation1D/`, matching the existing tutorial layout)

Each case has `Allrun`, `Allclean`, `constant/`, `system/`, `referenceData/` (copied from the Boileau supporting data) and a plotting script.

| Tutorial | Boileau § | Reference | Notes |
|---|---|---|---|
| `singlePulseInviscid` | 4.1 | linear theory | L = 10 m, r = 1 cm, E = 400 kPa, h = 1.5 mm, ρ = 1050, ζ = 9, μ = 0; absorbing outlet; $c_0 = 6.17$ m/s |
| `singlePulseViscous` | 4.1 | Eq. 15 | μ = 4 mPa·s; peak decays as $\exp(-(\zeta+2)\pi\mu x/(\rho c_0 A_0))$, about 0.068 m⁻¹ |
| `commonCarotid` | 4.2 | 3D (Xiao 2014) | Table III, ζ = 2, RCR |
| `upperThoracicAorta` | 4.3 | 3D | Table V, ζ = 9, RCR |
| `aorticBifurcation` | 4.4 | 3D | Table VII; 3-vessel junction; probes at aorta midpoint, aorta end, iliac midpoint |
| `thirtySevenArtery` | 4.5 | in vitro | Table IX plus per-segment data from the supporting files; ρ = 1050, μ = 2.5 mPa·s, E = 1.2 MPa, $P_d = 0$, $P_{out} = 432.6$ Pa; resistance terminals |
| `ADAN56` | 4.6 | other schemes | Tables XI–XIII; tapering; h from Eq. 17 (constants in cm units); $P_d = P_0 = 10$ kPa; initial state $(A_0, 0, P_0)$; **phase 3** |

A shared script `tutorials/circulation1D/common/boileauMetrics.py` should:

- interpolate solids4foam output onto the reference time points (1 kHz);
- compute $\mathcal E^{RMS}$, $\mathcal E^{MAX}$, $\mathcal E^{SYS}$ and $\mathcal E^{DIAS}$ exactly as in Boileau Eqs. 10–13. Pressure is normalised pointwise; $Q$, $\Delta P$ and $\Delta r$ are normalised by the maximum reference value;
- evaluate over the **last cycle** once periodic (stop when the cycle-to-cycle change in probe RMS is below 1e-4);
- reproduce Figs. 2–16 with solids4foam overlaid on the six published schemes, and Tables IV, VI, VIII and X with an added solids4foam column.

### 7.2 Extra verification (beyond Boileau)

1. **Order of accuracy.** Run the single pulse inviscid with limiter `none` and refine Δx at fixed CFL. Measure error against $p = Z_0 Q_{in}(t - x/c_0)$ with $Z_0 = \rho c_0/A_0$; the amplitude is about 20 Pa, so the problem is effectively linear. Expect a slope of about 2.
2. **Fluid at rest.** A tapered vessel with $Q = 0$ and $A = A_d(x)$ must stay at rest to machine precision. This gates phase 3.
3. **Junction reflection.** Send a linear pulse through a 1→2 junction and compare the reflected amplitude with $(Y_0 - \sum Y_i)/(Y_0 + \sum Y_i)$.
4. **Windkessel matching.** With $R_1 = Z_0$, a high-frequency pulse should leave with negligible reflection.
5. **Riemann problems** (optional): exact solutions from Toro & Siviglia, CiCP 2013.

### 7.3 Acceptance criteria

| Case | Criterion |
|---|---|
| Single pulse, inviscid | peak decay < 1% over 10 m at Boileau-like resolution |
| Single pulse, viscous | peak within 1% of Eq. 15 |
| Order test | observed order 1.9–2.1 |
| Carotid, aorta, bifurcation | all metrics within (or very close to) the range of the six schemes in Tables IV, VI and VIII |
| 37-artery | within the range of Table X |
| ADAN56 | waveforms visually coincident with the five consistent schemes (DCG, LCG, FEM, FVM, McC) |

### 7.4 Conservation diagnostic

Each step, check

$$\frac{d}{dt}\sum_{\text{vessels}}\int A\,dx = \sum_{\text{inlets}} Q - \sum_{\text{outlets}} Q$$

and report the cumulative drift. Junctions conserve mass exactly, so the drift should sit at round-off and solver-tolerance level.

---

## 8. API for future 3D–1D coupling (design now, use later)

`network1D` must expose:

- `advance(t0, t1)`: subcycles with CFL steps, with the final substep shortened to hit `t1` exactly.
- `saveState()` / `restoreState()`: a full snapshot of all vessel fields and node ODE states. **Required** for implicit coupling, where the 1D model is re-integrated from $t^n$ on every coupling iteration (as in the 3D–0D implicit I-variants).
- `coupledFlow` nodes accepting $Q(t)$ as a callable over $[t^n, t^{n+1}]$, so the polynomial interface reconstructions (E1–E4 / I1–I4) from the 3D–0D framework can be reused unchanged, each evaluated at every RK stage time.
- `interfacePressure(name)`: returns $P(A_b)$ at the end of the step.

The flow-to-pressure exchange mirrors the existing 3D–0D coupling: the 3D side supplies $Q$, and the 1D side returns $P$.

---

## 9. Phasing

| Phase | Scope | Tutorials |
|---|---|---|
| 1 | Tube law, single vessel, MUSCL + HLL + SSP-RK2; `flowInlet`, `absorbing`, `resistance`, `windkesselRCR`; `circulation1DFoam`; probes; metrics script | single pulse ×2, carotid, aorta, order test |
| 2 | Junctions; network input | aortic bifurcation, 37-artery, junction reflection test |
| 3 | Variable properties, well-balancing, Eq. 17 wall thickness | fluid-at-rest test, ADAN56 |
| 4 (later) | `coupledFlow` node wired into `circulationCoupling` | — |

After each phase, report back with results against the acceptance criteria before starting the next.
