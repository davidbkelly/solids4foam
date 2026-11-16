# Draft: the benchmark section of the 3D-1D coupling paper

Draft of 16 November 2026 (David Kelly with Claude Code), to be edited
into the manuscript. Numbers are those of `PAPER_RESULTS.md` at the
commit stated there; where a second-round item will change a number it is
marked **[R-n]** with the item that changes it. Figure files are named
where they exist; captions are drafts.

## 1 A ladder of test cases for a 3D-1D coupling

Verifying a 3D-1D coupling has two difficulties that the usual FSI or
CFD verification does not: the exact solution of the coupled problem is
rarely available, because one side is a reduced model whose error is not
the coupling's; and the quantity passed across the interface, a mean
pressure for a flow, is the only quantity whose convergence is at stake.
We therefore use a ladder of five cases, each isolating one ingredient,
with the reference stated for each (Table 1), and we specify them so that
another code can run them (Appendix / supplement: `BENCHMARKS3D1D.md` with
CSV references).

**Table 1.** The test cases.

| case | 3D side | 1D side | reference | what it isolates |
|---|---|---|---|---|
| B1 rigid segment | rigid pipe, slip walls, 8 cm, R = 1 cm | two 1 m vessels, c0 = 6.17 m/s, inviscid; incident Gaussian pulse | exact: the pipe is the series inertance L = rho l/A, Gamma = i w tau/(1 + i w tau), tau = l/(2c) = 6.48 ms; closed form in time | the coupling alone: plug flow is exact in space, so every error is the coupling's or the 1D scheme's |
| B2 Womersley segment | rigid pipe, no-slip, 10 cm, alpha = 3, 6, 12 at 1 Hz | as B1, incident sine | exact: Womersley's series impedance; periodic state and the transient from rest (FFT convolution) | the 3D viscous solution's spatial error and the exchange of a mean pressure when the profile is not uniform |
| B3 resonator (Dobroserdova et al. 2016) | rigid pipe, 5 cm, open end | one 5 cm vessel between a flow source and the pipe, c0 = 3.5, 7, 10.5 m/s | exact: delay equations of the lossless line with the inductive load | a published 3D-1D case; the time order of the coupling against an exact solution; the limiter behaviour on a step input |
| B4 compliant tube (Formaggia et al. 2007) | FSI tube of the solids4foam 3dTube tutorial, 5 cm, pressure pulse | 10 cm vessel matched in static compliance | numerical: the same FSI tube three times as long | the physical case: a 1D vessel continuing an FSI segment; what the 1D model misses |
| B6 network (Boileau et al. 2015) | rigid no-slip cone of the thoracic aorta II | the other 36 vessels of the 37-artery in vitro network | numerical: the all-1D network, its exact inviscid rigid limit, the measurements | a 3D segment in a validated physiological network: cost, the effect of the segment on the network, the time order in a network |

Two conventions turned out to be essential and are stated in the
specification:

1. a finite-volume code's cell averages are compared with cell averages
   of the closed form (against point values the comparison is second
   order, (dx/(c0 sigma))^2/24, whatever the scheme);
2. the closed form uses the cross-section the 3D mesh has (a 32-gon
   has 0.9936 pi R^2; a 5-degree wedge with a flat wall 0.9987 pi R^2):
   against pi R^2 the comparison floors at the area deficit times
   |Gamma| (6e-5 for the wedge, 8e-3 for the resonator over 21 periods).
   **[R4: if the meshes are scaled to pi R^2 this convention becomes a
   requirement on the mesh instead, and every number below changes in
   its third digit.]**

## 2 The coupling and its scalar model

The 3D side imposes the 1D vessel end's pressure on its port patch and
returns the patch flux; the 1D network sees each port's flow over the
3D step as a polynomial in time (degree p - 1 for coupling order p) and
is re-integrated from the start of the step at each outer iteration
(implicit) or once per step from extrapolated flows (explicit). Within
one step each vessel end answers the port flow as a resistor of its
characteristic impedance Zc = rho c/A, so a 3D domain of inertance L
between two ports sees Pin - Pout = 2 P_inc - 2 Zc Q, and the coupling
number is dt/tau, tau = L/(2 Zc). The model gives the explicit stability
limits (Table 2) and the implicit contraction dt/(a0 tau), both verified
to four digits on a 1D surrogate of the 3D side (`Test-circulation1DPorts`).

**Table 2.** Explicit stability limits dt/tau of the Dirichlet interface.

| explicit order p | 1 | 2 | 3 | 4 |
|---|---|---|---|---|
| Euler in the 3D | 2 | 0.667 | 0.286 | 0.133 |
| BDF2 | 4 | 1.333 | 0.571 | 0.267 |
| BDF3 | 4.5 | 2.222 | 0.952 | 0.444 |
| BDF4 | 1.643 | 2.25 | 1.524 | 0.711 |

**The impedance (Robin) interface [R7].** Imposing instead p = P(Q*) -/+
Zc (Q - Q*) implicitly in the pressure equation makes the iterate Q*
drop out of the model: the implicit iteration converges in one exchange
(surrogate: one exchange per step reproduces the converged Dirichlet
iteration to 3e-9 of the peak) and the explicit scheme loses its dt/tau
limit (surrogate: orders 1-4 stable to dt/tau = 16; 3D rigid segment:
explicit order 2 stable at dt/tau = 2.5 and 4.9 where the Dirichlet
interface diverges). In a segregated 3D solver the mixed condition
slows the pressure-velocity iteration itself, so it does not save outer
iterations there; the same response applied as an Aitken relaxation of
the Dirichlet value does **[R7: result pending]**.

## 3 Results

### 3.1 B1: fourth order against the closed form

Linearised vessels (the 1D equations about the initial state, so the
closed form is exact at any amplitude), kExact degree 3 + RK4 in the
vessels, BDF4 in the pipe, dx and dt halved together over five levels,
errors relative to the peaks, profiles as cell averages
(`rigidSegment3D/highOrder_linear.md`, Fig. `plots/highOrder3D1D_linear.png`):

| coupling | level (dx, dt) | 3D flux | order | pipe-end pressure | order | vessel profiles | order |
|---|---|---|---|---|---|---|---|
| implicit 4 | 2.5 mm, 2 ms | 1.43e-03 | | 1.10e-03 | | 1.75e-03 | |
| | 1.25 mm, 1 ms | 9.41e-05 | 3.93 | 7.30e-05 | 3.92 | 1.14e-04 | 3.95 |
| | 0.625 mm, 0.5 ms | 5.99e-06 | 3.97 | 4.62e-06 | 3.98 | 7.14e-06 | 3.99 |
| | 0.3125 mm, 0.25 ms | 3.77e-07 | 3.99 | 2.89e-07 | 4.00 | 4.46e-07 | 4.00 |
| | 0.15625 mm, 0.125 ms | 2.36e-08 | 4.00 | 1.80e-08 | 4.00 | 2.82e-08 | 3.98 |
| explicit 4 | five levels | 4.11e-03 to 7.05e-08 | 3.88 to 3.99 | 4.16e-03 to 6.91e-08 | 3.89 to 4.00 | 2.96e-03 to 4.65e-08 | 3.96 to 4.00 |
| implicit 2 | five levels | as implicit 4 | 3.93 to 4.00 | as implicit 4 | 3.92 to 4.00 | 2.69e-03 to 1.52e-05 | 1.61, 1.89, 1.97, 1.99 |
| implicit 1 | five levels | as implicit 4 | 3.93 to 4.00 | as implicit 4 | 3.92 to 4.00 | 6.14e-02 to 3.90e-03 | 0.98, 1.00, 1.00, 1.00 |

Fourth order with no floor to 2e-8; coupling of order p gives order p on
the waves the ports send into the vessels. At physiological amplitude
(nonlinear vessels) the self-convergence study floors at 2e-5 for the
low-order couplings: the erroneous reflected wave interacts with the
incident one at second order in amplitude, a nonlinearity of the 1D
equations and not of the coupling. **[R4 changes the numbers in the
third digit.]**

### 3.2 Time order against exact solutions

B3 with a smooth inflow (a 15 ms Gaussian pulse; the published sine and
step have jumps in dQ/dt that floor every scheme at the 1D's resolution
of the jump) **[R3]**, RMS error / RMS of the signal of the port flow:

| dt [ms] | Euler | backward (BDF2) | BDF3 | BDF4 |
|---|---|---|---|---|
| 2 | 1.59e-1 | 1.73e-2 | 1.59e-3 | 1.55e-4 |
| 1 | 1.55e-1 | 4.16e-3 | 1.94e-4 | 9.02e-6 |
| 0.5 | 8.14e-2 | 9.26e-4 | 2.20e-5 | 5.25e-7 |
| 0.25 | 4.12e-2 | 2.31e-4 | 2.74e-6 | 3.31e-8 |
| 0.125 | 2.07e-2 | 5.77e-5 | 3.44e-7 | 2.53e-9 |
| order (finest pairs) | 0.99 | 2.00 | 3.00 | 3.99, 3.71 |

B2 (alpha = 6, linearised vessels, exact transient): Euler at 0.99 to
the finest step, backward at 1.65 before the 50 x 100 wedge's spatial
floor of 5.6e-6, which BDF3 and BDF4 meet from the coarsest step
(self-convergence of the same runs 2.86 and 3.32).

B6 (the network, cylinder segment) **[R9]**, self-convergence of the
port flow / pressure:

| dt triple [ms] | Euler | backward | BDF3 | BDF4 |
|---|---|---|---|---|
| 2 / 1 / 0.5 | 0.95 / 0.87 | 1.99 / 2.01 | 2.99 / 2.99 | 3.83 / 3.72 |
| 1 / 0.5 / 0.25 | 0.98 / 0.94 | 2.00 / 2.01 | 3.02 / 3.02 | 3.82 / 3.91 |

On the cone (the tutorial's tapered segment) BDF3 and BDF4 floor at 3e-6
of the signal: pimpleFoam on a non-orthogonal mesh (its flux correction
and the non-orthogonal treatment) does not go beyond second order in
time, whatever the ddt scheme; the coupling takes the 3D scheme's order
when the 3D discretisation has it.

### 3.3 B2: spatial order of the 3D side and the exchange of a mean pressure

Periodic state, alpha = 6 and 12, against the pipe of the wedge's
cross-section: amplitude errors 3e-6, phases 0.0002 to 0.0004 degrees
(alpha 3: 1.2e-5). Order in space on the uniform mesh family **[R2:
table pending; the two coarser levels give 2.01 on the ports]**. The
coupling is exact for a 3D side whose profile is not uniform: the
exchange of a mean pressure and a flow loses nothing when the 3D flow
is a genuine viscous solution.

### 3.4 B3: the resonator as published, and the step input

The 3D run against the delay equations with the pipe's patch area: RMS
flow errors 7.7e-3 (sine) and 9.3e-3 (step) of the inflow over 21
resonant periods, peaks 1.596 against 1.596 and 2.263 against 2.265;
the remaining error is BDF2's phase error of the resonance. Step input:
no scheme overshoots the exact peak; the limiters show as a deficit of
total variation (14 % minmod, 5 % MC at 100 cells); the unlimited
k-exact reconstruction has a quarter to a sixth of the RMS error. The
jumps are first order in dx for every scheme.

### 3.5 B4: the compliant tube

Coupled 5 cm FSI tube + 1D vessel against the 15 cm FSI tube, static
compliance matched (c0 = 5.25 m/s), 3 ms pulse. The clamped cut end of
the wall puts the coupled tube 3 to 7 % of the pulse (RMS) from the full
tube at the sections; the free end 2 to 4 %; an elastic support with the
continuing shell's static edge stiffness (1.04e7 Pa/m, a sixth of the
hoop stiffness) behaves as the free end. The pulse returns through the
interface from the 1D side (closed far end) within 1-2 ms of the full
tube at every section, with the 1D model's amplitude (844-1124 Pa
against 475-680 Pa): the mismatch is the thick wall's dispersion, counted
twice. **[R1: which solids4foam build; R5: mesh refinement per end
condition; R6: the pulse-length sweep, R/(c0 T) from 0.3 to 0.04, both
pending.]**

### 3.6 B6: the network

The rigid segment changes the network by 4 % of the pulse against its
exact inviscid rigid limit; static against total pressure at the ports
moves the network by 0.06 to 0.2 % of the pulse and the segment's own
pressure by 0.3 % (0.9 % at peak flow); the energy-flux form of the
dynamic pressure is singular at flow reversal and not usable. Cost: the
1D network is never the cost (37 vessels at 2 mm cells, 0.5 ms per CFL
step); the implicit coupling's price is the outer iterations, 5 per step
here, 11 on the rigid segment at dt/tau = 0.077, 21 on the viscous
wedge where the segregated algorithm sets the count.

## 4 Limitations (one sentence to add)

The 3D solver is not tested along the direction of wave travel: in the
rigid and viscous segments the pressure wave crosses the 3D part
instantaneously (incompressible) and in the FSI segment the wall
carries it; a 3D segment that itself propagates the wave at the
vessels' speed (a compressible fluid of matched impedance, or a long
compliant segment resolved in 3D) is left for later work.

## 5 Reproducibility

All cases are tutorials of solids4foam (`tutorials/circulationModels/
coupling3D1D/`), each with an `Allrun`, an `analyse.py` writing its
`results.md`, and the study scripts named in `PAPER_RESULTS.md`; the
exact references are written as CSV by `benchmarks3D1D.py csv` and the
cases' `analyse.py --csv`, and the numerical references of the
physiological-amplitude variants are committed with their headers
stating the runs that made them. Code state: `feature/circulation1D` at
the commit stated in `PAPER_RESULTS.md`, OpenFOAM v2412, plus the BDF
time scheme of solids4foam PR #674 **[to be committed]**.
