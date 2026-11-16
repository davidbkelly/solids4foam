# Paper plan: 3D–1D coupling in finite volume (version 2)

13 November 2026. Supersedes the first plan. Sources: the technical note of 6 November, the 32-paper literature map, and the decisions made since.

## The paper in one paragraph

A 3D domain coupled to a 1D arterial network behaves, as seen from the network, as an impedance element in a transmission line. That picture gives exact, unforced reference solutions for coupled 3D–1D problems, and a scalar model that predicts the stability and order of partitioned coupling. We implement the coupling in open-source finite volume (solids4foam), show it reaches the predicted order up to four with 1D sub-cycling, and package five cases as a benchmark suite with code-independent references.

Working title: "Exact benchmarks and a fourth-order finite-volume method for coupling 3D and 1D blood-flow models". Target: CMAME, with Journal of Computational Physics as the alternative.

## Research gaps: what the literature has, what we have, what remains

"Literature" means the 32 papers plus a short web check.

| Gap | Where the literature stands | What we have done | Still to do |
| --- | --- | --- | --- |
| 1. No verified open-source finite-volume 3D–1D coupling for blood flow | 18 of the 22 coupling papers use finite elements. Papadakis 2009 is the only FV blood-flow coupling: in-house, first-order characteristics, one test, no convergence study, interface errors 2.1–10.1%. Zhang 2023 is OpenFOAM but water hammer. Katsoudas 2025 is one-way. | 1D solver and coupler in solids4foam, on pimpleFoam and on FSI, explicit and implicit, serial and parallel. Five coupled cases run. | Case specifications and data release. |
| 2. No coupled convergence above second order | Best coupled result is about second order (Dobroserdova & Olshanskii 2013; Dobroserdova et al. 2016), in FE, with a rigid 3D wall, space and time refined together. The LifeV papers are first order. Müller & Toro 2013 reach fifth order in 1D only. | Fourth order in space and time on the rigid segment, by self-convergence. Order min(p, q) for coupling order p and time-scheme order q. | Measure against the exact solution. Show a 3D spatial order on the Womersley case. |
| 3. Temporal order of partitioned coupling never measured, including sub-cycling | Staggered schemes lag interface data by a step (Formaggia 2007, Papadakis 2009, Zhang 2023). Sub-cycling is routine (ratios of 25 to 100) and Malossi 2011 note the conditions then hold "just at the outer time step". Nobody reports the order that survives. | Explicit and implicit coupling of orders 1 to 4 with 3 to 200 sub-steps. Orders 1 to 4 observed on four rigid-domain cases; order 2 on the FSI case. | Trace the 3 × 10⁻⁶ floor on the network case. |
| 4. No predictive stability and order model for the partitioned coupling | Energy estimates for the interface conditions (Formaggia 2007; Dobroserdova & Olshanskii 2013). Iteration counts are reported empirically (Malossi 2011, 2013; Blanco 2013). | Scalar model in Δt/τ. It predicts explicit stability limits and implicit iteration rates to three or four digits against 16 3D runs. | Nothing required. Optional: extend it to an impedance interface condition. |
| 5. No exact unforced reference for a coupled problem, and no standard benchmark | Two manufactured solutions with forcing terms. Otherwise each paper compares with a longer 3D run or a pure 1D run. The 2026 review: "The absence of standardised benchmarks for mixed-dimensional coupling is itself a gap in the field." | Three exact references: closed form (rigid segment), Womersley impedance (viscous segment), delay equation (resonator). | Mesh-independent specifications, reference data, results template. |
| 6. Interface reflections reported but not explained | Single numbers: about 2% (Blanco 2007), 4.88% upstream error (Papadakis 2009). Malossi 2013: refinement reduced the reflection for one wall-end condition and increased it for two. Passerini 2009 argue total reflection at a rigid pipe from a reflection coefficient. | Closed-form reflection of a rigid segment, Γ = iωτ/(1 + iωτ). Hard against soft coupling explained as a resonance being damped. FSI-tube error split into an impedance part and a wall-end part. | Elastic-support wall end and a refinement study on the FSI tube. |

## What we reproduce or borrow, and do not claim

| Element | Who did it | Our use of it |
| --- | --- | --- |
| Interface conditions (flow rate with mean pressure; total pressure) | Formaggia 2001, 2007; Blanco 2007, 2013; Malossi 2013 | Used as they stand. |
| High-order 1D finite volume | Müller & Toro 2013 (ADER, WENO, fifth order, handles jumps) | Ours is a different and simpler scheme: k-exact reconstruction with Runge–Kutta, fourth order, smooth flow only. A means, not a contribution. |
| 1D benchmark set | Boileau 2015 | Our 1D solver passes it (127 of 132 errors in range). One summary table. |
| Compliant tube test | Formaggia 2001, 2007; close to Papadakis 2009 and Malossi 2013 | Reproduced in FV as B4. |
| Resonator set-up, hard and soft coupling | Dobroserdova 2016 | Reproduced as B3, with a new exact reference. |
| Womersley's solution as a check | Malossi 2011 (3D–3D); Filonova 2020 (3D alone) | Used as a coupled 3D–1D reference in B2. |
| 3D region embedded in a 1D network | Urquiza 2006; Blanco 2007–2013; Passerini 2009; Malossi 2013 | Done on a network with in vitro data as B6. |

## Contributions, as they will be stated

1. Exact, unforced reference solutions for coupled 3D–1D flow, derived from the impedance of the 3D domain: inviscid, viscous and resonant cases.
2. A scalar model for the stability and order of partitioned 3D–1D coupling, confirmed against 3D runs.
3. An open-source finite-volume coupling with 1D sub-cycling that reaches order min(p, q) up to four in the 1D discretisation, at the interface and in time.
4. A quantitative account of interface reflections reported in the literature.
5. A five-case benchmark suite with specifications, reference data and expected results.

## The benchmark suite (fixed)

| ID | Case | Reference | New or not | Status |
| --- | --- | --- | --- | --- |
| B1 | Rigid inviscid segment between two vessels | Exact closed form | New reference and new use | Done; needs order against exact |
| B2 | Womersley segment, α = 3, 6, 12 | Exact | New as a coupled benchmark | Done; needs 3D mesh refinement |
| B3 | Open-ended rigid pipe resonator | Exact delay equation | Published set-up, new reference | Done; compress; document limiter behaviour on the step |
| B4 | Compliant FSI tube continued by 1D | Longer 3D FSI tube | Established test, first open-source FV | Done; upgrade wall end, returning wave, refinement |
| B6 | 3D segment in the 37-artery network | In vitro data; exact inertance limit | Demonstration | Done; add static against total pressure, cost |

B5 (wave-carrying segment with matched impedance) stays out of the fixed suite. See "Optional extensions".

## Remaining work

This is what the Claude Code prompt runs.

**Phase 1, required.** Fourth order against the exact solution on B1 through a linearised 1D mode. Time order against exact on B2 and B3. 3D spatial order on B2. The network floor. Static against total pressure on B6. Cost table. Benchmark specifications and reference data.

**Phase 2, strengthens.** B4: elastic-support wall end, returning-wave variant, refinement of the reflection. B3: limiter behaviour on the step input.

**Then** a LaTeX draft of the benchmark section from verified numbers only.

### Optional extensions, with the rule for deciding

- **Impedance interface condition.** Go if Phases 1 and 2 finish without a surprise that changes the claims. It adds a method contribution and uses the model we already have. The paper stands without it.
- **B5.** Go only if a referee-proof answer is wanted to "the 3D solver is never tested along the wave direction" and there is time for a compressible solver variant. Otherwise state the limitation in the discussion.

## What we will and will not claim

| We will claim | We will not claim |
| --- | --- |
| Fourth order in the 1D discretisation, at the interface and in time, with the 3D solution exact in B1 | "Fourth-order 3D–1D coupling" without that qualification |
| Second-order 3D spatial convergence, shown on B2 | High order in 3D space |
| Second order on the FSI case, set by the FSI solver | Orders three and four with a compliant 3D wall |
| A closed form for the reflection of a rigid segment | That we discovered that rigid segments reflect |
| Exact references for the linear regime; one nonlinear check on B6 | Verification in the nonlinear regime |
| New relative to the literature we surveyed | "First" without having read Formaggia, Quarteroni & Vergara 2013 and the Businaro thesis |

## Outline

| Section | Content | Pages |
| --- | --- | --- |
| 1 Introduction | Gaps 1 to 6, the impedance idea, contributions | 2 |
| 2 Models | 1D equations and transmission-line limit; 3D models; interface conditions | 3 |
| 3 Numerical method | 1D finite volume; BDF in 3D; coupler with sub-cycling | 4 |
| 4 Analysis | Scalar model: stability, iteration rate, order | 3 |
| 5 Exact benchmarks | B1, B2, B3 | 9 |
| 6 Compliant benchmark | B4 | 4 |
| 7 Physiological case | B6, cost | 3 |
| 8 Discussion | Reflections in the literature explained; limits | 2 |
| Appendices | Closed forms; moving-mesh BDF; 1D verification summary | 4 |
| Supplementary | Boileau results, full tables, cases and data with a DOI | — |

## Items only you can settle

1. Whether to state that the printed solution in Dobroserdova et al. 2016 is inconsistent, after checking the journal version.
2. Whether to invite an outside group to run B1 to B3 with a finite-element code, before or after submission.
3. Go or no-go on the two optional extensions once Phase 2 is in.
