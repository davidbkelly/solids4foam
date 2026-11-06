/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::network1D

Description
    A network of 1D vessels joined by nodes, read from
    constant/circulation1DProperties (SI throughout):

        fluid     { rho 1060; mu 4e-3; profileOrder 2; alpha 1; }
        numerics  { CFL 0.9; flux HLL; limiter minmod; timeScheme SSPRK2; }

        vesselDefaults  { youngsModulus 1.2e6; }    // optional: every vessel

        vessels
        {
            carotid { length 0.126; nCells 126; radius 3e-3; ... }
        }

        nodes
        {
            inlet   { type flowInlet; connection carotid.in;  ... }
            outlet  { type lumped0D;  connection carotid.out; circuit {...} }
        }

        probes
        {
            mid     { vessel carotid; x midpoint; }    // or in, out, metres
        }

    CFL defaults to 0.9, but a limited reconstruction (any limiter other than
    none) is TVD only for CFL <= 0.5 and warns above it.

    Higher order (uniform vessels): the k-exact reconstruction of degree k,
    for face values of order k + 1, with a time scheme of that order:

        numerics
        {
            reconstruction  kExact;     // default MUSCL (with the limiter)
            polynomialOrder 3;          // 1 to 3: order k + 1 in space
            timeScheme      RK4;        // SSPRK2, SSPRK3 or RK4
            CFL             0.5;
        }

    Entries of "vesselDefaults" apply to every vessel that does not give its
    own.

    The friction coefficient is KR = 2 pi (profileOrder + 2) mu/rho. alpha
    (default 1) must be 1 for now, as the nodes use the alpha = 1 Riemann
    invariants. Every vessel end must belong to exactly one node.

    advance(t0, t1) integrates with explicit Runge-Kutta steps of
    dt = CFL min(dx/max|lambda|) over the whole network, shortening the last
    steps to land on t1 exactly. Boundary states and node ODE states are
    updated in every stage.

    For coupling: saveState()/restoreState() keep and restore a full
    snapshot, so that a step can be repeated, and interfacePressure() is
    the pressure at the vessel end of a node. The coupledFlow nodes are the
    network's ports, with the interface of circulationNetwork's (portNames,
    portPatch, portInlet, setPortForcing, portPressure), and state() /
    setState() hold the whole state as one field, so that a coupler can keep
    the state at the start of a step and integrate from it again.

    The volume balance, V(t) - V(t0) - integral of the net inflow at the
    vessel ends, is integrated with the same stages and should stay at
    round-off.

Author
    David Kelly, UCD.

SourceFiles
    network1D.C

\*---------------------------------------------------------------------------*/

#ifndef network1D_H
#define network1D_H

#include "vessel1D.H"
#include "node1D.H"
#include "coupledFlow.H"
#include "rungeKutta1D.H"

namespace Foam
{

class network1D
{
    // Private Data

        //- Fluid
        scalar rho_;
        scalar mu_;
        scalar profileOrder_;
        scalar alpha_;

        //- Numerics
        scalar CFL_;
        autoPtr<slopeLimiter1D> limiter_;

        //- k-exact reconstruction, if selected (else MUSCL with limiter_)
        autoPtr<kExactReconstruction1D> kExact_;
        autoPtr<riemannFlux1D> flux_;
        rungeKutta1D rk_;

        PtrList<vessel1D> vessels_;
        PtrList<node1D> nodes_;

        //- Probes
        wordList probeNames_;
        labelList probeVessel_;
        scalarList probeX_;

        scalar t_;

        //- Volume balance: initial volume and the integral of the net inflow
        scalar V0_;
        scalar inflow_;
        scalar inflow0_;
        scalar dInflowdt_;

        //- Butcher-form Runge-Kutta: the stage derivatives of the inflow
        scalarList kInflow_;

        //- Step statistics
        label nSteps_;
        scalar lastDeltaT_;
        bool frictionWarned_;

        //- Snapshot kept by saveState()
        struct snapshot
        {
            scalar t;
            List<scalarField> A;
            List<scalarField> Q;
            List<scalarField> y;
            scalar inflow;
            label nSteps;
            bool valid = false;
        };

        snapshot saved_;


    // Private Member Functions

        //- Reconstruct every vessel (MUSCL or k-exact)
        void reconstructVessels();

        //- Right-hand side of a stage at time t: reconstruct, boundary
        //  states, vessel and node derivatives
        void evaluateStage(const scalar t);

        //- One Butcher-form Runge-Kutta step
        void butcherStep(const scalar dt);

        //- Stable time step of the current state
        scalar stableDeltaT();

        //- One Runge-Kutta step
        void step(const scalar dt);

        //- Boundary states at the current time, for output
        void updateBoundaries();

        label nodeLabel(const word& name) const;

        //- The coupledFlow node named port; fatal if there is none
        const coupledFlow& port(const word& name) const;
        coupledFlow& port(const word& name);


public:

    // Constructors

        //- Construct from circulation1DProperties, with time starting at t0
        network1D(const dictionary& dict, const scalar t0);

        network1D(const network1D&) = delete;
        void operator=(const network1D&) = delete;


    // Time stepping

        scalar time() const { return t_; }

        //- Integrate from t0 (the current time) to t1
        void advance(const scalar t0, const scalar t1);

        label nSteps() const { return nSteps_; }

        scalar lastDeltaT() const { return lastDeltaT_; }


    // Coupling

        //- Keep a full snapshot of the state (vessel fields, node states,
        //  time and volume balance)
        void saveState();

        //- Restore the snapshot kept by saveState()
        void restoreState();

        //- Pressure at the vessel end of a node at the current time [Pa]
        scalar interfacePressure(const word& node) const;

        //- Names of the coupledFlow nodes
        wordList portNames() const;

        //- CFD patch of a port (empty if none)
        const word& portPatch(const word& name) const;

        //- True if the port's patch is an inlet of the 3D domain
        bool portInlet(const word& name) const;

        //- Factor from the flow through the port's patch to the vessel flow
        scalar portFlowScale(const word& name) const;

        //- Flow of a port over [t0, t1], in alpha = (t - t0)/(t1 - t0)
        void setPortForcing
        (
            const word& name,
            const scalar t0,
            const scalar t1,
            const circulationForcing::coeffs& c
        );

        //- Pressure at a port at the current time [Pa]
        scalar portPressure(const word& name) const;

        //- Flow of a port's forcing at the current time [m3/s]
        scalar portFlow(const word& name) const;

        //- Size of state()
        label nEqns() const;

        //- The whole state: vessel A and Q, node states, net inflow
        scalarField state() const;

        //- Set the state and the time; the boundaries are updated with the
        //  forcing then set
        void setState(const scalar t, const scalarField& y);


    // Checkpoint

        //- The state as a dictionary [SI]. Values below VSMALL (1e-300) in
        //  magnitude are written as 0: they arise where the friction source
        //  keeps shrinking round-off ahead of a wave, and OpenFOAM cannot read
        //  numbers that small back.
        dictionary stateDict() const;

        //- A copy of f with the values below VSMALL in magnitude set to 0
        static scalarField flushDenormals(const scalarField& f);

        //- Restore from stateDict(); also sets the time
        void readState(const dictionary& d);


    // Access

        //- Blood density [kg/m3]
        scalar rho() const { return rho_; }

        const PtrList<vessel1D>& vessels() const { return vessels_; }

        PtrList<node1D>& nodes() { return nodes_; }

        const wordList& probeNames() const { return probeNames_; }


    // Results

        //- P, Q, A, U and the change in radius from the reference area at
        //  probe i
        void probeValues
        (
            const label i,
            scalar& P,
            scalar& Q,
            scalar& A,
            scalar& U,
            scalar& dr
        ) const;

        //- Pressures and flows at the two ends of vessel v
        void vesselEndValues
        (
            const label v,
            scalar& Pin,
            scalar& Pout,
            scalar& Qin,
            scalar& Qout
        ) const;

        //- Volume of the network [m3]
        scalar volume() const;

        //- Integral of the net inflow at the vessel ends since t0 [m3]
        scalar netInflow() const { return inflow_; }

        //- V - V0 - netInflow [m3]
        scalar volumeDrift() const { return volume() - V0_ - inflow_; }

        //- Initial volume [m3]
        scalar initialVolume() const { return V0_; }

        //- x A Q P U along vessel v: the in end, the cell centres, the out end
        void writeProfile(const label v, Ostream& os) const;


    // Info

        void summary(Ostream& os) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
