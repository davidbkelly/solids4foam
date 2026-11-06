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
    Foam::node1D

Description
    Base class of the entries of the "nodes" dictionary of a 1D network,
    selected at run time by their "type". A node joins one or more vessel
    ends, named "<vessel>.in" (x = 0) or "<vessel>.out" (x = L):

        inlet  { type flowInlet; connection carotid.in; ... }
        joint  { type ...;       connections (aorta.out iliac.in); ... }

    A node with one connection is a boundary condition. Every stage the node
    sets the boundary state (A, Q) of each of its vessel ends from the
    outgoing characteristic of the vessel, taken from the reconstructed face
    value of the cell next to the end, and its own condition.

    With alpha = 1 the Riemann invariants are W = U +- I(A), I the tube
    law's invariant integral. With s = endSign (+1 at an out end, -1 at an in
    end) the outgoing one is

        W = U + s I(A)

    so on the characteristic U = W - s I(A), and the flow into the node,
    Qn = s A U, has dQn/dA = s (U - s c).

    A node may carry ODE states (y), advanced with the same Runge-Kutta
    stages as the vessels from the derivatives it computes in derivatives().

Author
    David Kelly, UCD.

SourceFiles
    node1D.C

\*---------------------------------------------------------------------------*/

#ifndef node1D_H
#define node1D_H

#include "vessel1D.H"
#include "rungeKutta1D.H"
#include "PtrList.H"
#include "wordList.H"

namespace Foam
{

class node1D
{
public:

    //- A vessel end joined by the node
    struct vesselEnd
    {
        label vessel;
        label end;      //!< vessel1D::inEnd or outEnd
    };


protected:

    // Protected Data

        const word name_;

        const dictionary dict_;

        PtrList<vessel1D>& vessels_;

        List<vesselEnd> ends_;

        //- ODE states, at the start of the step and their derivatives
        scalarField y_;
        scalarField y0_;
        scalarField dydt_;

        //- Butcher-form Runge-Kutta: the stage derivatives of y
        List<scalarField> ky_;


    // Protected Member Functions

        //- Fatal error with the node name and its dictionary location
        [[noreturn]] void fail(const string& message) const;

        //- Fail unless the node joins exactly n vessel ends
        void requireConnections(const label n) const;

        vessel1D& vessel(const label i) { return vessels_[ends_[i].vessel]; }

        const vessel1D& vessel(const label i) const
        {
            return vessels_[ends_[i].vessel];
        }

        label end(const label i) const { return ends_[i].end; }

        //- endSign of connection i
        scalar sign(const label i) const
        {
            return vessel1D::endSign(ends_[i].end);
        }

        //- Outgoing invariant of connection i, from the reconstructed face
        //  state; also returns that face area as a starting guess
        scalar outgoingInvariant(const label i, scalar& Aface) const;

        //- Set the boundary state of connection i for a prescribed flow Qn
        //  into the node: A solves s A (W - s I(A)) = Qn by Newton
        void setFlowBoundary(const label i, const scalar Qn);

        //- Set the number of ODE states
        void setNStates(const label n);

        //- Newton solve of residual(A, dRdA) = 0 for the area, from A, to a
        //  relative change of 1e-12
        template<class Residual>
        scalar solveArea(const Residual& residual, scalar A) const;


public:

    TypeName("node1D");

    declareRunTimeSelectionTable
    (
        autoPtr,
        node1D,
        dictionary,
        (
            const word& name,
            const dictionary& dict,
            PtrList<vessel1D>& vessels,
            const scalar t0
        ),
        (name, dict, vessels, t0)
    );

    //- Construct reading the connections
    node1D
    (
        const word& name,
        const dictionary& dict,
        PtrList<vessel1D>& vessels
    );

    //- Select by dict.get<word>("type")
    static autoPtr<node1D> New
    (
        const word& name,
        const dictionary& dict,
        PtrList<vessel1D>& vessels,
        const scalar t0
    );

    virtual ~node1D() = default;


    // Member Functions

        const word& name() const { return name_; }

        const List<vesselEnd>& ends() const { return ends_; }

        //- Set the boundary states of the vessel ends at time t, from the
        //  current stage (vessels reconstructed, node states y)
        virtual void setBoundary(const scalar t) = 0;

        //- Derivatives of the ODE states at time t, from the boundary states
        //  just set
        virtual void derivatives(const scalar t) {}


    // ODE states

        label nStates() const { return y_.size(); }

        const scalarField& states() const { return y_; }

        void setStates(const scalarField& y);

        void storeOld() { y0_ = y_; }

        void stageUpdate(const scalar a, const scalar b, const scalar dt);

        //- Butcher form: keep dy/dt as stage derivative k
        void storeStageRhs(const label k, const label nStages);

        //- Butcher form: y = y^n + dt sum_{j<k} A_kj (dy/dt)_j; k = nStages
        //  gives the end of the step, with the weights B
        void butcherStage
        (
            const rungeKutta1D& rk,
            const label k,
            const scalar dt
        );


    // Output

        //- Names of extra output columns of the node, if any
        virtual wordList columns() const { return wordList(); }

        //- Values of those columns at time t (boundary states up to date)
        virtual void outputValues(const scalar t, scalarField& values) {}

        virtual void summary(Ostream& os) const;
};

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Templates  * * * * * * * * * * * * * //

template<class Residual>
Foam::scalar Foam::node1D::solveArea(const Residual& residual, scalar A) const
{
    const label maxIter = 100;

    for (label iter = 0; iter < maxIter; ++iter)
    {
        scalar dRdA = 0;
        const scalar R = residual(A, dRdA);

        if (R == 0)
        {
            return A;
        }

        if (dRdA == 0 || !std::isfinite(dRdA))
        {
            fail("singular Newton derivative at A = " + std::to_string(A));
        }

        scalar An = A - R/dRdA;

        // Stay positive
        if (An <= 0)
        {
            An = 0.5*A;
        }

        const bool converged = mag(An - A) <= 1e-12*An;
        A = An;

        if (converged)
        {
            return A;
        }
    }

    fail
    (
        "Newton solve for the boundary area did not converge in "
      + std::to_string(maxIter) + " iterations (A = " + std::to_string(A)
      + " m2)"
    );
}

#endif

// ************************************************************************* //
