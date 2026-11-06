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
    Foam::coupledFlow

Description
    A vessel end coupled to a 3D domain: the 3D side supplies the flow and
    the vessel answers with the pressure at its end, as a circulation0D
    externalPort does for a 0D circuit.

        to3D
        {
            type        coupledFlow;
            connection  aorta.out;
            patch       inlet;          // the CFD patch, for the coupler
            flowScale   4;              // optional: the patch is a quarter
        }

    flowScale (default 1) multiplies the flow sampled on the patch before it
    enters the vessel, for a CFD domain that is a symmetric part of the
    vessel's cross-section; the pressure returned is the same either way.

    The flow over a step [t0, t1] is a polynomial of up to third degree in
    alpha = (t - t0)/(t1 - t0), set by network1D::setPortForcing() as for a
    0D port, and evaluated at every Runge-Kutta stage; the end takes it
    through the outgoing characteristic, as flowInlet does. It is zero until
    a forcing is set.

    The port flow Q is the flow along the vessel at the end, positive from
    x = 0 to x = L. So at an out end the 3D domain lies downstream, the patch
    is an inlet of the 3D domain and Q enters it; at an in end the patch is
    an outlet of the 3D domain and Q leaves it. This is the sign convention
    of the coupler's ports (portInlet()).

Author
    David Kelly, UCD.

SourceFiles
    coupledFlow.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFlow_H
#define coupledFlow_H

#include "node1D.H"
#include "circulationForcing.H"

namespace Foam
{

class coupledFlow
:
    public node1D
{
    // Private Data

        //- CFD patch the port is coupled to (may be empty)
        word patch_;

        //- Patch flow to vessel flow
        scalar flowScale_;

        //- Forcing over [t0_, t1_]
        circulationForcing::coeffs coeffs_;
        scalar t0_;
        scalar t1_;


public:

    TypeName("coupledFlow");

    coupledFlow
    (
        const word& name,
        const dictionary& dict,
        PtrList<vessel1D>& vessels,
        const scalar t0
    );


    // Member Functions

        const word& patch() const { return patch_; }

        //- Factor from the flow through the patch to the flow in the vessel
        scalar flowScale() const { return flowScale_; }

        //- True if the patch is an inlet of the 3D domain (an out end)
        bool portInlet() const { return end(0) == vessel1D::outEnd; }

        //- Set the flow over [t0, t1]
        void setForcing
        (
            const scalar t0,
            const scalar t1,
            const circulationForcing::coeffs& c
        );

        //- Port flow at time t [m3/s]
        scalar flow(const scalar t) const;

        //- Pressure at the vessel end, boundary state up to date [Pa]
        scalar pressure() const;

        virtual void setBoundary(const scalar t);

        virtual wordList columns() const;

        virtual void outputValues(const scalar t, scalarField& values);

        virtual void summary(Ostream& os) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
