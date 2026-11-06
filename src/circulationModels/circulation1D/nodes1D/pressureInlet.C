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

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#include "pressureInlet.H"

namespace Foam
{
    defineTypeNameAndDebug(pressureInlet, 0);
    addToRunTimeSelectionTable(node1D, pressureInlet, dictionary);
}


Foam::pressureInlet::pressureInlet
(
    const word& name,
    const dictionary& dict,
    PtrList<vessel1D>& vessels,
    const scalar t0
)
:
    node1D(name, dict, vessels),
    pressure_(Function1<scalar>::New("pressure", dict))
{
    requireConnections(1);
}


void Foam::pressureInlet::setBoundary(const scalar t)
{
    const scalar s = sign(0);
    const tubeLaw1D& law = vessel(0).endLaw(end(0));

    scalar Af;
    const scalar W = outgoingInvariant(0, Af);

    const scalar P = pressure(t);
    if (P <= law.minPressure())
    {
        fail
        (
            "pressure " + std::to_string(P) + " Pa at t = " + std::to_string(t)
          + " is not above the tube law's minimum "
          + std::to_string(law.minPressure()) + " Pa"
        );
    }

    const scalar A = law.area(P);
    const scalar U = W - s*law.invariant(A);

    vessel(0).setBoundaryState(end(0), A, A*U);
}

// ************************************************************************* //
