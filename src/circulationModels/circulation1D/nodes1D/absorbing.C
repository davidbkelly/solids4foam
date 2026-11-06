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

#include "absorbing.H"

namespace Foam
{
    defineTypeNameAndDebug(absorbing, 0);
    addToRunTimeSelectionTable(node1D, absorbing, dictionary);
}


Foam::absorbing::absorbing
(
    const word& name,
    const dictionary& dict,
    PtrList<vessel1D>& vessels,
    const scalar t0
)
:
    node1D(name, dict, vessels),
    Win0_(0)
{
    requireConnections(1);

    // The initial state of the cell at the end (the vessel's initial
    // condition; the same on a restart, as the node is built before the
    // state is read)
    const vessel1D& v = vessel(0);
    const label i = end(0) == vessel1D::inEnd ? 0 : v.nCells() - 1;
    const scalar A0 = v.A()[i];
    const scalar Q0 = v.Q()[i];

    Win0_ = Q0/A0 - sign(0)*v.endLaw(end(0)).invariant(A0);
}


void Foam::absorbing::setBoundary(const scalar t)
{
    const scalar s = sign(0);

    scalar Af;
    const scalar W = outgoingInvariant(0, Af);

    const scalar A =
        vessel(0).endLaw(end(0)).areaFromInvariant(0.5*s*(W - Win0_));
    const scalar U = 0.5*(W + Win0_);

    vessel(0).setBoundaryState(end(0), A, A*U);
}

// ************************************************************************* //
