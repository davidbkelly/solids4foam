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

#include "resistance.H"

namespace Foam
{
    defineTypeNameAndDebug(resistance, 0);
    addToRunTimeSelectionTable(node1D, resistance, dictionary);
}


Foam::resistance::resistance
(
    const word& name,
    const dictionary& dict,
    PtrList<vessel1D>& vessels,
    const scalar t0
)
:
    node1D(name, dict, vessels),
    R_(dict.get<scalar>("R")),
    Pout_(dict.getOrDefault<scalar>("Pout", 0))
{
    requireConnections(1);

    if (R_ < 0)
    {
        fail("R must not be negative");
    }
}


void Foam::resistance::setBoundary(const scalar t)
{
    const scalar s = sign(0);
    const tubeLaw1D& law = vessel(0).endLaw(end(0));

    scalar Af;
    const scalar W = outgoingInvariant(0, Af);

    const scalar A = solveArea
    (
        [&](const scalar a, scalar& dRdA)
        {
            const scalar U = W - s*law.invariant(a);
            const scalar Qn = s*a*U;
            const scalar dQn = s*(U - s*law.waveSpeed(a));
            dRdA = law.dPdA(a) - R_*dQn;
            return law.pressure(a) - Pout_ - R_*Qn;
        },
        Af
    );

    vessel(0).setBoundaryState(end(0), A, A*(W - s*law.invariant(A)));
}


void Foam::resistance::summary(Ostream& os) const
{
    node1D::summary(os);
    os  << "    R " << R_ << " Pa s/m3, Pout " << Pout_ << " Pa" << nl;
}

// ************************************************************************* //
