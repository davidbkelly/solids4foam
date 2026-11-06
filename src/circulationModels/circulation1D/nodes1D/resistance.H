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
    Foam::resistance

Description
    A single terminal resistance R to an outflow pressure Pout (SI):

        P(A) - Pout - R Qn(A) = 0

    with Qn the flow from the vessel into the node on the outgoing
    characteristic, solved for A by Newton.

        outlet { type resistance; connection tube.out; R 1e8; Pout 0; }

    The same terminal written as a lumped0D circuit (a resistor to a ground
    or pressureSource) gives the same answer; this one has no 0D overhead.

Author
    David Kelly, UCD.

SourceFiles
    resistance.C

\*---------------------------------------------------------------------------*/

#ifndef resistance_H
#define resistance_H

#include "node1D.H"

namespace Foam
{

class resistance
:
    public node1D
{
    // Private Data

        //- Resistance [Pa s/m3]
        const scalar R_;

        //- Outflow pressure [Pa]
        const scalar Pout_;


public:

    TypeName("resistance");

    resistance
    (
        const word& name,
        const dictionary& dict,
        PtrList<vessel1D>& vessels,
        const scalar t0
    );


    // Member Functions

        virtual void setBoundary(const scalar t);

        virtual void summary(Ostream& os) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
