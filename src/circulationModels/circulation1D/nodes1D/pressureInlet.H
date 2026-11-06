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
    Foam::pressureInlet

Description
    Prescribed pressure at a vessel end, imposed through the outgoing
    characteristic: the area is the tube law's A(P(t)) and the velocity
    follows from the outgoing invariant, U = W - s I(A). The pressure is the
    tube law's pressure (external and reference pressures included) and is a
    Function1 of time [Pa], so for instance the pressure pulse of the
    Greenshields-Weller / Formaggia tube:

        inlet
        {
            type        pressureInlet;
            connection  tube.in;
            pressure
            {
                type        expression;
                expression  #{ (arg() <= 3e-3) ? 1333.2 : 0 #};
            }
        }

    The end reflects nothing of what reaches it in pressure (a pressure
    source has zero impedance): an incoming wave f returns as -f.

Author
    David Kelly, UCD.

SourceFiles
    pressureInlet.C

\*---------------------------------------------------------------------------*/

#ifndef pressureInlet_H
#define pressureInlet_H

#include "node1D.H"
#include "Function1.H"

namespace Foam
{

class pressureInlet
:
    public node1D
{
    // Private Data

        autoPtr<Function1<scalar>> pressure_;


public:

    TypeName("pressureInlet");

    pressureInlet
    (
        const word& name,
        const dictionary& dict,
        PtrList<vessel1D>& vessels,
        const scalar t0
    );


    // Member Functions

        //- Pressure at the vessel end at time t [Pa]
        scalar pressure(const scalar t) const
        {
            return pressure_->value(t);
        }

        virtual void setBoundary(const scalar t);
};

} // End namespace Foam

#endif

// ************************************************************************* //
