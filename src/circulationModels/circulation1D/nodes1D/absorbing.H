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
    Foam::absorbing

Description
    Non-reflecting vessel end: the incoming Riemann invariant is held at its
    value in the initial state of the end cell, Win0 = U0 - s I(A0). With the
    outgoing invariant W from the vessel,

        U = (W + Win0)/2,       I(A) = s (W - Win0)/2

    in closed form.

        outlet { type absorbing; connection tube.out; }

Author
    David Kelly, UCD.

SourceFiles
    absorbing.C

\*---------------------------------------------------------------------------*/

#ifndef absorbing_H
#define absorbing_H

#include "node1D.H"

namespace Foam
{

class absorbing
:
    public node1D
{
    // Private Data

        //- Incoming invariant [m/s]
        scalar Win0_;


public:

    TypeName("absorbing");

    absorbing
    (
        const word& name,
        const dictionary& dict,
        PtrList<vessel1D>& vessels,
        const scalar t0
    );


    // Member Functions

        virtual void setBoundary(const scalar t);
};

} // End namespace Foam

#endif

// ************************************************************************* //
