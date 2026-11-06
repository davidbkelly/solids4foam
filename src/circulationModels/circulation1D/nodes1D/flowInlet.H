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
    Foam::flowInlet

Description
    Prescribed flow rate into a vessel end (Boileau et al. 2015, Eq. 6),
    imposed through the outgoing characteristic: A solves

        A (W - s I(A)) = Qv,        Qv = -s Qin(t)

    by Newton, with Qin positive into the vessel. Qin is a Function1 of time
    [m3/s], so a constant, a table (periodic with "outOfBounds repeat"), a
    CSV file or an expression:

        inlet
        {
            type        flowInlet;
            connection  carotid.in;
            flowRate
            {
                type        expression;
                expression  #{ 1e-6*exp(-1e4*sqr(arg() - 0.05)) #};
            }
        }

Author
    David Kelly, UCD.

SourceFiles
    flowInlet.C

\*---------------------------------------------------------------------------*/

#ifndef flowInlet_H
#define flowInlet_H

#include "node1D.H"
#include "Function1.H"

namespace Foam
{

class flowInlet
:
    public node1D
{
    // Private Data

        autoPtr<Function1<scalar>> flowRate_;


public:

    TypeName("flowInlet");

    flowInlet
    (
        const word& name,
        const dictionary& dict,
        PtrList<vessel1D>& vessels,
        const scalar t0
    );


    // Member Functions

        //- Flow into the vessel at time t [m3/s]
        scalar flowRate(const scalar t) const
        {
            return flowRate_->value(t);
        }

        virtual void setBoundary(const scalar t);
};

} // End namespace Foam

#endif

// ************************************************************************* //
