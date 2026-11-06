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
    Foam::junction

Description
    A junction of N >= 2 vessel ends (Boileau et al. 2015, Eqs. 7-8): mass
    is conserved and the total pressure is continuous,

        sum_i s_i A_i U_i = 0
        P_1 + rho U_1^2/2 = P_i + rho U_i^2/2,      i = 2..N

    with U_i = W_i - s_i I_i(A_i) on the outgoing characteristic of each
    vessel (s_i = +1 at an out end, -1 at an in end). The N areas are solved
    by Newton with the analytic Jacobian, whose rows i >= 2 involve only A_1
    and A_i, so each iteration is solved in O(N).

        aortaEnd { type junction; connections (aorta.out iliacL.in iliacR.in); }

    A two-vessel junction joins two vessels of different properties, for
    example the pieces of a tapered vessel. For a vessel at rest (U = 0) the
    conditions reduce to equal pressures, so a network at uniform pressure
    stays at rest.

    In the linear limit the reflection coefficient of the pressure for a wave
    arriving in vessel 0 is (Y0 - sum Yi)/(Y0 + sum Yi), Y = A/(rho c).

Author
    David Kelly, UCD.

SourceFiles
    junction.C

\*---------------------------------------------------------------------------*/

#ifndef junction_H
#define junction_H

#include "node1D.H"

namespace Foam
{

class junction
:
    public node1D
{
public:

    TypeName("junction");

    junction
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
