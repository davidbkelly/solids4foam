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
    Foam::inertance

Description
    A rigid straight pipe between two vessel ends, or between one vessel end
    and an outlet pressure (an open-ended pipe), as plug flow: the 0D
    reduction of a rigid 3D segment. Incompressible flow in a rigid pipe has
    the same flow rate Q at both ends, and with slip walls the 3D solution is
    uniform velocity Q/Ap and linear pressure, so exactly

        L dQ/dt = P1 - P2 - R Q,        L = rho length/Ap

    with P1, P2 the pressures at the two vessel ends and R an optional series
    resistance. Q is positive from the first connection to the second. With
    one connection, P2 is the outlet pressure Pout (default 0): a rigid pipe
    open to a reservoir, which reflects an incident wave with
    Gamma = (i w L - Zc)/(i w L + Zc), of modulus 1.

        rigid
        {
            type         inertance;
            connections  (upstream.out downstream.in);
            length       0.05;          // [m]
            radius       0.01;          // [m], or area [m2]
            R            0;             // optional [Pa s/m3]
            pressure     static;        // static | total
            Q0           0;             // optional initial flow [m3/s]
        }

        open { type inertance; connection vessel.out; length 0.05;
               radius 0.01; Pout 0; }

    "pressure" is the quantity taken as continuous between each vessel end
    and the pipe: static (P) or total (P + rho U^2/2). For total, the
    pipe's own kinetic pressure is the same at both of its ends and cancels,
    so P1, P2 above become P + rho (Q/A)^2/2 with A the vessel end area.

    Each vessel end takes the flow Q through its outgoing characteristic, as
    flowInlet does; Q is an ODE state advanced in the same Runge-Kutta stages
    as the vessels. In the linear limit the pipe reflects an incident wave
    with Gamma = i w tau/(1 + i w tau), tau = L/(2 Zc) = length/(2c) for two
    equal vessels (see src/circulationModels/docs/BENCHMARKS3D1D.md).

    The explicit stages need dt 2c/length below about 2: vessel cells no
    longer than the pipe.

Author
    David Kelly, UCD.

SourceFiles
    inertance.C

\*---------------------------------------------------------------------------*/

#ifndef inertance_H
#define inertance_H

#include "node1D.H"

namespace Foam
{

class inertance
:
    public node1D
{
    // Private Data

        //- Pipe length [m] and area [m2]
        scalar length_;
        scalar area_;

        //- Series resistance [Pa s/m3]
        scalar R_;

        //- Total rather than static pressure continuity
        bool total_;

        //- Outlet pressure with one connection [Pa]
        scalar Pout_;


    // Private Member Functions

        //- rho length/area [kg/m4]
        scalar L() const { return vessel(0).rho()*length_/area_; }

        //- Pressure of connection i that the pipe sees, for flow Q
        scalar endPressure(const label i, const scalar Q) const;


public:

    TypeName("inertance");

    inertance
    (
        const word& name,
        const dictionary& dict,
        PtrList<vessel1D>& vessels,
        const scalar t0
    );


    // Member Functions

        virtual void setBoundary(const scalar t);

        virtual void derivatives(const scalar t);

        virtual wordList columns() const;

        virtual void outputValues(const scalar t, scalarField& values);

        virtual void summary(Ostream& os) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
