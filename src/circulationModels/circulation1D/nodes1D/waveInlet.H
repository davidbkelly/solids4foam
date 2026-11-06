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
    Foam::waveInlet

Description
    A non-reflecting end that sends a prescribed wave into the vessel: the
    incoming Riemann invariant is set to its initial value plus the wave,

        W_in(t) = W_in0 - 2 s p(t)/(rho c0),

    with p the pressure of the incident wave (a Function1 of time, Pa), s
    the end sign and c0 the wave speed of the initial state, so that in the
    linear regime a wave of pressure p and velocity p/(rho c0) enters the
    vessel and whatever reaches the end leaves it without reflection (as for
    the absorbing node, of which this is the forced version).

        inlet
        {
            type        waveInlet;
            connection  proximal.in;
            pressure
            {
                type        expression;
                expression  #{ 50*sin(6.283185307179586*arg()) #};
            }
        }

    A flowInlet or pressureInlet reflects the returning waves (as a source
    of zero or infinite impedance); this end does not, which makes a short
    vessel behave as a semi-infinite line: the single-frequency benchmarks
    use it.

Author
    David Kelly, UCD.

SourceFiles
    waveInlet.C

\*---------------------------------------------------------------------------*/

#ifndef waveInlet_H
#define waveInlet_H

#include "node1D.H"
#include "Function1.H"

namespace Foam
{

class waveInlet
:
    public node1D
{
    // Private Data

        //- Incoming invariant of the initial state [m/s]
        scalar Win0_;

        //- Pressure to velocity of the incident wave, -2 s/(rho c0) [m/s/Pa]
        scalar scale_;

        autoPtr<Function1<scalar>> pressure_;


public:

    TypeName("waveInlet");

    waveInlet
    (
        const word& name,
        const dictionary& dict,
        PtrList<vessel1D>& vessels,
        const scalar t0
    );


    // Member Functions

        //- Pressure of the incident wave at time t [Pa]
        scalar pressure(const scalar t) const
        {
            return pressure_->value(t);
        }

        virtual void setBoundary(const scalar t);
};

} // End namespace Foam

#endif

// ************************************************************************* //
