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
    Foam::linearisedTubeLaw

Description
    A tube law linearised about a reference state (Aref, Q = 0): the wall of
    a vessel solved with the linearised 1D equations ("linearised true" in
    the vessel dictionary). Wraps the vessel's own tube law and replaces it
    by its tangent at Aref,

        P(A) = P(Aref) + k (A - Aref),   k = dP/dA(Aref),
        c0^2 = Aref k/rho,

    with the quantities the linear system

        dA/dt + dQ/dx = 0
        dQ/dt + c0^2 dA/dx = -KR Q/Aref

    needs: constant wave speed c0, momentum flux c0^2 (A - Aref) and no
    convective term, velocity U = Q/Aref (velocityArea() returns Aref), and
    the linear Riemann invariants W = Q/Aref +- c0 (A - Aref)/Aref, so that
    invariant(A) = c0 (A - Aref)/Aref. nonlinear() is false, which tells the
    vessel and the nodes to drop the convective and dynamic-pressure terms.

    With it, the linear closed forms of the 3D-1D benchmarks (the rigid
    segment, the Womersley segment, the delay equation of the open-ended
    pipe) are the exact solutions of the 1D equations at any amplitude.

    Not selectable by "tubeLaw": the vessel builds it around the law it
    selected, with Aref its initial area.

Author
    David Kelly, UCD.

SourceFiles
    linearisedTubeLaw.C

\*---------------------------------------------------------------------------*/

#ifndef linearisedTubeLaw_H
#define linearisedTubeLaw_H

#include "tubeLaw1D.H"

namespace Foam
{

class linearisedTubeLaw
:
    public tubeLaw1D
{
    // Private Data

        //- The law linearised
        autoPtr<tubeLaw1D> law_;

        //- Reference area [m2] and pressure [Pa]
        const scalar Aref_;
        const scalar Pref_;

        //- k = dP/dA at Aref [Pa/m2]
        const scalar k_;

        //- Wave speed [m/s]
        const scalar c0_;


public:

    TypeName("linearised");

    //- Construct about the area Aref, taking the law
    linearisedTubeLaw
    (
        const word& vesselName,
        const dictionary& dict,
        const scalar rho,
        autoPtr<tubeLaw1D>&& law,
        const scalar Aref
    );


    // Member Functions

        virtual bool uniform() const { return law_->uniform(); }

        virtual scalar referenceArea() const { return law_->referenceArea(); }

        virtual scalar pressure(const scalar A) const;

        virtual scalar minPressure() const;

        virtual scalar dPdA(const scalar A) const { return k_; }

        virtual scalar waveSpeed(const scalar A) const { return c0_; }

        virtual scalar area(const scalar P) const;

        virtual scalar fluxPressure(const scalar A) const;

        virtual scalar invariant(const scalar A) const;

        virtual scalar areaFromInvariant(const scalar I) const;

        virtual bool nonlinear() const { return false; }

        virtual scalar velocityArea(const scalar A) const { return Aref_; }

        virtual void summary(Ostream& os) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
