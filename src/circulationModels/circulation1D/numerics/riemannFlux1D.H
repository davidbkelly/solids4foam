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
    Foam::riemannFlux1D

Description
    Approximate Riemann solver for the interior faces of a vessel: the face
    flux of (A, Q) from the reconstructed states either side of it. Selected
    by the "flux" entry of "numerics":

        HLL         (default) Harten-Lax-van Leer, with
                    SL = min(lambda-(L), lambda-(R)),
                    SR = max(lambda+(L), lambda+(R))
        Rusanov     local Lax-Friedrichs, with the largest |lambda| of L, R

    The physical flux and the eigenvalues are those of the vessel.

Author
    David Kelly, UCD.

SourceFiles
    riemannFlux1D.C

\*---------------------------------------------------------------------------*/

#ifndef riemannFlux1D_H
#define riemannFlux1D_H

#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"
#include "typeInfo.H"

namespace Foam
{

class vessel1D;

class riemannFlux1D
{
public:

    TypeName("riemannFlux1D");

    declareRunTimeSelectionTable
    (
        autoPtr,
        riemannFlux1D,
        dictionary,
        (const dictionary& dict),
        (dict)
    );

    riemannFlux1D() = default;

    //- Select by dict.getOrDefault<word>("flux", "HLL")
    static autoPtr<riemannFlux1D> New(const dictionary& dict);

    virtual ~riemannFlux1D() = default;


    // Member Functions

        //- Face flux (FA, FQ) for the left and right states
        virtual void flux
        (
            const vessel1D& v,
            const scalar AL,
            const scalar QL,
            const scalar AR,
            const scalar QR,
            scalar& FA,
            scalar& FQ
        ) const = 0;
};


class HLLRiemannFlux1D
:
    public riemannFlux1D
{
public:

    TypeName("HLL");

    explicit HLLRiemannFlux1D(const dictionary&) {}

    virtual void flux
    (
        const vessel1D& v,
        const scalar AL,
        const scalar QL,
        const scalar AR,
        const scalar QR,
        scalar& FA,
        scalar& FQ
    ) const;
};


class RusanovRiemannFlux1D
:
    public riemannFlux1D
{
public:

    TypeName("Rusanov");

    explicit RusanovRiemannFlux1D(const dictionary&) {}

    virtual void flux
    (
        const vessel1D& v,
        const scalar AL,
        const scalar QL,
        const scalar AR,
        const scalar QR,
        scalar& FA,
        scalar& FQ
    ) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
