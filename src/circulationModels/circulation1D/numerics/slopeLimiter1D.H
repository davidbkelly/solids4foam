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
    Foam::slopeLimiter1D

Description
    Limited slope of the MUSCL reconstruction in a cell, from the differences
    to its left and right neighbours, dL = u_i - u_(i-1), dR = u_(i+1) - u_i.
    The cell's face values are u_i -+ slope/2. Selected by the "limiter"
    entry of "numerics":

        minmod      minmod(dL, dR)
        vanLeer     2 dL dR/(dL + dR) if dL dR > 0, else 0
        MC          minmod(2 dL, 2 dR, (dL + dR)/2)
        none        (dL + dR)/2, unlimited central slope (for convergence
                    tests on smooth problems: limiters clip at extrema)

Author
    David Kelly, UCD.

SourceFiles
    slopeLimiter1D.C

\*---------------------------------------------------------------------------*/

#ifndef slopeLimiter1D_H
#define slopeLimiter1D_H

#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"
#include "typeInfo.H"

namespace Foam
{

class slopeLimiter1D
{
public:

    TypeName("slopeLimiter1D");

    declareRunTimeSelectionTable
    (
        autoPtr,
        slopeLimiter1D,
        dictionary,
        (const dictionary& dict),
        (dict)
    );

    slopeLimiter1D() = default;

    //- Select by dict.getOrDefault<word>("limiter", "minmod")
    static autoPtr<slopeLimiter1D> New(const dictionary& dict);

    virtual ~slopeLimiter1D() = default;


    // Member Functions

        //- Limited slope from the left and right differences
        virtual scalar slope(const scalar dL, const scalar dR) const = 0;
};


//- minmod(dL, dR)
class minmodSlopeLimiter1D
:
    public slopeLimiter1D
{
public:

    TypeName("minmod");

    explicit minmodSlopeLimiter1D(const dictionary&) {}

    virtual scalar slope(const scalar dL, const scalar dR) const;
};


//- van Leer's harmonic mean
class vanLeerSlopeLimiter1D
:
    public slopeLimiter1D
{
public:

    TypeName("vanLeer");

    explicit vanLeerSlopeLimiter1D(const dictionary&) {}

    virtual scalar slope(const scalar dL, const scalar dR) const;
};


//- Monotonised central
class MCSlopeLimiter1D
:
    public slopeLimiter1D
{
public:

    TypeName("MC");

    explicit MCSlopeLimiter1D(const dictionary&) {}

    virtual scalar slope(const scalar dL, const scalar dR) const;
};


//- Unlimited central slope
class noneSlopeLimiter1D
:
    public slopeLimiter1D
{
public:

    TypeName("none");

    explicit noneSlopeLimiter1D(const dictionary&) {}

    virtual scalar slope(const scalar dL, const scalar dR) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
