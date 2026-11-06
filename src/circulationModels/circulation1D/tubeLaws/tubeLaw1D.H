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
    Foam::tubeLaw1D

Description
    The pressure-area relation P(A) of a 1D vessel, selected at run time by
    the "tubeLaw" entry of the vessel (default sqrtElastic).

    Besides P(A), its derivative and its inverse, a tube law supplies what the
    hyperbolic system needs:

        c(A)          wave speed, c^2 = (A/rho) dP/dA
        fluxPressure  the pressure part of the momentum flux,
                      the integral of (a/rho) dP/da from 0 to A
        invariant     I(A), the integral of c/a da, so that the Riemann
                      invariants are W = U +- I(A) (for alpha = 1), and its
                      inverse

    dI/dA = c/A for every tube law.

    A tube law describes the wall at one position s = x/L along the vessel.
    A vessel whose properties vary along it (a tapered vessel, uniform() is
    false) holds one tube law per cell, built at the cell centre; a uniform
    vessel holds one. Everything is SI.

Author
    David Kelly, UCD.

SourceFiles
    tubeLaw1D.C

\*---------------------------------------------------------------------------*/

#ifndef tubeLaw1D_H
#define tubeLaw1D_H

#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"
#include "typeInfo.H"

namespace Foam
{

class tubeLaw1D
{
protected:

    //- Name of the vessel, for messages
    const word vesselName_;

    //- Blood density [kg/m3]
    const scalar rho_;


public:

    TypeName("tubeLaw1D");

    declareRunTimeSelectionTable
    (
        autoPtr,
        tubeLaw1D,
        dictionary,
        (
            const word& vesselName,
            const dictionary& dict,
            const scalar rho,
            const scalar s
        ),
        (vesselName, dict, rho, s)
    );

    //- Construct from the vessel dictionary
    tubeLaw1D
    (
        const word& vesselName,
        const dictionary& dict,
        const scalar rho
    );

    //- Select by dict.getOrDefault<word>("tubeLaw", "sqrtElastic"), for the
    //  wall at s = x/L along the vessel
    static autoPtr<tubeLaw1D> New
    (
        const word& vesselName,
        const dictionary& dict,
        const scalar rho,
        const scalar s = 0.5
    );

    virtual ~tubeLaw1D() = default;


    // Member Functions

        //- True if the wall is the same all along the vessel
        virtual bool uniform() const = 0;

        //- Reference (diastolic) area Ad [m2]
        virtual scalar referenceArea() const = 0;

        //- Pressure, external and reference pressures included [Pa]
        virtual scalar pressure(const scalar A) const = 0;

        //- Pressure as the area goes to zero; area(P) needs P above it [Pa]
        virtual scalar minPressure() const = 0;

        //- dP/dA [Pa/m2]
        virtual scalar dPdA(const scalar A) const = 0;

        //- Wave speed c [m/s]
        virtual scalar waveSpeed(const scalar A) const = 0;

        //- Area for pressure P, the inverse of pressure() [m2]
        virtual scalar area(const scalar P) const = 0;

        //- Pressure part of the momentum flux [m4/s2]
        virtual scalar fluxPressure(const scalar A) const = 0;

        //- Invariant integral I(A) [m/s]
        virtual scalar invariant(const scalar A) const = 0;

        //- Area for invariant integral I, the inverse of invariant() [m2]
        virtual scalar areaFromInvariant(const scalar I) const = 0;

        //- Print the parameters
        virtual void summary(Ostream& os) const = 0;
};

} // End namespace Foam

#endif

// ************************************************************************* //
