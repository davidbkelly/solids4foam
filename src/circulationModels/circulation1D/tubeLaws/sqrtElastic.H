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
    Foam::sqrtElastic

Description
    Elastic tube law of Boileau et al. (2015), Eq. 3:

        P = Pext + Pd + (beta/Ad) (sqrt(A) - sqrt(Ad)),
        beta = 4/3 sqrt(pi) E h

    so that, with k = beta/Ad,

        c  = sqrt(k/(2 rho)) A^(1/4)
        I  = 4 c
        fluxPressure = k/(3 rho) A^(3/2)

    Entries of the vessel dictionary (SI):

        radius              3e-3;       // at the diastolic pressure Pd
        wallThickness       0.3e-3;
        youngsModulus       700e3;
        diastolicPressure   10933;      // Pd, default 0
        externalPressure    0;          // Pext, default 0

    A tapered vessel gives radiusProximal (at x = 0) and radiusDistal (at
    x = L) instead of radius; the radius is linear in x between them.

    The wall thickness may instead follow the radius by the empirical law of
    Blanco et al. (2015), Boileau Eq. 17,

        wallThicknessModel  ADAN56;
        h = R (a exp(b R) + c exp(d R)),
        a = 0.2802, b = -5.053 /cm, c = 0.1324, d = -0.1114 /cm

    with R the radius at the position.

Author
    David Kelly, UCD.

SourceFiles
    sqrtElastic.C

\*---------------------------------------------------------------------------*/

#ifndef sqrtElastic_H
#define sqrtElastic_H

#include "tubeLaw1D.H"

namespace Foam
{

class sqrtElastic
:
    public tubeLaw1D
{
    // Private Data

        //- True if the radius is the same along the vessel
        bool uniform_;

        //- Diastolic radius at this position [m]
        scalar rd_;

        //- Wall thickness [m]
        scalar h_;

        //- Young's modulus [Pa]
        const scalar E_;

        //- Diastolic and external pressure [Pa]
        const scalar Pd_;
        const scalar Pext_;

        //- Diastolic area [m2]
        scalar Ad_;

        //- beta [Pa m]
        scalar beta_;

        //- k = beta/Ad [Pa/m]
        scalar k_;

        //- sqrt(k/(2 rho)), c = cCoeff_ A^(1/4)
        scalar cCoeff_;


public:

    TypeName("sqrtElastic");

    //- Construct for the wall at s = x/L along the vessel
    sqrtElastic
    (
        const word& vesselName,
        const dictionary& dict,
        const scalar rho,
        const scalar s
    );


    // Static Functions

        //- Wall thickness of Boileau Eq. 17 for radius R [m]
        static scalar ADAN56WallThickness(const scalar R);


    // Member Functions

        virtual bool uniform() const { return uniform_; }

        virtual scalar referenceArea() const { return Ad_; }

        virtual scalar pressure(const scalar A) const;

        virtual scalar minPressure() const;

        virtual scalar dPdA(const scalar A) const;

        virtual scalar waveSpeed(const scalar A) const;

        virtual scalar area(const scalar P) const;

        virtual scalar fluxPressure(const scalar A) const;

        virtual scalar invariant(const scalar A) const;

        virtual scalar areaFromInvariant(const scalar I) const;

        virtual void summary(Ostream& os) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
