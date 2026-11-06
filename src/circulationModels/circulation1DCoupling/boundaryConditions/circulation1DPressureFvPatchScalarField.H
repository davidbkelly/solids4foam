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
    Foam::circulation1DPressureFvPatchScalarField

Description
    Kinematic pressure p/rho on a patch coupled to a port (coupledFlow node) of
    a 1D arterial network. Whether the exchange is explicit or implicit, and of
    what order, is chosen in the "coupling" block of
    constant/circulation1DProperties, not by the BC type.

        inlet
        {
            type    circulation1DPressure;
            value   uniform 0;
            // network circulation1DProperties;   // optional: which description
        }

    The patch is bound to its port by the port's "patch" entry. Each time the
    boundary is evaluated the coupler is asked to sample the flux and, on a
    new exchange, advance the network; the pressure returned is the one the
    network's vessel end has, divided by rho.

    A copy of circulationPressure (circulationCoupling) for circulation1D,
    until the two share one coupler.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#ifndef circulation1DPressureFvPatchScalarField_H
#define circulation1DPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

namespace Foam
{

class circulation1DPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    word propertiesName_;

public:

    TypeName("circulation1DPressure");

    circulation1DPressureFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    circulation1DPressureFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    circulation1DPressureFvPatchScalarField
    (
        const circulation1DPressureFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    circulation1DPressureFvPatchScalarField
    (
        const circulation1DPressureFvPatchScalarField&
    );

    circulation1DPressureFvPatchScalarField
    (
        const circulation1DPressureFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
            new circulation1DPressureFvPatchScalarField(*this)
        );
    }

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new circulation1DPressureFvPatchScalarField(*this, iF)
        );
    }

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
