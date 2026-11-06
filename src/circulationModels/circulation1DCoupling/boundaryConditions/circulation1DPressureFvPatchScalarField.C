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

Description
    Sets the patch pressure of the circulation1DPressure boundary condition
    from the coupler, and writes the boundary condition.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#include "circulation1DPressureFvPatchScalarField.H"
#include "circulation1DCoupler.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::circulation1DPressureFvPatchScalarField::
circulation1DPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    propertiesName_("circulation1DProperties")
{}


Foam::circulation1DPressureFvPatchScalarField::
circulation1DPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    propertiesName_
    (
        dict.getOrDefault<word>("network", "circulation1DProperties")
    )
{}


Foam::circulation1DPressureFvPatchScalarField::
circulation1DPressureFvPatchScalarField
(
    const circulation1DPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    propertiesName_(ptf.propertiesName_)
{}


Foam::circulation1DPressureFvPatchScalarField::
circulation1DPressureFvPatchScalarField
(
    const circulation1DPressureFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    propertiesName_(ptf.propertiesName_)
{}


Foam::circulation1DPressureFvPatchScalarField::
circulation1DPressureFvPatchScalarField
(
    const circulation1DPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    propertiesName_(ptf.propertiesName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::circulation1DPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    circulation1DCoupler& coupler =
        circulation1DCoupler::lookupOrCreate(mesh, propertiesName_);

    // Sample the fluxes and, on a new exchange, advance the network. Repeated
    // evaluations within one exchange reuse the result.
    coupler.sampleAndAdvance(mesh);

    // The network works in Pa; the incompressible solvers expect p/rho
    operator==
    (
        coupler.patchPressure(coupler.portOfPatch(patch().index()))
       /coupler.rho()
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::circulation1DPressureFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>
    (
        "network", "circulation1DProperties", propertiesName_
    );
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        circulation1DPressureFvPatchScalarField
    );
}

// ************************************************************************* //
