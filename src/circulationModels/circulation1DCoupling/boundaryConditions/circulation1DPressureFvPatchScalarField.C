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
    from the coupler (a fixed value, or the impedance condition), and writes
    the boundary condition.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#include "circulation1DPressureFvPatchScalarField.H"
#include "circulation1DCoupler.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::circulation1DPressureFvPatchScalarField::
circulation1DPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    propertiesName_("circulation1DProperties"),
    snGradSet_(false)
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 1;
}


Foam::circulation1DPressureFvPatchScalarField::
circulation1DPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    propertiesName_
    (
        dict.getOrDefault<word>("network", "circulation1DProperties")
    ),
    snGradSet_(false)
{
    // A fixed value unless the mixed entries were written (a restart of a
    // Robin run)
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
    }
    else
    {
        fvPatchScalarField::operator=(Zero);
    }

    if (dict.found("refValue"))
    {
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        refValue() = *this;
        refGrad() = Zero;
        valueFraction() = 1;
    }
}


Foam::circulation1DPressureFvPatchScalarField::
circulation1DPressureFvPatchScalarField
(
    const circulation1DPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    propertiesName_(ptf.propertiesName_),
    snGradSet_(false)
{}


Foam::circulation1DPressureFvPatchScalarField::
circulation1DPressureFvPatchScalarField
(
    const circulation1DPressureFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    propertiesName_(ptf.propertiesName_),
    snGradCoeff_(ptf.snGradCoeff_),
    snGradSet_(ptf.snGradSet_),
    rAUf_(ptf.rAUf_)
{}


Foam::circulation1DPressureFvPatchScalarField::
circulation1DPressureFvPatchScalarField
(
    const circulation1DPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    propertiesName_(ptf.propertiesName_),
    snGradCoeff_(ptf.snGradCoeff_),
    snGradSet_(ptf.snGradSet_),
    rAUf_(ptf.rAUf_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::circulation1DPressureFvPatchScalarField::updateSnGrad
(
    const scalarField& snGradp
)
{
    snGradCoeff_ = snGradp;
    snGradSet_ = true;
}


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

    const word& port = coupler.portOfPatch(patch().index());

    // The network works in Pa; the incompressible solvers expect p/rho
    const scalar Pstar = coupler.patchPressure(port)/coupler.rho();

    if (!coupler.robin())
    {
        // The Dirichlet interface: a fixed value
        refValue() = Pstar;
        refGrad() = Zero;
        valueFraction() = 1;

        mixedFvPatchScalarField::updateCoeffs();
        return;
    }

    // ---- The impedance (Robin) interface -----------------------------------

    const scalarField& magSf = patch().magSf();
    const scalar zu = coupler.robinCoefficient(port);
    const scalarField& phiStar = coupler.iterateFlux(port);

    const auto* phiHbyAPtr = db().cfindObject<surfaceScalarField>("phiHbyA");
    if (!phiHbyAPtr)
    {
        FatalErrorInFunction
            << "The Robin interface needs the flux predictor phiHbyA of a"
            << " PIMPLE/PISO pressure equation in the registry; none found."
            << " Use  interface dirichlet;  in the coupling dictionary."
            << exit(FatalError);
    }
    const scalarField& phiHbyA = phiHbyAPtr->boundaryField()[patch().index()];

    // rAU_f of the faces: from the solver's updateSnGrad call, which passes
    // (phiHbyA - Sf & U)/(|Sf| rAU_f); the time step until it has been seen
    if (rAUf_.size() != patch().size())
    {
        rAUf_ = scalarField(patch().size(), mesh.time().deltaTValue());
    }
    if (snGradSet_ && snGradCoeff_.size() == patch().size())
    {
        const fvPatchVectorField& U =
            patch().lookupPatchField<volVectorField, vector>("U");
        const scalarField num((phiHbyA - (patch().Sf() & U))/magSf);
        forAll(num, i)
        {
            if (mag(snGradCoeff_[i]) > VSMALL)
            {
                const scalar r = num[i]/snGradCoeff_[i];
                if (std::isfinite(r) && r > 0)
                {
                    rAUf_[i] = r;
                }
            }
        }
        snGradSet_ = false;
    }

    const scalarField k(zu*rAUf_*patch().deltaCoeffs());

    valueFraction() = 1/(1 + k);
    refGrad() = Zero;
    refValue() = Pstar + zu*(phiHbyA - phiStar)/magSf;

    if (debug)
    {
        Pout<< "circulation1DPressure Robin " << patch().name()
            << " t " << mesh.time().value()
            << " exchange " << coupler.exchangeCount()
            << " k " << gAverage(k) << " rAUf " << gAverage(rAUf_)
            << " zu " << zu << " Pstar " << Pstar
            << " max|phiHbyA - phi*| " << gMax(mag(phiHbyA - phiStar))
            << " max|phiHbyA - phi| "
            << gMax(mag(phiHbyA - db().lookupObject<surfaceScalarField>("phi").boundaryField()[patch().index()]))
            << " max|phi*| " << gMax(mag(phiStar))
            << " refValue " << gAverage(refValue())
            << " pb " << gAverage(*this)
            << " pP " << gAverage(patchInternalField()) << endl;
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::circulation1DPressureFvPatchScalarField::write(Ostream& os) const
{
    mixedFvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>
    (
        "network", "circulation1DProperties", propertiesName_
    );
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
