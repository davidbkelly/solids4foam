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

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#include "riemannFlux1D.H"
#include "vessel1D.H"

namespace Foam
{
    defineTypeNameAndDebug(riemannFlux1D, 0);
    defineRunTimeSelectionTable(riemannFlux1D, dictionary);

    defineTypeNameAndDebug(HLLRiemannFlux1D, 0);
    addToRunTimeSelectionTable(riemannFlux1D, HLLRiemannFlux1D, dictionary);

    defineTypeNameAndDebug(RusanovRiemannFlux1D, 0);
    addToRunTimeSelectionTable
    (
        riemannFlux1D,
        RusanovRiemannFlux1D,
        dictionary
    );
}


Foam::autoPtr<Foam::riemannFlux1D> Foam::riemannFlux1D::New
(
    const dictionary& dict
)
{
    const word type = dict.getOrDefault<word>("flux", "HLL");

    auto* ctorPtr = dictionaryConstructorTable(type);

    if (!ctorPtr)
    {
        FatalIOErrorInFunction(dict)
            << "Unknown flux '" << type << "'. Fluxes:" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalIOError);
    }

    return autoPtr<riemannFlux1D>(ctorPtr(dict));
}


void Foam::HLLRiemannFlux1D::flux
(
    const vessel1D& v,
    const scalar AL,
    const scalar QL,
    const scalar AR,
    const scalar QR,
    scalar& FA,
    scalar& FQ
) const
{
    scalar lmL, lpL, lmR, lpR;
    v.eigenvalues(AL, QL, lmL, lpL);
    v.eigenvalues(AR, QR, lmR, lpR);

    const scalar SL = Foam::min(lmL, lmR);
    const scalar SR = Foam::max(lpL, lpR);

    scalar FAL, FQL, FAR, FQR;

    if (SL >= 0)
    {
        v.physicalFlux(AL, QL, FA, FQ);
        return;
    }
    if (SR <= 0)
    {
        v.physicalFlux(AR, QR, FA, FQ);
        return;
    }

    v.physicalFlux(AL, QL, FAL, FQL);
    v.physicalFlux(AR, QR, FAR, FQR);

    const scalar r = 1/(SR - SL);
    FA = (SR*FAL - SL*FAR + SL*SR*(AR - AL))*r;
    FQ = (SR*FQL - SL*FQR + SL*SR*(QR - QL))*r;
}


void Foam::RusanovRiemannFlux1D::flux
(
    const vessel1D& v,
    const scalar AL,
    const scalar QL,
    const scalar AR,
    const scalar QR,
    scalar& FA,
    scalar& FQ
) const
{
    scalar lmL, lpL, lmR, lpR;
    v.eigenvalues(AL, QL, lmL, lpL);
    v.eigenvalues(AR, QR, lmR, lpR);

    const scalar s =
        Foam::max
        (
            Foam::max(mag(lmL), mag(lpL)),
            Foam::max(mag(lmR), mag(lpR))
        );

    scalar FAL, FQL, FAR, FQR;
    v.physicalFlux(AL, QL, FAL, FQL);
    v.physicalFlux(AR, QR, FAR, FQR);

    FA = 0.5*(FAL + FAR) - 0.5*s*(AR - AL);
    FQ = 0.5*(FQL + FQR) - 0.5*s*(QR - QL);
}

// ************************************************************************* //
