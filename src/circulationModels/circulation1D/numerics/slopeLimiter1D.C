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

#include "slopeLimiter1D.H"

namespace Foam
{
    defineTypeNameAndDebug(slopeLimiter1D, 0);
    defineRunTimeSelectionTable(slopeLimiter1D, dictionary);

    defineTypeNameAndDebug(minmodSlopeLimiter1D, 0);
    addToRunTimeSelectionTable
    (
        slopeLimiter1D,
        minmodSlopeLimiter1D,
        dictionary
    );

    defineTypeNameAndDebug(vanLeerSlopeLimiter1D, 0);
    addToRunTimeSelectionTable
    (
        slopeLimiter1D,
        vanLeerSlopeLimiter1D,
        dictionary
    );

    defineTypeNameAndDebug(MCSlopeLimiter1D, 0);
    addToRunTimeSelectionTable(slopeLimiter1D, MCSlopeLimiter1D, dictionary);

    defineTypeNameAndDebug(noneSlopeLimiter1D, 0);
    addToRunTimeSelectionTable(slopeLimiter1D, noneSlopeLimiter1D, dictionary);
}


namespace
{
    inline Foam::scalar minmod(const Foam::scalar a, const Foam::scalar b)
    {
        if (a*b <= 0)
        {
            return 0;
        }
        return a > 0 ? Foam::min(a, b) : Foam::max(a, b);
    }
}


Foam::autoPtr<Foam::slopeLimiter1D> Foam::slopeLimiter1D::New
(
    const dictionary& dict
)
{
    const word type = dict.getOrDefault<word>("limiter", "minmod");

    auto* ctorPtr = dictionaryConstructorTable(type);

    if (!ctorPtr)
    {
        FatalIOErrorInFunction(dict)
            << "Unknown limiter '" << type << "'. Limiters:" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalIOError);
    }

    return autoPtr<slopeLimiter1D>(ctorPtr(dict));
}


Foam::scalar Foam::minmodSlopeLimiter1D::slope
(
    const scalar dL,
    const scalar dR
) const
{
    return minmod(dL, dR);
}


Foam::scalar Foam::vanLeerSlopeLimiter1D::slope
(
    const scalar dL,
    const scalar dR
) const
{
    const scalar p = dL*dR;
    return p > 0 ? 2*p/(dL + dR) : 0;
}


Foam::scalar Foam::MCSlopeLimiter1D::slope
(
    const scalar dL,
    const scalar dR
) const
{
    return minmod(0.5*(dL + dR), minmod(2*dL, 2*dR));
}


Foam::scalar Foam::noneSlopeLimiter1D::slope
(
    const scalar dL,
    const scalar dR
) const
{
    return 0.5*(dL + dR);
}

// ************************************************************************* //
