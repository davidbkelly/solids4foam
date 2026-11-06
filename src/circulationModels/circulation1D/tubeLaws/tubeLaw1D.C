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

#include "tubeLaw1D.H"

namespace Foam
{
    defineTypeNameAndDebug(tubeLaw1D, 0);
    defineRunTimeSelectionTable(tubeLaw1D, dictionary);
}


Foam::tubeLaw1D::tubeLaw1D
(
    const word& vesselName,
    const dictionary& dict,
    const scalar rho
)
:
    vesselName_(vesselName),
    rho_(rho)
{}


Foam::autoPtr<Foam::tubeLaw1D> Foam::tubeLaw1D::New
(
    const word& vesselName,
    const dictionary& dict,
    const scalar rho,
    const scalar s
)
{
    const word type = dict.getOrDefault<word>("tubeLaw", "sqrtElastic");

    auto* ctorPtr = dictionaryConstructorTable(type);

    if (!ctorPtr)
    {
        FatalIOErrorInFunction(dict)
            << "Unknown tubeLaw '" << type << "' for vessel '" << vesselName
            << "'. Tube laws:" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalIOError);
    }

    return autoPtr<tubeLaw1D>(ctorPtr(vesselName, dict, rho, s));
}

// ************************************************************************* //
