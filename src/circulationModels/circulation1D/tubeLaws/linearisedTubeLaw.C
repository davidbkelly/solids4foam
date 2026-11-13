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

#include "linearisedTubeLaw.H"

namespace Foam
{
    defineTypeNameAndDebug(linearisedTubeLaw, 0);
}


Foam::linearisedTubeLaw::linearisedTubeLaw
(
    const word& vesselName,
    const dictionary& dict,
    const scalar rho,
    autoPtr<tubeLaw1D>&& law,
    const scalar Aref
)
:
    tubeLaw1D(vesselName, dict, rho),
    law_(std::move(law)),
    Aref_(Aref),
    Pref_(law_->pressure(Aref)),
    k_(law_->dPdA(Aref)),
    c0_(law_->waveSpeed(Aref))
{
    if (Aref_ <= 0 || k_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "Vessel " << vesselName_ << ": cannot linearise the tube law"
            << " about A = " << Aref_ << " m2 (dP/dA = " << k_ << " Pa/m2)"
            << exit(FatalIOError);
    }
}


Foam::scalar Foam::linearisedTubeLaw::pressure(const scalar A) const
{
    return Pref_ + k_*(A - Aref_);
}


Foam::scalar Foam::linearisedTubeLaw::minPressure() const
{
    return Pref_ - k_*Aref_;
}


Foam::scalar Foam::linearisedTubeLaw::area(const scalar P) const
{
    return Aref_ + (P - Pref_)/k_;
}


Foam::scalar Foam::linearisedTubeLaw::fluxPressure(const scalar A) const
{
    // (Aref/rho) (P - Pref) = c0^2 (A - Aref)
    return c0_*c0_*(A - Aref_);
}


Foam::scalar Foam::linearisedTubeLaw::invariant(const scalar A) const
{
    return c0_*(A - Aref_)/Aref_;
}


Foam::scalar Foam::linearisedTubeLaw::areaFromInvariant(const scalar I) const
{
    return Aref_*(1 + I/c0_);
}


void Foam::linearisedTubeLaw::summary(Ostream& os) const
{
    os  << "linearised about A " << Aref_ << " m2, P " << Pref_
        << " Pa: dP/dA " << k_ << " Pa/m2, c0 " << c0_ << " m/s, Zc "
        << rho_*c0_/Aref_ << " Pa s/m3" << nl
        << "    of ";
    law_->summary(os);
}

// ************************************************************************* //
