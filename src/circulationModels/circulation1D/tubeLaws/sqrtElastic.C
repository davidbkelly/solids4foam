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

#include "sqrtElastic.H"
#include "mathematicalConstants.H"

namespace Foam
{
    defineTypeNameAndDebug(sqrtElastic, 0);
    addToRunTimeSelectionTable(tubeLaw1D, sqrtElastic, dictionary);
}


Foam::sqrtElastic::sqrtElastic
(
    const word& vesselName,
    const dictionary& dict,
    const scalar rho,
    const scalar s
)
:
    tubeLaw1D(vesselName, dict, rho),
    uniform_(true),
    rd_(0),
    h_(0),
    E_(dict.get<scalar>("youngsModulus")),
    Pd_(dict.getOrDefault<scalar>("diastolicPressure", 0)),
    Pext_(dict.getOrDefault<scalar>("externalPressure", 0)),
    Ad_(0),
    beta_(0),
    k_(0),
    cCoeff_(0)
{
    if (dict.found("radius"))
    {
        if (dict.found("radiusProximal") || dict.found("radiusDistal"))
        {
            FatalIOErrorInFunction(dict)
                << "Vessel " << vesselName_ << ": give either radius, or"
                << " radiusProximal and radiusDistal" << exit(FatalIOError);
        }
        rd_ = dict.get<scalar>("radius");
    }
    else
    {
        const scalar r0 = dict.get<scalar>("radiusProximal");
        const scalar r1 = dict.getOrDefault<scalar>("radiusDistal", r0);

        uniform_ = mag(r1 - r0) <= SMALL*r0;
        rd_ = uniform_ ? r0 : r0 + (r1 - r0)*s;
    }

    if (dict.found("wallThicknessModel"))
    {
        const word model = dict.get<word>("wallThicknessModel");

        if (dict.found("wallThickness") || model != "ADAN56")
        {
            FatalIOErrorInFunction(dict)
                << "Vessel " << vesselName_ << ": give either wallThickness"
                << " or wallThicknessModel ADAN56" << exit(FatalIOError);
        }
        h_ = ADAN56WallThickness(rd_);
    }
    else
    {
        h_ = dict.get<scalar>("wallThickness");
    }

    if (rd_ <= 0 || h_ <= 0 || E_ <= 0 || rho_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "Vessel " << vesselName_ << ": radius, wallThickness,"
            << " youngsModulus and rho must be positive"
            << exit(FatalIOError);
    }

    Ad_ = constant::mathematical::pi*sqr(rd_);
    beta_ = 4.0/3.0*Foam::sqrt(constant::mathematical::pi)*E_*h_;
    k_ = beta_/Ad_;
    cCoeff_ = Foam::sqrt(k_/(2*rho_));
}


Foam::scalar Foam::sqrtElastic::ADAN56WallThickness(const scalar R)
{
    // Boileau Eq. 17 with b and d in 1/m
    const scalar a = 0.2802, b = -505.3, c = 0.1324, d = -11.14;
    return R*(a*Foam::exp(b*R) + c*Foam::exp(d*R));
}


Foam::scalar Foam::sqrtElastic::pressure(const scalar A) const
{
    return Pext_ + Pd_ + k_*(Foam::sqrt(A) - Foam::sqrt(Ad_));
}


Foam::scalar Foam::sqrtElastic::minPressure() const
{
    return Pext_ + Pd_ - k_*Foam::sqrt(Ad_);
}


Foam::scalar Foam::sqrtElastic::dPdA(const scalar A) const
{
    return 0.5*k_/Foam::sqrt(A);
}


Foam::scalar Foam::sqrtElastic::waveSpeed(const scalar A) const
{
    return cCoeff_*Foam::sqrt(Foam::sqrt(A));
}


Foam::scalar Foam::sqrtElastic::area(const scalar P) const
{
    const scalar sqrtA = Foam::sqrt(Ad_) + (P - Pext_ - Pd_)/k_;

    if (sqrtA <= 0)
    {
        FatalErrorInFunction
            << "Vessel " << vesselName_ << ": no positive area for pressure "
            << P << " Pa" << exit(FatalError);
    }

    return sqr(sqrtA);
}


Foam::scalar Foam::sqrtElastic::fluxPressure(const scalar A) const
{
    return k_/(3*rho_)*A*Foam::sqrt(A);
}


Foam::scalar Foam::sqrtElastic::invariant(const scalar A) const
{
    return 4*waveSpeed(A);
}


Foam::scalar Foam::sqrtElastic::areaFromInvariant(const scalar I) const
{
    if (I <= 0)
    {
        FatalErrorInFunction
            << "Vessel " << vesselName_ << ": non-positive wave speed "
            << I/4 << " m/s from the Riemann invariants" << exit(FatalError);
    }

    return sqr(sqr(0.25*I/cCoeff_));
}


void Foam::sqrtElastic::summary(Ostream& os) const
{
    // Area at zero pressure, Boileau Eq. 16 (with Pext)
    const scalar sqrtA0 = Foam::sqrt(Ad_) - (Pd_ + Pext_)/k_;

    os  << "sqrtElastic" << (uniform_ ? "" : " (here)") << ": rd " << rd_
        << " m, h " << h_ << " m, E " << E_
        << " Pa, Pd " << Pd_ << " Pa, Pext " << Pext_ << " Pa" << nl
        << "    Ad " << Ad_ << " m2, beta " << beta_ << " Pa m, c(Ad) "
        << waveSpeed(Ad_) << " m/s";

    if (sqrtA0 > 0)
    {
        os  << ", A(P=0) " << sqr(sqrtA0) << " m2, c(A(P=0)) "
            << waveSpeed(sqr(sqrtA0)) << " m/s";
    }
    os  << nl;
}

// ************************************************************************* //
