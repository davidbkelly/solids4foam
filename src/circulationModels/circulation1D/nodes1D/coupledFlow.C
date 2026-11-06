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

#include "coupledFlow.H"

namespace Foam
{
    defineTypeNameAndDebug(coupledFlow, 0);
    addToRunTimeSelectionTable(node1D, coupledFlow, dictionary);
}


Foam::coupledFlow::coupledFlow
(
    const word& name,
    const dictionary& dict,
    PtrList<vessel1D>& vessels,
    const scalar t0
)
:
    node1D(name, dict, vessels),
    patch_(dict.getOrDefault<word>("patch", word::null)),
    flowScale_(dict.getOrDefault<scalar>("flowScale", 1)),
    coeffs_(Zero),
    t0_(t0),
    t1_(t0)
{
    requireConnections(1);

    if (flowScale_ <= 0)
    {
        fail("flowScale must be positive");
    }
}


void Foam::coupledFlow::setForcing
(
    const scalar t0,
    const scalar t1,
    const circulationForcing::coeffs& c
)
{
    if (t1 < t0)
    {
        fail
        (
            "forcing interval [" + std::to_string(t0) + ", "
          + std::to_string(t1) + "] is reversed"
        );
    }

    t0_ = t0;
    t1_ = t1;
    coeffs_ = c;
}


Foam::scalar Foam::coupledFlow::flow(const scalar t) const
{
    const scalar alpha = t1_ > t0_ ? (t - t0_)/(t1_ - t0_) : 0;

    return circulationForcing::evaluate(coeffs_, alpha);
}


Foam::scalar Foam::coupledFlow::pressure() const
{
    return vessel(0).endLaw(end(0)).pressure(vessel(0).boundaryArea(end(0)));
}


void Foam::coupledFlow::setBoundary(const scalar t)
{
    // Q along the vessel; the flow into the node is s Q
    setFlowBoundary(0, sign(0)*flow(t));
}


Foam::wordList Foam::coupledFlow::columns() const
{
    return wordList({"Q", "P"});
}


void Foam::coupledFlow::outputValues(const scalar t, scalarField& values)
{
    values.setSize(2);
    values[0] = vessel(0).boundaryFlow(end(0));
    values[1] = pressure();
}


void Foam::coupledFlow::summary(Ostream& os) const
{
    node1D::summary(os);
    os  << "    coupled to patch "
        << (patch_.empty() ? word("(none)") : patch_) << ", a 3D "
        << (portInlet() ? "inlet" : "outlet");
    if (flowScale_ != 1)
    {
        os  << ", patch flow x " << flowScale_;
    }
    os  << nl;
}

// ************************************************************************* //
