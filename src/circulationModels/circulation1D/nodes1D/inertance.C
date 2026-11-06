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

#include "inertance.H"
#include "mathematicalConstants.H"

namespace Foam
{
    defineTypeNameAndDebug(inertance, 0);
    addToRunTimeSelectionTable(node1D, inertance, dictionary);
}


Foam::inertance::inertance
(
    const word& name,
    const dictionary& dict,
    PtrList<vessel1D>& vessels,
    const scalar t0
)
:
    node1D(name, dict, vessels),
    length_(dict.get<scalar>("length")),
    area_(0),
    R_(dict.getOrDefault<scalar>("R", 0)),
    total_(false),
    Pout_(dict.getOrDefault<scalar>("Pout", 0))
{
    if (ends_.size() != 1 && ends_.size() != 2)
    {
        fail("an inertance joins one or two vessel ends");
    }

    if (dict.found("radius") == dict.found("area"))
    {
        fail("give the pipe's radius or its area, not both or neither");
    }

    area_ =
        dict.found("area")
      ? dict.get<scalar>("area")
      : constant::mathematical::pi*sqr(dict.get<scalar>("radius"));

    if (length_ <= 0 || area_ <= 0 || R_ < 0)
    {
        fail("length and the radius or area must be positive, R not negative");
    }

    const word p = dict.getOrDefault<word>("pressure", "static");
    if (p == "total")
    {
        total_ = true;
    }
    else if (p != "static")
    {
        fail("pressure must be static or total, not '" + p + "'");
    }

    setNStates(1);
    y_[0] = dict.getOrDefault<scalar>("Q0", 0);
}


Foam::scalar Foam::inertance::endPressure
(
    const label i,
    const scalar Q
) const
{
    if (i >= ends_.size())
    {
        return Pout_;
    }

    const vessel1D& v = vessel(i);
    const scalar A = v.boundaryArea(end(i));
    const scalar P = v.endLaw(end(i)).pressure(A);

    return total_ ? P + 0.5*v.rho()*sqr(Q/A) : P;
}


void Foam::inertance::setBoundary(const scalar t)
{
    // Q leaves the first vessel into the pipe and enters the second
    setFlowBoundary(0, y_[0]);
    if (ends_.size() == 2)
    {
        setFlowBoundary(1, -y_[0]);
    }
}


void Foam::inertance::derivatives(const scalar t)
{
    const scalar Q = y_[0];

    dydt_[0] = (endPressure(0, Q) - endPressure(1, Q) - R_*Q)/L();
}


Foam::wordList Foam::inertance::columns() const
{
    return wordList({"Q", "P1", "P2"});
}


void Foam::inertance::outputValues(const scalar t, scalarField& values)
{
    const scalar Q = y_[0];

    values.setSize(3);
    values[0] = Q;
    values[1] = endPressure(0, Q);
    values[2] = endPressure(1, Q);
}


void Foam::inertance::summary(Ostream& os) const
{
    node1D::summary(os);
    os  << "    rigid pipe: length " << length_ << " m, area " << area_
        << " m2, L " << L() << " kg/m4, R " << R_ << " Pa s/m3, "
        << (total_ ? "total" : "static") << " pressure continuity";
    if (ends_.size() == 1)
    {
        os  << ", open to Pout " << Pout_ << " Pa";
    }
    os  << nl;
}

// ************************************************************************* //
