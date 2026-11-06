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

#include "node1D.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(node1D, 0);
    defineRunTimeSelectionTable(node1D, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::node1D::node1D
(
    const word& name,
    const dictionary& dict,
    PtrList<vessel1D>& vessels
)
:
    name_(name),
    dict_(dict),
    vessels_(vessels)
{
    wordList terminals;

    if (dict.found("connection"))
    {
        terminals = wordList(1, dict.get<word>("connection"));
    }
    else if (dict.found("connections"))
    {
        terminals = dict.get<wordList>("connections");
    }
    else
    {
        fail("needs a \"connection\" (or \"connections\"), e.g. carotid.in");
    }

    ends_.setSize(terminals.size());

    forAll(terminals, i)
    {
        const word& t = terminals[i];
        const auto dot = t.rfind('.');

        if (dot == std::string::npos)
        {
            fail("connection '" + t + "' is not <vessel>.in or <vessel>.out");
        }

        const word vesselName = t.substr(0, dot);
        const word endName = t.substr(dot + 1);

        label v = -1;
        forAll(vessels_, j)
        {
            if (vessels_[j].name() == vesselName)
            {
                v = j;
                break;
            }
        }

        if (v < 0)
        {
            fail("connection '" + t + "': no vessel '" + vesselName + "'");
        }

        if (endName == "in")
        {
            ends_[i] = {v, vessel1D::inEnd};
        }
        else if (endName == "out")
        {
            ends_[i] = {v, vessel1D::outEnd};
        }
        else
        {
            fail
            (
                "connection '" + t + "': the end must be 'in' or 'out', not '"
              + endName + "'"
            );
        }
    }
}


Foam::autoPtr<Foam::node1D> Foam::node1D::New
(
    const word& name,
    const dictionary& dict,
    PtrList<vessel1D>& vessels,
    const scalar t0
)
{
    const word type(dict.get<word>("type"));

    auto* ctorPtr = dictionaryConstructorTable(type);

    if (!ctorPtr)
    {
        FatalIOErrorInFunction(dict)
            << "Unknown node type '" << type << "' for node '" << name
            << "'. Node types:" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalIOError);
    }

    return autoPtr<node1D>(ctorPtr(name, dict, vessels, t0));
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::node1D::fail(const string& message) const
{
    FatalIOErrorInFunction(dict_)
        << "Node '" << name_ << "': " << message << nl
        << exit(FatalIOError);

    std::abort(); // not reached; satisfies [[noreturn]]
}


void Foam::node1D::requireConnections(const label n) const
{
    if (ends_.size() != n)
    {
        fail
        (
            "a " + std::string(dict_.get<word>("type")) + " node joins "
          + std::to_string(n) + " vessel end(s), not "
          + std::to_string(ends_.size())
        );
    }
}


Foam::scalar Foam::node1D::outgoingInvariant
(
    const label i,
    scalar& Aface
) const
{
    scalar Qface;
    vessel(i).faceState(end(i), Aface, Qface);

    return Qface/Aface + sign(i)*vessel(i).endLaw(end(i)).invariant(Aface);
}


void Foam::node1D::setFlowBoundary(const label i, const scalar Qn)
{
    const scalar s = sign(i);
    const tubeLaw1D& law = vessel(i).endLaw(end(i));

    scalar Af;
    const scalar W = outgoingInvariant(i, Af);

    // Flow along the vessel at the end
    const scalar Qv = s*Qn;

    const scalar A = solveArea
    (
        [&](const scalar a, scalar& dRdA)
        {
            const scalar U = W - s*law.invariant(a);
            dRdA = U - s*law.waveSpeed(a);
            return a*U - Qv;
        },
        Af
    );

    vessel(i).setBoundaryState(end(i), A, Qv);
}


void Foam::node1D::setNStates(const label n)
{
    y_.setSize(n, 0);
    y0_.setSize(n, 0);
    dydt_.setSize(n, 0);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::node1D::setStates(const scalarField& y)
{
    if (y.size() != y_.size())
    {
        fail
        (
            "has " + std::to_string(y_.size()) + " states, given "
          + std::to_string(y.size())
        );
    }
    y_ = y;
}


void Foam::node1D::stageUpdate
(
    const scalar a,
    const scalar b,
    const scalar dt
)
{
    forAll(y_, i)
    {
        y_[i] = a*y0_[i] + b*(y_[i] + dt*dydt_[i]);
    }
}


void Foam::node1D::storeStageRhs(const label k, const label nStages)
{
    if (ky_.size() != nStages)
    {
        ky_.setSize(nStages, scalarField(y_.size(), 0));
    }
    ky_[k] = dydt_;
}


void Foam::node1D::butcherStage
(
    const rungeKutta1D& rk,
    const label k,
    const scalar dt
)
{
    if (y_.empty())
    {
        return;
    }

    const bool last = (k == rk.nStages());

    y_ = y0_;
    for (label j = 0; j < k; ++j)
    {
        const scalar w = dt*(last ? rk.B(j) : rk.A(k, j));
        if (w != 0)
        {
            y_ += w*ky_[j];
        }
    }
}


void Foam::node1D::summary(Ostream& os) const
{
    os  << "Node " << name_ << " (" << dict_.get<word>("type") << "):";
    for (const vesselEnd& e : ends_)
    {
        os  << ' ' << vessels_[e.vessel].name() << '.'
            << vessel1D::endName(e.end);
    }
    os  << nl;
}

// ************************************************************************* //
