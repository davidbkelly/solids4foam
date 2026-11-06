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

#include "lumped0D.H"
#include "OStringStream.H"

namespace Foam
{
    defineTypeNameAndDebug(lumped0D, 0);
    addToRunTimeSelectionTable(node1D, lumped0D, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lumped0D::lumped0D
(
    const word& name,
    const dictionary& dict,
    PtrList<vessel1D>& vessels,
    const scalar t0
)
:
    node1D(name, dict, vessels),
    net_(new circulationNetwork(dict.subDict("circuit"), t0)),
    port_(),
    portInlet_(false)
{
    requireConnections(1);

    const wordList ports = net_->portNames();

    if (dict.found("port"))
    {
        port_ = dict.get<word>("port");
    }
    else if (ports.size() == 1)
    {
        port_ = ports[0];
    }
    else
    {
        OStringStream os;
        os  << flatOutput(ports);
        fail
        (
            "the circuit must have one externalPort for the vessel end, or"
            " \"port\" must name one. Ports: " + os.str()
        );
    }

    if (!ports.found(port_))
    {
        fail("the circuit has no externalPort '" + port_ + "'");
    }

    if (!net_->portExternal(port_))
    {
        fail("port '" + port_ + "' is table driven; it must be external");
    }

    portInlet_ = net_->portInlet(port_);

    setNStates(net_->nEqns());
    y_ = net_->state();
}


// * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * * //

Foam::scalar Foam::lumped0D::portPressure(const scalar t, const scalar Qp)
{
    externalPort::forcingCoeffs c(Zero);
    c[0] = Qp;
    net_->setPortForcing(port_, t, t, c);

    return net_->portPressure(port_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lumped0D::setBoundary(const scalar t)
{
    const scalar s = sign(0);
    const tubeLaw1D& law = vessel(0).endLaw(end(0));

    net_->setState(t, y_);

    scalar Af;
    const scalar W = outgoingInvariant(0, Af);

    // Flow scale for the finite-difference step of dPport/dQ
    const scalar Qscale = Af*law.waveSpeed(Af);

    const scalar A = solveArea
    (
        [&](const scalar a, scalar& dRdA)
        {
            const scalar U = W - s*law.invariant(a);
            const scalar Qn = s*a*U;
            const scalar dQn = s*(U - s*law.waveSpeed(a));

            const scalar Qp = portFlow(Qn);
            const scalar h = 1e-7*(mag(Qp) + Qscale);
            const scalar P = portPressure(t, Qp);
            const scalar dPdQp = (portPressure(t, Qp + h) - P)/h;

            // dQp/dQn = +-1
            dRdA = law.dPdA(a) - dPdQp*portFlow(dQn);
            return law.pressure(a) - P;
        },
        Af
    );

    const scalar Q = A*(W - s*law.invariant(A));

    // Leave the circuit forced with the converged flow
    portPressure(t, portFlow(s*Q));

    vessel(0).setBoundaryState(end(0), A, Q);
}


void Foam::lumped0D::derivatives(const scalar t)
{
    const scalar Qp = portFlow(sign(0)*vessel(0).boundaryFlow(end(0)));

    externalPort::forcingCoeffs c(Zero);
    c[0] = Qp;
    net_->setPortForcing(port_, t, t, c);

    net_->derivatives(t, y_, dydt_);
}


Foam::wordList Foam::lumped0D::columns() const
{
    return net_->columns();
}


void Foam::lumped0D::outputValues(const scalar t, scalarField& values)
{
    const scalar Qp = portFlow(sign(0)*vessel(0).boundaryFlow(end(0)));

    net_->setState(t, y_);

    externalPort::forcingCoeffs c(Zero);
    c[0] = Qp;
    net_->setPortForcing(port_, t, t, c);

    net_->outputValues(values);
}


void Foam::lumped0D::summary(Ostream& os) const
{
    node1D::summary(os);
    os  << "    port " << port_ << " ("
        << (portInlet_ ? "inlet" : "outlet") << "), ";
    net_->summary(os);
}

// ************************************************************************* //
