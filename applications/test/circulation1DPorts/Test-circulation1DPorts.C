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

Application
    Test-circulation1DPorts

Description
    Drives the coupledFlow ports of a 1D network as the 3D-1D coupler will,
    with a surrogate in place of the CFD: a rigid pipe in plug flow, the
    exact reduction of a rigid 3D segment,

        L dQ/dt = Pin - Pout,       L = rho length/area,

    integrated with the 3D time step deltaT by a BDF formula of order 1 to 4
    (the CFD's "Euler", "backward", "BDF 3", "BDF 4"). The network has two
    coupledFlow ports, one at a
    3D inlet (the pipe's upstream end) and one at a 3D outlet; both carry the
    pipe flow Q.

    The forcing polynomials are those of circulationCoupler:

        explicit  order p: Q(t0) and p - 1 earlier flows, extrapolated
        implicit  order 1: the iterate, held constant
                  order 2+: Q(t0) at alpha 0, the iterate at alpha 1, and
                            p - 2 earlier flows

    Implicit coupling is a plain fixed-point iteration (as PIMPLE outer
    correctors), to a change of the flow below "tolerance" times the largest
    flow so far, or "nIter" iterations. Its contraction rate (the third
    change over the second) is written for each step: for this problem it
    is dt/(a0 tau), a0 the leading BDF coefficient (1, 3/2, 11/6, 25/12) and
    tau = L/(2 Zc). The system starts at rest, so the history is seeded with
    zeros ("seed zero") and every order is used from the first step.

    In constant/circulation1DProperties, besides the network:

        surrogate3D
        {
            length      0.08;
            radius      0.01;
            ddtScheme   backward;       // Euler | backward | BDF3 | BDF4
        }

    and in system/controlDict

        coupling { scheme implicit; order 2; nIter 50; tolerance 1e-10; }

    deltaT of controlDict is the 3D time step. Writes
    postProcessing/surrogate3D.dat: t Q Pin Pout iterations rate. Stops with a
    message (exit status 2) if the flow grows beyond "maxFlow" (default 1).

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IOdictionary.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "mathematicalConstants.H"
#include "network1D.H"

#include <limits>

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addNote("3D-1D port test with a rigid-pipe surrogate for the CFD");
    argList::noParallel();

    #include "setRootCase.H"
    #include "createTime.H"

    IOdictionary props
    (
        IOobject
        (
            "circulation1DProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    network1D network(props, runTime.value());

    // ---- Surrogate and coupling settings -----------------------------------

    const dictionary& sd = props.subDict("surrogate3D");
    const scalar length = sd.get<scalar>("length");
    const scalar area =
        constant::mathematical::pi*sqr(sd.get<scalar>("radius"));
    // BDF coefficients: ddt(Q)^(n+1) = (1/dt) sum_k a_k Q^(n+1-k)
    const word ddt = sd.get<word>("ddtScheme");
    scalarList a;
    if (ddt == "Euler")
    {
        a = scalarList({1, -1});
    }
    else if (ddt == "backward")
    {
        a = scalarList({1.5, -2, 0.5});
    }
    else if (ddt == "BDF3")
    {
        a = scalarList({11.0/6, -3, 1.5, -1.0/3});
    }
    else if (ddt == "BDF4")
    {
        a = scalarList({25.0/12, -4, 3, -4.0/3, 0.25});
    }
    else
    {
        FatalIOErrorInFunction(sd)
            << "ddtScheme must be Euler, backward, BDF3 or BDF4"
            << exit(FatalIOError);
    }

    const scalar rho = props.subDict("fluid").get<scalar>("rho");
    const scalar L = rho*length/area;

    const dictionary& cd = runTime.controlDict().subDict("coupling");
    const bool imp = (cd.get<word>("scheme") == "implicit");
    const label order = cd.get<label>("order");
    const label nIter = cd.getOrDefault<label>("nIter", 50);
    const scalar tol = cd.getOrDefault<scalar>("tolerance", 1e-10);
    const scalar maxFlow = cd.getOrDefault<scalar>("maxFlow", 1);

    if (order < 1 || order > 4)
    {
        FatalIOErrorInFunction(cd) << "order 1..4" << exit(FatalIOError);
    }

    // The two ports: the 3D inlet and the 3D outlet
    word inPort, outPort;
    for (const word& p : network.portNames())
    {
        (network.portInlet(p) ? inPort : outPort) = p;
    }
    if (inPort.empty() || outPort.empty() || network.portNames().size() != 2)
    {
        FatalErrorInFunction
            << "Need one coupledFlow port at a vessel out end (3D inlet) and"
            << " one at an in end (3D outlet); ports: "
            << network.portNames() << exit(FatalError);
    }

    network.summary(Info);
    Info<< "Surrogate: rigid pipe, L = " << L << " kg/m4, " << ddt
        << "; coupling " << (imp ? "implicit" : "explicit") << " order "
        << order << ", deltaT " << runTime.deltaTValue() << nl << endl;

    // ---- Output ------------------------------------------------------------

    mkDir(runTime.globalPath()/"postProcessing");
    OFstream os(runTime.globalPath()/"postProcessing"/"surrogate3D.dat");
    os.precision(std::numeric_limits<scalar>::max_digits10);
    os  << "# t Q Pin Pout iterations rate" << nl;

    // ---- Time loop ---------------------------------------------------------

    // Accepted pipe flows: Q at the step start, and three earlier (zero: at
    // rest)
    scalar Qn = 0;
    FixedList<scalar, 3> history(0.0);
    // Pipe flows Q^n, Q^(n-1), Q^(n-2), Q^(n-3) for the BDF formula
    FixedList<scalar, 4> pipe(0.0);

    label totalIter = 0;

    // Largest accepted flow, the scale of the iteration tolerance
    scalar Qscale = VSMALL;

    while (runTime.run())
    {
        const scalar t0 = runTime.value();
        ++runTime;
        const scalar t1 = runTime.value();
        const scalar dt = t1 - t0;

        const scalarField y0 = network.state();

        auto forcing = [&](const scalar Qiter)
        {
            scalarList alpha, Q;
            if (!imp)
            {
                alpha.setSize(order);
                Q.setSize(order);
                alpha[0] = 0;
                Q[0] = Qn;
                for (label j = 1; j < order; ++j)
                {
                    alpha[j] = -j;
                    Q[j] = history[j-1];
                }
            }
            else if (order == 1)
            {
                alpha.setSize(1, scalar(1));
                Q.setSize(1, Qiter);
            }
            else
            {
                alpha.setSize(order);
                Q.setSize(order);
                alpha[0] = 1;
                Q[0] = Qiter;
                alpha[1] = 0;
                Q[1] = Qn;
                for (label j = 2; j < order; ++j)
                {
                    alpha[j] = 1 - j;
                    Q[j] = history[j-2];
                }
            }
            return circulationForcing::polynomialThrough(alpha, Q);
        };

        // One exchange: integrate the network over the step with the
        // forcing of the iterate, then advance the pipe with the pressures
        scalar Pin = 0, Pout = 0;
        auto exchange = [&](const scalar Qiter)
        {
            const circulationForcing::coeffs c = forcing(Qiter);
            network.setState(t0, y0);
            network.setPortForcing(inPort, t0, t1, c);
            network.setPortForcing(outPort, t0, t1, c);
            network.advance(t0, t1);

            Pin = network.portPressure(inPort);
            Pout = network.portPressure(outPort);

            // L ddt(Q) = Pin - Pout at t1
            scalar rhs = dt*(Pin - Pout)/L;
            for (label k = 1; k < a.size(); ++k)
            {
                rhs -= a[k]*pipe[k-1];
            }
            return rhs/a[0];
        };

        scalar Q = Qn;
        label iter = 0;
        scalarList changes;

        if (!imp)
        {
            Q = exchange(Qn);
            iter = 1;
        }
        else
        {
            for (iter = 1; iter <= nIter; ++iter)
            {
                const scalar Qnew = exchange(Q);
                const scalar change = mag(Qnew - Q);
                changes.append(change);
                Q = Qnew;
                if
                (
                    !std::isfinite(Q)
                 || change <= tol*Foam::max(mag(Q), Qscale)
                )
                {
                    break;
                }
            }
            iter = Foam::min(iter, nIter);
        }

        totalIter += iter;
        Qscale = Foam::max(Qscale, mag(Q));

        const scalar rate =
            changes.size() > 2 && changes[1] > 0 ? changes[2]/changes[1] : 0;

        os  << t1 << ' ' << Q << ' ' << Pin << ' ' << Pout << ' ' << iter
            << ' ' << rate << nl;

        if (!std::isfinite(Q) || mag(Q) > maxFlow)
        {
            Info<< "UNSTABLE at t = " << t1 << ": Q = " << Q << endl;
            os.flush();
            return 2;
        }

        for (label k = 3; k > 0; --k)
        {
            pipe[k] = pipe[k-1];
        }
        pipe[0] = Q;
        history[2] = history[1];
        history[1] = history[0];
        history[0] = Qn;
        Qn = Q;
    }

    Info<< "Steps " << runTime.timeIndex() << ", exchanges " << totalIter
        << ", 1D steps " << network.nSteps() << nl
        << "volume drift " << network.volumeDrift() << " m3" << nl
        << nl << "End" << nl << endl;

    return 0;
}

// ************************************************************************* //
