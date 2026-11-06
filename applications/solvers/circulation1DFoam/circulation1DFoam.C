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
    circulation1DFoam

Description
    Solves a network of 1D arterial vessels described in
    constant/circulation1DProperties: finite volumes with a MUSCL (second
    order) or k-exact (up to fourth order) reconstruction, an approximate
    Riemann solver and an explicit Runge-Kutta scheme (SSPRK2, SSPRK3, RK4)
    in each vessel, joined by nodes (inflows, terminals, lumped 0D circuits,
    junctions).
    No mesh is needed.

    controlDict supplies startTime, endTime, deltaT and the write times.
    deltaT is the output interval; the network takes its own CFL-limited
    steps inside it and lands on every output time.

    Every deltaT a row is appended to, in postProcessing/circulation1D/
    <startTime>/,

        probes/<probe>.dat      t P Q A U deltaR at each probe
        vessels/<vessel>.dat    t Pin Pout dP Qin Qout volume of each vessel
        nodes/<node>.dat        the circuit of each lumped0D node, as
                                circulationFoam writes it
        balance.dat             t volume netInflow drift: the volume balance

    all SI. At write times the profiles along every vessel are written to
    profiles/<time>/<vessel>.dat, and the state to <time>/circulation1DState,
    which is read from the start time if present, so a run can be continued.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IOdictionary.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "network1D.H"

#include <limits>

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addNote("1D arterial network solver");
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

    // Continue from a stored state if there is one at the start time
    {
        IOobject stateIO
        (
            "circulation1DState",
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (stateIO.typeHeaderOk<IOdictionary>(true))
        {
            const IOdictionary state(stateIO);
            network.readState(state);
            Info<< "Restarted from " << stateIO.objectPath() << nl << endl;
        }
    }

    network.summary(Info);
    Info<< endl;

    const fileName outDir
    (
        runTime.globalPath()/"postProcessing"/"circulation1D"
       /runTime.timeName()
    );
    mkDir(outDir);

    const int precision = std::numeric_limits<scalar>::max_digits10;

    // Probes, vessels and nodes each in their own directory, as their names
    // may coincide
    auto openFile = [&](const fileName& dir, const word& name)
    {
        mkDir(outDir/dir);
        OFstream* os = new OFstream(outDir/dir/(name + ".dat"));
        os->precision(precision);
        return os;
    };

    // Probes
    const wordList& probes = network.probeNames();
    PtrList<OFstream> probeFiles(probes.size());
    forAll(probes, i)
    {
        probeFiles.set(i, openFile("probes", probes[i]));
        probeFiles[i]
            << "# circulation1DFoam probe " << probes[i] << nl
            << "# t P Q A U deltaR" << nl
            << "# s Pa m3/s m2 m/s m" << nl;
    }

    // Vessels
    const PtrList<vessel1D>& vessels = network.vessels();
    PtrList<OFstream> vesselFiles(vessels.size());
    forAll(vessels, v)
    {
        vesselFiles.set(v, openFile("vessels", vessels[v].name()));
        vesselFiles[v]
            << "# circulation1DFoam vessel " << vessels[v].name() << nl
            << "# t Pin Pout dP Qin Qout volume" << nl
            << "# s Pa Pa Pa m3/s m3/s m3" << nl;
    }

    // Nodes with output columns (lumped0D circuits)
    PtrList<node1D>& nodes = network.nodes();
    PtrList<OFstream> nodeFiles(nodes.size());
    forAll(nodes, n)
    {
        const wordList cols = nodes[n].columns();

        if (cols.empty())
        {
            continue;
        }

        nodeFiles.set(n, openFile("nodes", nodes[n].name()));
        nodeFiles[n]
            << "# circulation1DFoam node " << nodes[n].name() << nl << "# t";
        for (const word& c : cols)
        {
            nodeFiles[n] << ' ' << c;
        }
        nodeFiles[n] << nl;
    }

    autoPtr<OFstream> balanceFile(openFile(".", "balance"));
    balanceFile()
        << "# circulation1DFoam volume balance: drift = volume - V0 - netInflow"
        << nl << "# t volume netInflow drift" << nl << "# s m3 m3 m3" << nl;

    scalarField values;

    auto writeRows = [&]()
    {
        const scalar t = network.time();
        scalar P, Q, A, U, dr;

        forAll(probeFiles, i)
        {
            network.probeValues(i, P, Q, A, U, dr);
            probeFiles[i]
                << t << ' ' << P << ' ' << Q << ' ' << A << ' ' << U << ' '
                << dr << nl;
        }

        scalar Pin, Pout, Qin, Qout;
        forAll(vesselFiles, v)
        {
            network.vesselEndValues(v, Pin, Pout, Qin, Qout);
            vesselFiles[v]
                << t << ' ' << Pin << ' ' << Pout << ' ' << Pin - Pout << ' '
                << Qin << ' ' << Qout << ' ' << vessels[v].volume() << nl;
        }

        forAll(nodeFiles, n)
        {
            if (nodeFiles.set(n))
            {
                nodes[n].outputValues(t, values);
                nodeFiles[n] << t;
                for (const scalar x : values)
                {
                    nodeFiles[n] << ' ' << x;
                }
                nodeFiles[n] << nl;
            }
        }

        balanceFile()
            << t << ' ' << network.volume() << ' ' << network.netInflow()
            << ' ' << network.volumeDrift() << nl;
    };

    auto writeProfiles = [&]()
    {
        const fileName dir(outDir/"profiles"/runTime.timeName());
        mkDir(dir);

        forAll(vessels, v)
        {
            OFstream os(dir/(vessels[v].name() + ".dat"));
            os.precision(precision);
            network.writeProfile(v, os);
        }
    };

    auto writeState = [&]()
    {
        // The restart state at full double precision whatever writePrecision
        // the case uses. The values are converted to text when the
        // dictionary is built, so the precision is raised before that.
        const unsigned int oldPrecision =
            IOstream::defaultPrecision(precision);

        IOdictionary state
        (
            IOobject
            (
                "circulation1DState",
                runTime.timeName(),
                runTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            network.stateDict()
        );
        state.regIOobject::writeObject
        (
            IOstreamOption(IOstreamOption::ASCII),
            true
        );

        IOstream::defaultPrecision(oldPrecision);
    };

    writeRows();
    writeProfiles();

    Info<< "Starting time loop\n" << endl;

    while (runTime.run())
    {
        const scalar t0 = runTime.value();
        runTime++;

        network.advance(t0, runTime.value());
        writeRows();

        if (runTime.writeTime())
        {
            writeProfiles();
            writeState();

            Info<< "Time = " << runTime.timeName()
                << "  steps " << network.nSteps()
                << "  dt " << network.lastDeltaT()
                << "  volume drift " << network.volumeDrift() << " m3 ("
                << network.volumeDrift()/network.initialVolume()
                << " of V0)" << endl;

            forAll(probeFiles, i)
            {
                probeFiles[i].flush();
            }
        }
    }

    Info<< nl << "Steps " << network.nSteps() << ", volume drift "
        << network.volumeDrift() << " m3 ("
        << network.volumeDrift()/network.initialVolume() << " of V0";
    if (mag(network.netInflow()) > SMALL*network.initialVolume())
    {
        Info<< ", " << network.volumeDrift()/network.netInflow()
            << " of the net inflow";
    }
    Info<< ")" << nl
        << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << nl
        << nl << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
