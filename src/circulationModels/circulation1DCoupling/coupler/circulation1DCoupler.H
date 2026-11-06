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

Class
    Foam::circulation1DCoupler

Description
    Couples the patches of a CFD case to a 1D arterial network (network1D),
    through the network's coupledFlow ports: the CFD supplies the flow at each
    coupled patch, the network the pressure to impose there.

    A copy of circulationCoupler (circulationCoupling), with network1D in
    place of the 0D circulationNetwork; the coupling schemes are the same.
    Temporary: the two are to share one coupler when circulation1D is merged.

    In constant/circulation1DProperties, besides the network:

        coupling
        {
            scheme      implicit;   // explicit | implicit
            order       2;          // 1..4: degree of the forcing + 1
            maxPressure 1e8;        // optional [Pa]
            seed        zero;       // optional, as for circulationCoupler
        }

    The density is the network's (fluid { rho ...; }); the CFD solves p/rho.

    The network sees each port's flow over a step [t0, t1] as a polynomial in
    alpha = (t - t0)/(t1 - t0), through samples of the flow:

    \table
        scheme    | order | samples (alpha)
        explicit  | 1     | Q(t0) held constant
        explicit  | 2     | Q(t0), Q(t0-dt)             extrapolated to t1
        explicit  | 3, 4  | Q(t0) and two / three earlier flows, extrapolated
        implicit  | 1     | the current outer iterate, held constant
        implicit  | 2     | Q(t0) at alpha 0, the iterate at alpha 1
        implicit  | 3, 4  | the above and one / two earlier flows
    \endtable

    Explicit coupling integrates the network once per time step. Implicit
    coupling re-integrates from the start of the step each time the CFD flux
    is updated, which makes the outer correctors of PIMPLE a fixed-point
    iteration on the interface flow.

    Within one step each vessel end answers the port flow as a resistor of
    its characteristic impedance Zc = rho c/A. For a 3D domain of inertance L
    between two such ports the coupling number is dt/tau, tau = L/(2 Zc);
    explicit coupling of order p is stable for dt/tau below about
    4/(2^p - 1) with BDF2 in the CFD (2/(2^p - 1) with Euler), and the
    implicit iteration contracts by dt/(1.5 tau) (dt/tau). See
    src/circulationModels/docs/BENCHMARKS3D1D.md.

    Ports are bound to patches by the "patch" entry of each coupledFlow node.
    Sign convention: the port flow is the flow along the vessel, so a port at
    a vessel's out end is a 3D inlet (Q enters the domain) and one at an in
    end a 3D outlet (Q leaves it).

    Results: postProcessing/circulation1DCoupling/<start>/coupling.dat, one row
    per accepted step: t, updates, then Q and P of each port (Q the forcing
    at the end of the step), P Q of each probe of the network, and the flux
    the CFD delivered at each port (<port>.Qsampled: for explicit coupling
    the flux at the start of the step, sampled once); and at every write time of the CFD the profiles
    along every vessel (cell averages, as circulation1DFoam writes them) in
    profiles/<time>/<vessel>.dat.

    Restart: the network state and the flow history are written with the
    fields to <time>/circulation1DState-<properties name> and read back when
    the case starts from a time that has one.

Author
    David Kelly, UCD.

SourceFiles
    circulation1DCoupler.C

\*---------------------------------------------------------------------------*/

#ifndef circulation1DCoupler_H
#define circulation1DCoupler_H

#include "regIOobject.H"
#include "network1D.H"
#include "autoPtr.H"
#include "OFstream.H"
#include "localIOdictionary.H"
#include "IOstreamOption.H"

namespace Foam
{

class fvMesh;
class Time;

//- The restart checkpoint: a dictionary kept up to date at every exchange and
//  written automatically, at full precision, with every output time
class circulation1DCheckpoint
:
    public localIOdictionary
{
public:

    circulation1DCheckpoint(const fvMesh& mesh, const word& name);

    //- True if a checkpoint was read at construction
    bool resumed() const { return resumed_; }

    virtual bool writeObject
    (
        IOstreamOption streamOpt,
        const bool writeOnProc
    ) const;

private:

    bool resumed_;
};


class circulation1DCoupler
:
    public regIOobject
{
public:

    enum class scheme { explicitScheme, implicitScheme };

private:

    // Private Data

        const fvMesh& mesh_;

        scheme scheme_;
        label order_;

        //- Pressure at a port beyond which the coupling is taken to have
        //  diverged [Pa]
        scalar maxPressure_;

        autoPtr<network1D> network_;

        //- One coupled patch
        struct portData
        {
            word name;
            label patchID;
            bool inlet;
            scalar flowScale;              //!< patch flow to vessel flow
            scalar Q;                      //!< flow most recently sampled
            scalar stepStartQ;             //!< accepted flow at the step start
            FixedList<scalar, 3> history;  //!< earlier flows, most recent first
            scalar P;                      //!< pressure to impose [Pa]
        };

        List<portData> ports_;

        //- Interval belonging to each history entry, and how many are valid
        FixedList<scalar, 3> historyDt_;
        label historyDepth_;

        //- Accepted network state and its time; and at the start of this step
        scalarField y_;
        scalar lastT_;
        scalarField stepStartY_;
        scalar stepStartT_;

        label activeTimeIndex_;
        label lastConsumedFluxEvent_;
        label stepUpdateCount_;
        label totalUpdateCount_;

        //- Restart checkpoint, <time>/circulation1DState-<name>
        autoPtr<circulation1DCheckpoint> checkpoint_;

        autoPtr<OFstream> log_;

        //- Profiles are written when the step now running is accepted (it is
        //  a write time of the CFD)
        bool profilesPending_;


    // Private Member Functions

        void beginTimeIndex(const label timeIndex);

        void advance(const Time& runTime);

        void writeRow();

        //- Profiles along every vessel at the accepted time
        void writeProfiles();

        //- Store what a restart from the end of the step just accepted needs
        void updateCheckpoint();

        //- Resume from a checkpoint found at the start time
        void restore();

        circulation1DCoupler(const circulation1DCoupler&) = delete;
        void operator=(const circulation1DCoupler&) = delete;


public:

    TypeName("circulation1DCoupler");

    //- Construct reading the network description "propertiesName" from
    //  constant/
    circulation1DCoupler
    (
        const word& propertiesName,
        const fvMesh& mesh
    );

    virtual ~circulation1DCoupler();

    //- The coupler registered with the mesh, constructed on first use
    static circulation1DCoupler& lookupOrCreate
    (
        const fvMesh& mesh,
        const word& propertiesName = "circulation1DProperties"
    );

    //- Name of the port coupled to a patch; fatal if there is none
    const word& portOfPatch(const label patchID) const;

    //- Sample the flows and, on a new exchange, advance the network
    void sampleAndAdvance(const fvMesh& mesh);

    //- Pressure to impose at a port [Pa]
    scalar patchPressure(const word& port) const;

    scalar rho() const { return network_().rho(); }

    const network1D& network() const { return network_(); }

    virtual bool writeData(Ostream&) const { return true; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
