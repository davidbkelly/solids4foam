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

#include "network1D.H"
#include "mathematicalConstants.H"
#include "ITstream.H"
#include "SubField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::network1D::network1D(const dictionary& dict, const scalar t0)
:
    rho_(dict.subDict("fluid").get<scalar>("rho")),
    mu_(dict.subDict("fluid").get<scalar>("mu")),
    profileOrder_(dict.subDict("fluid").get<scalar>("profileOrder")),
    alpha_(dict.subDict("fluid").getOrDefault<scalar>("alpha", 1)),
    CFL_(dict.subOrEmptyDict("numerics").getOrDefault<scalar>("CFL", 0.9)),
    limiter_(slopeLimiter1D::New(dict.subOrEmptyDict("numerics"))),
    flux_(riemannFlux1D::New(dict.subOrEmptyDict("numerics"))),
    rk_
    (
        dict.subOrEmptyDict("numerics").getOrDefault<word>
        (
            "timeScheme",
            "SSPRK2"
        )
    ),
    t_(t0),
    V0_(0),
    inflow_(0),
    inflow0_(0),
    dInflowdt_(0),
    nSteps_(0),
    lastDeltaT_(0),
    frictionWarned_(false)
{
    const dictionary& fluidDict = dict.subDict("fluid");

    if (rho_ <= 0 || mu_ < 0)
    {
        FatalIOErrorInFunction(fluidDict)
            << "rho must be positive and mu non-negative"
            << exit(FatalIOError);
    }

    if (mag(alpha_ - 1) > SMALL)
    {
        FatalIOErrorInFunction(fluidDict)
            << "alpha = " << alpha_ << " is not supported yet: the boundary"
            << " nodes use the Riemann invariants U +- I(A), which hold for"
            << " alpha = 1 only." << exit(FatalIOError);
    }

    if (CFL_ <= 0 || CFL_ > 1)
    {
        FatalIOErrorInFunction(dict)
            << "numerics CFL must be in (0, 1], not " << CFL_
            << exit(FatalIOError);
    }

    {
        const dictionary& nd = dict.subOrEmptyDict("numerics");
        const word rec = nd.getOrDefault<word>("reconstruction", "MUSCL");

        if (rec == "kExact")
        {
            kExact_.reset
            (
                new kExactReconstruction1D(nd.get<label>("polynomialOrder"))
            );
        }
        else if (rec != "MUSCL")
        {
            FatalIOErrorInFunction(nd)
                << "reconstruction must be MUSCL or kExact, not '" << rec
                << "'" << exit(FatalIOError);
        }
    }

    // A limited MUSCL reconstruction is TVD under forward Euler for
    // CFL <= 1/2, and so under the SSP schemes (SSP coefficient 1); beyond
    // that the limiter's clipping of extrema accumulates step by step
    if (!kExact_.valid() && limiter_->type() != "none" && CFL_ > 0.5)
    {
        WarningInFunction
            << "CFL " << CFL_ << " with limiter " << limiter_->type()
            << " is above the TVD bound of 0.5 of the limited MUSCL"
            << " reconstruction: smooth peaks are eroded, increasingly so as"
            << " the mesh is refined. Use CFL <= 0.5." << endl;
    }

    const scalar KR =
        2*constant::mathematical::pi*(profileOrder_ + 2)*mu_/rho_;

    // ---- Vessels -----------------------------------------------------------

    const dictionary& vesselsDict = dict.subDict("vessels");
    const dictionary& vesselDefaults = dict.subOrEmptyDict("vesselDefaults");

    for (const entry& e : vesselsDict)
    {
        if (!e.isDict())
        {
            FatalIOErrorInFunction(vesselsDict)
                << "Vessel '" << e.keyword() << "' must be a dictionary"
                << exit(FatalIOError);
        }

        if (e.keyword().find('.') != std::string::npos)
        {
            FatalIOErrorInFunction(vesselsDict)
                << "Vessel name '" << e.keyword() << "' must not contain '.'"
                << exit(FatalIOError);
        }

        // The vessel's own entries over the shared defaults
        dictionary vesselDict(vesselDefaults);
        vesselDict.merge(e.dict());

        vessels_.append
        (
            new vessel1D(e.keyword(), vesselDict, rho_, alpha_, KR)
        );
    }

    if (vessels_.empty())
    {
        FatalIOErrorInFunction(vesselsDict) << "No vessels defined"
            << exit(FatalIOError);
    }

    if (kExact_.valid())
    {
        for (const vessel1D& v : vessels_)
        {
            if (!v.uniform() || v.nCells() < kExact_->width())
            {
                FatalIOErrorInFunction(vesselsDict)
                    << "Vessel " << v.name() << ": the kExact reconstruction"
                    << " of polynomialOrder " << kExact_->polynomialOrder()
                    << " needs a uniform vessel of at least "
                    << kExact_->width() << " cells" << exit(FatalIOError);
            }
        }
    }

    // ---- Nodes -------------------------------------------------------------

    const dictionary& nodesDict = dict.subDict("nodes");

    for (const entry& e : nodesDict)
    {
        if (!e.isDict())
        {
            FatalIOErrorInFunction(nodesDict)
                << "Node '" << e.keyword() << "' must be a dictionary"
                << exit(FatalIOError);
        }

        nodes_.append
        (
            node1D::New(e.keyword(), e.dict(), vessels_, t0).ptr()
        );
    }

    // Every vessel end in exactly one node
    {
        labelList owner(2*vessels_.size(), label(-1));

        forAll(nodes_, n)
        {
            for (const node1D::vesselEnd& ve : nodes_[n].ends())
            {
                label& o = owner[2*ve.vessel + ve.end];

                if (o >= 0)
                {
                    FatalIOErrorInFunction(nodesDict)
                        << "Vessel end " << vessels_[ve.vessel].name() << '.'
                        << vessel1D::endName(ve.end) << " belongs to nodes "
                        << nodes_[o].name() << " and " << nodes_[n].name()
                        << exit(FatalIOError);
                }
                o = n;
            }
        }

        forAll(owner, i)
        {
            if (owner[i] < 0)
            {
                FatalIOErrorInFunction(nodesDict)
                    << "Vessel end " << vessels_[i/2].name() << '.'
                    << vessel1D::endName(i % 2) << " belongs to no node"
                    << exit(FatalIOError);
            }
        }
    }

    // ---- Probes ------------------------------------------------------------

    const dictionary& probesDict = dict.subOrEmptyDict("probes");

    for (const entry& e : probesDict)
    {
        if (!e.isDict())
        {
            FatalIOErrorInFunction(probesDict)
                << "Probe '" << e.keyword() << "' must be a dictionary"
                << exit(FatalIOError);
        }

        const dictionary& pd = e.dict();
        const word vesselName = pd.get<word>("vessel");

        label v = -1;
        forAll(vessels_, j)
        {
            if (vessels_[j].name() == vesselName)
            {
                v = j;
            }
        }

        if (v < 0)
        {
            FatalIOErrorInFunction(pd)
                << "Probe '" << e.keyword() << "': no vessel '" << vesselName
                << "'" << exit(FatalIOError);
        }

        const scalar L = vessels_[v].length();
        scalar x = 0;

        ITstream& is = pd.lookup("x");
        const token tok(is);

        if (tok.isWord())
        {
            if (tok.wordToken() == "midpoint")
            {
                x = 0.5*L;
            }
            else if (tok.wordToken() == "in")
            {
                x = 0;
            }
            else if (tok.wordToken() == "out")
            {
                x = L;
            }
            else
            {
                FatalIOErrorInFunction(pd)
                    << "Probe '" << e.keyword() << "': x must be a distance"
                    << " [m], midpoint, in or out" << exit(FatalIOError);
            }
        }
        else if (tok.isNumber())
        {
            x = tok.number();
        }
        else
        {
            FatalIOErrorInFunction(pd)
                << "Probe '" << e.keyword() << "': x must be a distance"
                << " [m], midpoint, in or out" << exit(FatalIOError);
        }

        if (x < 0 || x > L)
        {
            FatalIOErrorInFunction(pd)
                << "Probe '" << e.keyword() << "': x = " << x
                << " is outside the vessel [0, " << L << "]"
                << exit(FatalIOError);
        }

        probeNames_.append(e.keyword());
        probeVessel_.append(v);
        probeX_.append(x);
    }

    V0_ = volume();

    updateBoundaries();
}


// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void Foam::network1D::reconstructVessels()
{
    for (vessel1D& v : vessels_)
    {
        if (kExact_.valid())
        {
            v.reconstruct(kExact_());
        }
        else
        {
            v.reconstruct(limiter_());
        }
    }
}


void Foam::network1D::evaluateStage(const scalar t)
{
    reconstructVessels();

    for (node1D& n : nodes_)
    {
        n.setBoundary(t);
    }

    dInflowdt_ = 0;

    for (vessel1D& v : vessels_)
    {
        v.evaluateRhs(flux_());

        dInflowdt_ +=
            v.boundaryFlow(vessel1D::inEnd) - v.boundaryFlow(vessel1D::outEnd);
    }

    for (node1D& n : nodes_)
    {
        n.derivatives(t);
    }
}


Foam::scalar Foam::network1D::stableDeltaT()
{
    scalar dt = GREAT;
    scalar frictionNumber = 0;

    for (const vessel1D& v : vessels_)
    {
        const scalar dtv = CFL_*v.dx()/v.maxWaveSpeed();
        dt = Foam::min(dt, dtv);
    }

    for (const vessel1D& v : vessels_)
    {
        frictionNumber =
            Foam::max(frictionNumber, v.frictionCoeff()*dt/v.minArea());
    }

    if (frictionNumber > 0.1 && !frictionWarned_)
    {
        WarningInFunction
            << "The friction number KR dt/A = " << frictionNumber
            << " exceeds 0.1 at t = " << t_ << ": the explicit friction"
            << " source may be inaccurate or unstable. Reduce the CFL number"
            << " or dx." << nl << "    (warned once)" << endl;
        frictionWarned_ = true;
    }

    return dt;
}


void Foam::network1D::butcherStep(const scalar dt)
{
    const label ns = rk_.nStages();
    kInflow_.setSize(ns, 0);

    auto checkArea = [&]()
    {
        for (const vessel1D& v : vessels_)
        {
            if (v.minArea() <= 0)
            {
                FatalErrorInFunction
                    << "Non-positive area in vessel " << v.name()
                    << " at t = " << t_ << " (step " << dt << ")"
                    << exit(FatalError);
            }
        }
    };

    for (label k = 0; k < ns; ++k)
    {
        if (k > 0)
        {
            for (vessel1D& v : vessels_)
            {
                v.butcherStage(rk_, k, dt);
            }
            for (node1D& n : nodes_)
            {
                n.butcherStage(rk_, k, dt);
            }
            inflow_ = inflow0_;
            for (label j = 0; j < k; ++j)
            {
                inflow_ += dt*rk_.A(k, j)*kInflow_[j];
            }
            checkArea();
        }

        evaluateStage(t_ + rk_.c(k)*dt);

        for (vessel1D& v : vessels_)
        {
            v.storeStageRhs(k, ns);
        }
        for (node1D& n : nodes_)
        {
            n.storeStageRhs(k, ns);
        }
        kInflow_[k] = dInflowdt_;
    }

    // u^(n+1) = u^n + dt sum_k B_k L_k
    for (vessel1D& v : vessels_)
    {
        v.butcherStage(rk_, ns, dt);
    }
    for (node1D& n : nodes_)
    {
        n.butcherStage(rk_, ns, dt);
    }
    inflow_ = inflow0_;
    for (label j = 0; j < ns; ++j)
    {
        inflow_ += dt*rk_.B(j)*kInflow_[j];
    }
    checkArea();
}


void Foam::network1D::step(const scalar dt)
{
    for (vessel1D& v : vessels_)
    {
        v.storeOld();
    }
    for (node1D& n : nodes_)
    {
        n.storeOld();
    }
    inflow0_ = inflow_;

    if (rk_.butcher())
    {
        butcherStep(dt);

        t_ += dt;
        ++nSteps_;
        lastDeltaT_ = dt;
        return;
    }

    for (label k = 0; k < rk_.nStages(); ++k)
    {
        const scalar a = rk_.a(k);
        const scalar b = rk_.b(k);

        evaluateStage(t_ + rk_.c(k)*dt);

        for (vessel1D& v : vessels_)
        {
            v.stageUpdate(a, b, dt);

            if (v.minArea() <= 0)
            {
                FatalErrorInFunction
                    << "Non-positive area in vessel " << v.name()
                    << " at t = " << t_ << " (step " << dt << ")"
                    << exit(FatalError);
            }
        }
        for (node1D& n : nodes_)
        {
            n.stageUpdate(a, b, dt);
        }
        inflow_ = a*inflow0_ + b*(inflow_ + dt*dInflowdt_);
    }

    t_ += dt;
    ++nSteps_;
    lastDeltaT_ = dt;
}


void Foam::network1D::updateBoundaries()
{
    reconstructVessels();

    for (node1D& n : nodes_)
    {
        n.setBoundary(t_);
    }
}


Foam::label Foam::network1D::nodeLabel(const word& name) const
{
    forAll(nodes_, n)
    {
        if (nodes_[n].name() == name)
        {
            return n;
        }
    }

    FatalErrorInFunction
        << "No node '" << name << "'" << exit(FatalError);

    return -1;
}


const Foam::coupledFlow& Foam::network1D::port(const word& name) const
{
    const coupledFlow* p = isA<coupledFlow>(nodes_[nodeLabel(name)]);

    if (!p)
    {
        FatalErrorInFunction
            << "Node '" << name << "' is not a coupledFlow port"
            << exit(FatalError);
    }

    return *p;
}


Foam::coupledFlow& Foam::network1D::port(const word& name)
{
    return const_cast<coupledFlow&>
    (
        static_cast<const network1D&>(*this).port(name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::network1D::advance(const scalar t0, const scalar t1)
{
    if (mag(t0 - t_) > 1e-12*Foam::max(scalar(1), mag(t_)))
    {
        FatalErrorInFunction
            << "advance from t0 = " << t0 << " but the network is at t = "
            << t_ << exit(FatalError);
    }

    t_ = t0;

    while (t_ < t1)
    {
        scalar dt = stableDeltaT();
        const scalar remaining = t1 - t_;
        bool last = false;

        if (remaining <= dt*(1 + 1e-12))
        {
            dt = remaining;
            last = true;
        }
        else if (remaining < 2*dt)
        {
            // Two equal steps rather than a full one and a sliver
            dt = 0.5*remaining;
        }

        step(dt);

        if (last)
        {
            t_ = t1;
        }
    }

    updateBoundaries();
}


void Foam::network1D::saveState()
{
    saved_.t = t_;
    saved_.A.setSize(vessels_.size());
    saved_.Q.setSize(vessels_.size());
    forAll(vessels_, v)
    {
        saved_.A[v] = vessels_[v].A();
        saved_.Q[v] = vessels_[v].Q();
    }
    saved_.y.setSize(nodes_.size());
    forAll(nodes_, n)
    {
        saved_.y[n] = nodes_[n].states();
    }
    saved_.inflow = inflow_;
    saved_.nSteps = nSteps_;
    saved_.valid = true;
}


void Foam::network1D::restoreState()
{
    if (!saved_.valid)
    {
        FatalErrorInFunction
            << "restoreState() without a saveState()" << exit(FatalError);
    }

    t_ = saved_.t;
    forAll(vessels_, v)
    {
        vessels_[v].setFields(saved_.A[v], saved_.Q[v]);
    }
    forAll(nodes_, n)
    {
        nodes_[n].setStates(saved_.y[n]);
    }
    inflow_ = saved_.inflow;
    nSteps_ = saved_.nSteps;

    updateBoundaries();
}


Foam::scalar Foam::network1D::interfacePressure(const word& node) const
{
    const node1D::vesselEnd& ve = nodes_[nodeLabel(node)].ends()[0];
    const vessel1D& v = vessels_[ve.vessel];

    return v.endLaw(ve.end).pressure(v.boundaryArea(ve.end));
}


Foam::wordList Foam::network1D::portNames() const
{
    wordList names;
    for (const node1D& n : nodes_)
    {
        if (isA<coupledFlow>(n))
        {
            names.append(n.name());
        }
    }
    return names;
}


const Foam::word& Foam::network1D::portPatch(const word& name) const
{
    return port(name).patch();
}


bool Foam::network1D::portInlet(const word& name) const
{
    return port(name).portInlet();
}


Foam::scalar Foam::network1D::portFlowScale(const word& name) const
{
    return port(name).flowScale();
}


void Foam::network1D::setPortForcing
(
    const word& name,
    const scalar t0,
    const scalar t1,
    const circulationForcing::coeffs& c
)
{
    port(name).setForcing(t0, t1, c);
}


Foam::scalar Foam::network1D::portPressure(const word& name) const
{
    return port(name).pressure();
}


Foam::scalar Foam::network1D::portFlow(const word& name) const
{
    return port(name).flow(t_);
}


Foam::label Foam::network1D::nEqns() const
{
    label n = 1;
    for (const vessel1D& v : vessels_)
    {
        n += 2*v.nCells();
    }
    for (const node1D& nd : nodes_)
    {
        n += nd.nStates();
    }
    return n;
}


Foam::scalarField Foam::network1D::state() const
{
    scalarField y(nEqns());
    label k = 0;

    for (const vessel1D& v : vessels_)
    {
        for (const scalar a : v.A())
        {
            y[k++] = a;
        }
        for (const scalar q : v.Q())
        {
            y[k++] = q;
        }
    }
    for (const node1D& n : nodes_)
    {
        for (const scalar s : n.states())
        {
            y[k++] = s;
        }
    }
    y[k] = inflow_;

    return y;
}


void Foam::network1D::setState(const scalar t, const scalarField& y)
{
    if (y.size() != nEqns())
    {
        FatalErrorInFunction
            << "State of size " << y.size() << " for a network of "
            << nEqns() << exit(FatalError);
    }

    label k = 0;

    for (vessel1D& v : vessels_)
    {
        const label n = v.nCells();
        v.setFields
        (
            scalarField(SubField<scalar>(y, n, k)),
            scalarField(SubField<scalar>(y, n, k + n))
        );
        k += 2*n;
    }
    for (node1D& nd : nodes_)
    {
        const label n = nd.nStates();
        if (n)
        {
            nd.setStates(scalarField(SubField<scalar>(y, n, k)));
            k += n;
        }
    }
    inflow_ = y[k];

    t_ = t;

    updateBoundaries();
}


Foam::scalarField Foam::network1D::flushDenormals(const scalarField& f)
{
    scalarField g(f);
    for (scalar& v : g)
    {
        if (mag(v) < VSMALL)
        {
            v = 0;
        }
    }
    return g;
}


Foam::dictionary Foam::network1D::stateDict() const
{
    dictionary d;
    d.add("time", t_);
    d.add("initialVolume", V0_);
    d.add("netInflow", inflow_);
    d.add("nSteps", nSteps_);

    dictionary vd;
    for (const vessel1D& v : vessels_)
    {
        dictionary s;
        s.add("A", flushDenormals(v.A()));
        s.add("Q", flushDenormals(v.Q()));
        vd.add(v.name(), s);
    }
    d.add("vessels", vd);

    dictionary nd;
    for (const node1D& n : nodes_)
    {
        if (n.nStates())
        {
            dictionary s;
            s.add("y", flushDenormals(n.states()));
            nd.add(n.name(), s);
        }
    }
    d.add("nodes", nd);

    return d;
}


void Foam::network1D::readState(const dictionary& d)
{
    t_ = d.get<scalar>("time");
    V0_ = d.get<scalar>("initialVolume");
    inflow_ = d.get<scalar>("netInflow");
    nSteps_ = d.getOrDefault<label>("nSteps", 0);

    const dictionary& vd = d.subDict("vessels");
    for (vessel1D& v : vessels_)
    {
        const dictionary& s = vd.subDict(v.name());
        v.setFields(s.get<scalarField>("A"), s.get<scalarField>("Q"));
    }

    const dictionary& nd = d.subDict("nodes");
    for (node1D& n : nodes_)
    {
        if (n.nStates())
        {
            n.setStates(nd.subDict(n.name()).get<scalarField>("y"));
        }
    }

    updateBoundaries();
}


void Foam::network1D::probeValues
(
    const label i,
    scalar& P,
    scalar& Q,
    scalar& A,
    scalar& U,
    scalar& dr
) const
{
    const vessel1D& v = vessels_[probeVessel_[i]];
    scalar Ad;
    v.sample(probeX_[i], P, Q, A, Ad);

    U = Q/A;

    const scalar rpi = 1/constant::mathematical::pi;
    dr = Foam::sqrt(A*rpi) - Foam::sqrt(Ad*rpi);
}


void Foam::network1D::vesselEndValues
(
    const label vi,
    scalar& Pin,
    scalar& Pout,
    scalar& Qin,
    scalar& Qout
) const
{
    const vessel1D& v = vessels_[vi];

    Pin = v.endLaw(vessel1D::inEnd).pressure(v.boundaryArea(vessel1D::inEnd));
    Pout =
        v.endLaw(vessel1D::outEnd).pressure(v.boundaryArea(vessel1D::outEnd));
    Qin = v.boundaryFlow(vessel1D::inEnd);
    Qout = v.boundaryFlow(vessel1D::outEnd);
}


Foam::scalar Foam::network1D::volume() const
{
    scalar V = 0;
    for (const vessel1D& v : vessels_)
    {
        V += v.volume();
    }
    return V;
}


void Foam::network1D::writeProfile(const label vi, Ostream& os) const
{
    const vessel1D& v = vessels_[vi];

    os  << "# vessel " << v.name() << " at t = " << t_ << nl
        << "# x A Q P U Ad" << nl
        << "# m m2 m3/s Pa m/s m2" << nl;

    auto row =
        [&](const scalar x, const scalar A, const scalar Q, const label i)
    {
        const tubeLaw1D& law = v.tubeLaw(i);
        os  << x << ' ' << A << ' ' << Q << ' ' << law.pressure(A) << ' '
            << Q/A << ' ' << law.referenceArea() << nl;
    };

    row
    (
        0,
        v.boundaryArea(vessel1D::inEnd),
        v.boundaryFlow(vessel1D::inEnd),
        0
    );

    forAll(v.A(), i)
    {
        row((i + 0.5)*v.dx(), v.A()[i], v.Q()[i], i);
    }

    row
    (
        v.length(),
        v.boundaryArea(vessel1D::outEnd),
        v.boundaryFlow(vessel1D::outEnd),
        v.nCells() - 1
    );
}


void Foam::network1D::summary(Ostream& os) const
{
    os  << "network1D: " << vessels_.size() << " vessel(s), " << nodes_.size()
        << " node(s), " << probeNames_.size() << " probe(s)" << nl
        << "  fluid: rho " << rho_ << " kg/m3, mu " << mu_
        << " Pa s, profileOrder " << profileOrder_ << ", alpha " << alpha_
        << nl
        << "  numerics: CFL " << CFL_ << ", flux " << flux_->type()
        << ", ";
    if (kExact_.valid())
    {
        os  << "reconstruction kExact, polynomialOrder "
            << kExact_->polynomialOrder() << " (order "
            << kExact_->order() << ")";
    }
    else
    {
        os  << "MUSCL, limiter " << limiter_->type();
    }
    os  << ", timeScheme " << rk_.name() << nl;

    for (const vessel1D& v : vessels_)
    {
        os  << "  ";
        v.summary(os);
    }

    for (const node1D& n : nodes_)
    {
        os  << "  ";
        n.summary(os);
    }

    forAll(probeNames_, i)
    {
        os  << "  Probe " << probeNames_[i] << ": "
            << vessels_[probeVessel_[i]].name() << " at x = " << probeX_[i]
            << " m" << nl;
    }

    os  << "  Initial volume " << V0_ << " m3" << nl;
}

// ************************************************************************* //
