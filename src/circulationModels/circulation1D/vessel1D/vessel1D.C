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

#include "vessel1D.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vessel1D::vessel1D
(
    const word& name,
    const dictionary& dict,
    const scalar rho,
    const scalar alpha,
    const scalar KR
)
:
    name_(name),
    length_(dict.get<scalar>("length")),
    nCells_(0),
    dx_(0),
    rho_(rho),
    alpha_(alpha),
    KR_(KR),
    laws_(),
    uniform_(true),
    gaussSource_(false)
{
    if (length_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "Vessel " << name_ << ": length must be positive"
            << exit(FatalIOError);
    }

    if (dict.found("nCells"))
    {
        nCells_ = dict.get<label>("nCells");
    }
    else
    {
        nCells_ = label(std::round(length_/dict.get<scalar>("dx")));
    }

    if (nCells_ < 2)
    {
        FatalIOErrorInFunction(dict)
            << "Vessel " << name_ << ": at least 2 cells are needed"
            << exit(FatalIOError);
    }

    dx_ = length_/nCells_;

    // The wall: one tube law, or one per cell centre if it varies
    {
        autoPtr<tubeLaw1D> law(tubeLaw1D::New(name, dict, rho));
        uniform_ = law->uniform();

        if (uniform_)
        {
            laws_.setSize(1);
            laws_.set(0, law.ptr());
        }
        else
        {
            laws_.setSize(nCells_);
            forAll(laws_, i)
            {
                laws_.set
                (
                    i,
                    tubeLaw1D::New(name, dict, rho, (i + 0.5)/nCells_).ptr()
                );
            }
        }
    }

    const scalar P0 = dict.getOrDefault<scalar>("initialPressure", 0);
    const scalar Q0 = dict.getOrDefault<scalar>("initialFlowRate", 0);

    A_.setSize(nCells_);
    forAll(A_, i)
    {
        A_[i] = tubeLaw(i).area(P0);
    }
    Q_.setSize(nCells_, Q0);
    A0_ = A_;
    Q0_ = Q_;
    dAdt_.setSize(nCells_, 0);
    dQdt_.setSize(nCells_, 0);
    AfL_ = A_;
    AfR_ = A_;
    QfL_ = Q_;
    QfR_ = Q_;
    s1_.setSize(nCells_, 0);
    s2_.setSize(nCells_, 0);
    Ag1_ = A_;
    Ag2_ = A_;
    Qg1_ = Q_;
    Qg2_ = Q_;
    Ab_[inEnd] = A_[0];
    Ab_[outEnd] = A_[nCells_ - 1];
    Qb_ = Q0;
}


// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void Foam::vessel1D::coupledFlux
(
    const tubeLaw1D& lawL,
    const scalar AL,
    const scalar QL,
    const tubeLaw1D& lawR,
    const scalar AR,
    const scalar QR,
    scalar& FA,
    scalar& FQL,
    scalar& FQR
) const
{
    // Outgoing invariants: the left state is the out end of the left cell
    // (s = +1), the right state the in end of the right cell (s = -1)
    const scalar WL = QL/AL + lawL.invariant(AL);
    const scalar WR = QR/AR - lawR.invariant(AR);

    scalar a = AL, b = AR, UL = 0, UR = 0;

    const label maxIter = 50;
    label iter = 0;
    scalar lastChange = GREAT;

    for (; iter < maxIter; ++iter)
    {
        const scalar cL = lawL.waveSpeed(a), cR = lawR.waveSpeed(b);
        UL = WL - lawL.invariant(a);
        UR = WR + lawR.invariant(b);

        // Mass, R0 = a UL - b UR, and total pressure, R1 = HL - HR, with
        // their derivatives in a and b
        const scalar R0 = a*UL - b*UR;
        const scalar R1 =
            lawL.pressure(a) + 0.5*rho_*UL*UL
          - lawR.pressure(b) - 0.5*rho_*UR*UR;

        const scalar dQL = UL - cL;
        const scalar dQR = -(UR + cR);
        const scalar dHL = lawL.dPdA(a) - rho_*UL*cL/a;
        const scalar dHR = lawR.dPdA(b) + rho_*UR*cR/b;

        // dQL da + dQR db = -R0, dHL da - dHR db = -R1
        const scalar da = (-R0 - dQR*R1/dHR)/(dQL + dQR*dHL/dHR);
        const scalar db = (R1 + dHL*da)/dHR;

        // Converged when the correction is at round-off: below 1e-15, or no
        // longer shrinking quadratically once small. It is then not applied,
        // so that a state at rest stays exactly at rest. (A looser test leaves
        // errors that build up step after step for small perturbations of
        // the area.)
        const scalar change = Foam::max(mag(da)/a, mag(db)/b);

        if (change <= 1e-15 || (change < 1e-11 && change > 0.5*lastChange))
        {
            break;
        }

        a = (a + da > 0) ? a + da : 0.5*a;
        b = (b + db > 0) ? b + db : 0.5*b;
        lastChange = change;
    }

    if (iter == maxIter)
    {
        FatalErrorInFunction
            << "Vessel " << name_ << ": the coupling of two cells of"
            << " different walls did not converge" << exit(FatalError);
    }

    FA = a*UL;
    FQL = alpha_*FA*FA/a + lawL.fluxPressure(a);
    FQR = alpha_*FA*FA/b + lawR.fluxPressure(b);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vessel1D::setFields(const scalarField& A, const scalarField& Q)
{
    if (A.size() != nCells_ || Q.size() != nCells_)
    {
        FatalErrorInFunction
            << "Vessel " << name_ << " has " << nCells_ << " cells, given "
            << A.size() << " values of A and " << Q.size() << " of Q"
            << exit(FatalError);
    }
    A_ = A;
    Q_ = Q;
}


void Foam::vessel1D::physicalFlux
(
    const scalar A,
    const scalar Q,
    scalar& FA,
    scalar& FQ,
    const label cell
) const
{
    FA = Q;
    FQ = alpha_*Q*Q/A + tubeLaw(cell).fluxPressure(A);
}


void Foam::vessel1D::eigenvalues
(
    const scalar A,
    const scalar Q,
    scalar& lm,
    scalar& lp,
    const label cell
) const
{
    const scalar U = Q/A;
    const scalar c = tubeLaw(cell).waveSpeed(A);
    const scalar d = Foam::sqrt(c*c + alpha_*(alpha_ - 1)*U*U);

    lm = alpha_*U - d;
    lp = alpha_*U + d;
}


Foam::scalar Foam::vessel1D::maxWaveSpeed() const
{
    scalar s = 0;
    scalar lm, lp;

    forAll(A_, i)
    {
        eigenvalues(A_[i], Q_[i], lm, lp, i);
        s = Foam::max(s, Foam::max(mag(lm), mag(lp)));
    }

    for (label e = 0; e < 2; ++e)
    {
        eigenvalues(Ab_[e], Qb_[e], lm, lp, e == inEnd ? 0 : nCells_ - 1);
        s = Foam::max(s, Foam::max(mag(lm), mag(lp)));
    }

    return s;
}


Foam::scalar Foam::vessel1D::minArea() const
{
    return Foam::min(A_);
}


Foam::scalar Foam::vessel1D::volume() const
{
    return sum(A_)*dx_;
}


void Foam::vessel1D::reconstruct(const slopeLimiter1D& limiter)
{
    const label n = nCells_;
    gaussSource_ = false;

    // Uniform vessel: (A, Q); varying wall: (P, Q), uniform at rest
    scalarField& u = s1_;
    scalarField& su = AfL_;     // the slopes, before the face states
    scalarField& sq = s2_;

    if (uniform_)
    {
        u = A_;
    }
    else
    {
        forAll(A_, i)
        {
            u[i] = tubeLaw(i).pressure(A_[i]);
        }
    }

    // End cells, one-sided
    su[0] = u[1] - u[0];
    sq[0] = Q_[1] - Q_[0];
    su[n-1] = u[n-1] - u[n-2];
    sq[n-1] = Q_[n-1] - Q_[n-2];

    for (label i = 1; i < n - 1; ++i)
    {
        su[i] = limiter.slope(u[i] - u[i-1], u[i+1] - u[i]);
        sq[i] = limiter.slope(Q_[i] - Q_[i-1], Q_[i+1] - Q_[i]);
    }

    forAll(A_, i)
    {
        scalar s = su[i];

        if (uniform_)
        {
            // Keep the face areas positive
            if (A_[i] - 0.5*mag(s) <= 0)
            {
                s = 0;
            }
            AfL_[i] = A_[i] - 0.5*s;
            AfR_[i] = A_[i] + 0.5*s;
        }
        else
        {
            const tubeLaw1D& law = tubeLaw(i);

            // Keep the face pressures above that of zero area
            if (u[i] - 0.5*mag(s) <= law.minPressure())
            {
                s = 0;
            }
            AfL_[i] = law.area(u[i] - 0.5*s);
            AfR_[i] = law.area(u[i] + 0.5*s);
        }

        QfL_[i] = Q_[i] - 0.5*sq[i];
        QfR_[i] = Q_[i] + 0.5*sq[i];
    }
}


void Foam::vessel1D::reconstruct(const kExactReconstruction1D& kExact)
{
    if (!uniform_)
    {
        FatalErrorInFunction
            << "Vessel " << name_ << ": the kExact reconstruction is for"
            << " uniform vessels only (this one is tapered)"
            << exit(FatalError);
    }

    const label n = nCells_;
    typedef kExactReconstruction1D kE;

    label start;
    forAll(A_, i)
    {
        const kE::weights& w = kExact.stencil(n, i, start);

        scalar a[4] = {0, 0, 0, 0}, q[4] = {0, 0, 0, 0};
        forAll(w[0], j)
        {
            const scalar Aj = A_[start + j], Qj = Q_[start + j];
            for (label o = 0; o < 4; ++o)
            {
                a[o] += w[o][j]*Aj;
                q[o] += w[o][j]*Qj;
            }
        }

        // Keep the areas positive: first order in a cell where they are not
        if (a[0] <= 0 || a[1] <= 0 || a[2] <= 0 || a[3] <= 0)
        {
            for (label o = 0; o < 4; ++o)
            {
                a[o] = A_[i];
                q[o] = Q_[i];
            }
        }

        AfL_[i] = a[kE::faceL];
        AfR_[i] = a[kE::faceR];
        QfL_[i] = q[kE::faceL];
        QfR_[i] = q[kE::faceR];
        Ag1_[i] = a[kE::gauss1];
        Ag2_[i] = a[kE::gauss2];
        Qg1_[i] = q[kE::gauss1];
        Qg2_[i] = q[kE::gauss2];
    }

    gaussSource_ = true;
}


void Foam::vessel1D::faceState
(
    const label end,
    scalar& A,
    scalar& Q
) const
{
    if (end == inEnd)
    {
        A = AfL_[0];
        Q = QfL_[0];
    }
    else
    {
        A = AfR_[nCells_ - 1];
        Q = QfR_[nCells_ - 1];
    }
}


void Foam::vessel1D::setBoundaryState
(
    const label end,
    const scalar A,
    const scalar Q
)
{
    Ab_[end] = A;
    Qb_[end] = Q;
}


void Foam::vessel1D::evaluateRhs(const riemannFlux1D& flux)
{
    const label n = nCells_;

    dAdt_ = 0;
    dQdt_ = 0;

    scalar FA, FQ;

    // Interior faces: face f lies between cells f-1 and f
    if (uniform_)
    {
        for (label f = 1; f < n; ++f)
        {
            const label l = f - 1;

            flux.flux(*this, AfR_[l], QfR_[l], AfL_[f], QfL_[f], FA, FQ);

            dAdt_[l] -= FA;
            dQdt_[l] -= FQ;
            dAdt_[f] += FA;
            dQdt_[f] += FQ;
        }
    }
    else
    {
        scalar FQR;

        for (label f = 1; f < n; ++f)
        {
            const label l = f - 1;

            coupledFlux
            (
                tubeLaw(l), AfR_[l], QfR_[l],
                tubeLaw(f), AfL_[f], QfL_[f],
                FA, FQ, FQR
            );

            dAdt_[l] -= FA;
            dQdt_[l] -= FQ;
            dAdt_[f] += FA;
            dQdt_[f] += FQR;
        }
    }

    // End faces, from the boundary states
    physicalFlux(Ab_[inEnd], Qb_[inEnd], FA, FQ, 0);
    dAdt_[0] += FA;
    dQdt_[0] += FQ;

    physicalFlux(Ab_[outEnd], Qb_[outEnd], FA, FQ, n - 1);
    dAdt_[n-1] -= FA;
    dQdt_[n-1] -= FQ;

    const scalar rdx = 1/dx_;

    if (gaussSource_)
    {
        // The cell average of the friction by two-point Gauss quadrature
        forAll(A_, i)
        {
            dAdt_[i] *= rdx;
            dQdt_[i] =
                dQdt_[i]*rdx
              - 0.5*KR_*(Qg1_[i]/Ag1_[i] + Qg2_[i]/Ag2_[i]);
        }
    }
    else
    {
        forAll(A_, i)
        {
            dAdt_[i] *= rdx;
            dQdt_[i] = dQdt_[i]*rdx - KR_*Q_[i]/A_[i];
        }
    }
}


void Foam::vessel1D::storeOld()
{
    A0_ = A_;
    Q0_ = Q_;
}


void Foam::vessel1D::stageUpdate
(
    const scalar a,
    const scalar b,
    const scalar dt
)
{
    forAll(A_, i)
    {
        A_[i] = a*A0_[i] + b*(A_[i] + dt*dAdt_[i]);
        Q_[i] = a*Q0_[i] + b*(Q_[i] + dt*dQdt_[i]);
    }
}


void Foam::vessel1D::storeStageRhs(const label k, const label nStages)
{
    if (kA_.size() != nStages)
    {
        kA_.setSize(nStages);
        kQ_.setSize(nStages);
        for (label j = 0; j < nStages; ++j)
        {
            kA_.set(j, new scalarField(nCells_, 0));
            kQ_.set(j, new scalarField(nCells_, 0));
        }
    }
    kA_[k] = dAdt_;
    kQ_[k] = dQdt_;
}


void Foam::vessel1D::butcherStage
(
    const rungeKutta1D& rk,
    const label k,
    const scalar dt
)
{
    const bool last = (k == rk.nStages());

    A_ = A0_;
    Q_ = Q0_;
    for (label j = 0; j < k; ++j)
    {
        const scalar w = dt*(last ? rk.B(j) : rk.A(k, j));
        if (w != 0)
        {
            A_ += w*kA_[j];
            Q_ += w*kQ_[j];
        }
    }
}


void Foam::vessel1D::sample
(
    const scalar x,
    scalar& P,
    scalar& Q,
    scalar& A,
    scalar& Ad
) const
{
    // The two points to interpolate between: index -1 is the in-end
    // boundary state, nCells the out-end one, others cell centres
    const scalar h = 0.5*dx_;
    label i0, i1;
    scalar w;

    if (x <= h)
    {
        i0 = -1;
        i1 = 0;
        w = Foam::max(x, scalar(0))/h;
    }
    else if (x >= length_ - h)
    {
        i0 = nCells_ - 1;
        i1 = nCells_;
        w = 1 - Foam::max(length_ - x, scalar(0))/h;
    }
    else
    {
        const scalar s = x/dx_ - 0.5;
        i0 = Foam::min(label(s), nCells_ - 2);
        i1 = i0 + 1;
        w = s - i0;
    }

    auto point = [&](const label i, scalar& Pi, scalar& Qi, scalar& Ai)
    {
        const label cell = Foam::min(Foam::max(i, label(0)), nCells_ - 1);
        Ai = i < 0 ? Ab_[inEnd] : i >= nCells_ ? Ab_[outEnd] : A_[i];
        Qi = i < 0 ? Qb_[inEnd] : i >= nCells_ ? Qb_[outEnd] : Q_[i];
        Pi = tubeLaw(cell).pressure(Ai);
        return cell;
    };

    scalar P0, Q0, A0, P1, Q1, A1;
    const label c0 = point(i0, P0, Q0, A0);
    const label c1 = point(i1, P1, Q1, A1);

    A = (1 - w)*A0 + w*A1;
    Q = (1 - w)*Q0 + w*Q1;
    Ad =
        (1 - w)*tubeLaw(c0).referenceArea()
      + w*tubeLaw(c1).referenceArea();

    // A uniform vessel takes the pressure of the interpolated area; a
    // varying one interpolates the pressure, which is uniform at rest
    P = uniform_ ? tubeLaw(0).pressure(A) : (1 - w)*P0 + w*P1;
}


void Foam::vessel1D::summary(Ostream& os) const
{
    const scalar A = A_[0];
    const tubeLaw1D& law = tubeLaw(0);
    const scalar c = law.waveSpeed(A);

    os  << "Vessel " << name_ << ": L " << length_ << " m, " << nCells_
        << " cells, dx " << dx_ << " m, KR " << KR_ << " m2/s, alpha "
        << alpha_ << (uniform_ ? "" : ", tapered") << nl
        << "    ";
    law.summary(os);
    if (!uniform_)
    {
        os  << "    to ";
        tubeLaw(nCells_ - 1).summary(os);
    }
    os  << "    initial state: A " << A << " m2, P " << law.pressure(A)
        << " Pa, c " << c << " m/s, Z = rho c/A " << rho_*c/A
        << " Pa s/m3" << nl;
}

// ************************************************************************* //
