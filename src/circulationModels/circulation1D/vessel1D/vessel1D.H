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
    Foam::vessel1D

Description
    One vessel of a 1D network: a uniform array of finite-volume cells on
    0 <= x <= L holding the cell averages of A and Q, solved for

        dA/dt + dQ/dx = 0
        dQ/dt + d/dx(alpha Q^2/A) + (A/rho) dP/dx = -KR Q/A

    with KR = 2 pi (zeta + 2) mu/rho.

    Uniform vessel: the momentum equation is in conservation form, with the
    flux alpha Q^2/A + fluxPressure(A). The MUSCL reconstruction is of (A, Q)
    and interior faces take a Riemann flux (HLL or Rusanov).

    Vessel with properties varying along it (a tapered vessel): each cell has
    the wall at its centre, uniform within the cell, so the cells are short
    uniform tubes. A face between cells of different walls is a two-vessel
    junction (mass conserved, total pressure P + rho U^2/2 continuous, on the
    outgoing characteristics of both sides): the mass flux is common and the
    momentum flux of each side is its own. The MUSCL reconstruction is of
    (P, Q), the variables that are uniform at rest, so a vessel at rest
    (uniform pressure, no flow) stays at rest exactly, whatever the taper: the
    scheme is well-balanced. No geometric source term is needed.

    The cells at the ends are reconstructed one-sided. The two end faces take
    the physical flux of the boundary state set by the node at that end
    (setBoundaryState()).

    Uniform vessels may instead use the k-exact reconstruction of degree k
    (kExactReconstruction1D), for face values of order k + 1; the friction
    source is then the average over the cell by two-point Gauss quadrature of
    the reconstruction, so that it keeps that order.

    The ends are "in" (x = 0) and "out" (x = L); flow is positive from in to
    out, so a flow Q at an end enters the node there as endSign*Q.

    Entries of a vessel dictionary, besides those of its tube law:

        length              0.126;
        nCells              126;        // or dx 1e-3;
        tubeLaw             sqrtElastic;
        initialPressure     0;          // uniform; default 0
        initialFlowRate     0;          // uniform; default 0

Author
    David Kelly, UCD.

SourceFiles
    vessel1D.C

\*---------------------------------------------------------------------------*/

#ifndef vessel1D_H
#define vessel1D_H

#include "tubeLaw1D.H"
#include "slopeLimiter1D.H"
#include "riemannFlux1D.H"
#include "kExactReconstruction1D.H"
#include "rungeKutta1D.H"
#include "scalarField.H"
#include "PtrList.H"
#include "FixedList.H"

namespace Foam
{

class vessel1D
{
public:

    //- The two ends of a vessel
    enum endType
    {
        inEnd = 0,      //!< x = 0
        outEnd = 1      //!< x = L
    };

    //- +1 at the out end, -1 at the in end: the flow from the vessel into
    //  the node at that end is endSign*Q
    static scalar endSign(const label end)
    {
        return end == outEnd ? 1 : -1;
    }

    //- "in" or "out"
    static word endName(const label end)
    {
        return end == outEnd ? "out" : "in";
    }


private:

    // Private Data

        const word name_;

        const scalar length_;

        label nCells_;

        scalar dx_;

        //- Blood density [kg/m3]
        const scalar rho_;

        //- Momentum-flux coefficient
        const scalar alpha_;

        //- Friction coefficient KR [m2/s]
        const scalar KR_;

        //- The wall: one tube law if the vessel is uniform, else one per
        //  cell, at the cell centre
        PtrList<tubeLaw1D> laws_;

        bool uniform_;

        //- Cell averages, current stage
        scalarField A_;
        scalarField Q_;

        //- Cell averages at the start of the step
        scalarField A0_;
        scalarField Q0_;

        //- Right-hand side of the current stage
        scalarField dAdt_;
        scalarField dQdt_;

        //- Reconstructed states at the left (x - dx/2) and right (x + dx/2)
        //  face of each cell
        scalarField AfL_;
        scalarField QfL_;
        scalarField AfR_;
        scalarField QfR_;

        //- Work space for the slopes
        scalarField s1_;
        scalarField s2_;

        //- k-exact reconstruction: A and Q at the two Gauss points of each
        //  cell, for the friction source; used if gaussSource_
        scalarField Ag1_;
        scalarField Ag2_;
        scalarField Qg1_;
        scalarField Qg2_;
        bool gaussSource_;

        //- Butcher-form Runge-Kutta: the stage derivatives
        PtrList<scalarField> kA_;
        PtrList<scalarField> kQ_;

        //- Boundary states at the in and out ends
        FixedList<scalar, 2> Ab_;
        FixedList<scalar, 2> Qb_;


    // Private Member Functions

        //- Flux across a face between cells of different walls: a two-vessel
        //  junction of the left state (out end of the left cell) and the right
        //  state (in end of the right cell). Returns the common mass flux FA
        //  and the momentum flux of each side.
        void coupledFlux
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
        ) const;


public:

    // Constructors

        vessel1D
        (
            const word& name,
            const dictionary& dict,
            const scalar rho,
            const scalar alpha,
            const scalar KR
        );

        vessel1D(const vessel1D&) = delete;
        void operator=(const vessel1D&) = delete;


    // Access

        const word& name() const { return name_; }

        scalar length() const { return length_; }

        label nCells() const { return nCells_; }

        scalar dx() const { return dx_; }

        scalar rho() const { return rho_; }

        scalar alpha() const { return alpha_; }

        scalar frictionCoeff() const { return KR_; }

        //- True if the wall is the same all along the vessel
        bool uniform() const { return uniform_; }

        //- The tube law of a cell
        const tubeLaw1D& tubeLaw(const label cell) const
        {
            return laws_[uniform_ ? 0 : cell];
        }

        //- The tube law of the cell at an end
        const tubeLaw1D& endLaw(const label end) const
        {
            return tubeLaw(end == inEnd ? 0 : nCells_ - 1);
        }

        const scalarField& A() const { return A_; }

        const scalarField& Q() const { return Q_; }

        //- Replace the cell averages, e.g. when restoring a state
        void setFields(const scalarField& A, const scalarField& Q);

        //- Boundary state at an end
        scalar boundaryArea(const label end) const { return Ab_[end]; }
        scalar boundaryFlow(const label end) const { return Qb_[end]; }


    // Physics

        //- Physical flux of (A, Q) for a tube law (the first by default)
        void physicalFlux
        (
            const scalar A,
            const scalar Q,
            scalar& FA,
            scalar& FQ,
            const label cell = 0
        ) const;

        //- Eigenvalues lambda- <= lambda+ (for the tube law of a cell)
        void eigenvalues
        (
            const scalar A,
            const scalar Q,
            scalar& lm,
            scalar& lp,
            const label cell = 0
        ) const;

        //- Largest |lambda| over the cells and boundary states [m/s]
        scalar maxWaveSpeed() const;

        //- Smallest cell area [m2]
        scalar minArea() const;

        //- Volume, the integral of A [m3]
        scalar volume() const;


    // Time stepping

        //- Limited slopes and face states of the current stage
        void reconstruct(const slopeLimiter1D& limiter);

        //- Face states (and Gauss-point states) of the current stage by the
        //  k-exact reconstruction; uniform vessels only
        void reconstruct(const kExactReconstruction1D& kExact);

        //- Reconstructed state of the cell next to an end, at that end
        void faceState(const label end, scalar& A, scalar& Q) const;

        //- Set the boundary state at an end
        void setBoundaryState(const label end, const scalar A, const scalar Q);

        //- Right-hand side of the current stage (after reconstruct() and the
        //  boundary states)
        void evaluateRhs(const riemannFlux1D& flux);

        //- Keep the fields at the start of the step
        void storeOld();

        //- Stage update, u = a u^n + b (u + dt L)
        void stageUpdate(const scalar a, const scalar b, const scalar dt);

        //- Butcher form: keep the right-hand side as stage derivative k
        void storeStageRhs(const label k, const label nStages);

        //- Butcher form: u = u^n + dt sum_{j<k} A_kj L_j, the state of stage
        //  k; k = nStages gives the end of the step, with the weights B
        void butcherStage
        (
            const rungeKutta1D& rk,
            const label k,
            const scalar dt
        );


    // Sampling

        //- State at x (0 <= x <= L), linear between the cell centres and the
        //  boundary states at the ends: pressure, flow, area and the
        //  reference area there
        void sample
        (
            const scalar x,
            scalar& P,
            scalar& Q,
            scalar& A,
            scalar& Ad
        ) const;


    // Info

        void summary(Ostream& os) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
