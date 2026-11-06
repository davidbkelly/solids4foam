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
    Foam::rungeKutta1D

Description
    Explicit Runge-Kutta schemes, selected by the "timeScheme" entry of
    "numerics":

        SSPRK2      (default) Heun, second order
        SSPRK3      Shu-Osher, third order
        RK4         the classical fourth-order scheme (not SSP)

    The SSP schemes are in Shu-Osher form: stage k (k = 1..n) of a step from
    u^n with u^(0) = u^n is

        u^(k) = a_k u^n + b_k (u^(k-1) + dt L(t^n + c_k dt, u^(k-1)))

    and u^(n+1) = u^(n). RK4 is in Butcher form (butcher() is true): stage k
    evaluates L_k = L(t^n + c_k dt, u^n + dt sum_{j<k} A_kj L_j), and
    u^(n+1) = u^n + dt sum_k B_k L_k, which needs the stage derivatives kept.
    Every state of the system (vessel fields and the ODE states of the nodes)
    is advanced with the same stages.

Author
    David Kelly, UCD.

SourceFiles
    rungeKutta1D.C

\*---------------------------------------------------------------------------*/

#ifndef rungeKutta1D_H
#define rungeKutta1D_H

#include "scalarList.H"
#include "word.H"
#include "List.H"

namespace Foam
{

class rungeKutta1D
{
    // Private Data

        word name_;

        scalarList a_;
        scalarList b_;
        scalarList c_;

        label order_;

        //- Butcher form: A (lower triangular, by stage) and B
        bool butcher_;
        List<scalarList> A_;
        scalarList B_;


public:

    //- Construct by name
    explicit rungeKutta1D(const word& name);


    // Member Functions

        const word& name() const { return name_; }

        label nStages() const { return c_.size(); }

        label order() const { return order_; }

        //- Weight of u^n in stage k
        scalar a(const label k) const { return a_[k]; }

        //- Weight of the previous stage and its update in stage k
        scalar b(const label k) const { return b_[k]; }

        //- Time of the right-hand side evaluated in stage k, as a fraction
        //  of the step
        scalar c(const label k) const { return c_[k]; }

        //- True for a scheme in Butcher form
        bool butcher() const { return butcher_; }

        //- Butcher coefficient of stage derivative j in stage k (j < k)
        scalar A(const label k, const label j) const { return A_[k][j]; }

        //- Butcher weight of stage derivative k in the update
        scalar B(const label k) const { return B_[k]; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
