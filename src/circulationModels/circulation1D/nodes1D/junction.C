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

#include "junction.H"

namespace Foam
{
    defineTypeNameAndDebug(junction, 0);
    addToRunTimeSelectionTable(node1D, junction, dictionary);
}


Foam::junction::junction
(
    const word& name,
    const dictionary& dict,
    PtrList<vessel1D>& vessels,
    const scalar t0
)
:
    node1D(name, dict, vessels)
{
    if (ends_.size() < 2)
    {
        fail("a junction joins at least 2 vessel ends");
    }
}


void Foam::junction::setBoundary(const scalar t)
{
    const label n = ends_.size();

    scalarField W(n), A(n), s(n);
    forAll(ends_, i)
    {
        W[i] = outgoingInvariant(i, A[i]);
        s[i] = sign(i);
    }

    // Per vessel: velocity, total pressure H = P + rho U^2/2, and the
    // derivatives of H and of the flow into the junction, Qn = s A U
    scalarField U(n), H(n), dH(n), Qn(n), dQn(n);

    auto evaluate = [&]()
    {
        forAll(ends_, i)
        {
            const tubeLaw1D& law = vessel(i).endLaw(end(i));
            const scalar rho = vessel(i).rho();
            const scalar c = law.waveSpeed(A[i]);

            U[i] = W[i] - s[i]*law.invariant(A[i]);
            H[i] = law.pressure(A[i]) + 0.5*rho*sqr(U[i]);
            dH[i] = law.dPdA(A[i]) - rho*U[i]*s[i]*c/A[i];
            Qn[i] = s[i]*A[i]*U[i];
            dQn[i] = s[i]*(U[i] - s[i]*c);
        }
    };

    const label maxIter = 100;
    label iter = 0;

    for (; iter < maxIter; ++iter)
    {
        evaluate();

        // Residuals: R0 = sum Qn, Ri = H0 - Hi. The Newton step solves
        //     sum_i dQn_i dA_i = -R0
        //     dH_0 dA_0 - dH_i dA_i = -Ri         (i >= 1)
        // so dA_i = (Ri + dH_0 dA_0)/dH_i, and dA_0 from the first row.
        const scalar R0 = sum(Qn);

        scalar num = -R0;
        scalar den = dQn[0];
        for (label i = 1; i < n; ++i)
        {
            const scalar Ri = H[0] - H[i];
            num -= dQn[i]*Ri/dH[i];
            den += dQn[i]*dH[0]/dH[i];
        }

        if (den == 0 || !std::isfinite(den))
        {
            fail("singular junction Jacobian");
        }

        scalarField dA(n);
        dA[0] = num/den;
        for (label i = 1; i < n; ++i)
        {
            dA[i] = (H[0] - H[i] + dH[0]*dA[0])/dH[i];
        }

        scalar change = 0;
        forAll(A, i)
        {
            scalar An = A[i] + dA[i];
            if (An <= 0)
            {
                An = 0.5*A[i];
            }
            change = Foam::max(change, mag(An - A[i])/An);
            A[i] = An;
        }

        if (change <= 1e-12)
        {
            break;
        }
    }

    if (iter == maxIter)
    {
        fail
        (
            "Newton solve of the junction did not converge in "
          + std::to_string(maxIter) + " iterations"
        );
    }

    evaluate();

    forAll(ends_, i)
    {
        vessel(i).setBoundaryState(end(i), A[i], A[i]*U[i]);
    }
}

// ************************************************************************* //
