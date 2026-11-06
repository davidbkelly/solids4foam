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

#include "kExactReconstruction1D.H"
#include "error.H"

#include <cmath>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

using Foam::scalar;
using Foam::label;

//- Average of xi^m over the cell at offset d (xi in cell widths)
scalar moment(const label m, const label d)
{
    return
    (
        std::pow(d + 0.5, m + 1) - std::pow(d - 0.5, m + 1)
    )/(m + 1);
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kExactReconstruction1D::kExactReconstruction1D(const label k)
:
    k_(k),
    W_(0),
    w_()
{
    // Degree 4 (fifth order) was tried: with its one-sided stencils at the
    // vessel ends it went unstable on fine meshes (RK4, CFL 0.5)
    if (k < 1 || k > 3)
    {
        FatalErrorInFunction
            << "polynomialOrder must be 1, 2 or 3 (order 2 to 4), not " << k
            << exit(FatalError);
    }

    W_ = 2*((k + 1)/2) + 1;
    w_.setSize(W_);

    const scalar g = 0.5/std::sqrt(3.0);
    const scalar xi[4] = {-0.5, 0.5, -g, g};

    // For each position p of the cell in its stencil
    for (label p = 0; p < W_; ++p)
    {
        // Rows: the other cells of the stencil, offsets d = q - p
        List<label> d;
        for (label q = 0; q < W_; ++q)
        {
            if (q != p)
            {
                d.append(q - p);
            }
        }
        const label nr = d.size();

        // M (nr x k): M[r][m-1] = <xi^m>_d - <xi^m>_0
        List<scalarList> M(nr, scalarList(k, 0));
        forAll(d, r)
        {
            for (label m = 1; m <= k; ++m)
            {
                M[r][m-1] = moment(m, d[r]) - moment(m, 0);
            }
        }

        // Normal matrix N = M^T M, and X = N^-1 M^T (k x nr) by Gauss-Jordan
        // elimination with partial pivoting on [N | M^T]
        List<scalarList> aug(k, scalarList(k + nr, 0));
        for (label a = 0; a < k; ++a)
        {
            for (label b = 0; b < k; ++b)
            {
                scalar s = 0;
                for (label r = 0; r < nr; ++r)
                {
                    s += M[r][a]*M[r][b];
                }
                aug[a][b] = s;
            }
            for (label r = 0; r < nr; ++r)
            {
                aug[a][k + r] = M[r][a];
            }
        }

        for (label col = 0; col < k; ++col)
        {
            label piv = col;
            for (label row = col + 1; row < k; ++row)
            {
                if (std::abs(aug[row][col]) > std::abs(aug[piv][col]))
                {
                    piv = row;
                }
            }
            std::swap(aug[col], aug[piv]);

            const scalar diag = aug[col][col];
            if (std::abs(diag) < 1e-14)
            {
                FatalErrorInFunction
                    << "Singular least-squares system for polynomialOrder "
                    << k << exit(FatalError);
            }
            for (label c = 0; c < k + nr; ++c)
            {
                aug[col][c] /= diag;
            }
            for (label row = 0; row < k; ++row)
            {
                if (row != col)
                {
                    const scalar f = aug[row][col];
                    for (label c = 0; c < k + nr; ++c)
                    {
                        aug[row][c] -= f*aug[col][c];
                    }
                }
            }
        }

        // Value at xi_o: u_p + e^T X (u_d - u_p), e_m = xi_o^m - <xi^m>_0
        for (label o = 0; o < 4; ++o)
        {
            scalarList& wo = w_[p][o];
            wo.setSize(W_, 0);

            scalar sumG = 0;
            forAll(d, r)
            {
                scalar gr = 0;
                for (label m = 1; m <= k; ++m)
                {
                    gr += (std::pow(xi[o], m) - moment(m, 0))*aug[m-1][k + r];
                }
                wo[p + d[r]] = gr;
                sumG += gr;
            }
            wo[p] = 1 - sumG;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::kExactReconstruction1D::weights&
Foam::kExactReconstruction1D::stencil
(
    const label n,
    const label i,
    label& start
) const
{
    const label h = W_/2;
    start = Foam::min(Foam::max(i - h, label(0)), n - W_);
    return w_[i - start];
}

// ************************************************************************* //
