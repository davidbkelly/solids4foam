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
    Foam::kExactReconstruction1D

Description
    High-order reconstruction of cell averages on a uniform 1D mesh: the 1D
    form of the k-exact least-squares reconstruction of solids4foam's
    higher-order finite volumes (higherOrderHelpers/kExactLeastSquares).

    In cell i, with xi = (x - x_i)/dx, the reconstruction is the polynomial
    of degree k (polynomialOrder)

        u(xi) = u_i + sum_{m=1..k} c_m (xi^m - <xi^m>_i)

    whose average over the cell is the cell average u_i, <.>_i the cell
    average. Its coefficients fit, by least squares, the averages of the
    other cells of a stencil of W = 2 ceil(k/2) + 1 cells:

        sum_m c_m (<xi^m>_j - <xi^m>_i) = u_j - u_i.

    k = 1, 2 use 3 cells (k = 2 exactly), k = 3 uses 5. (k = 4, the 5-cell
    exact fit, was unstable with its one-sided stencils at the vessel ends
    on fine meshes, and is not offered.)
    Reconstructions of degree k reproduce polynomials of degree k exactly, so
    the face values are of order k + 1: k = 2 gives the third-order face
    values (-u_{i-1} + 5 u_i + 2 u_{i+1})/6, k = 3 fourth-order ones.

    Near the ends of a vessel the stencil is shifted to stay inside it, so
    the reconstruction there is one-sided, of the same degree. On a uniform
    mesh every value is a fixed linear combination of the stencil's averages;
    the weights are computed once for each of the W positions of a cell in
    its stencil, for four points: the left and right faces (xi = -1/2, 1/2)
    and the two Gauss points (xi = -+1/(2 sqrt 3)), with which the cell
    average of a source term is integrated to fourth order.

        numerics { reconstruction kExact; polynomialOrder 3; ... }

    The reconstruction is unlimited: for smooth solutions.

Author
    David Kelly, UCD.

SourceFiles
    kExactReconstruction1D.C

\*---------------------------------------------------------------------------*/

#ifndef kExactReconstruction1D_H
#define kExactReconstruction1D_H

#include "scalarList.H"
#include "FixedList.H"
#include "List.H"

namespace Foam
{

class kExactReconstruction1D
{
public:

    //- The points reconstructed in a cell
    enum point
    {
        faceL = 0,      //!< xi = -1/2
        faceR = 1,      //!< xi = 1/2
        gauss1 = 2,     //!< xi = -1/(2 sqrt 3)
        gauss2 = 3      //!< xi = 1/(2 sqrt 3)
    };

    typedef FixedList<scalarList, 4> weights;


private:

    // Private Data

        //- Degree of the polynomial
        label k_;

        //- Stencil width
        label W_;

        //- Weights by the position of the cell in its stencil
        List<weights> w_;


public:

    //- Construct for polynomial degree k (1 to 3)
    explicit kExactReconstruction1D(const label k);


    // Member Functions

        label polynomialOrder() const { return k_; }

        //- Order of the face values, k + 1
        label order() const { return k_ + 1; }

        //- Cells in a stencil
        label width() const { return W_; }

        //- First cell of the stencil of cell i of n, and the weights of its
        //  cells
        const weights& stencil(const label n, const label i, label& start)
            const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
