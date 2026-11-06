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

#include "rungeKutta1D.H"
#include "error.H"

Foam::rungeKutta1D::rungeKutta1D(const word& name)
:
    name_(name),
    order_(0),
    butcher_(false)
{
    if (name == "SSPRK2")
    {
        a_ = scalarList({0, 0.5});
        b_ = scalarList({1, 0.5});
        c_ = scalarList({0, 1});
        order_ = 2;
    }
    else if (name == "SSPRK3")
    {
        a_ = scalarList({0, 0.75, 1.0/3.0});
        b_ = scalarList({1, 0.25, 2.0/3.0});
        c_ = scalarList({0, 1, 0.5});
        order_ = 3;
    }
    else if (name == "RK4")
    {
        butcher_ = true;
        c_ = scalarList({0, 0.5, 0.5, 1});
        A_.setSize(4);
        A_[0] = scalarList();
        A_[1] = scalarList({0.5});
        A_[2] = scalarList({0, 0.5});
        A_[3] = scalarList({0, 0, 1});
        B_ = scalarList({1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0});
        order_ = 4;
    }
    else
    {
        FatalErrorInFunction
            << "Unknown timeScheme '" << name << "'. Use SSPRK2, SSPRK3 or"
            << " RK4" << exit(FatalError);
    }
}

// ************************************************************************* //
