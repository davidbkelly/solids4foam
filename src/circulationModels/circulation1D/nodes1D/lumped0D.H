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
    Foam::lumped0D

Description
    A vessel end coupled to a lumped-parameter (0D) circuit of the
    circulation0D library, as a 3D domain is coupled to one: the vessel
    supplies the flow through an externalPort of the circuit and the circuit
    answers with the pressure at the port.

    The "circuit" sub-dictionary is a circulationProperties of its own
    (components, connections, units, modules, presets), with one externalPort
    for the vessel end. A three-element Windkessel (Boileau et al. 2015,
    Eq. 9) with the windkessel3 preset:

        outlet
        {
            type        lumped0D;
            connection  carotid.out;
            circuit
            {
                units       SI;
                components
                {
                    port    { type externalPort; direction outlet; }
                    wk      { type windkessel3;
                              Rch 2.4875e8; C 1.7529e-10; R 1.8697e9; }
                }
                connections ( (port wk) );
            }
        }

    The port's flow is positive in its "direction": outlet is flow leaving
    the vessel into the circuit, inlet flow from the circuit into the vessel.
    "port" names the port if the circuit has more than one.

    Every stage, with the circuit's states y of that stage, A at the vessel
    end solves

        P(A) - Pport(Q(A); t, y) = 0

    by Newton, Q(A) the port flow on the outgoing characteristic and
    dPport/dQ by a finite difference (exact for linear circuits). The
    circuit's states are then advanced in the same Runge-Kutta stage as the
    vessels, from circulationNetwork::derivatives(); its own ODE solver is not
    used, so the coupled system keeps the order of the scheme.

    Its output columns are those of circulationFoam for the circuit.

    circulationNetwork needs at least one state (a compliance, chamber or
    inductor), so a purely resistive terminal is the resistance node instead.

Author
    David Kelly, UCD.

SourceFiles
    lumped0D.C

\*---------------------------------------------------------------------------*/

#ifndef lumped0D_H
#define lumped0D_H

#include "node1D.H"
#include "circulationNetwork.H"

namespace Foam
{

class lumped0D
:
    public node1D
{
    // Private Data

        autoPtr<circulationNetwork> net_;

        word port_;

        //- True if the port is an inlet (flow positive into the vessel)
        bool portInlet_;


    // Private Member Functions

        //- Port pressure for port flow Qp at time t with the current states
        //  (the circuit's state set beforehand)
        scalar portPressure(const scalar t, const scalar Qp);

        //- Port flow for the flow Qn from the vessel into the node
        scalar portFlow(const scalar Qn) const
        {
            return portInlet_ ? -Qn : Qn;
        }


public:

    TypeName("lumped0D");

    lumped0D
    (
        const word& name,
        const dictionary& dict,
        PtrList<vessel1D>& vessels,
        const scalar t0
    );


    // Member Functions

        virtual void setBoundary(const scalar t);

        virtual void derivatives(const scalar t);

        virtual wordList columns() const;

        virtual void outputValues(const scalar t, scalarField& values);

        virtual void summary(Ostream& os) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
