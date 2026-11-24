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
    Foam::circulation1DPressureFvPatchScalarField

Description
    Kinematic pressure p/rho on a patch coupled to a port (coupledFlow node) of
    a 1D arterial network, or the pressure in Pa when the field has pressure
    dimensions (a compressible solver; Dirichlet or Aitken interface only).
    Whether the exchange is explicit or implicit, and of
    what order, is chosen in the "coupling" block of
    constant/circulation1DProperties, not by the BC type.

        inlet
        {
            type    circulation1DPressure;
            value   uniform 0;
            // network circulation1DProperties;   // optional: which description
        }

    The patch is bound to its port by the port's "patch" entry. Each time the
    boundary is evaluated the coupler is asked to sample the flux and, on a
    new exchange, advance the network.

    With "interface dirichlet" (the coupler's default) the patch gets the
    uniform pressure the network's vessel end has, divided by rho: a fixed
    value, as before.

    With "interface robin" the patch gets the impedance (Robin) condition

        p = P* / rho + zu (u_n - u_n*)

    per face: P* the vessel end's pressure for the flow the network was last
    integrated with (the iterate, or the extrapolated flow of explicit
    coupling), u_n the face's outward normal velocity as the pressure
    equation is about to set it, u_n* the same velocity of the flow the
    network saw, and zu = Zc flowScale A3D / rho the vessel end's
    characteristic impedance Zc = rho c/A per unit of mean velocity of the
    patch. Within one step a vessel end answers a change of the port flow
    with exactly this slope (its reflections take a round trip to return),
    so the condition gives the pressure equation the 1D side's response to
    the flux it is computing, instead of a value lagged by one exchange: the
    outer iteration of implicit coupling converges in one or two passes
    where the fixed point needed dt/(a0 tau) per pass, and explicit coupling
    loses its dt/tau stability limit for linearised vessels (the lagged part
    is only the incoming wave). At the fixed point u_n = u_n* on every face
    and the patch pressure is P*, uniform, as with the Dirichlet interface:
    the converged solution is the same.

    The condition is a mixed one, p = f R + (1 - f) p_P, with f = 1/(1 + k),
    k = zu rAU_f deltaCoeffs the implicit flux-pressure coefficient of the
    face (the pressure equation sets the face flux to phiHbyA - rAU_f |Sf|
    snGrad p), and R = P* / rho + zu (phiHbyA - phi*) / |Sf| with phi* the
    face flux the network was integrated with: the flux sampled at the
    exchange for implicit coupling, the sampled flux shifted to the
    extrapolated flow for explicit coupling. rAU_f is recovered from the
    solver's call to updateSnGrad (constrainPressure passes
    (phiHbyA - Sf & U)/(|Sf| rAU_f)) and phiHbyA, so the relation is the
    pressure equation's own; until the first such call the time step is
    used. Referring the velocity difference to the previous solve's
    boundary values instead (p_b* - p_P*) would make the fixed point
    independent of k but relaxes it by k/(1 + k) per pass, which is what
    the direct form avoids.

    A copy of circulationPressure (circulationCoupling) for circulation1D,
    until the two share one coupler.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#ifndef circulation1DPressureFvPatchScalarField_H
#define circulation1DPressureFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "updateableSnGrad.H"

namespace Foam
{

class circulation1DPressureFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public updateablePatchTypes::updateableSnGrad
{
    word propertiesName_;

    //- (phiHbyA - Sf & U)/(|Sf| rAU_f) from the solver, when it has
    //  called updateSnGrad since the last update
    scalarField snGradCoeff_;
    bool snGradSet_;

    //- rAU_f of each face as last recovered (deltaT until then)
    scalarField rAUf_;

public:

    TypeName("circulation1DPressure");

    circulation1DPressureFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    circulation1DPressureFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    circulation1DPressureFvPatchScalarField
    (
        const circulation1DPressureFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    circulation1DPressureFvPatchScalarField
    (
        const circulation1DPressureFvPatchScalarField&
    );

    circulation1DPressureFvPatchScalarField
    (
        const circulation1DPressureFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
            new circulation1DPressureFvPatchScalarField(*this)
        );
    }

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new circulation1DPressureFvPatchScalarField(*this, iF)
        );
    }

    //- From constrainPressure: (phiHbyA - Sf & U)/(|Sf| rAU_f)
    virtual void updateSnGrad(const scalarField& snGradp);

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
