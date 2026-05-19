/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "electronWallFluxCorrectedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
electronWallFluxCorrectedFvPatchScalarField::electronWallFluxCorrectedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    electricFieldName_("E"),
    speciesName_(iF.name()),
    uthPrefactor_(0.0),
    reflectionFraction_(0.0)
{
    refValue()      = 0.0;
    refGrad()       = 0.0;
    valueFraction() = 0.0;
}


electronWallFluxCorrectedFvPatchScalarField::electronWallFluxCorrectedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField
    (
        p,
        iF,
        dict,
        IOobjectOption::NO_READ
    ),
    electricFieldName_(dict.getOrDefault<word>("electricField", "E")),
    speciesName_(iF.name()),
    uthPrefactor_(8.0*constant::physicoChemical::k.value()/(constant::mathematical::pi*constant::atomic::me.value())),
    reflectionFraction_(dict.getOrDefault<scalar>("reflectionFraction", 0.0))
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;

    // Initialisation of patch field
    if (!this->readValueEntry(dict))
    {
        // Initialise patch values equal to internal values
        fvPatchField<scalar>::operator=(this->patchInternalField());
    }
}


electronWallFluxCorrectedFvPatchScalarField::electronWallFluxCorrectedFvPatchScalarField
(
    const electronWallFluxCorrectedFvPatchScalarField& rhs,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& m
)
:
    mixedFvPatchScalarField(rhs, p, iF, m),
    electricFieldName_(rhs.electricFieldName_),
    speciesName_(rhs.speciesName_),
    uthPrefactor_(rhs.uthPrefactor_),
    reflectionFraction_(rhs.reflectionFraction_)
{}


electronWallFluxCorrectedFvPatchScalarField::electronWallFluxCorrectedFvPatchScalarField
(
    const electronWallFluxCorrectedFvPatchScalarField& rhs
)
:
    mixedFvPatchScalarField(rhs),
    electricFieldName_(rhs.electricFieldName_),
    speciesName_(rhs.speciesName_),
    uthPrefactor_(rhs.uthPrefactor_),
    reflectionFraction_(rhs.reflectionFraction_)
{}


electronWallFluxCorrectedFvPatchScalarField::electronWallFluxCorrectedFvPatchScalarField
(
    const electronWallFluxCorrectedFvPatchScalarField& rhs,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(rhs, iF),
    electricFieldName_(rhs.electricFieldName_),
    speciesName_(rhs.speciesName_),
    uthPrefactor_(rhs.uthPrefactor_),
    reflectionFraction_(rhs.reflectionFraction_)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void electronWallFluxCorrectedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Lookup from the registered fields
    const scalarField& D = patch().lookupPatchField<volScalarField, scalar>("diffusionCoeffSpecies_e");
    const scalarField& mu = patch().lookupPatchField<volScalarField, scalar>("mobilityCoeffSpecies_e");
    const vectorField& E = patch().lookupPatchField<volVectorField, vector>(electricFieldName_);
    const scalarField& T = patch().lookupPatchField<volScalarField, scalar>("TEle_K");
    const surfaceScalarField& sumPosFluxField = this->db().lookupObject<surfaceScalarField>("sumPosFlux");
    const scalarField sumPosFlux(sumPosFluxField.boundaryField()[patch().index()]);

    // Check flux direction
    const vectorField nHat = patch().nf();
    const scalarField EdotN = E & nHat;
    const scalarField alpha = pos(-EdotN);

    scalarField C1 = -(1.0-reflectionFraction_)/(1.0+reflectionFraction_)*0.5*sqrt(uthPrefactor_*T) - ((1.0-alpha*(1.0-reflectionFraction_))/(1.0+reflectionFraction_))*2.0*mu*EdotN;
    scalarField C2 = -D;
    scalarField C3 = -(2.0/(1.0+reflectionFraction_))*sumPosFlux;

    // Get mixed boundary condition references
    this->valueFraction() = C1/(C1 + C2*patch().deltaCoeffs());
    this->refValue() = 0.0;
    this->refGrad() = C3/(C2 + SMALL);

    mixedFvPatchScalarField::updateCoeffs();
}


void electronWallFluxCorrectedFvPatchScalarField::write(Ostream& os) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("electricField")
        << electricFieldName_ << token::END_STATEMENT << nl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        electronWallFluxCorrectedFvPatchScalarField
    );
}

// ************************************************************************* //
