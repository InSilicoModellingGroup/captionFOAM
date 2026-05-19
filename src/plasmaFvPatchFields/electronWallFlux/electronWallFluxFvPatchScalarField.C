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

#include "electronWallFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
electronWallFluxFvPatchScalarField::electronWallFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureFieldName_("TGas"),
    electricFieldName_("E"),
    ionFluxFieldName_("ionFlux"),
    seec_(0.0),
    speciesName_(iF.name()),
    mass_(0.0),
    chargeNumber_(0),
    uthPrefactor_(0.0),
    initializedScalars_(false)
{
    refValue()      = 0.0;
    refGrad()       = 0.0;
    valueFraction() = 0.0;
}


electronWallFluxFvPatchScalarField::electronWallFluxFvPatchScalarField
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
    temperatureFieldName_(dict.getOrDefault<word>("tempField", "TGas")),
    electricFieldName_(dict.getOrDefault<word>("electricField", "E")),
    ionFluxFieldName_(dict.getOrDefault<word>("ionFluxField", "ionFlux")),
    seec_(dict.getOrDefault<scalar>("seec", 0.0)),
    speciesName_(iF.name()),
    mass_(0.0),
    chargeNumber_(0),
    uthPrefactor_(0.0),
    initializedScalars_(false)
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


electronWallFluxFvPatchScalarField::electronWallFluxFvPatchScalarField
(
    const electronWallFluxFvPatchScalarField& rhs,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& m
)
:
    mixedFvPatchScalarField(rhs, p, iF, m),
    temperatureFieldName_(rhs.temperatureFieldName_),
    electricFieldName_(rhs.electricFieldName_),
    ionFluxFieldName_(rhs.ionFluxFieldName_),
    seec_(rhs.seec_),
    speciesName_(rhs.speciesName_),
    mass_(rhs.mass_),
    chargeNumber_(rhs.chargeNumber_),
    uthPrefactor_(rhs.uthPrefactor_),
    initializedScalars_(rhs.initializedScalars_)
{}


electronWallFluxFvPatchScalarField::electronWallFluxFvPatchScalarField
(
    const electronWallFluxFvPatchScalarField& rhs
)
:
    mixedFvPatchScalarField(rhs),
    temperatureFieldName_(rhs.temperatureFieldName_),
    electricFieldName_(rhs.electricFieldName_),
    ionFluxFieldName_(rhs.ionFluxFieldName_),
    seec_(rhs.seec_),
    speciesName_(rhs.speciesName_),
    mass_(rhs.mass_),
    chargeNumber_(rhs.chargeNumber_),
    uthPrefactor_(rhs.uthPrefactor_),
    initializedScalars_(rhs.initializedScalars_)
{}


electronWallFluxFvPatchScalarField::electronWallFluxFvPatchScalarField
(
    const electronWallFluxFvPatchScalarField& rhs,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(rhs, iF),
    temperatureFieldName_(rhs.temperatureFieldName_),
    electricFieldName_(rhs.electricFieldName_),
    ionFluxFieldName_(rhs.ionFluxFieldName_),
    seec_(rhs.seec_),
    speciesName_(rhs.speciesName_),
    mass_(rhs.mass_),
    chargeNumber_(rhs.chargeNumber_),
    uthPrefactor_(rhs.uthPrefactor_),
    initializedScalars_(rhs.initializedScalars_)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void electronWallFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Initialise these scalar variables only once.
    // This cannot be done in the constructor because the dictionaries are not
    // available at that point from the solver.
    if (!initializedScalars_)
    {
        const IOdictionary& molarMassDict = this->db().lookupObject<IOdictionary>("molarMass");
        const IOdictionary& chargeNumberDict = this->db().lookupObject<IOdictionary>("chargeNumber");
        const scalar molarMass = molarMassDict.get<scalar>(speciesName_);

        mass_ = molarMass / constant::physicoChemical::NA.value();
        chargeNumber_ = chargeNumberDict.get<int>(speciesName_);
        uthPrefactor_ = 8.0*constant::physicoChemical::k.value()/(constant::mathematical::pi*mass_);

        initializedScalars_ = true;
    }

    // Lookup from the registered fields
    const scalarField& D = patch().lookupPatchField<volScalarField, scalar>("diffusionCoeffSpecies_" + speciesName_);
    const scalarField& mu = patch().lookupPatchField<volScalarField, scalar>("mobilityCoeffSpecies_" + speciesName_);
    const vectorField& E = patch().lookupPatchField<volVectorField, vector>(electricFieldName_);
    const scalarField& T = patch().lookupPatchField<volScalarField, scalar>(temperatureFieldName_);
    const vectorField& ionFlux = patch().lookupPatchField<volVectorField, vector>(ionFluxFieldName_);

    // Check flux direction
    const vectorField nHat = patch().nf();
    const scalarField EdotN = E & nHat;
    const scalarField ionFluxDotN = ionFlux & nHat;
    const scalarField positiveIonFlux = pos(ionFluxDotN)*ionFluxDotN;
    const scalarField alpha = pos(chargeNumber_*EdotN);

    scalarField C2 = D;
    scalarField C1 = 0.25*sqrt(uthPrefactor_*T) - (1.0-alpha)*chargeNumber_*mu*EdotN;
    const scalarField gammaSe = seec_*positiveIonFlux;

    // Get mixed boundary condition references
    this->valueFraction() = C1/(C1 + C2*patch().deltaCoeffs());
    this->refValue() = 0.0;
    this->refGrad() = gammaSe/(C2 + SMALL);

    mixedFvPatchScalarField::updateCoeffs();
}


void electronWallFluxFvPatchScalarField::write(Ostream& os) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("tempField")
        << temperatureFieldName_ << token::END_STATEMENT << nl;
    os.writeKeyword("electricField")
        << electricFieldName_ << token::END_STATEMENT << nl;
    os.writeKeyword("ionFluxField")
        << ionFluxFieldName_ << token::END_STATEMENT << nl;
    os.writeKeyword("seec")
        << seec_ << token::END_STATEMENT << nl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        electronWallFluxFvPatchScalarField
    );
}

// ************************************************************************* //
