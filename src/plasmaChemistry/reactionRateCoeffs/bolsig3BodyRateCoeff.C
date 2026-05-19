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

#include "bolsig3BodyRateCoeff.H"
#include "addToRunTimeSelectionTable.H"
#include "plasmaChemistryModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(bolsig3BodyRateCoeff, 0);
addToRunTimeSelectionTable(reactionRateCoeffsBase, bolsig3BodyRateCoeff, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
bolsig3BodyRateCoeff::bolsig3BodyRateCoeff
(
    const dictionary& dict,
    const word& name,
    plasmaChemistryModel& chemistry
)
:
    Foam::bolsigRateCoeff(dict, name, chemistry),
    NGas_(chemistry.mesh().lookupObjectRef<volScalarField>("NGas"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void bolsig3BodyRateCoeff::calculate(plasmaChemistryModel& chemistry, const label j) const
{
    if (!localCoeffApproxPtr_)
    {
        if (localCoeffApproxType_ == "LFA")
        {
            localCoeffApproxPtr_ =
                &chemistry.mesh().lookupObject<volScalarField>("EN");
        }
        else if (localCoeffApproxType_ == "LMEA")
        {
            localCoeffApproxPtr_ =
                &chemistry.mesh().lookupObject<volScalarField>("em");
        }
    }

    const volScalarField& localCoeffApprox = *localCoeffApproxPtr_;
    volScalarField& kj = chemistry.k()[j];

    if (backgndAtmosphereType_ == "singleGas")
    {
        // Internal field
        forAll(kj, cellI)
        {
            scalar twoBodyRate = rates1DTable_[j](localCoeffApprox[cellI]);
            kj[cellI] = twoBodyRate * NGas_[cellI];
        }

        // Boundary fields
        forAll(kj.boundaryFieldRef(), patchI)
        {
            fvPatchScalarField& kjPatch = kj.boundaryFieldRef()[patchI];
            const fvPatchScalarField& localCoeffApproxPatch = localCoeffApprox.boundaryField()[patchI];
            const fvPatchScalarField& NGasPatch = NGas_.boundaryField()[patchI];

            forAll(kjPatch, faceI)
            {
                scalar twoBodyRate = rates1DTable_[j](localCoeffApproxPatch[faceI]);
                kjPatch[faceI] = twoBodyRate * NGasPatch[faceI];
            }
        }
    }
    else if (backgndAtmosphereType_ == "twoGasMixture")
    {
        if (!workingGasPtr_)
        {
            if (chemistry.species().find(workingGasName_) < 0)
            {
                FatalIOErrorInFunction(bolsigProperties_)
                    << "Working gas '" << workingGasName_
                    << "' was not found in reactionsDict species list: "
                    << chemistry.species() << nl
                    << "The twoGasMixture BolSig model expects workingGas "
                    << "to name a species number-density field."
                    << exit(FatalIOError);
            }

            workingGasPtr_ =
                &chemistry.mesh().lookupObject<volScalarField>(workingGasName_);
        }

        const volScalarField& workingGas = *workingGasPtr_;

        if (localCoeffApproxType_ == "LFA")
        {
            // Internal field
            forAll(kj, cellI)
            {
                scalar twoBodyRate = rates2DTable_[j](workingGas[cellI], localCoeffApprox[cellI]);
                kj[cellI] = twoBodyRate * NGas_[cellI];
            }

            // Boundary fields
            forAll(kj.boundaryFieldRef(), patchI)
            {
                fvPatchScalarField& kjPatch = kj.boundaryFieldRef()[patchI];
                const fvPatchScalarField& workingGasPatch = workingGas.boundaryField()[patchI];
                const fvPatchScalarField& localCoeffApproxPatch = localCoeffApprox.boundaryField()[patchI];
                const fvPatchScalarField& NGasPatch = NGas_.boundaryField()[patchI];

                forAll(kjPatch, faceI)
                {
                    scalar twoBodyRate = rates2DTable_[j](workingGasPatch[faceI], localCoeffApproxPatch[faceI]);
                    kjPatch[faceI] = twoBodyRate * NGasPatch[faceI];
                }
            }
        }
        else if (localCoeffApproxType_ == "LMEA")
        {
            // Internal field
            forAll(kj, cellI)
            {
                const scalar gasValue = workingGas[cellI];
                const scalar approxValue = localCoeffApprox[cellI];

                label i0, i1;
                scalar w;
                lowerUpperListValues(workingGasTableValues_, gasValue, i0, i1, w);

                const scalar k0 = ratesLMEATables_[j][i0](approxValue);
                const scalar k1 = ratesLMEATables_[j][i1](approxValue);
                const scalar twoBodyRate = (i0 == i1 ? k0 : (1.0 - w)*k0 + w*k1);

                kj[cellI] = twoBodyRate * NGas_[cellI];
            }

            // Boundary fields
            forAll(kj.boundaryFieldRef(), patchI)
            {
                fvPatchScalarField& kjPatch = kj.boundaryFieldRef()[patchI];
                const fvPatchScalarField& workingGasPatch = workingGas.boundaryField()[patchI];
                const fvPatchScalarField& localCoeffApproxPatch = localCoeffApprox.boundaryField()[patchI];
                const fvPatchScalarField& NGasPatch = NGas_.boundaryField()[patchI];

                forAll(kjPatch, faceI)
                {
                    const scalar gasValue = workingGasPatch[faceI];
                    const scalar approxValue = localCoeffApproxPatch[faceI];

                    label i0, i1;
                    scalar w;
                    lowerUpperListValues(workingGasTableValues_, gasValue, i0, i1, w);

                    const scalar k0 = ratesLMEATables_[j][i0](approxValue);
                    const scalar k1 = ratesLMEATables_[j][i1](approxValue);
                    const scalar twoBodyRate = (i0 == i1 ? k0 : (1.0 - w)*k0 + w*k1);

                    kjPatch[faceI] = twoBodyRate * NGasPatch[faceI];
                }
            }
        }
    }

    // Correct boundary conditions (useful for mpi interfaces)
    kj.correctBoundaryConditions();
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
