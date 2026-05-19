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

#include "bolsigRateCoeff.H"
#include "plasmaChemistryModel.H"
#include "volFields.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(bolsigRateCoeff, 0);
addToRunTimeSelectionTable(reactionRateCoeffsBase, bolsigRateCoeff, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
bolsigRateCoeff::bolsigRateCoeff
(
    const dictionary& dict, 
    const word& name, 
    plasmaChemistryModel& chemistry
)
:
    reactionIndex_(),
    bolsigProperties_
    (
        IOobject
        (
            "bolsigProperties",
            chemistry.mesh().time().constant(),
            chemistry.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    rates1DTable_(),
    rates2DTable_(),
    ratesLMEATables_(),
    workingGasTableValues_(),
    localCoeffApproxType_(),
    backgndAtmosphereType_(),
    workingGasName_(),
    workingGasPtr_(nullptr),
    localCoeffApproxName_(),
    localCoeffApproxPtr_(nullptr)
{
    // Read physicalProperties dictionary from the constant directory
    Info<< "Reading physicalProperties\n" << endl;
    IOdictionary physicalProperties_
    (
        IOobject
        (
            "physicalProperties",
            chemistry.mesh().time().constant(),
            chemistry.mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    const dictionary& plasmaCoeffApprox = physicalProperties_.subDict("plasmaCoeffApprox");
    const dictionary& localCoeffApprox = plasmaCoeffApprox.subDict("localCoeffApprox");
    const dictionary& backgndAtmosphere = plasmaCoeffApprox.subDict("backgndAtmosphere");
    localCoeffApproxType_ = localCoeffApprox.get<word>("type");
    backgndAtmosphereType_ = backgndAtmosphere.get<word>("type");

    // Read for LFA or LMEA
    if (localCoeffApproxType_ == "LFA")
    {
        localCoeffApproxName_ = "EN";
    }
    else if (localCoeffApproxType_ == "LMEA")
    {
        localCoeffApproxName_ = "meanEnergy";
    }
    else
    {
        FatalIOErrorInFunction(physicalProperties_)
            << "Unknown localCoeffApprox type: " << localCoeffApproxType_ << nl
            << ". Valid options are LFA and LMEA"
            << exit(FatalIOError);
    }

    scalarList localCoeffApproxList
    (
        bolsigProperties_.lookup(localCoeffApproxName_)
    );
    
    List<List<scalar>> rateCoefficientsList_
    (
        bolsigProperties_.lookup("rateCoefficients")
    );

    if (rateCoefficientsList_.empty())
    {
        FatalIOErrorInFunction(bolsigProperties_)
            << "rateCoefficients list is empty"
            << exit(FatalIOError);
    }

    // Read for background atmosphere single or two gas mixture
    if (backgndAtmosphereType_ == "singleGas")
    {
        // Safety guard for size mismatch in 1D table
        if (rateCoefficientsList_.size() != localCoeffApproxList.size())
        {
            FatalIOErrorInFunction(bolsigProperties_)
                << "Size mismatch in 1D table:" << nl
                << "  localCoeffApprox size = " << localCoeffApproxList.size() << nl
                << "  rateCoefficients size = " << rateCoefficientsList_.size()
                << exit(FatalIOError);
        }

        // Initialise pointer list size
        rates1DTable_.setSize(rateCoefficientsList_[0].size());

        forAll(rateCoefficientsList_[0], r)
        {
            // Create lookup table for each reaction
            List<Tuple2<scalar, scalar>> rTable(localCoeffApproxList.size());

            forAll(localCoeffApproxList, i)
            {
                rTable[i] =
                    Tuple2<scalar, scalar>(localCoeffApproxList[i], rateCoefficientsList_[i][r]);
            }

            rates1DTable_.set
            (
                r,
                new interpolationTable<scalar>
                (
                    rTable,
                    bounds::repeatableBounding(bounds::repeatableBounding::CLAMP),
                    fileName("rates1DTable")
                )
            );
        }
    }
    else if (backgndAtmosphereType_ == "twoGasMixture")
    {
        workingGasName_ = backgndAtmosphere.get<word>("workingGas");
        workingGasTableValues_ = scalarList(bolsigProperties_.lookup(workingGasName_));

        if (localCoeffApproxType_ == "LFA")
        {
            // Safety guard for size mismatch in 2D table
            if (rateCoefficientsList_.size() != workingGasTableValues_.size()*localCoeffApproxList.size())
            {
                FatalIOErrorInFunction(bolsigProperties_)
                    << "Size mismatch in 2D LFA table:" << nl
                    << "  workingGas size = " << workingGasTableValues_.size() << nl
                    << "  localCoeffApprox size = " << localCoeffApproxList.size() << nl
                    << "  rateCoefficients size = " << rateCoefficientsList_.size() << nl
                    << "  expected size = "
                    << workingGasTableValues_.size()*localCoeffApproxList.size()
                    << exit(FatalIOError);
            }

            // Initialise pointer list size
            rates2DTable_.setSize(rateCoefficientsList_[0].size());

            forAll(rateCoefficientsList_[0], r)
            {
                // Create lookup table for each reaction
                List<Tuple2<scalar, List<Tuple2<scalar, scalar>>>> rTable(workingGasTableValues_.size());

                forAll(workingGasTableValues_, i)
                {
                    List<Tuple2<scalar, scalar>> rRow(localCoeffApproxList.size());

                    forAll(localCoeffApproxList, j)
                    {
                        const label index = i * localCoeffApproxList.size() + j;
                        rRow[j] =
                            Tuple2<scalar, scalar>(localCoeffApproxList[j], rateCoefficientsList_[index][r]);
                    }
                    rTable[i] = Tuple2<scalar, List<Tuple2<scalar, scalar>>>(workingGasTableValues_[i], rRow);
                }

                rates2DTable_.set
                (
                    r,
                    new interpolation2DTable<scalar>
                    (
                        rTable,
                        bounds::normalBounding(bounds::normalBounding::CLAMP),
                        fileName("rates2DTable")
                    )
                );
            }
        }
        else if (localCoeffApproxType_ == "LMEA")
        {
            if (workingGasTableValues_.empty())
            {
                FatalIOErrorInFunction(bolsigProperties_)
                    << "Working-gas list is empty for LMEA table"
                    << exit(FatalIOError);
            }

            if (rateCoefficientsList_.size() != localCoeffApproxList.size())
            {
                FatalIOErrorInFunction(bolsigProperties_)
                    << "Size mismatch in LMEA table:" << nl
                    << "  meanEnergy size = " << localCoeffApproxList.size() << nl
                    << "  rateCoefficients size = " << rateCoefficientsList_.size()
                    << exit(FatalIOError);
            }

            if (localCoeffApproxList.size() % workingGasTableValues_.size() != 0)
            {
                FatalIOErrorInFunction(bolsigProperties_)
                    << "Invalid flattened LMEA table:" << nl
                    << "  meanEnergy size = " << localCoeffApproxList.size() << nl
                    << "  workingGas size = " << workingGasTableValues_.size() << nl
                    << "  meanEnergy size must be divisible by workingGas size"
                    << exit(FatalIOError);
            }

            const label nLocalCoeff = localCoeffApproxList.size()/workingGasTableValues_.size();

            // Initialise pointer list size
            ratesLMEATables_.setSize(rateCoefficientsList_[0].size());

            forAll(rateCoefficientsList_[0], r)
            {
                PtrList<interpolationTable<scalar>>* rTablePtr =
                    new PtrList<interpolationTable<scalar>>(workingGasTableValues_.size());

                forAll(workingGasTableValues_, i)
                {
                    List<Tuple2<scalar, scalar>> rRow(nLocalCoeff);

                    for (label j = 0; j < nLocalCoeff; ++j)
                    {
                        const label index = i*nLocalCoeff + j;
                        rRow[j] =
                            Tuple2<scalar, scalar>(localCoeffApproxList[index], rateCoefficientsList_[index][r]);
                    }

                    rTablePtr->set
                    (
                        i,
                        new interpolationTable<scalar>
                        (
                            rRow,
                            bounds::repeatableBounding(bounds::repeatableBounding::CLAMP),
                            fileName("ratesLMEATable")
                        )
                    );
                }

                ratesLMEATables_.set(r, rTablePtr);
            }
        }
    }
    else
    {
        FatalIOErrorInFunction(physicalProperties_)
            << "Unknown backgndAtmosphere type: "
            << backgndAtmosphereType_ << nl
            << "Valid options: singleGas, twoGasMixture"
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void bolsigRateCoeff::calculate(plasmaChemistryModel& chemistry, const label j) const
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
            kj[cellI] = rates1DTable_[j](localCoeffApprox[cellI]);
        }

        // Boundary fields
        forAll(kj.boundaryFieldRef(), patchI)
        {
            fvPatchScalarField& kjPatch = kj.boundaryFieldRef()[patchI];
            const fvPatchScalarField& localCoeffApproxPatch = localCoeffApprox.boundaryField()[patchI];

            forAll(kjPatch, faceI)
            {
                kjPatch[faceI] = rates1DTable_[j](localCoeffApproxPatch[faceI]);
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
                kj[cellI] = rates2DTable_[j](workingGas[cellI], localCoeffApprox[cellI]);
            }

            // Boundary fields
            forAll(kj.boundaryFieldRef(), patchI)
            {
                fvPatchScalarField& kjPatch = kj.boundaryFieldRef()[patchI];
                const fvPatchScalarField& workingGasPatch = workingGas.boundaryField()[patchI];
                const fvPatchScalarField& localCoeffApproxPatch = localCoeffApprox.boundaryField()[patchI];

                forAll(kjPatch, faceI)
                {
                    kjPatch[faceI] = rates2DTable_[j](workingGasPatch[faceI], localCoeffApproxPatch[faceI]);
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

                kj[cellI] = (i0 == i1 ? k0 : (1.0 - w)*k0 + w*k1);
            }

            // Boundary fields
            forAll(kj.boundaryFieldRef(), patchI)
            {
                fvPatchScalarField& kjPatch = kj.boundaryFieldRef()[patchI];
                const fvPatchScalarField& workingGasPatch = workingGas.boundaryField()[patchI];
                const fvPatchScalarField& localCoeffApproxPatch = localCoeffApprox.boundaryField()[patchI];

                forAll(kjPatch, faceI)
                {
                    const scalar gasValue = workingGasPatch[faceI];
                    const scalar approxValue = localCoeffApproxPatch[faceI];

                    label i0, i1;
                    scalar w;
                    lowerUpperListValues(workingGasTableValues_, gasValue, i0, i1, w);

                    const scalar k0 = ratesLMEATables_[j][i0](approxValue);
                    const scalar k1 = ratesLMEATables_[j][i1](approxValue);

                    kjPatch[faceI] = (i0 == i1 ? k0 : (1.0 - w)*k0 + w*k1);
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
