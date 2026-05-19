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

#include "mobilityCoeffBolsig.H" 
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(mobilityCoeffBolsig, 0);
addToRunTimeSelectionTable(transportCoeffsBase, mobilityCoeffBolsig, mobility);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
mobilityCoeffBolsig::mobilityCoeffBolsig(const fvMesh& mesh, const dictionary& dict, volScalarField& transportCoeff):
    transportCoeffsBase(mesh, dict, transportCoeff),
    mobility1DTable_(),
    mobility2DTable_(),
    workingGasTableValues_(),
    mobilityLMEATables_(),
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
            mesh.time().constant(),
            mesh,
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
        localCoeffApproxPtr_ = &mesh_.lookupObject<volScalarField>("EN");
    }
    else if (localCoeffApproxType_ == "LMEA")
    {
        localCoeffApproxName_ = "meanEnergy";
        localCoeffApproxPtr_ = &mesh_.lookupObject<volScalarField>("em");
    }
    else
    {
        FatalIOErrorInFunction(physicalProperties_)
            << "Unknown localCoeffApprox type: " << localCoeffApproxType_ << nl
            << ". Valid options are LFA and LMEA"
            << exit(FatalIOError);
    }


    // Read the bolsigProperties dictionary from the constant directory
    Info<< "Reading bolsigProperties\n" << endl;
    IOdictionary bolsigProperties_
    (
        IOobject
        (
            "bolsigProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    scalarList localCoeffApproxList
    (
        bolsigProperties_.lookup(localCoeffApproxName_)
    );

    scalarList eMobList
    (
        bolsigProperties_.lookup("electronMobility")
    );

    // Read for background atmosphere single or two gas mixture
    if (backgndAtmosphereType_ == "singleGas")
    {
        // Safety guard for size mismatch in 1D table
        if (eMobList.size() != localCoeffApproxList.size())
        {
            FatalIOErrorInFunction(bolsigProperties_)
                << "Size mismatch in 1D table:" << nl
                << "  localCoeffApprox size = " << localCoeffApproxList.size() << nl
                << "  electronMobility size = " << eMobList.size()
                << exit(FatalIOError);
        }

        // Create lookup tables for mobility
        List<Tuple2<scalar, scalar>> table(localCoeffApproxList.size());

        forAll(localCoeffApproxList, i)
        {
            table[i] = Tuple2<scalar, scalar>(localCoeffApproxList[i],eMobList[i]);
        }

        mobility1DTable_ = interpolationTable<scalar>
        (
            table,
            bounds::repeatableBounding(bounds::repeatableBounding::CLAMP),
            fileName("mobility1DTable")
        );

    }
    else if (backgndAtmosphereType_ == "twoGasMixture")
    {
        workingGasName_ = backgndAtmosphere.get<word>("workingGas");
        workingGasTableValues_ = scalarList(bolsigProperties_.lookup(workingGasName_));

        if (localCoeffApproxType_ == "LFA")
        {
            // Safety guard for size mismatch in 2D table
            if (eMobList.size() != workingGasTableValues_.size()*localCoeffApproxList.size())
            {
                FatalIOErrorInFunction(bolsigProperties_)
                    << "Size mismatch in 2D LFA table:" << nl
                    << "  workingGas size = " << workingGasTableValues_.size() << nl
                    << "  localCoeffApprox size = " << localCoeffApproxList.size() << nl
                    << "  electronMobility size = " << eMobList.size() << nl
                    << "  expected size = "
                    << workingGasTableValues_.size()*localCoeffApproxList.size()
                    << exit(FatalIOError);
            }

            // Create lookup table for LFA mobility
            List<Tuple2<scalar, List<Tuple2<scalar, scalar>>>> mobTable_
            (
                workingGasTableValues_.size()
            );

            forAll(workingGasTableValues_, i)
            {
                List<Tuple2<scalar, scalar>> mobRow_(localCoeffApproxList.size());

                forAll(localCoeffApproxList, j)
                {
                    const label index = i * localCoeffApproxList.size() + j;
                    mobRow_[j] =
                        Tuple2<scalar, scalar>(localCoeffApproxList[j], eMobList[index]);
                }

                mobTable_[i] =
                    Tuple2<scalar, List<Tuple2<scalar, scalar>>>
                    (
                        workingGasTableValues_[i],
                        mobRow_
                    );
            }

            mobility2DTable_ = interpolation2DTable<scalar>
            (
                mobTable_,
                bounds::normalBounding(bounds::normalBounding::CLAMP),
                fileName("mobility2DTable")
            );
        }
        else
        {
            if (workingGasTableValues_.empty())
            {
                FatalIOErrorInFunction(bolsigProperties_)
                    << "Working-gas list is empty for LMEA table"
                    << exit(FatalIOError);
            }

            if (eMobList.size() != localCoeffApproxList.size())
            {
                FatalIOErrorInFunction(bolsigProperties_)
                    << "Size mismatch in LMEA table:" << nl
                    << "  meanEnergy size = " << localCoeffApproxList.size() << nl
                    << "  electronMobility size = " << eMobList.size()
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

            mobilityLMEATables_.setSize(workingGasTableValues_.size());

            forAll(workingGasTableValues_, i)
            {
                List<Tuple2<scalar, scalar>> mobRow_(nLocalCoeff);

                for (label j = 0; j < nLocalCoeff; ++j)
                {
                    const label index = i*nLocalCoeff + j;
                    mobRow_[j] = Tuple2<scalar, scalar>(localCoeffApproxList[index], eMobList[index]);
                }

                mobilityLMEATables_.set
                (
                    i,
                    new interpolationTable<scalar>
                    (
                        mobRow_,
                        bounds::repeatableBounding(bounds::repeatableBounding::CLAMP),
                        fileName("mobilityLMEATable_row_" + name(i))
                    )
                );
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


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void mobilityCoeffBolsig::calcMobilityCoeffs()
{
    const volScalarField& localCoeffApprox = *localCoeffApproxPtr_;

    if (backgndAtmosphereType_ == "singleGas")
    {
        forAll(transportCoeff_, cellI)
        {
            transportCoeff_[cellI] =  mobility1DTable_(localCoeffApprox[cellI]);
        }

        forAll(transportCoeff_.boundaryFieldRef(), patchI)
        {
            fvPatchScalarField& transportCoeffPatch  = transportCoeff_.boundaryFieldRef()[patchI];
            const fvPatchScalarField& localCoeffApproxPatch = localCoeffApprox.boundaryField()[patchI];

            forAll(transportCoeffPatch, faceI)
            {
                transportCoeffPatch[faceI] = mobility1DTable_(localCoeffApproxPatch[faceI]);
            }
        }
    }
    else if (backgndAtmosphereType_ == "twoGasMixture")
    {
        if (!workingGasPtr_)
        {
            workingGasPtr_ = &mesh_.lookupObject<volScalarField>(workingGasName_);
        }

        const volScalarField& workingGas = *workingGasPtr_;

        if (localCoeffApproxType_ == "LFA")
        {
            // Internal field
            forAll(transportCoeff_, cellI)
            {
                transportCoeff_[cellI] =
                    mobility2DTable_(workingGas[cellI], localCoeffApprox[cellI]);
            }

            // Boundary fields
            forAll(transportCoeff_.boundaryFieldRef(), patchI)
            {
                fvPatchScalarField& transportCoeffPatch  = transportCoeff_.boundaryFieldRef()[patchI];
                const fvPatchScalarField& workingGasPatch = workingGas.boundaryField()[patchI];
                const fvPatchScalarField& localCoeffApproxPatch = localCoeffApprox.boundaryField()[patchI];

                forAll(transportCoeffPatch, faceI)
                {
                    transportCoeffPatch[faceI] =
                        mobility2DTable_(workingGasPatch[faceI], localCoeffApproxPatch[faceI]);
                }
            }
        }
        else if (localCoeffApproxType_ == "LMEA")
        {
            // Internal field
            forAll(transportCoeff_, cellI)
            {
                const scalar gasValue = workingGas[cellI];
                const scalar approxValue = localCoeffApprox[cellI];

                label i0, i1;
                scalar w;
                lowerUpperListValues(workingGasTableValues_, gasValue, i0, i1, w);

                const scalar mu0 = mobilityLMEATables_[i0](approxValue);
                const scalar mu1 = mobilityLMEATables_[i1](approxValue);

                transportCoeff_[cellI] = (i0 == i1 ? mu0 : (1.0 - w)*mu0 + w*mu1);
            }

            // Boundary fields
            forAll(transportCoeff_.boundaryFieldRef(), patchI)
            {
                fvPatchScalarField& transportCoeffPatch  = transportCoeff_.boundaryFieldRef()[patchI];
                const fvPatchScalarField& workingGasPatch = workingGas.boundaryField()[patchI];
                const fvPatchScalarField& localCoeffApproxPatch = localCoeffApprox.boundaryField()[patchI];

                forAll(transportCoeffPatch, faceI)
                {
                    const scalar gasValue = workingGasPatch[faceI];
                    const scalar approxValue = localCoeffApproxPatch[faceI];

                    label i0, i1;
                    scalar w;
                    lowerUpperListValues(workingGasTableValues_, gasValue, i0, i1, w);

                    const scalar mu0 = mobilityLMEATables_[i0](approxValue);
                    const scalar mu1 = mobilityLMEATables_[i1](approxValue);

                    transportCoeffPatch[faceI] = (i0 == i1 ? mu0 : (1.0 - w)*mu0 + w*mu1);
                }
            }
        }
    }

    // Correct boundary conditions (useful for mpi interfaces)
    transportCoeff_.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

