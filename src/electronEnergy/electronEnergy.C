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

#include "electronEnergy.H" 
#include "volFields.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
electronEnergy::electronEnergy
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    energyMobility1DTable_(),
    energyDiffusivity1DTable_(),
    powerLoss1DTable_(),
    energyMobility2DTable_(),
    energyDiffusivity2DTable_(),
    powerLoss2DTable_(),
    workingGasTableValues_(),
    energyMobilityLMEATables_(),
    energyDiffusivityLMEATables_(),
    powerLossLMEATables_(),
    EN_(mesh_.lookupObjectRef<volScalarField>("EN")),
    em_(mesh_.lookupObjectRef<volScalarField>("em")),
    localCoeffApproxType_(),
    backgndAtmosphereType_(),
    workingGasName_(),
    workingGasPtr_(nullptr),
    localCoeffApproxName_(),
    localCoeffApproxPtr_(nullptr),
    muE_(mesh_.lookupObjectRef<volScalarField>("muE")),
    DE_(mesh_.lookupObjectRef<volScalarField>("DE")),
    theta_
    (
        IOobject
        (
            "theta",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("temp", dimensionSet(1,2,-3,0,0,0,0), scalar(0))
    ),
    isConstant_(false)
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
        localCoeffApproxPtr_ = &EN_;
    }
    else if (localCoeffApproxType_ == "LMEA")
    {
        localCoeffApproxName_ = "meanEnergy";
        localCoeffApproxPtr_ = &em_;
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

    scalarList energyMobList
    (
        bolsigProperties_.lookup("energyMobility")
    );

    scalarList energyDiffList
    (
        bolsigProperties_.lookup("energyDiffusivity")
    );

    scalarList powerLossList
    (
        bolsigProperties_.lookup("powerLoss")
    );

    // Read for background atmosphere single or two gas mixture
    if (backgndAtmosphereType_ == "singleGas")
    {
        // Safety guard for size mismatch in 1D table
        if (energyMobList.size() != localCoeffApproxList.size())
        {
            FatalIOErrorInFunction(bolsigProperties_)
                << "Size mismatch in 1D table:" << nl
                << "  localCoeffApprox size = " << localCoeffApproxList.size() << nl
                << "  energyMobility size = " << energyMobList.size()
                << exit(FatalIOError);
        }

        if (energyDiffList.size() != localCoeffApproxList.size())
        {
            FatalIOErrorInFunction(bolsigProperties_)
                << "Size mismatch in 1D table:" << nl
                << "  localCoeffApprox size = " << localCoeffApproxList.size() << nl
                << "  energyDiffusivity size = " << energyDiffList.size()
                << exit(FatalIOError);
        }

        if (powerLossList.size() != localCoeffApproxList.size())
        {
            FatalIOErrorInFunction(bolsigProperties_)
                << "Size mismatch in 1D table:" << nl
                << "  localCoeffApprox size = " << localCoeffApproxList.size() << nl
                << "  powerLoss size = " << powerLossList.size()
                << exit(FatalIOError);
        }

        // Create lookup tables for energy transport coefficients
        List<Tuple2<scalar, scalar>> energyMobTable(localCoeffApproxList.size());
        List<Tuple2<scalar, scalar>> energyDiffTable(localCoeffApproxList.size());
        List<Tuple2<scalar, scalar>> powerLossTable(localCoeffApproxList.size());

        forAll(localCoeffApproxList, i)
        {
            energyMobTable[i] = Tuple2<scalar, scalar>(localCoeffApproxList[i],energyMobList[i]);
            energyDiffTable[i] = Tuple2<scalar, scalar>(localCoeffApproxList[i],energyDiffList[i]);
            powerLossTable[i] = Tuple2<scalar, scalar>(localCoeffApproxList[i],powerLossList[i]);
        }

        energyMobility1DTable_ = interpolationTable<scalar>
        (
            energyMobTable,
            bounds::repeatableBounding(bounds::repeatableBounding::CLAMP),
            fileName("energyMobility1DTable")
        );

        energyDiffusivity1DTable_ = interpolationTable<scalar>
        (
            energyDiffTable,
            bounds::repeatableBounding(bounds::repeatableBounding::CLAMP),
            fileName("energyDiffusivity1DTable")
        );

        powerLoss1DTable_ = interpolationTable<scalar>
        (
            powerLossTable,
            bounds::repeatableBounding(bounds::repeatableBounding::CLAMP),
            fileName("powerLoss1DTable")
        );

    }
    else if (backgndAtmosphereType_ == "twoGasMixture")
    {
        workingGasName_ = backgndAtmosphere.get<word>("workingGas");
        workingGasPtr_ = &mesh_.lookupObjectRef<volScalarField>(workingGasName_);
        workingGasTableValues_ = scalarList(bolsigProperties_.lookup(workingGasName_));

        if (localCoeffApproxType_ == "LFA")
        {
            // Safety guard for size mismatch in 2D table
            if (energyMobList.size() != workingGasTableValues_.size()*localCoeffApproxList.size())
            {
                FatalIOErrorInFunction(bolsigProperties_)
                    << "Size mismatch in 2D LFA table:" << nl
                    << "  workingGas size = " << workingGasTableValues_.size() << nl
                    << "  localCoeffApprox size = " << localCoeffApproxList.size() << nl
                    << "  energyMobility size = " << energyMobList.size() << nl
                    << "  expected size = "
                    << workingGasTableValues_.size()*localCoeffApproxList.size()
                    << exit(FatalIOError);
            }

            if (energyDiffList.size() != workingGasTableValues_.size()*localCoeffApproxList.size())
            {
                FatalIOErrorInFunction(bolsigProperties_)
                    << "Size mismatch in 2D LFA table:" << nl
                    << "  workingGas size = " << workingGasTableValues_.size() << nl
                    << "  localCoeffApprox size = " << localCoeffApproxList.size() << nl
                    << "  energyDiffusivity size = " << energyDiffList.size() << nl
                    << "  expected size = "
                    << workingGasTableValues_.size()*localCoeffApproxList.size()
                    << exit(FatalIOError);
            }

            if (powerLossList.size() != workingGasTableValues_.size()*localCoeffApproxList.size())
            {
                FatalIOErrorInFunction(bolsigProperties_)
                    << "Size mismatch in 2D LFA table:" << nl
                    << "  workingGas size = " << workingGasTableValues_.size() << nl
                    << "  localCoeffApprox size = " << localCoeffApproxList.size() << nl
                    << "  powerLoss size = " << powerLossList.size() << nl
                    << "  expected size = "
                    << workingGasTableValues_.size()*localCoeffApproxList.size()
                    << exit(FatalIOError);
            }

            // Create lookup tables for LFA energy transport coefficients
            List<Tuple2<scalar, List<Tuple2<scalar, scalar>>>> energyMobTable_
            (
                workingGasTableValues_.size()
            );
            List<Tuple2<scalar, List<Tuple2<scalar, scalar>>>> energyDiffTable_
            (
                workingGasTableValues_.size()
            );
            List<Tuple2<scalar, List<Tuple2<scalar, scalar>>>> powerLossTable_
            (
                workingGasTableValues_.size()
            );

            forAll(workingGasTableValues_, i)
            {
                List<Tuple2<scalar, scalar>> energyMobRow_(localCoeffApproxList.size());
                List<Tuple2<scalar, scalar>> energyDiffRow_(localCoeffApproxList.size());
                List<Tuple2<scalar, scalar>> powerLossRow_(localCoeffApproxList.size());

                forAll(localCoeffApproxList, j)
                {
                    const label index = i * localCoeffApproxList.size() + j;
                    energyMobRow_[j] =
                        Tuple2<scalar, scalar>(localCoeffApproxList[j], energyMobList[index]);
                    energyDiffRow_[j] =
                        Tuple2<scalar, scalar>(localCoeffApproxList[j], energyDiffList[index]);
                    powerLossRow_[j] =
                        Tuple2<scalar, scalar>(localCoeffApproxList[j], powerLossList[index]);
                }

                energyMobTable_[i] =
                    Tuple2<scalar, List<Tuple2<scalar, scalar>>>
                    (
                        workingGasTableValues_[i],
                        energyMobRow_
                    );

                energyDiffTable_[i] =
                    Tuple2<scalar, List<Tuple2<scalar, scalar>>>
                    (
                        workingGasTableValues_[i],
                        energyDiffRow_
                    );

                powerLossTable_[i] =
                    Tuple2<scalar, List<Tuple2<scalar, scalar>>>
                    (
                        workingGasTableValues_[i],
                        powerLossRow_
                    );
            }

            energyMobility2DTable_ = interpolation2DTable<scalar>
            (
                energyMobTable_,
                bounds::normalBounding(bounds::normalBounding::CLAMP),
                fileName("energyMobility2DTable")
            );

            energyDiffusivity2DTable_ = interpolation2DTable<scalar>
            (
                energyDiffTable_,
                bounds::normalBounding(bounds::normalBounding::CLAMP),
                fileName("energyDiffusivity2DTable")
            );

            powerLoss2DTable_ = interpolation2DTable<scalar>
            (
                powerLossTable_,
                bounds::normalBounding(bounds::normalBounding::CLAMP),
                fileName("powerLoss2DTable")
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

            if (energyMobList.size() != localCoeffApproxList.size())
            {
                FatalIOErrorInFunction(bolsigProperties_)
                    << "Size mismatch in LMEA table:" << nl
                    << "  meanEnergy size = " << localCoeffApproxList.size() << nl
                    << "  energyMobility size = " << energyMobList.size()
                    << exit(FatalIOError);
            }

            if (energyDiffList.size() != localCoeffApproxList.size())
            {
                FatalIOErrorInFunction(bolsigProperties_)
                    << "Size mismatch in LMEA table:" << nl
                    << "  meanEnergy size = " << localCoeffApproxList.size() << nl
                    << "  energyDiffusivity size = " << energyDiffList.size()
                    << exit(FatalIOError);
            }

            if (powerLossList.size() != localCoeffApproxList.size())
            {
                FatalIOErrorInFunction(bolsigProperties_)
                    << "Size mismatch in LMEA table:" << nl
                    << "  meanEnergy size = " << localCoeffApproxList.size() << nl
                    << "  powerLoss size = " << powerLossList.size()
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

            energyMobilityLMEATables_.setSize(workingGasTableValues_.size());
            energyDiffusivityLMEATables_.setSize(workingGasTableValues_.size());
            powerLossLMEATables_.setSize(workingGasTableValues_.size());

            forAll(workingGasTableValues_, i)
            {
                List<Tuple2<scalar, scalar>> energyMobRow_(nLocalCoeff);
                List<Tuple2<scalar, scalar>> energyDiffRow_(nLocalCoeff);
                List<Tuple2<scalar, scalar>> powerLossRow_(nLocalCoeff);

                for (label j = 0; j < nLocalCoeff; ++j)
                {
                    const label index = i*nLocalCoeff + j;
                    energyMobRow_[j] = Tuple2<scalar, scalar>(localCoeffApproxList[index], energyMobList[index]);
                    energyDiffRow_[j] = Tuple2<scalar, scalar>(localCoeffApproxList[index], energyDiffList[index]);
                    powerLossRow_[j] = Tuple2<scalar, scalar>(localCoeffApproxList[index], powerLossList[index]);
                }

                energyMobilityLMEATables_.set
                (
                    i,
                    new interpolationTable<scalar>
                    (
                        energyMobRow_,
                        bounds::repeatableBounding(bounds::repeatableBounding::CLAMP),
                        fileName("energyMobilityLMEATable_row_" + name(i))
                    )
                );

                energyDiffusivityLMEATables_.set
                (
                    i,
                    new interpolationTable<scalar>
                    (
                        energyDiffRow_,
                        bounds::repeatableBounding(bounds::repeatableBounding::CLAMP),
                        fileName("energyDiffusivityLMEATable_row_" + name(i))
                    )
                );

                powerLossLMEATables_.set
                (
                    i,
                    new interpolationTable<scalar>
                    (
                        powerLossRow_,
                        bounds::repeatableBounding(bounds::repeatableBounding::CLAMP),
                        fileName("powerLossLMEATable_row_" + name(i))
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
void electronEnergy::calculateTransportCoeffs()
{
    if (isConstant_)
    {
        return;
    }

    const volScalarField& localCoeffApprox = *localCoeffApproxPtr_;

    if (backgndAtmosphereType_ == "singleGas")
    {
        forAll(muE_, cellI)
        {
            muE_[cellI] =  energyMobility1DTable_(localCoeffApprox[cellI]);
            DE_[cellI] =  energyDiffusivity1DTable_(localCoeffApprox[cellI]);
        }

        forAll(muE_.boundaryFieldRef(), patchI)
        {
            fvPatchScalarField& muEPatch  = muE_.boundaryFieldRef()[patchI];
            fvPatchScalarField& DEPatch  = DE_.boundaryFieldRef()[patchI];
            const fvPatchScalarField& localCoeffApproxPatch = localCoeffApprox.boundaryField()[patchI];

            forAll(muEPatch, faceI)
            {
                muEPatch[faceI] = energyMobility1DTable_(localCoeffApproxPatch[faceI]);
                DEPatch[faceI] = energyDiffusivity1DTable_(localCoeffApproxPatch[faceI]);
            }
        }
    }
    else if (backgndAtmosphereType_ == "twoGasMixture")
    {
        const volScalarField& workingGas = *workingGasPtr_;

        if (localCoeffApproxType_ == "LFA")
        {
            // Internal field
            forAll(muE_, cellI)
            {
                muE_[cellI] =
                    energyMobility2DTable_(workingGas[cellI], localCoeffApprox[cellI]);
                DE_[cellI] =
                    energyDiffusivity2DTable_(workingGas[cellI], localCoeffApprox[cellI]);
            }

            // Boundary fields
            forAll(muE_.boundaryFieldRef(), patchI)
            {
                fvPatchScalarField& muEPatch  = muE_.boundaryFieldRef()[patchI];
                fvPatchScalarField& DEPatch  = DE_.boundaryFieldRef()[patchI];
                const fvPatchScalarField& workingGasPatch = workingGas.boundaryField()[patchI];
                const fvPatchScalarField& localCoeffApproxPatch = localCoeffApprox.boundaryField()[patchI];

                forAll(muEPatch, faceI)
                {
                    muEPatch[faceI] =
                        energyMobility2DTable_(workingGasPatch[faceI], localCoeffApproxPatch[faceI]);
                    DEPatch[faceI] =
                        energyDiffusivity2DTable_(workingGasPatch[faceI], localCoeffApproxPatch[faceI]);
                }
            }
        }
        else if (localCoeffApproxType_ == "LMEA")
        {
            // Internal field
            forAll(muE_, cellI)
            {
                const scalar gasValue = workingGas[cellI];
                const scalar approxValue = localCoeffApprox[cellI];

                label i0, i1;
                scalar w;
                lowerUpperListValues(workingGasTableValues_, gasValue, i0, i1, w);

                const scalar mu0 = energyMobilityLMEATables_[i0](approxValue);
                const scalar mu1 = energyMobilityLMEATables_[i1](approxValue);

                muE_[cellI] = (i0 == i1 ? mu0 : (1.0 - w)*mu0 + w*mu1);

                const scalar D0 = energyDiffusivityLMEATables_[i0](approxValue);
                const scalar D1 = energyDiffusivityLMEATables_[i1](approxValue);

                DE_[cellI] = (i0 == i1 ? D0 : (1.0 - w)*D0 + w*D1);
            }

            // Boundary fields
            forAll(muE_.boundaryFieldRef(), patchI)
            {
                fvPatchScalarField& muEPatch  = muE_.boundaryFieldRef()[patchI];
                fvPatchScalarField& DEPatch  = DE_.boundaryFieldRef()[patchI];
                const fvPatchScalarField& workingGasPatch = workingGas.boundaryField()[patchI];
                const fvPatchScalarField& localCoeffApproxPatch = localCoeffApprox.boundaryField()[patchI];

                forAll(muEPatch, faceI)
                {
                    const scalar gasValue = workingGasPatch[faceI];
                    const scalar approxValue = localCoeffApproxPatch[faceI];

                    label i0, i1;
                    scalar w;
                    lowerUpperListValues(workingGasTableValues_, gasValue, i0, i1, w);

                    const scalar mu0 = energyMobilityLMEATables_[i0](approxValue);
                    const scalar mu1 = energyMobilityLMEATables_[i1](approxValue);

                    muEPatch[faceI] = (i0 == i1 ? mu0 : (1.0 - w)*mu0 + w*mu1);

                    const scalar D0 = energyDiffusivityLMEATables_[i0](approxValue);
                    const scalar D1 = energyDiffusivityLMEATables_[i1](approxValue);

                    DEPatch[faceI] = (i0 == i1 ? D0 : (1.0 - w)*D0 + w*D1);
                }
            }
        }
    }

    // Correct boundary conditions (useful for mpi interfaces)
    muE_.correctBoundaryConditions();
    DE_.correctBoundaryConditions();
}

const volScalarField& electronEnergy::ThetaE()
{
    const volScalarField& localCoeffApprox = *localCoeffApproxPtr_;

    if (backgndAtmosphereType_ == "singleGas")
    {
        forAll(theta_, cellI)
        {
            theta_[cellI] =  powerLoss1DTable_(localCoeffApprox[cellI]);
        }

        forAll(theta_.boundaryFieldRef(), patchI)
        {
            fvPatchScalarField& thetaPatch = theta_.boundaryFieldRef()[patchI];
            const fvPatchScalarField& localCoeffApproxPatch = localCoeffApprox.boundaryField()[patchI];

            forAll(thetaPatch, faceI)
            {
                thetaPatch[faceI] = powerLoss1DTable_(localCoeffApproxPatch[faceI]);
            }
        }
    }
    else if (backgndAtmosphereType_ == "twoGasMixture")
    {
        const volScalarField& workingGas = *workingGasPtr_;

        if (localCoeffApproxType_ == "LFA")
        {
            // Internal field
            forAll(theta_, cellI)
            {
                theta_[cellI] =
                    powerLoss2DTable_(workingGas[cellI], localCoeffApprox[cellI]);
            }

            // Boundary fields
            forAll(theta_.boundaryFieldRef(), patchI)
            {
                fvPatchScalarField& thetaPatch = theta_.boundaryFieldRef()[patchI];
                const fvPatchScalarField& workingGasPatch = workingGas.boundaryField()[patchI];
                const fvPatchScalarField& localCoeffApproxPatch = localCoeffApprox.boundaryField()[patchI];

                forAll(thetaPatch, faceI)
                {
                    thetaPatch[faceI] =
                        powerLoss2DTable_(workingGasPatch[faceI], localCoeffApproxPatch[faceI]);
                }
            }
        }
        else if (localCoeffApproxType_ == "LMEA")
        {
            // Internal field
            forAll(theta_, cellI)
            {
                const scalar gasValue = workingGas[cellI];
                const scalar approxValue = localCoeffApprox[cellI];

                label i0, i1;
                scalar w;
                lowerUpperListValues(workingGasTableValues_, gasValue, i0, i1, w);

                const scalar PL0 = powerLossLMEATables_[i0](approxValue);
                const scalar PL1 = powerLossLMEATables_[i1](approxValue);

                theta_[cellI] = (i0 == i1 ? PL0 : (1.0 - w)*PL0 + w*PL1);
            }

            // Boundary fields
            forAll(theta_.boundaryFieldRef(), patchI)
            {
                fvPatchScalarField& thetaPatch = theta_.boundaryFieldRef()[patchI];
                const fvPatchScalarField& workingGasPatch = workingGas.boundaryField()[patchI];
                const fvPatchScalarField& localCoeffApproxPatch = localCoeffApprox.boundaryField()[patchI];

                forAll(thetaPatch, faceI)
                {
                    const scalar gasValue = workingGasPatch[faceI];
                    const scalar approxValue = localCoeffApproxPatch[faceI];

                    label i0, i1;
                    scalar w;
                    lowerUpperListValues(workingGasTableValues_, gasValue, i0, i1, w);

                    const scalar PL0 = powerLossLMEATables_[i0](approxValue);
                    const scalar PL1 = powerLossLMEATables_[i1](approxValue);

                    thetaPatch[faceI] = (i0 == i1 ? PL0 : (1.0 - w)*PL0 + w*PL1);
                }
            }
        }
    }

    // Correct boundary conditions (useful for mpi interfaces)
    theta_.correctBoundaryConditions();
    return theta_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
