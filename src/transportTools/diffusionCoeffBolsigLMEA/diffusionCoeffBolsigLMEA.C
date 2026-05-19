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

#include "diffusionCoeffBolsigLMEA.H" 
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "SortableList.H"
#include "Tuple2.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(diffusionCoeffBolsigLMEA, 0);
addToRunTimeSelectionTable(transportCoeffsBase, diffusionCoeffBolsigLMEA, diffusion);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
diffusionCoeffBolsigLMEA::diffusionCoeffBolsigLMEA
(
    const fvMesh& mesh, 
    const dictionary& dict, 
    volScalarField& transportCoeff
)
:
    transportCoeffsBase(mesh, dict, transportCoeff),
    He_(mesh_.lookupObjectRef<volScalarField>("He")),
    em_(mesh_.lookupObjectRef<volScalarField>("em"))
{
    // Read the Bolsig properties from the constant directory
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

    bolsigProperties_.lookup("He")  >> HeList_;

    scalarList emList_
    (
        bolsigProperties_.lookup("meanEnergy")
    );

    scalarList difList_
    (
        bolsigProperties_.lookup("electronDiffusivity")
    );

    nHe_ = HeList_.size();
    nEN_ = emList_.size()/nHe_;

    diffusivityTable_.setSize(nHe_);
    for (label i = 0; i < nHe_; ++i)
    {
        // Build (ε, μ) pairs for row i using GLOBAL index idx = i*nEN + j
        List<Tuple2<scalar,scalar>> pairs(nEN_);
        for (label j = 0; j < nEN_; ++j)
        {
            const label idx = i*nEN_ + j;
            // x = ε (mean energy), y = μ (or μ·N, depending on your data)
            pairs[j] = Tuple2<scalar,scalar>(emList_[idx], difList_[idx]);
        }

        // Build the 1-D table 
        diffusivityTable_.set
        (
            i,
            new interpolationTable<scalar>
            (
                pairs,
                bounds::repeatableBounding(bounds::repeatableBounding::CLAMP),
                fileName("diffusivityTable_row_"+name(i))
            )
        );
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void diffusionCoeffBolsigLMEA::calcDiffusionCoeffs()
{
    // Internal field
    forAll(transportCoeff_, cellI)
    {
        const scalar he  = He_[cellI];   // local He*
        const scalar eps = em_[cellI];   // local mean energy (eV)

        // Find bracketing He indices and linear weight
        label i0, i1; 
        scalar w;
        lowerUpperListValues(HeList_, he, i0, i1, w);

        // Evaluate diffusivity on the two per-He tables and blend across He
        const scalar D0 = diffusivityTable_[i0](eps);
        const scalar D1 = diffusivityTable_[i1](eps);

        transportCoeff_[cellI] = (i0 == i1 ? D0 : (1.0 - w)*D0 + w*D1);
    }

    // Boundary fields
    forAll(transportCoeff_.boundaryFieldRef(), patchI)
    {
        fvPatchScalarField& transportCoeffPatch = transportCoeff_.boundaryFieldRef()[patchI];
        const fvPatchScalarField& HePatch = He_.boundaryField()[patchI];
        const fvPatchScalarField& emPatch = em_.boundaryField()[patchI];

        forAll(transportCoeffPatch, faceI)
        {
            const scalar he  = HePatch[faceI]; 
            const scalar eps = emPatch[faceI];

            // Find bracketing indices in HeList_
            label i0, i1;
            scalar w;
            lowerUpperListValues(HeList_, he, i0, i1, w);

            // Evaluate interpolated D
            const scalar D0 = diffusivityTable_[i0](eps);
            const scalar D1 = diffusivityTable_[i1](eps);

            transportCoeffPatch[faceI] = (i0 == i1 ? D0 : (1.0 - w)*D0 + w*D1);
        }
    }

    // Correct boundary conditions (useful for mpi interfaces)
    transportCoeff_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



