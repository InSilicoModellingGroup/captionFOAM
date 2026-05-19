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

#include "bolsigLMEARateCoeff.H"
#include "plasmaChemistryModel.H"
#include "volFields.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"
#include "Tuple2.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(bolsigLMEARateCoeff, 0);
addToRunTimeSelectionTable(reactionRateCoeffsBase, bolsigLMEARateCoeff, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
bolsigLMEARateCoeff::bolsigLMEARateCoeff
(
    const dictionary& dict, 
    const word& name, 
    plasmaChemistryModel& chemistry
)
:
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
    em_(chemistry.mesh().lookupObjectRef<volScalarField>("em")),
    HeList_(bolsigProperties_.lookup("He"))
{
    HePtr_ = nullptr;
    
    // Read mean energy list
    scalarList emList_
    (
        bolsigProperties_.lookup("meanEnergy")
    );

    List<List<scalar>> rateCoefficientsList_
    (
        bolsigProperties_.lookup("rateCoefficients")
    );

    // List dimensions
    nHe_ = HeList_.size();
    nEN_ = emList_.size()/nHe_;
    
    // Initialise pointer list size
    ratesTable_.setSize(rateCoefficientsList_[0].size());

    forAll(rateCoefficientsList_[0], r)
    {
        // Create lookup table for each reaction
        PtrList<interpolationTable<scalar>>* rTablePtr =
            new PtrList<interpolationTable<scalar>>(nHe_);

        forAll(HeList_, i)
        {
            List<Tuple2<scalar, scalar>> rRow(nEN_);

            forAll(rRow, j)
            {
                label index = i * nEN_ + j;
                rRow[j] = Tuple2<scalar, scalar>(emList_[index], rateCoefficientsList_[index][r]);
            }

            rTablePtr->set
            (
                i,
                new interpolationTable<scalar>
                (
                    rRow,
                    bounds::repeatableBounding(bounds::repeatableBounding::CLAMP),
                    fileName("ratesTable")
                )
            );
        }

        ratesTable_.set(r, rTablePtr);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void bolsigLMEARateCoeff::calculate(plasmaChemistryModel& chemistry, const label j) const
{
    // Only look up He once
    if (!HePtr_)
    {
        HePtr_ = &chemistry.mesh().lookupObjectRef<volScalarField>("He");
    }

    const volScalarField& He = *HePtr_;

    volScalarField& kj = chemistry.k()[j];

    // Internal field
    forAll(kj, cellI)
    {
        const scalar he  = He[cellI];
        const scalar eps = em_[cellI];

        // Find bracketing He indices and linear weight
        label i0, i1;
        scalar w;
        lowerUpperListValues(HeList_, he, i0, i1, w);

        // Evaluate rate on the two per-He tables and blend across He
        const scalar k0 = ratesTable_[j][i0](eps);
        const scalar k1 = ratesTable_[j][i1](eps);

        kj[cellI] = (i0 == i1 ? k0 : (1.0 - w)*k0 + w*k1);
    }

    // Boundary fields
    forAll(kj.boundaryFieldRef(), patchI)
    {
        fvPatchScalarField& kjPatch = kj.boundaryFieldRef()[patchI];
        const fvPatchScalarField& HePatch = He.boundaryField()[patchI];
        const fvPatchScalarField& emPatch = em_.boundaryField()[patchI];

        forAll(kjPatch, faceI)
        {
            const scalar he  = HePatch[faceI];
            const scalar eps = emPatch[faceI];

            // Find bracketing He indices and linear weight
            label i0, i1;
            scalar w;
            lowerUpperListValues(HeList_, he, i0, i1, w);

            // Evaluate rate on the two per-He tables and blend across He
            const scalar k0 = ratesTable_[j][i0](eps);
            const scalar k1 = ratesTable_[j][i1](eps);

            kjPatch[faceI] = (i0 == i1 ? k0 : (1.0 - w)*k0 + w*k1);
        }
    }

    // Correct boundary conditions (useful for mpi interfaces)
    kj.correctBoundaryConditions();
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

