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

#include "mobilityCoeffIonHePlusHe.H" 
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(mobilityCoeffIonHePlusHe, 0);
addToRunTimeSelectionTable(transportCoeffsBase, mobilityCoeffIonHePlusHe, mobility);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
mobilityCoeffIonHePlusHe::mobilityCoeffIonHePlusHe(const fvMesh& mesh, const dictionary& dict, volScalarField& transportCoeff):
    transportCoeffsBase(mesh, dict, transportCoeff),
    EN_(mesh_.lookupObjectRef<volScalarField>("EN"))
{
    // Give values to the polynomial coefficients
    polyCoeff_.setSize(5);
    polyCoeff_[0] = 1.15875342e-03;   // a0
    polyCoeff_[1] =-4.01605282e-06;   // a1
    polyCoeff_[2] = 1.02480441e-08;   // a2
    polyCoeff_[3] =-1.34329667e-11;   // a3
    polyCoeff_[4] = 6.79261506e-15;   // a4

    // Calculate upon construction
    calcMobilityCoeffs();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void mobilityCoeffIonHePlusHe::calcMobilityCoeffs()
{
    // Calculate coeff
    const volScalarField& ENnonDim_ = nonDim(EN_);
    transportCoeff_ = dim(transportCoeff_)*((((polyCoeff_[4]*ENnonDim_ + polyCoeff_[3])*ENnonDim_ + polyCoeff_[2])*ENnonDim_ + polyCoeff_[1])*ENnonDim_ + polyCoeff_[0]);

    // Bound control to avoid unphysical values
    transportCoeff_.clamp_max(0.00113120515448821);
    transportCoeff_.clamp_min(0.000392077903060477);

    // Correct boundary conditions (useful for mpi interfaces)
    transportCoeff_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

