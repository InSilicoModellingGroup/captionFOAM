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

Class
    Foam::diffusionCoeffBolsig

Description
    This class reads the mobility and diffusion coefficients from a table
    and interpolates them based on the electric field and He concentration.

SourceFiles
    diffusionCoeffBolsig.C

\*---------------------------------------------------------------------------*/

#ifndef diffusionCoeffBolsigLMEA_H
#define diffusionCoeffBolsigLMEA_H

#include "transportCoeffsBase.H"
#include "interpolationTable.H"
#include "tableBounds.H"
#include "PtrList.H"

namespace Foam
{

class diffusionCoeffBolsigLMEA : public transportCoeffsBase
{
protected:
    // Interpolation tables for mobility and diffusion coefficients
    List<scalar> HeList_;
    PtrList<interpolationTable<scalar>> diffusivityTable_;

    // Fields from the mesh
    const volScalarField& He_;
    const volScalarField& em_;

    // Dimensions
    label nHe_;
    label nEN_;

    // Helper: find immediate upper and lower values of He
    static inline void lowerUpperListValues
    (
        const scalarList& HeList, 
        const scalar HeValue,
        label& i0, 
        label& i1, 
        scalar& w
    )
    {
        const label n = HeList.size();

        if (HeValue <= HeList[0])      
        {
            i0 = i1 = 0;       
            w = 0.0; 
            return; 
        }

        if (HeValue >= HeList[n-1])    
        { 
            i0 = i1 = n-1;     
            w = 0.0; 
            return; 
        }

        i0 = findLower(HeList, HeValue); // index of immediate lower value
        i1 = i0 + 1; // index of immediate upper value

        const scalar d = HeList[i1] - HeList[i0];
        w = (d > SMALL) ? (HeValue - HeList[i0]) / d : 0.0; // if d is larger than small, calculate the weight, else set w to 0
    }

public:
    //- Runtime type information
    TypeName("bolsigLMEA");

    // Constructor
    diffusionCoeffBolsigLMEA
    (
        const fvMesh& mesh,
        const dictionary& dict,
        volScalarField& transportCoeff
    );

    // Destructor
    // virtual ~transportCoeffsBolsig() {}
    ~diffusionCoeffBolsigLMEA() {}

    // Member functions
    //- Calculate diffusion coefficients
    void calcDiffusionCoeffs();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
