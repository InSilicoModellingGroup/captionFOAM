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
    Foam::electronEnergy

Description
    This is a class that provides the calculation of the transport
    coefficients (mobility and diffusion) for the ions.

SourceFiles
    electronEnergy.C

\*---------------------------------------------------------------------------*/

#ifndef electronEnergy_H
#define electronEnergy_H

#include "volFields.H"
#include "constants.H"
#include "IOstreams.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "interpolationTable.H"
#include "interpolation2DTable.H"
#include "tableBounds.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class electronEnergy Declaration
\*---------------------------------------------------------------------------*/

class electronEnergy
{
protected:
    // fvMesh reference
    const fvMesh& mesh_;

    // Interpolation tables for energy transport coefficients
    interpolationTable<scalar> energyMobility1DTable_;
    interpolationTable<scalar> energyDiffusivity1DTable_;
    interpolationTable<scalar> powerLoss1DTable_;
    interpolation2DTable<scalar> energyMobility2DTable_;
    interpolation2DTable<scalar> energyDiffusivity2DTable_;
    interpolation2DTable<scalar> powerLoss2DTable_;
    scalarList workingGasTableValues_;
    PtrList<interpolationTable<scalar>> energyMobilityLMEATables_;
    PtrList<interpolationTable<scalar>> energyDiffusivityLMEATables_;
    PtrList<interpolationTable<scalar>> powerLossLMEATables_;

    // Tabulation variables
    const volScalarField& EN_;
    const volScalarField& em_;
    word localCoeffApproxType_;
    word backgndAtmosphereType_;
    word workingGasName_;
    const volScalarField* workingGasPtr_;
    word localCoeffApproxName_;
    const volScalarField* localCoeffApproxPtr_;

    // Fields from the mesh
    volScalarField& muE_;
    volScalarField& DE_;
    volScalarField theta_;

    // Flag indicating constant coefficient throughout the simulation
    bool isConstant_;

    // Helper: find immediate lower/upper values and interpolation weight
    static inline void lowerUpperListValues
    (
        const scalarList& values,
        const scalar value,
        label& i0,
        label& i1,
        scalar& w
    )
    {
        const label n = values.size();

        if (value <= values[0])
        {
            i0 = i1 = 0;
            w = 0.0;
            return;
        }

        if (value >= values[n - 1])
        {
            i0 = i1 = n - 1;
            w = 0.0;
            return;
        }

        i0 = 0;
        while ((i0 < n - 1) && (values[i0 + 1] <= value))
        {
            ++i0;
        }

        i1 = i0 + 1;

        const scalar d = values[i1] - values[i0];
        w = (d > SMALL) ? (value - values[i0])/d : 0.0;
    }

public:

    // Constructor
    electronEnergy(const fvMesh& mesh);

    // Destructor
    virtual ~electronEnergy() {}

    // Member functions
    void calculateTransportCoeffs();

    const volScalarField& ThetaE();

    //- inline - return true to ignore calculation of coefficients
    bool isConstant() const 
    {
        return isConstant_;
    }
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
