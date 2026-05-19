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
    and interpolates them based on the local approximation variable and
    working-gas
    concentration.

SourceFiles
    diffusionCoeffBolsig.C

\*---------------------------------------------------------------------------*/

#ifndef diffusionCoeffBolsig_H
#define diffusionCoeffBolsig_H

#include "transportCoeffsBase.H"
#include "PtrList.H"
#include "interpolationTable.H"
#include "interpolation2DTable.H"
#include "tableBounds.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class diffusionCoeffBolsig Declaration
\*---------------------------------------------------------------------------*/

class diffusionCoeffBolsig : public transportCoeffsBase
{
protected:
    // Interpolation tables for mobility and diffusion coefficients
    interpolationTable<scalar> diffusivity1DTable_;
    interpolation2DTable<scalar> diffusivity2DTable_;
    scalarList workingGasTableValues_;
    PtrList<interpolationTable<scalar>> diffusivityLMEATables_;

    // Tabulation variables
    word localCoeffApproxType_;
    word backgndAtmosphereType_;
    word workingGasName_;
    const volScalarField* workingGasPtr_;
    word localCoeffApproxName_;
    const volScalarField* localCoeffApproxPtr_;

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
    //- Runtime type information
    TypeName("bolsig");

    // Constructor
    diffusionCoeffBolsig
    (
        const fvMesh& mesh,
        const dictionary& dict,
        volScalarField& transportCoeff
    );

    // Destructor
    // virtual ~transportCoeffsBolsig() {}
    ~diffusionCoeffBolsig() {}

    // Member functions
    //- Calculate diffusion coefficients
    void calcDiffusionCoeffs();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
