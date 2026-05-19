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
    Foam::bolsigRateCoeff

Description
    This class provides the reads and defines the plasma chemistry.

SourceFiles
    bolsigRateCoeff.C

\*---------------------------------------------------------------------------*/

#ifndef bolsigRateCoeff_H
#define bolsigRateCoeff_H

#include "reactionRateCoeffsBase.H"
#include "PtrList.H"
#include "interpolationTable.H"
#include "interpolation2DTable.H"
#include "tableBounds.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class bolsigRateCoeff Declaration
\*---------------------------------------------------------------------------*/

class bolsigRateCoeff : public reactionRateCoeffsBase
{
protected:
    label reactionIndex_;
    IOdictionary bolsigProperties_;

    PtrList<interpolationTable<scalar>> rates1DTable_;
    PtrList<interpolation2DTable<scalar>> rates2DTable_;
    PtrList< PtrList<interpolationTable<scalar>> > ratesLMEATables_;
    scalarList workingGasTableValues_;

    // Tabulation variables
    word localCoeffApproxType_;
    word backgndAtmosphereType_;
    word workingGasName_;
    mutable const volScalarField* workingGasPtr_;
    word localCoeffApproxName_;
    mutable const volScalarField* localCoeffApproxPtr_;

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
    TypeName("bolsig");

    // Constructor
    bolsigRateCoeff
    (
        const dictionary& dict, 
        const word& name, 
        plasmaChemistryModel& chemistry
    );

    // Member functions
    // virtual void read(const dictionary& reactionsSubDict) override;
    virtual void calculate(plasmaChemistryModel& chemistry, const label reactionIndex) const override;
    
    // Virtual destructor
    virtual ~bolsigRateCoeff() {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
