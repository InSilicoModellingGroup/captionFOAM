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
    Foam::bolsigLMEARateCoeff

Description
    This class provides the reads and defines the plasma chemistry.

SourceFiles
    bolsigLMEARateCoeff.C

\*---------------------------------------------------------------------------*/

#ifndef bolsigLMEARateCoeff_H
#define bolsigLMEARateCoeff_H

#include "reactionRateCoeffsBase.H"
#include "interpolationTable.H"
#include "tableBounds.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class bolsigLMEARateCoeff Declaration
\*---------------------------------------------------------------------------*/

class bolsigLMEARateCoeff : public reactionRateCoeffsBase
{
protected:
    label reactionIndex_;
    IOdictionary bolsigProperties_;

    // For each reaction: list of 1D tables (meanEnergy) for each He row
    PtrList< PtrList<interpolationTable<scalar>> > ratesTable_;
    // PtrList<interpolation2DTable<scalar>> ratesTable_;
    // static interpolation2DTable<scalar>* bolsigTablePtr_;
    // static bool initialized_;

    // Tabulation variables
    const volScalarField& em_;
    mutable const volScalarField* HePtr_;

    List<scalar> HeList_;

    // List dimensions
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
    TypeName("bolsigLMEA");

    // Constructor
    bolsigLMEARateCoeff
    (
        const dictionary& dict, 
        const word& name, 
        plasmaChemistryModel& chemistry
    );

    // Member functions
    // virtual void read(const dictionary& reactionsSubDict) override;
    virtual void calculate(plasmaChemistryModel& chemistry, const label reactionIndex) const override;
    
    // Virtual destructor
    virtual ~bolsigLMEARateCoeff() {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
