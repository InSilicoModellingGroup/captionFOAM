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
    Foam::TownsendRateCoeff

Description
    The ionisation rate coefficient using the Townsend formula.

SourceFiles
    TownsendRateCoeff.C

\*---------------------------------------------------------------------------*/

#ifndef TownsendRateCoeff_H
#define TownsendRateCoeff_H

#include "reactionRateCoeffsBase.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class TownsendRateCoeff Declaration
\*---------------------------------------------------------------------------*/

class TownsendRateCoeff : public reactionRateCoeffsBase
{
protected:
    // Parameters (A[m2], B[Td}])
    scalar A_, B_; 
    dimensionedScalar Adim_;

    // Number density
    volScalarField& N_;
    
    // Name of the field
    word fieldName_;

    // Reduced field
    volScalarField& EN_;

    // Electric field
    volVectorField& E_;

    // Range of E/N [Td]
    List<scalar> rangeValues_;

    // Dimensioned scalar for E/N
    List<dimensionedScalar> rangeEN_;


public:
    TypeName("TownsendFormula");

    // Constructor
    TownsendRateCoeff(const dictionary& dict, const word& name, plasmaChemistryModel& chemistry);

    // Member functions
    // virtual void read(const dictionary& reactionsSubDict) override;
    virtual void calculate(plasmaChemistryModel& chemistry, const label reactionIndex) const override;
    
    // Virtual destructor
    virtual ~TownsendRateCoeff() {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
