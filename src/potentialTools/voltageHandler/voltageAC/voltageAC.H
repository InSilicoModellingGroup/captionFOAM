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
    Foam::voltageAC

Description
    Calculate the electric potential that corresponds to an AC external
    electric field.

SourceFiles
    voltageAC.C

\*---------------------------------------------------------------------------*/

#ifndef voltageAC_H
#define voltageAC_H

#include "voltageHandlerBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class voltageAC Declaration
\*---------------------------------------------------------------------------*/

class voltageAC : public voltageHandlerBase
{
public:
    TypeName("AC");

    // Constructor
    voltageAC
    (
        const fvMesh& mesh, 
        const dictionary& dict, 
        volScalarField& voltExt, 
        const volScalarField& voltExtAmp
    );

    // destructor
    ~voltageAC() {}

    // member functions
    //- calculate voltage field
    volScalarField calcVoltage();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
