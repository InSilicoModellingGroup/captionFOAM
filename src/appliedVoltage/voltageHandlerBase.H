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
    Foam::voltageHandlerBase

Description
    This is the base class for the calculation of the electric 
    potential that corresponds to the external electric field, in
    the current time step.

SourceFiles
    voltageHandlerBase.C

\*---------------------------------------------------------------------------*/

#ifndef voltageHandlerBase_H
#define voltageHandlerBase_H

#include "fvMesh.H"
#include "volFields.H"
#include <stdexcept>
#include "IOstreams.H" 
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class voltageHandlerBase Declaration
\*---------------------------------------------------------------------------*/

class voltageHandlerBase
{
protected:
    // fvMesh reference
    const fvMesh& mesh_;

    // Dictionary reference
    const dictionary& dict_;

    // Potential of the external electric field reference
    volScalarField& voltExt_;

    // Potential of the external electric field amplitude reference 
    const volScalarField& voltExtAmp_;

public:
    //- Runtime type information
    TypeName("voltageHandlerBase");

    // Declare run-time selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        voltageHandlerBase,
        voltageWaveform,
        (
            const fvMesh& mesh,
            const dictionary& dict,
            volScalarField& voltExt,
            const volScalarField& voltExtAmp
        ),
        (mesh, dict, voltExt, voltExtAmp)
    );

    // Selectors
    static autoPtr<voltageHandlerBase> NewVoltageWaveform
    (
        const fvMesh& mesh,
        const dictionary& dict,
        volScalarField& voltExt,
        const volScalarField& voltExtAmp
    );

    // Constructor
    voltageHandlerBase(const fvMesh& mesh, const dictionary& dict, volScalarField& voltExt, const volScalarField& voltExtAmp);

    // Destructor
    virtual ~voltageHandlerBase() {}

    // Member operators
    //- Equality operator
    voltageHandlerBase operator = (const voltageHandlerBase& voltageHandlerBase);
    // voltageHandler& operator=(const voltageHandler& other);

    // Member functions
    //- calculate voltage field
    virtual volScalarField calcVoltage() { return voltExt_; }
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
