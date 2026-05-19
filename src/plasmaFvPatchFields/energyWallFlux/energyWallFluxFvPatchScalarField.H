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
    Foam::energyWallFluxFvPatchScalarField

Description
    Mixed wall boundary condition for charged-species number density.
    The current simple model applies wall absorption from drift-to-wall
    (and optional thermal) flux without secondary emission.

SourceFiles
    energyWallFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef energyWallFluxFvPatchScalarField_H
#define energyWallFluxFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "volFields.H"
#include "fundamentalConstants.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class energyWallFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class energyWallFluxFvPatchScalarField
:
    public mixedFvPatchScalarField
{
protected:
    //- Name of temperature field
    word temperatureFieldName_;

    //- Name of electric field
    word electricFieldName_;

    // Species-specific constants
    word speciesName_;
    scalar mass_;
    int chargeNumber_;
    scalar uthPrefactor_;
    bool initializedScalars_;

public:

    //- Runtime type information
    TypeName("energyWallFlux");

    // Constructors
    //- Construct from patch and internal field
    energyWallFluxFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    energyWallFluxFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping given energyWallFluxFvPatchScalarField onto a new patch
    energyWallFluxFvPatchScalarField
    (
        const energyWallFluxFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Construct as copy
    energyWallFluxFvPatchScalarField
    (
        const energyWallFluxFvPatchScalarField&
    );

    //- Construct as copy setting internal field reference
    energyWallFluxFvPatchScalarField
    (
        const energyWallFluxFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    // Member functions
    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();

    //- Write
    virtual void write(Ostream&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
