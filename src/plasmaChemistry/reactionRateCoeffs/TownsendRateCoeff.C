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

#include "TownsendRateCoeff.H"
#include "plasmaChemistryModel.H"
#include "volFields.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"
#include "Swap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(TownsendRateCoeff, 0);
    addToRunTimeSelectionTable(reactionRateCoeffsBase, TownsendRateCoeff, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
TownsendRateCoeff::TownsendRateCoeff(const dictionary& dict, const word& name, plasmaChemistryModel& chemistry)
:
    // reactionRateCoeffsBase(dict, name, chemistry),
    A_(dict.lookupOrDefault<scalar>("A", 1.0)),
    B_(dict.lookupOrDefault<scalar>("B", 0.0)),
    Adim_("Adim", dimensionSet(0,2,0,0,0,0,0), A_),
    N_(chemistry.mesh().lookupObjectRef<volScalarField>("NGas")),
    fieldName_(dict.lookupOrDefault<word>("field", "EN")),
    EN_(chemistry.mesh().lookupObjectRef<volScalarField>(fieldName_)),
    E_(chemistry.mesh().lookupObjectRef<volVectorField>("E")),
    rangeValues_(dict.lookup("range"))
{
    
    // Throw error if range is not a list of two scalars
    if (rangeValues_.size() != 2)
    {
        FatalIOErrorInFunction(dict)
            << "Entry 'range' must be a list of two scalars: (lower upper). "
            << "Got " << rangeValues_ << exit(FatalIOError);
    }

    // Swap if given reversed by user
    if (rangeValues_[1] < rangeValues_[0])
    {
        Info<< "TownsendRateCoeff: swapping reversed range (" << rangeValues_[0]
            << ", " << rangeValues_[1] << ")\n";
        Swap(rangeValues_[0], rangeValues_[1]);
    }

    if (mag(rangeValues_[1] - rangeValues_[0]) < SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "range upper and lower limits are identical (" << rangeValues_[0]
            << ")." << exit(FatalIOError);
    }

    // // Populate dimensioned scalar range
    // // Build the list with the right dimensions from the start
    // rangeEN_ = List<dimensionedScalar>
    // (
    //     2,
    //     dimensionedScalar("ENtmp", EN_.dimensions(), 0.0)  // sets correct dims
    // );

    // // Now only change the numeric values (dimensions remain intact)
    // rangeEN_[0].value() = rangeValues_[0];  // ENmin
    // rangeEN_[1].value() = rangeValues_[1];  // ENmax

    // Give to the A coefficient units
    // Adim_ = dimensionedScalar("A", dimensionSet(0,2,0,0,0,0,0), A_);

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void TownsendRateCoeff::calculate(plasmaChemistryModel& chemistry, const label j) const
{
    volScalarField& kj = chemistry.k()[j];

    const dimensionedScalar ENmin("ENmin", EN_.dimensions(), rangeValues_[0]);
    const dimensionedScalar ENmax("ENmax", EN_.dimensions(), rangeValues_[1]);

    // Calculate clamped reduced electric field
    const tmp<volScalarField> clampedEN = min(max(EN_, ENmin), ENmax);

    // Direct lookup of the mobility coefficient (guaranteed to exist)
    const volScalarField& mu = chemistry.mesh().lookupObject<volScalarField>("mobilityCoeffSpecies_e");

    // ki = mue*E*(a/N)
    // kj = dim(kj)*mu/dim(mu);
    kj = mu * mag(E_) *  (Adim_*exp(B_/(clampedEN()/dim(clampedEN()))));

    // Info<< "TOWNSEND " << mu << endl;
}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

