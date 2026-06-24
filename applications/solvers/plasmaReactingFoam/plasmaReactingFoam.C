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

Application
    plasmaReactingFoam

Description
    Solver for plasma discharge with detailed chemistry.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "regionProperties.H"
#include "voltageHandlerBase.H"
#include "transportCoeffsBase.H"
#include "plasmaChemistryModel.H"
#include "electronEnergy.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for plasma discharge with detailed chemistry."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"

    // Calculate amplitude of external field
    #include "EExtAmpEqn.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;
    while( runTime.loop() )
    {
        // ++runTime;
        if (runTime.timeIndex() % printScreenResults == 0 || runTime.timeIndex() == 1)
        {
            Info<< "\nTime = " << runTime.timeName() << "  Time step = " << runTime.timeIndex() << nl << endl;
            solverPerformance::debug = 1;
        }

        // Calculate the induced electric field
        #include "EIndEqn.H"
        
        // Calculate the local electric field
        #include "calcE.H"

        // Control time step according to Co number
        #include "CourantNo.H"
        #include "setDeltaT.H" 

        // Equations for species continuity
        #include "NEqn.H"

        // Equation for mean electron energy density
        #include "nEEqn.H"

        if (runTime.timeIndex() % printScreenResults == 0 || runTime.timeIndex() == 1)
        {
            forAll(N, i)
            {
                Info<< "Max Species " << composition.species()[i] << " = " << gMax(N[i]) << " and Min = " << gMin(N[i]) << endl;
            }
            Info<< "Max EInd = " << mag(gMax(EInd)) << endl;
            Info<< "Max E = " << mag(gMax(E)) << endl;
            Info<< "Max E/N = " << gMax(mag(EN)()) << endl;
            Info<< "Max volt = " << gMax(mag(volt)()) << endl;
            Info<< "Max uth = " << gMax(uth) << endl;
            Info<< "Max nE = " << gMax(nE) << endl;
            Info<< "Max em = " << gMax(em) << endl;
            Info<< "Max TEle (eV) = " << gMax(TEle) << endl;
            Info<< "I_cond = " << Icond.value() << " A" << endl;

            // Info<< "At elecrode patch:" << endl;
            // label patchID = mesh.boundaryMesh().findPatchID("electrode");
            // Info<< "Max De patch = " << gMax(diffusionCoeffSpecies[0].boundaryFieldRef()[patchID]) << endl;
            // Info<< "Max ne = " << gMax(N[0].boundaryFieldRef()[patchID]) << endl;
            // Info<< "Max uth patch = " << gMax(uth.boundaryFieldRef()[patchID]) << endl;
            // Info<< "Max E = " << mag(gMax(E.boundaryFieldRef()[patchID])) << endl;
            // Info<< "Max em (eV) = " << gMax(em.boundaryFieldRef()[patchID]) << endl;
            // Info<< "Max TEle (eV) = " << gMax(TEle.boundaryFieldRef()[patchID]) << endl;
            // Info<< "Max TEle (K) = " << gMax(TEle_K.boundaryFieldRef()[patchID]) << endl;
        }

        solverPerformance::debug = 0;
        if (runTime.timeIndex() % printScreenResults == 0 || runTime.timeIndex() == 1)
        {
            runTime.printExecutionTime(Info);
        }
        runTime.write();
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
