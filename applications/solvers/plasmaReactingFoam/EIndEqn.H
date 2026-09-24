// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Calculate for the time varying External Electric Field
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Calculate the potential of the external electric field
voltExt = voltHandler->calcVoltage();
forAll(solidRegions, i)
{
    #include "setRegionSolidFields.H"
    voltExtSolid[i] = voltHandlerSolid[i].calcVoltage();
}

// Calculate the external electric field
EExt = -fvc::grad(voltExt);
forAll(solidRegions, i)
{
    #include "setRegionSolidFields.H"
    EExt = -fvc::grad(voltExt);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Equations for the time varying Induced Electric Field 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Coefficients for semi implicit integration
volScalarField coeffSemiImplicit
(
    "coeffSemiImplicit",
    ones + runTime.deltaT()*sigmaCond/constant::electromagnetic::epsilon0
);

volScalarField sourceSemiImplicit
(
    "sourceSemiImplicit",
    (rhoq + runTime.deltaT()*divDiffFlux + runTime.deltaT()*fvc::laplacian(sigmaCond, voltExt))/constant::electromagnetic::epsilon0
);

// Reset counters
regionLoopCounter = 0;
if (runTime.timeIndex() % printScreenResults == 0 || runTime.timeIndex() == 1)
{
    solverPerformance::debug = 1;
    Info<< "Current Iteration = " << 1 << endl;
}

voltEqnIter.clear();
voltEqnIter.push_back(1);

// Loop through different regions to solve the Poisson equations for voltInd
while (*std::max_element(voltEqnIter.begin(), voltEqnIter.end()) && regionLoopCounter<maxRegionLoops)
{
    voltEqnIter.clear();

    // Solve gas region
    Foam::solverPerformance solvPerfVoltEqn = solve
    (
        fvm::laplacian(coeffSemiImplicit, voltInd) + sourceSemiImplicit
    );
    voltEqnIter.push_back(solvPerfVoltEqn.nIterations());

    // Solve solid regions
    forAll(solidRegions, i)
    {
        #include "setRegionSolidFields.H"
        Foam::solverPerformance solvPerfVoltEqn = solve 
        (
            fvm::laplacian(epsilon, voltInd)
        );
        voltEqnIter.push_back(solvPerfVoltEqn.nIterations());
    }

    // Print performance at custom iteration intervals
    regionLoopCounter++;
    solverPerformance::debug = 0;
    int printPerformance = 500;
    if (regionLoopCounter % printPerformance == 0 && (runTime.timeIndex() % printScreenResults == 0 || runTime.timeIndex() == 1))
    {
        Info<< "Current Iteration = " << regionLoopCounter << endl;
        solverPerformance::debug = 1;
    }
}

if (runTime.timeIndex() % printScreenResults == 0 || runTime.timeIndex() == 1)
{
    Info<< "Region Inner Loops = " << regionLoopCounter << endl;
    solverPerformance::debug = 1;
}

// Calculate the induced electric field
EInd = -fvc::grad(voltInd);
forAll(solidRegions, i)
{
    #include "setRegionSolidFields.H"
    EInd = -fvc::grad(voltInd);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Calculate local electric field
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
volt = voltExt + voltInd;
forAll(solidRegions, i)
{
    #include "setRegionSolidFields.H"
    volt = voltExt + voltInd;
}

E = -fvc::grad(volt);
forAll(solidRegions, i)
{
    #include "setRegionSolidFields.H"
    E = -fvc::grad(volt);
}

// Calculate reduced electric field in Td
EN = 1.0e21*mag(E/NGas);
