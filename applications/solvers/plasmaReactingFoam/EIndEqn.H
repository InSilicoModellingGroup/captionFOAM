// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Equations for the time varying Induced Electric Field 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
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
        fvm::laplacian(voltInd) + rhoq/constant::electromagnetic::epsilon0
    );
    voltEqnIter.push_back(solvPerfVoltEqn.nIterations());

    // Solve solid regions
    forAll(solidRegions, i)
    {
        #include "setRegionSolidFields.H"
        Foam::solverPerformance solvPerfVoltEqn = solve 
        (
            fvm::laplacian(voltInd)
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
