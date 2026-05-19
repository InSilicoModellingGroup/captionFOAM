// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Mean electron energy equation
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
if (solveEleEnergyEqn)
{
    // Calculate coefficients
    label eID = composition.species().find("e");
    energy.calculateTransportCoeffs();
    surfaceScalarField EFlux = linearInterpolate(muE)*mesh.magSf()*fvc::snGrad(volt);
    volVectorField je = -mobilityCoeffSpecies[eID]*N[eID]*E - diffusionCoeffSpecies[eID]*fvc::grad(N[eID]);
    volScalarField SE = -constant::electromagnetic::e*(E & je) - N[eID]*energy.ThetaE();

    // Solve electron energy equation
    fvScalarMatrix nEEqn
    (
        fvm::ddt(nE) + fvm::div(EFlux, nE) - fvm::laplacian(DE, nE)
        ==
        SE
    );
    nEEqn.solve();

    // Mean electron energy
    em_J = nE/N[eID];
    em = J2eV*em_J;

    // Calculate electron temperature
    TEle = (2.0/3.0)*em;
    TEle_K = TEle*eV2K;
}
