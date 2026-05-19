// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Mean electron energy equation
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
{
    label eID = composition.species().find("e");
    // energy.calculateTransportCoeffs();
    surfaceScalarField EFlux = (5/3)*linearInterpolate(mobilityCoeffSpecies[eID])*mesh.magSf()*fvc::snGrad(volt);
    volVectorField je = -mobilityCoeffSpecies[eID]*N[eID]*E - diffusionCoeffSpecies[eID]*fvc::grad(N[eID]);
    volScalarField SE = -constant::electromagnetic::e*(E & je) - N[eID]*energy.ThetaE();
    
    // Solve electron energy equation
    fvScalarMatrix nEEqn
    (
        fvm::ddt(nE) + fvm::div(EFlux, nE) - fvm::laplacian((5/3)*diffusionCoeffSpecies[eID], nE) 
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
