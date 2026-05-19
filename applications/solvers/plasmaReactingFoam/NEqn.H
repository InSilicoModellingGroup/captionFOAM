{
    volScalarField Nt(0.0*N[0]);
    surfaceScalarField jt(0.0*j[0]);
    volVectorField Gammat(0.0*Gamma[0]);
    rhoq = 0.0*rhoq;
    sumPosFlux = 0.0*sumPosFlux;
    
    // Update only the non constant reaction rate coefficients at selected time steps
    // Always update at first time step though (if we want BOLSIG coeffs from the laplacian solution only)
    if ((customUpdatePlasmaCoeffs && runTime.timeIndex() % plasmaCoeffsUpdateInterval == 0) || runTime.timeIndex() == 1)
    {
        composition.updateReactionRateCoeffs();
    }

    // Calculate species fluxes and solve for species number densities
    forAll(N, i)
    {
        if (i != inertIndex)
        {
            if (runTime.timeIndex() % plasmaCoeffsUpdateInterval == 0)
            {
                // Update mobility and diffusion coefficients
                if (!mobilityCoeffHandler[i].isConstant())
                {
                    mobilityCoeffHandler[i].calcMobilityCoeffs();
                }

                if (!diffusionCoeffHandler[i].isConstant())
                {
                    diffusionCoeffHandler[i].calcDiffusionCoeffs();
                }
            }

            volScalarField& Ni = N[i];
            volScalarField& Di = diffusionCoeffSpecies[i];
            volScalarField& mui = mobilityCoeffSpecies[i];
            surfaceScalarField& NiFlux = NFlux[i];
            surfaceScalarField& ji = j[i];
            volVectorField& Gammai = Gamma[i];
            const scalar& zi = chargeNumber.get<scalar>(composition.species()[i]);

            NiFlux = -zi*linearInterpolate(mobilityCoeffSpecies[i])*mesh.magSf()*fvc::snGrad(volt);

            fvScalarMatrix NiEqn
            (
                fvm::ddt(Ni)
                + fvm::div(NiFlux, Ni)
                - fvm::laplacian(Di, Ni)
              ==
                composition.R(i,N)
            );   
            
            NiEqn.relax();
            
            NiEqn.solve("Ni");
            
            Ni.clamp_min(1.0e9);
            Nt += Ni;
            Ni.correctBoundaryConditions();

            // Update species integrated face particle flux with latest Ni values.
            ji = zi*constant::electromagnetic::e*(-zi*linearInterpolate(mui*Ni)*fvc::snGrad(volt) - linearInterpolate(Di)*fvc::snGrad(Ni));
            Gammai = -zi*mui*Ni*fvc::grad(volt) - Di*fvc::grad(Ni);
            jt += ji;
            Gammat += Gammai;
            rhoq += zi*constant::electromagnetic::e*N[i];

            if (zi > 0)
            {
                sumPosFlux += gammaSE*ji/(zi*constant::electromagnetic::e);
            }
        }

        N[composition.species().find("e")].clamp_max(1.0e25);
        // if (composition.species().find("He_plus") >= 0)
        // {
        //     N[composition.species().find("He_plus")].clamp_max(1.0e25);
        // }
    }

    N[inertIndex] = NInert - Nt;
    N[inertIndex].clamp_min(0);
    N[composition.species().find("e")].clamp_min(1.0e9);

    // Calculate surface charge on air-dielecrtic interface
    forAll(mesh.boundary(), patchI)
    {
        const fvPatch& p = mesh.boundary()[patchI];

        if (p.type() == "mappedWall")
        {
            sigma.boundaryFieldRef()[patchI] = sigma.oldTime().boundaryField()[patchI] + runTime.deltaTValue()*jt.boundaryField()[patchI];
        }
    }

    // Calculate electric current
    Icond = constant::electromagnetic::e*fvc::domainIntegrate(Gammat & EExt)/1.0e-3;
}



    // // Sum all charged species contributions into wall-normal charge flux density.
    // sigmaFlux = dimensionedScalar("zero", sigmaFlux.dimensions(), scalar(0));
    // forAll(N, i)
    // {
    //     if (i != inertIndex)
    //     {
    //         const scalar& chargeNumberValue = chargeNumber.get<scalar>(composition.species()[i]);
    //         sigmaFlux += chargeNumberValue*constant::electromagnetic::e*(j[i]/mesh.magSf());
    //     }
    // }

    // // Integrate surface charge only on wall patches: sigma^{n+1} = sigma^n + dt*J_n.
    // sigma.primitiveFieldRef() = 0.0;

    // scalar maxWallSigma = 0.0;
    // forAll(mesh.boundary(), patchI)
    // {
    //     const fvPatch& p = mesh.boundary()[patchI];

    //     if (p.type() == "wall" || p.type() == "mappedWall")
    //     {
    //         sigma.boundaryFieldRef()[patchI] +=
    //             runTime.deltaTValue()*sigmaFlux.boundaryField()[patchI];

    //         maxWallSigma =
    //             Foam::max
    //             (
    //                 maxWallSigma,
    //                 max(mag(sigma.boundaryField()[patchI]))
    //             );
    //     }
    // }
    // maxWallSigma = returnReduce(maxWallSigma, maxOp<scalar>());

    // if (runTime.timeIndex() % printScreenResults == 0 || runTime.timeIndex() == 1)
    // {
    //     Info<< "Max |sigma| on walls [C/m2] = " << maxWallSigma << endl;
    // }
// }
