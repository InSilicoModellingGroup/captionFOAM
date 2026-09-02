// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Define transport coefficient related objects
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
Info<< "Reading species properties\n" << endl;
IOdictionary speciesProperties
(
    IOobject
    (
        "speciesProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

IOdictionary molarMass
(
    IOobject
    (
        "molarMass",
        runTime.constant(),   // or runTime.timeName()
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    )
);

IOdictionary chargeNumber
(
    IOobject
    (
        "chargeNumber",
        runTime.constant(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    )
);

// Define pointer lists
PtrList<volScalarField> mobilityCoeffSpecies(composition.species().size());
PtrList<volScalarField> diffusionCoeffSpecies(composition.species().size());
PtrList<transportCoeffsBase> mobilityCoeffHandler(composition.species().size());
PtrList<transportCoeffsBase> diffusionCoeffHandler(composition.species().size());

// Loop over species
forAll(composition.species(), i)
{
    const word& specie = composition.species()[i];
    const dictionary& specieDict(speciesProperties.subDict(specie));

    chargeNumber.add
    (
        specie,
        specieDict.get<int>("chargeNumber")
    );

    molarMass.add
    (
        specie,
        1.0e-3 * specieDict.get<scalar>("molarMass")
    );

    // Info<< "Species: " << specie;

    // --- Mobility
    const dictionary& mobDict = specieDict.subDict("mobility");

    mobilityCoeffSpecies.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "mobilityCoeffSpecies_" + specie,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimensionSet(-1, 0, 2, 0, 0, 1, 0), scalar(0))
        )
    );

    autoPtr<transportCoeffsBase> model = transportCoeffsBase::NewMobility
    (
        mesh, 
        mobDict, 
        mobilityCoeffSpecies[i]
    );

    mobilityCoeffHandler.set
    (
        i, 
        model.release()
    );  

    // Create diffusion coefficient
    const dictionary& diffDict = specieDict.subDict("diffusion");

    diffusionCoeffSpecies.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "diffusionCoeffSpecies_" + specie,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimViscosity, scalar(0))
        )
    );

    autoPtr<transportCoeffsBase> diffModel = transportCoeffsBase::NewDiffusion
    (
        mesh, 
        diffDict, 
        diffusionCoeffSpecies[i]
    );

    diffusionCoeffHandler.set
    (
        i, 
        diffModel.release()
    );
}
