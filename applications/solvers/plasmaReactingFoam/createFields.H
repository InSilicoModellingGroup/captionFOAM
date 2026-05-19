// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Create general fields
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
Info<< "Reading physicalProperties\n" << endl;
IOdictionary physicalProperties
(
    IOobject
    (
        "physicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Access background gas properties sub-dictionary
const dictionary& atmosphericConditions = physicalProperties.subDict("atmosphericConditions");

// Gas pressure value
const dimensionedScalar pGasValue
(
    "pGasValue",
    dimensionSet(1, -1, -2, 0, 0, 0, 0),
    atmosphericConditions
);

// Electron temperature value
const dimensionedScalar TEleInitValue
(
    "TEleInitValue",
    dimensionSet(1,2,-2,0,0,0,0),
    atmosphericConditions
);

// Gas temperature value
const dimensionedScalar TGasValue
(
    "TGasValue",
    dimensionSet(0, 0, 0, 1, 0, 0, 0),
    atmosphericConditions
);

// Gas density value
const dimensionedScalar rhoGasValue
(
    "rhoGasValue",
    dimDensity,
    atmosphericConditions
);

// Secondary electron emission coefficient
const dimensionedScalar gammaSE
(
    "gammaSE",
    dimless,
    physicalProperties
);

// Convert eV to Kelvin
const dimensionedScalar eV2K
(   
    "eV2K",
    dimensionSet(-1, -2, 2, 1, 0, 0, 0),
    scalar(11606)
);

// Convert Joule to eV
const dimensionedScalar J2eV 
(   
    "J2eV",
    dimless,
    scalar(6.241509e+18)
);

// Convert eV to Joule
const dimensionedScalar eV2J
(   
    "eV2J",
    dimless,
    scalar(1.602176634e-19)
);

// Gas number density field
volScalarField NGas
(
    IOobject
    (
        "NGas",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    pGasValue/(constant::physicoChemical::k * TGasValue)
);

// Ambient pressure field
volScalarField pGas
(
    IOobject
    (
        "pGas",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    pGasValue
);

// Electron temperature field (eV)
volScalarField TEle
(
    IOobject
    (
        "TEle",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    TEleInitValue
);

// Electron temperature field (Kelvin)
volScalarField TEle_K
(
    IOobject
    (
        "TEle_K",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    TEleInitValue*eV2K
);

// Gas temperature field
volScalarField TGas
(
    IOobject
    (
        "TGas",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    TGasValue
);

volScalarField ones
(
    IOobject
    (
        "ones",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("temp",dimensionSet(0,0,0,0,0,0,0),scalar(1))
);

// Electrical current - conduction
dimensionedScalar Icond
(
    "Icond", 
    dimCurrent, 
    scalar(0.0)
);

// Temporary field - voltage amplitude
dimensionedScalar voltAmpTemp
(
    "voltAmpTemp",
    dimensionSet(1,2,-3,0,0,-1,0),
    scalar(1.0e3)
);

// Electron thermal velocity field
volScalarField uth
(
    IOobject
    (
        "uth",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    sqrt(8.0*constant::physicoChemical::k*TEle/(constant::mathematical::pi*constant::atomic::me))
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Create electrical fields
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "createSolidFields.H"

volScalarField epsilon
(
    IOobject
    (
        "epsilon",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("temp",dimless,scalar(1))
);

Info<< "Reading field voltExtAmp\n" << endl;
volScalarField voltExtAmp
(
    IOobject
    (
        "voltExtAmp",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field voltExt\n" << endl;
volScalarField voltExt
(
    IOobject
    (
        "voltExt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("temp", dimensionSet(1,2,-3,0,0,-1,0), scalar(0))
);

Info<< "Reading field voltInd\n" << endl;
volScalarField voltInd
(
    IOobject
    (
        "voltInd",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field volt\n" << endl;
volScalarField volt
(
    IOobject
    (
        "volt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("temp", dimensionSet(1,2,-3,0,0,-1,0), scalar(0))
);

Info<< "Calculating field EExtAmp\n" << endl;
volVectorField EExtAmp
(
    IOobject
    (
        "EExtAmp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltExtAmp)
);

Info<< "Calculating field EExt\n" << endl;
volVectorField EExt
(
    IOobject
    (
        "EExt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltExt)
);

Info<< "Calculating field EInd\n" << endl;
volVectorField EInd
(
    IOobject
    (
        "EInd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(voltInd)
);

Info<< "Calculating field E\n" << endl;
volVectorField E
(
    IOobject
    (
        "E",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(volt)
);

Info<< "Calculating field EN\n" << endl;
volScalarField EN
(
    IOobject
    (
        "EN",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    1.0e21*mag(E/NGas)
);

dimensionedScalar maxVal("maxVal", EInd.dimensions(), 5.0e6);

// Create voltage waveform handler
autoPtr<voltageHandlerBase> voltHandler;
voltHandler = voltageHandlerBase::NewVoltageWaveform
(
    mesh, 
    physicalProperties, 
    voltExt, 
    voltExtAmp
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Create electron mean energy fields
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Check if the electron energy equation will be solved
const dictionary& localCoeffApprox = physicalProperties.subDict("plasmaCoeffApprox").subDict("localCoeffApprox");
bool solveEleEnergyEqn = localCoeffApprox.lookupOrDefault("solveEleEnergyEqn", true);
if (localCoeffApprox.get<word>("type") == "LMEA")
{
    solveEleEnergyEqn = true;
}

Info<< "Reading field nE\n" << endl;
volScalarField nE
(
    IOobject
    (
        "nE",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Mean electron energy field (eV)
Info<< "Creating field em\n" << endl;
volScalarField em
(
    IOobject
    (
        "em",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("temp",dimensionSet(1,2,-2,0,0,0,0),scalar(0.04))
);

// Mean electron energy field (Joule)
volScalarField em_J
(
    IOobject
    (
        "em_J",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("temp",dimensionSet(1,2,-2,0,0,0,0),scalar(6.408707e-21))
);

Info<< "Creating field muE\n" << endl;
volScalarField muE
(
    IOobject
    (
        "muE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("temp", dimensionSet(-1,0,2,0,0,1,0), scalar(1))
);

Info<< "Creating field DE\n" << endl;
volScalarField DE
(
    IOobject
    (
        "DE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("temp", dimensionSet(0,2,-1,0,0,0,0), scalar(1))
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Create species
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Instantiate plasma chemistry model
plasmaChemistryModel composition(mesh);
int inertIndex = composition.inertIndex();

PtrList<volScalarField> N(composition.species().size());

// Create a temporary fallback field
tmp<volScalarField> tNdefault;

// Define name for default field
word NdefaultName("Ndefault");

// Create an IOobject to check if "nDefault" exists in the case
IOobject NdefaultIO
(
    NdefaultName,
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::NO_WRITE
);

// Check if nDefault exists
if (NdefaultIO.typeHeaderOk<volScalarField>(true))
{    
    tNdefault = new volScalarField
    (
        NdefaultIO, 
        mesh
    );
}
else
{
    // If not found, create a default empty field
    tNdefault = new volScalarField
    (
        IOobject
        (
            NdefaultName, 
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("temp", dimensionSet(0, -3, 0, 0, 0, 0, 0), scalar(0))
    );
}

forAll(composition.species(), i)
{
    IOobject header
    (
        composition.species()[i],
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    
    if (header.typeHeaderOk<volScalarField>(true))
    {
        N.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    composition.species()[i],
                    runTime.timeName(),
                    mesh, 
                    IOobject::MUST_READ, 
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }
    else
    {
        N.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    composition.species()[i], 
                    runTime.timeName(), 
                    mesh, 
                    IOobject::NO_READ, 
                    IOobject::AUTO_WRITE
                ),
                tNdefault()
            )
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Create background gas mixing and other fields
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Define and initialize species transport coefficients
#include "createTransportCoeffs.H"

// Initalise all reaction rate coefficients
composition.defineReactionRateCoeffs();

// Initialize electron energy model
electronEnergy energy(mesh);

// Sum of posivite ion fluxes for secondary electron emission wall BC
surfaceScalarField sumPosFlux
(
    IOobject
    (
        "sumPosFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("temp", dimensionSet(0, -2, -1, 0, 0, 0, 0), scalar(0))
);

PtrList<surfaceScalarField> NFlux(composition.species().size());
PtrList<surfaceScalarField> j(composition.species().size());
PtrList<volVectorField> Gamma(composition.species().size());

// Define species flux
forAll(composition.species(), i)
{
    NFlux.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "NFlux",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            -chargeNumber.get<scalar>(composition.species()[i])*linearInterpolate(mobilityCoeffSpecies[i])*mesh.magSf()*fvc::snGrad(volt)
        )
    );

    j.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "j",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            constant::electromagnetic::e*(-chargeNumber.get<scalar>(composition.species()[i])*linearInterpolate(mobilityCoeffSpecies[i]*N[i])*fvc::snGrad(volt) - linearInterpolate(diffusionCoeffSpecies[i])*fvc::snGrad(N[i]))
        )
    );

    Gamma.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "Gamma",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            -chargeNumber.get<scalar>(composition.species()[i])*mobilityCoeffSpecies[i]*N[i]*fvc::grad(volt) - diffusionCoeffSpecies[i]*fvc::grad(N[i])
        )
    );
}

// Info<< "Creating total current density\n" << endl;
// surfaceScalarField jTot
// (
//     IOobject
//     (
//         "jTot",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("zero", dimCurrent/dimArea, 0.0)
// );

// Info<< "Creating field sigma\n" << endl;
// volScalarField sigma
// (
//     IOobject
//     (
//         "sigma",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("zero", dimensionSet(0, -2, 1, 0, 0, 1, 0), scalar(0))
// );

Info<< "Creating field sigma\n" << endl;
volScalarField sigma
(
    IOobject
    (
        "sigma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, -2, 1, 0, 0, 1, 0), scalar(0))
);

Info<< "Creating field sigmaFlux\n" << endl;
surfaceScalarField sigmaFlux
(
    IOobject
    (
        "sigmaFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, -2, 0, 0, 0, 1, 0), scalar(0))
);

// // Read mass fractions of neutral species if exist
// if (composition.species().find("He") >= 0)
// {
//     Info<< "Reading field YHe\n" << endl;
//     volScalarField YHe
//     (
//         IOobject
//         (
//             "YHe",
//             runTime.timeName(),
//             mesh,
//             IOobject::MUST_READ,
//             IOobject::NO_WRITE
//         ),
//         mesh
//     );
//     N[composition.species().find("He")] = YHe*NGas;
// }

if (composition.species().find("N2") >= 0)
{
    Info<< "Reading field YN2\n" << endl;
    volScalarField YN2
    (
        IOobject
        (
            "YN2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    N[composition.species().find("N2")] = YN2*NGas;
}

if (composition.species().find("O2") >= 0)
{
    Info<< "Reading field YO2\n" << endl;
    volScalarField YO2
    (
        IOobject
        (
            "YO2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    N[composition.species().find("O2")] = YO2*NGas;
}

volScalarField NInert = N[inertIndex];

composition.updateReactionRateCoeffs();

Info<< "Reading field rhoq\n" << endl;
volScalarField rhoq
(
    IOobject
    (
        "rhoq",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("temp", dimensionSet(0,-3,1,0,0,1,0), scalar(0))
);
