#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

LSimmiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

mixture.correct();
//volScalarField& alpha1(mixture.alpha1());
//volScalarField& alpha2(mixture.alpha2());
volScalarField& d(mixture.d());
volScalarField d0
(
    IOobject
    (
        "d0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    d//mesh,
        //dimensionedScalar("psi0",dimless, 0.0)
);


volScalarField& H(mixture.H());
const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    rho2+(rho1-rho2)*H
//    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();

// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);
/*
// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);

*/
#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
//mesh.setFluxRequired(alpha1.name());

#include "createMRF.H"

//set prescribed U

const dictionary& UDict(mesh.solverDict(U.name()));

bool prescribedU
(
    UDict.lookupOrDefault<bool>("prescribedU",false)
);

//set periodic or reverse
const scalar period(UDict.lookupOrDefault<scalar>("period",0));
scalar reverseTime(UDict.lookupOrDefault<>("reverseTime",0));


const volVectorField U0(U);
const surfaceScalarField phi0(phi);
