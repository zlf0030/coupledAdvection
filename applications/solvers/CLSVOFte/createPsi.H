// Create level set field psi
volScalarField psi
(
    IOobject
    (
        "psi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField psi0
(
    IOobject
    (
        "psi0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    psi
);

volScalarField H
(
    IOobject
    (
        "H",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("H",dimensionSet(0,0,0,0,0,0,0),scalar(0))
);

volScalarField delta
(
    IOobject
    (
        "delta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("delta",dimless,0.0)
);

volScalarField Curve
(
    IOobject
    (
        "Curve",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Curve",dimless/dimLength,0.0)
);

// create Normal vector
Info << "get nHat" <<endl;
volVectorField nHat
(
    IOobject
    (
        "nHat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//Info<<"epsilon = "<<epsilon<<endl;
//Info <<"psi_max="<<max(psi).value()<<tab<<"psi_min="<<min(psi).value()<<endl;

//#include "calcPsiFromVOF.H"
#include "reinitialization.H"
// Instantiating CLSCorrection object for advancing alpha1
CLSCorrection corrector(alpha0,psi);
Info <<"calculate normal vector" <<endl;
#include "calcNormalVector.H"

#include "createCLSAdvection.H"