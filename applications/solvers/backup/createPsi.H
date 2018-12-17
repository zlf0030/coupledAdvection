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

// mapping alpha value to psi0
psi == (double(2.0)*alpha1 - double(1.0))*epsilon;

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

volVectorField nHat
(
    IOobject
    (
        "nHat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
#include "reinitialization.H"

#include "createCLSAdvection.H"