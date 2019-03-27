// Create level set field psi
volScalarField psi
(
    IOobject
    (
        "psi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("psi",dimensionSet(0,0,0,0,0,0,0),scalar(0))
);


volScalarField band
(
    IOobject
    (
        "band",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Band",dimensionSet(0,0,0,0,0,0,0),scalar(4))
);

volScalarField band0 = band;

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
        IOobject::NO_WRITE
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
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("nHat",dimless,vector(0,0,0))
);

//Info<<"epsilon = "<<epsilon<<endl;
//Info <<"psi_max="<<max(psi).value()<<tab<<"psi_min="<<min(psi).value()<<endl;

//#include "calcPsiFromVOF.H"
psi==Foam::asin(double(2.0)*alpha1 - double(1.0))*epsilon/M_PI;
#include "reinitialization.H"
// Instantiating CLSCorrection object for advancing alpha1

//#include "calcHeaviside.H"
forAll (band,celli)
{
//    H[celli] = double(1.0)/double(2.0)*(double(1.0)+psi[celli]/epsilon+Foam::sin(M_PI*psi[celli]/epsilon)/M_PI);
    delta[celli] = double(1.0)/(double(2.0)*epsilon)*(double(1.0)+Foam::cos(M_PI*psi[celli]/epsilon));
};

//rho == limitedH*rho1 + (1.0 - limitedH)*rho2;
//const_cast<volScalarField&>(mixture.nu()()) = limitedH*nu1 + (1.0 - limitedH)*nu2;
//volScalarField& nuTemp = const_cast<volScalarField&>(mixture.nu()());
//nuTemp == limitedH*nu1 + (1.0 - limitedH)*nu2;
#include "makeBand.H"
Info <<"calculate normal vector" <<endl;
#include "calcNormalVector.H"
#include "createIsoAdvection.H"