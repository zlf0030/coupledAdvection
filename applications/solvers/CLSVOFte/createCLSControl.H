IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dictionary LevelSetParameter
(
    transportProperties.subDict("LevelSetParameter")
);

const scalar ENUM
(
    readScalar(LevelSetParameter.lookup("interface_thickness"))
);

const scalar DX
(
    2.0/(max(mesh.deltaCoeffs()).value())
);

const scalar epsilon(ENUM*DX);

const scalar CFL
(
    readScalar(LevelSetParameter.lookup("cfl_number"))
);

dimensionedScalar sigma
(
    transportProperties.lookup("sigma")
);