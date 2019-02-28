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
/*
const dimensionedScalar deltaN_
(
    "deltaN_",
    1e-10/pow(average(mesh.V()),1.0/3.0)
);*/
word phase1Name
(
    wordList(transportProperties.lookup("phases"))[0]
);

word phase2Name
(
    wordList(transportProperties.lookup("phases"))[1]
);



dimensionedScalar nu1
(
    transportProperties.subDict(phase1Name).lookup("nu")
);

dimensionedScalar nu2
(
    transportProperties.subDict(phase2Name).lookup("nu")
);
