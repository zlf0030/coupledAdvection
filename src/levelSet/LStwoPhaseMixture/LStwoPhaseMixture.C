/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "LStwoPhaseMixture.H"
#include "fvc.H"
#include "FUNCTIONS.H"
#include "surfaceFields.H"
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LStwoPhaseMixture::calculateH()
{
    H_ = CALC_HEAVISIDE(d_,EPS);
//    Info<<"max(H)="<<max(H_).value()<<endl;
}

void Foam::LStwoPhaseMixture::calculateGrad()
{
    G_ = CALC_GRADIENT(d_,EPS);
//    Info<<"max(G)="<<max(G_).value()<<endl;
}

void Foam::LStwoPhaseMixture::reinitialization()
{
    Info <<"Reinitialization" <<endl;
    const fvMesh& mesh = d_.mesh();
    Info <<"Before Godunov  max= "<<max(d_).value()<<tab<<"min="<<min(d_).value()<<endl;
    dimensionedScalar one = dimensionedScalar("one", dimensionSet(0,1,0,0,0,0,0), 1.0);
    scalar Dtau  = CFL*DX;
    scalar NITER = ENUM/CFL;
//    const scalar& Thic = EPS;
//    volScalarField d0 = d_;
//    scalar ERR(1.0);
//    #include "Reinitialization.H"
    volScalarField & psi = d_;
    #include"Godunov.H"
    Info <<"After Godunov  max="<<max(d_).value()<<tab<<"min="<<min(d_).value()<<endl;

}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LStwoPhaseMixture::LStwoPhaseMixture
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    phase1Name_(wordList(dict.lookup("phases"))[0]),
    phase2Name_(wordList(dict.lookup("phases"))[1]),
   
    alpha1_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase1Name_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    alpha2_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase2Name_),
            mesh.time().timeName(),
            mesh
        ),
        1.0 - alpha1_
    ),

    d_
    (
        IOobject
        (
            "d",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    
    H_
    (
         IOobject
        (
            "H",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("H",dimensionSet(0,0,0,0,0,0,0),scalar(0))
    ),

    G_
    (
         IOobject
        (
            "G",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("G",dimensionSet(0,0,0,0,0,0,0),scalar(0))
    ),
    LevelSetParameter(dict.subDict("LevelSetParameter")),
    ENUM(readScalar(LevelSetParameter.lookup("interface_thickness"))),
    DX(2.0/(max(mesh.deltaCoeffs()).value())),
    EPS(ENUM*DX),
    CFL(readScalar(LevelSetParameter.lookup("cfl_number")))
{
    Info<<"Initialization" <<tab<<"DX="<<DX<<"           EPS="<<EPS<<endl;
    
    d_ == (double(2.0)*alpha1_ - double(1.0))*EPS;
//    Info <<"max= "<<max(d_).value()<<tab<<"min="<<min(d_).value()<<endl;
/*
    reinitialization();
    calculateH();
    calculateGrad();
*/
}

void Foam::LStwoPhaseMixture::correct()
{
//    Info <<"correct at LStwoPhaseMixture"<<endl;
    reinitialization();
    calculateH();
    calculateGrad();
}

// ************************************************************************* //
