/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    along with OpenFOAM.  If Not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "LSincompressibleTwoPhaseMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LSincompressibleTwoPhaseMixture, 0);
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void Foam::LSincompressibleTwoPhaseMixture::calcNu()
{
    nuModel1_->correct();
    nuModel2_->correct();

/*    const volScalarField limitedH
    (
        "limitedH",
        min(max(H_, scalar(0)), scalar(1))
    );
*/
  
    auto n1 = H_*rho1_ + (scalar(1) - H_)*rho2_;
    auto n3 = mu();
    auto n2 = n3/n1;
    nu_ = n2;
    //    nu_ = mu()/(rho2_ + (rho1_ - rho2_)*H_);
    // Average kinematic viscosity calculated from dynamic viscosity
    
    //    nu_ = mu()/(limitedAlpha1*rho1_ + (scalar(1) - limitedAlpha1)*rho2_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LSincompressibleTwoPhaseMixture::LSincompressibleTwoPhaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    LStwoPhaseMixture(U.mesh(),*this),
    nuModel1_
    (
        viscosityModel::New
        (
            "nu1",
            subDict(phase1Name_),
            U,
            phi
        )
    ),
    nuModel2_
    (
        viscosityModel::New
        (
            "nu2",
            subDict(phase2Name_),
            U,
            phi
        )
    ),

    rho1_("rho", dimDensity, nuModel1_->viscosityProperties()),
    rho2_("rho", dimDensity, nuModel2_->viscosityProperties()),

    U_(U),
    phi_(phi),

    nu_
    (
        IOobject
        (
            "nu",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar("nu", dimViscosity, 0)
//        calculatedFvPatchScalarField::typeName
    )
{
    calcNu();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::LSincompressibleTwoPhaseMixture::mu() const
{
/*
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );
*/
    return tmp<volScalarField>
    (
        new volScalarField
        (
            "mu",
            rho2_*nuModel2_->nu()       + 
            rho1_*nuModel1_->nu()*H_    - 
            rho2_*nuModel2_->nu()*H_  
//            limitedAlpha1*rho1_*nuModel1_->nu()
//          + (scalar(1) - limitedAlpha1)*rho2_*nuModel2_->nu()
        )
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::LSincompressibleTwoPhaseMixture::muf() const
{
/*    const surfaceScalarField alpha1f
    (
        min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1))
    );
*/
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "muf",
//            alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
//          + (scalar(1) - alpha1f)*rho2_*fvc::interpolate(nuModel2_->nu())
            rho2_*fvc::interpolate(nuModel2_->nu()) +
            (rho1_*fvc::interpolate(nuModel1_->nu())-
             rho2_*fvc::interpolate(nuModel2_->nu()))*fvc::interpolate(H_)
        )
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::LSincompressibleTwoPhaseMixture::nuf() const
{
/*    const surfaceScalarField alpha1f
    (
        min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1))
    );
*/
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "nuf",
//            (
//                alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
//              + (scalar(1) - alpha1f)*rho2_*fvc::interpolate(nuModel2_->nu())
//            )/(alpha1f*rho1_ + (scalar(1) - alpha1f)*rho2_)
            (
            (rho2_*fvc::interpolate(nuModel2_->nu()) +
            (rho1_*fvc::interpolate(nuModel1_->nu())-
             rho2_*fvc::interpolate(nuModel2_->nu()))*fvc::interpolate(H_)           
            )/(rho2_+(rho1_-rho2_)*fvc::interpolate(H_))
            )   
        )
    );
}


bool Foam::LSincompressibleTwoPhaseMixture::read()
{
    if (regIOobject::read())
    {
        if
        (
            nuModel1_().read
            (
                subDict(phase1Name_ == "1" ? "phase1": phase1Name_)
            )
         && nuModel2_().read
            (
                subDict(phase2Name_ == "2" ? "phase2": phase2Name_)
            )
        )
        {
            nuModel1_->viscosityProperties().lookup("rho") >> rho1_;
            nuModel2_->viscosityProperties().lookup("rho") >> rho2_;

            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
