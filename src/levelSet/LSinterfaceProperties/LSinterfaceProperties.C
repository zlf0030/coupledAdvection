/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "LSinterfaceProperties.H"
#include "LSContactAngleFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const Foam::scalar Foam::LSinterfaceProperties::convertToRad =
    Foam::constant::mathematical::pi/180.0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.


void Foam::LSinterfaceProperties::correctContactAngle
(
    surfaceVectorField::Boundary& nHatb,
    const surfaceVectorField::Boundary& gradLSf
) const
{
    const fvMesh& mesh = d_.mesh();
    const volScalarField::Boundary& abf = d_.boundaryField();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<LSContactAngleFvPatchScalarField>(abf[patchi]))
        {
            LSContactAngleFvPatchScalarField& acap =
                const_cast<LSContactAngleFvPatchScalarField&>
                (
                    refCast<const LSContactAngleFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                convertToRad*acap.theta(U_.boundaryField()[patchi], nHatp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle

            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

            acap.gradient() = (nf & nHatp)*mag(gradLSf[patchi]);
            acap.evaluate();
        }
    }
}

void Foam::LSinterfaceProperties::calculateK()
{
    const fvMesh& mesh = d_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();

    // Cell gradient of leveset 
    const volVectorField gradLS(fvc::grad(d_, "nHat"));

    // Interpolated face-gradient of leveset
    surfaceVectorField gradLSf(fvc::interpolate(gradLS));

    //gradLSf -=
    //    (mesh.Sf()/mesh.magSf())
    //   *(fvc::snGrad(d_) - (mesh.Sf() & gradLSf)/mesh.magSf());

    // Face unit interface normal
    surfaceVectorField nHatfv(gradLSf/(mag(gradLSf) + deltaN_));
    // surfaceVectorField nHatfv
    // (
    //     (gradLSf + deltaN_*vector(0, 0, 1)
    //    *sign(gradLSf.component(vector::Z)))/(mag(gradLSf) + deltaN_)
    // );
    correctContactAngle(nHatfv.boundaryFieldRef(), gradLSf.boundaryField());

    // Face unit interface normal flux
    nHatf_ = nHatfv & Sf;

    // Simple expression for curvature
    K_ = -fvc::div(nHatf_);

    // Complex expression for curvature.
    // Correction is formally zero but numerically non-zero.
    /*
    volVectorField nHat(gradLS/(mag(gradLS) + deltaN_));
    forAll(nHat.boundaryField(), patchi)
    {
        nHat.boundaryField()[patchi] = nHatfv.boundaryField()[patchi];
    }

    K_ = -fvc::div(nHatf_) + (nHat & fvc::grad(nHatfv) & nHat);
    */
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LSinterfaceProperties::LSinterfaceProperties
(
//    const volScalarField& alpha1,
    const volScalarField& d,
    const volScalarField& H,
    const volVectorField& U,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),
/*   
    cAlpha_
    (
        readScalar
        (
            alpha1.mesh().solverDict(alpha1.name()).lookup("cAlpha")
        )
    ),

    alpha1_(alpha1),
*/
    U_(U),

    d_(d),

    H_(H),
    
    deltaN_
    (
        "deltaN",
        1e-8/pow(average(d.mesh().V()), 1.0/3.0)
    ),

    sigmaPtr_(surfaceTensionModel::New(dict, d.mesh())),

    nHatf_
    (
        IOobject
        (
            "nHatf",
            d_.time().timeName(),
            d_.mesh()
        ),
        d_.mesh(),
        dimensionedScalar("nHatf", dimArea, 0.0)
    ),
    
    K_
    (
        IOobject
        (
            "LSinterfaceProperties:K",
            d_.time().timeName(),
            d_.mesh()
        ),
        d_.mesh(),
        dimensionedScalar("K", dimless/dimLength, 0.0)
    )
{
    calculateK();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::LSinterfaceProperties::sigmaK() const
{
    return sigmaPtr_->sigma()*K_;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::LSinterfaceProperties::surfaceTensionForce() const
{
    return fvc::interpolate(sigmaK())*fvc::snGrad(H_);
}

/*
Foam::tmp<Foam::volScalarField>
Foam::LSinterfaceProperties::nearInterface() const
{
    return pos0(alpha1_ - 0.01)*pos0(0.99 - alpha1_);
}
*/

void Foam::LSinterfaceProperties::correct()
{
//    Info <<"Now correct at LSinterfaceProperties!"<<endl;
    calculateK();
}


bool Foam::LSinterfaceProperties::read()
{
    sigmaPtr_->readDict(transportPropertiesDict_);
    return true;
}


// ************************************************************************* //
