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

Class
    Foam::LSinterfaceProperties

Description
    Contains the interface properties.

    Properties to aid interFoam:
    -# Correct the levelset boundary condition for dynamic contact angle.
    -# Calculate interface curvature with levelset method.

SourceFiles
    interfaceProperties.C

\*---------------------------------------------------------------------------*/

#ifndef LSinterfaceProperties_H
#define LSinterfaceProperties_H

#include "IOdictionary.H"
#include "surfaceTensionModel.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class interfaceProperties Declaration
\*---------------------------------------------------------------------------*/

class LSinterfaceProperties
{
    //Private data
        //- Keep a reference to the transportProperties dictionary
        const dictionary& transportPropertiesDict_;
/*        
        //- Compression coeffient
        scalar cAlpha_;
*/        
        //- Surface tension
        autoPtr<surfaceTensionModel> sigmaPtr_;

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;

        //- levelset Functions
        const volScalarField& d_;

        //- Heaveside Functions
        const volScalarField& H_;
        
        //- Phase fraction
//        const volScalarField& alpha1_;
        
        const volVectorField& U_;
        
        //- Normal vector
        surfaceScalarField nHatf_;
        
        //- curvature
        volScalarField K_;
    

    //Private Member Functions
        
        //-Disallow default bitwise copy construct and assignment
        LSinterfaceProperties(const LSinterfaceProperties&);
        void operator=(const LSinterfaceProperties&);

        //- Correction for the boundary condition on the unit normal nHat on
        //  walls to produce the correct contact dynamic angle
        //  calculated from the component of U parallel to the wall
        void correctContactAngle
        (
            surfaceVectorField::Boundary& nHat,
            const surfaceVectorField::Boundary& gradLSf
        ) const;

//        void calculateH();
        //- Re-calculate the interface curvature
        void calculateK();


public:

    //- Conversion factor for degrees into radians
    static const scalar convertToRad;


    // Constructors

        //- Construct from volume fraction field gamma and IOdictionary
        LSinterfaceProperties
        (
//            const volScalarField& alpha1, 
            const volScalarField& d,
            const volScalarField& H,
            const volVectorField& U,
            const IOdictionary&
        );


    // Member Functions
/*
        scalar cAlpha() const
        {
            return cAlpha_;
        }
*/
        const dimensionedScalar& deltaN() const
        {
            return deltaN_;
        }

        const surfaceScalarField& nHatf() const
        {
            return nHatf_;
        }

        tmp<volScalarField> sigmaK() const;

        tmp<surfaceScalarField> surfaceTensionForce() const;
/*
        //- Indicator of the proximity of the interface
        //  Field values are 1 near and 0 away for the interface.
        tmp<volScalarField> nearInterface() const;
*/
        void correct();

        //- Read transportProperties dictionary
        bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
