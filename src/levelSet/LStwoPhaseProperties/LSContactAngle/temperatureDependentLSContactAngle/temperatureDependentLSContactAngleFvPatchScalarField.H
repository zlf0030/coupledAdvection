/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 OpenFOAM Foundation
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
    Foam::temperatureDependentLSContactAngleFvPatchScalarField

Description
    Temperature-dependent LSContactAngle scalar boundary condition.

    The contact angle is evaluated from the specified Foam::Function1 for the
    temperature field looked-up from the mesh database the name of which
    may optionally be provided.

Usage
    \table
        Property     | Description             | Required    | Default value
        T            | Temperature field name  | no          | T
        theta0       | Contact angle function  | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            temperatureDependentLSContactAngle;
        theta0          constant 60;
    }
    \endverbatim

See also
    Foam::LSContactAngleFvPatchScalarField
    Foam::constantLSContactAngleFvPatchScalarField
    Foam::Function1

SourceFiles
    temperatureDependentLSContactAngleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef temperatureDependentLSContactAngleFvPatchScalarField_H
#define temperatureDependentLSContactAngleFvPatchScalarField_H

#include "LSContactAngleFvPatchScalarField.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class temperatureDependentLSContactAngleFvPatch Declaration
\*---------------------------------------------------------------------------*/

class temperatureDependentLSContactAngleFvPatchScalarField
:
    public LSContactAngleFvPatchScalarField
{
    // Private data

        //- Name of temperature field, default = "T"
        word TName_;

        //- Equilibrium contact angle function
        autoPtr<Function1<scalar>> theta0_;


public:

    //- Runtime type information
    TypeName("temperatureDependentLSContactAngle");


    // Constructors

        //- Construct from patch and internal field
        temperatureDependentLSContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        temperatureDependentLSContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  temperatureDependentLSContactAngleFvPatchScalarField
        //  onto a new patch
        temperatureDependentLSContactAngleFvPatchScalarField
        (
            const temperatureDependentLSContactAngleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        temperatureDependentLSContactAngleFvPatchScalarField
        (
            const temperatureDependentLSContactAngleFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new temperatureDependentLSContactAngleFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        temperatureDependentLSContactAngleFvPatchScalarField
        (
            const temperatureDependentLSContactAngleFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new temperatureDependentLSContactAngleFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Return the equilibrium contact-angle
        virtual tmp<scalarField> theta
        (
            const fvPatchVectorField& Up,
            const fvsPatchVectorField& nHat
        ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
