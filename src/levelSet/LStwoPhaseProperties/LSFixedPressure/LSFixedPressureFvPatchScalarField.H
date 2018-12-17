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

Class
    Foam::LSFixedPressureFvPatchScalarField

Description
    A fixed-pressure LSContactAngle boundary

SourceFiles
    LSFixedPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef LSFixedPressureFvPatchScalarField_H
#define LSFixedPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class LSFixedPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class LSFixedPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Fixed pressure
        scalarField p_;


public:

    //- Runtime type information
    TypeName("LSFixedPressure");


    // Constructors

        //- Construct from patch and internal field
        LSFixedPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        LSFixedPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given LSFixedPressureFvPatchScalarField
        //  onto a new patch
        LSFixedPressureFvPatchScalarField
        (
            const LSFixedPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        LSFixedPressureFvPatchScalarField
        (
            const LSFixedPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new LSFixedPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        LSFixedPressureFvPatchScalarField
        (
            const LSFixedPressureFvPatchScalarField&,
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
                new LSFixedPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the LSFixed pressure
            const scalarField& p() const
            {
                return p_;
            }

            //- Return reference to the LSFixed pressure to allow adjustment
            scalarField& p()
            {
                return p_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
