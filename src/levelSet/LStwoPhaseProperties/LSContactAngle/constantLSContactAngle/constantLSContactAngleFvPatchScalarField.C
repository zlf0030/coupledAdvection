/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "constantLSContactAngleFvPatchScalarField.H"
#include "volMesh.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantLSContactAngleFvPatchScalarField::
constantLSContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    LSContactAngleFvPatchScalarField(p, iF),
    theta0_(0.0)
{}


Foam::constantLSContactAngleFvPatchScalarField::
constantLSContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    LSContactAngleFvPatchScalarField(p, iF, dict),
    theta0_(readScalar(dict.lookup("theta0")))
{
    evaluate();
}


Foam::constantLSContactAngleFvPatchScalarField::
constantLSContactAngleFvPatchScalarField
(
    const constantLSContactAngleFvPatchScalarField& gcpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    LSContactAngleFvPatchScalarField(gcpsf, p, iF, mapper),
    theta0_(gcpsf.theta0_)
{}


Foam::constantLSContactAngleFvPatchScalarField::
constantLSContactAngleFvPatchScalarField
(
    const constantLSContactAngleFvPatchScalarField& gcpsf
)
:
    LSContactAngleFvPatchScalarField(gcpsf),
    theta0_(gcpsf.theta0_)
{}


Foam::constantLSContactAngleFvPatchScalarField::
constantLSContactAngleFvPatchScalarField
(
    const constantLSContactAngleFvPatchScalarField& gcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    LSContactAngleFvPatchScalarField(gcpsf, iF),
    theta0_(gcpsf.theta0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::constantLSContactAngleFvPatchScalarField::theta
(
    const fvPatchVectorField&,
    const fvsPatchVectorField&
) const
{
    return tmp<scalarField>(new scalarField(size(), theta0_));
}


void Foam::constantLSContactAngleFvPatchScalarField::write
(
    Ostream& os
) const
{
    LSContactAngleFvPatchScalarField::write(os);
    os.writeKeyword("theta0") << theta0_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        constantLSContactAngleFvPatchScalarField
    );
}

// ************************************************************************* //
