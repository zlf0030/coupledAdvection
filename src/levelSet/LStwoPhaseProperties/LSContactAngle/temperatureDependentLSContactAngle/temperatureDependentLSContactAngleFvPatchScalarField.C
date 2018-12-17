/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenFOAM Foundation
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

#include "temperatureDependentLSContactAngleFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperatureDependentLSContactAngleFvPatchScalarField::
temperatureDependentLSContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    LSContactAngleFvPatchScalarField(p, iF),
    TName_("T"),
    theta0_()
{}


Foam::temperatureDependentLSContactAngleFvPatchScalarField::
temperatureDependentLSContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    LSContactAngleFvPatchScalarField(p, iF, dict),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    theta0_(Function1<scalar>::New("theta0", dict))
{
    evaluate();
}


Foam::temperatureDependentLSContactAngleFvPatchScalarField::
temperatureDependentLSContactAngleFvPatchScalarField
(
    const temperatureDependentLSContactAngleFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    LSContactAngleFvPatchScalarField(psf, p, iF, mapper),
    TName_(psf.TName_),
    theta0_(psf.theta0_, false)
{}


Foam::temperatureDependentLSContactAngleFvPatchScalarField::
temperatureDependentLSContactAngleFvPatchScalarField
(
    const temperatureDependentLSContactAngleFvPatchScalarField& psf
)
:
    LSContactAngleFvPatchScalarField(psf),
    TName_(psf.TName_),
    theta0_(psf.theta0_, false)
{}


Foam::temperatureDependentLSContactAngleFvPatchScalarField::
temperatureDependentLSContactAngleFvPatchScalarField
(
    const temperatureDependentLSContactAngleFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    LSContactAngleFvPatchScalarField(psf, iF),
    TName_(psf.TName_),
    theta0_(psf.theta0_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::temperatureDependentLSContactAngleFvPatchScalarField::theta
(
    const fvPatchVectorField&,
    const fvsPatchVectorField&
) const
{
    return theta0_->value
    (
        patch().lookupPatchField<volScalarField, scalar>(TName_)
    );
}


void Foam::temperatureDependentLSContactAngleFvPatchScalarField::write
(
    Ostream& os
) const
{
    LSContactAngleFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "T", "T", TName_);
    theta0_->writeData(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        temperatureDependentLSContactAngleFvPatchScalarField
    );
}

// ************************************************************************* //
