/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "thermalContactResistanceFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermalContactResistanceFvPatchScalarField::
thermalContactResistanceFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), "undefined", "undefined-K"),
    neighbourFieldName_("undefined-neighbourFieldName")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


thermalContactResistanceFvPatchScalarField::
thermalContactResistanceFvPatchScalarField
(
    const thermalContactResistanceFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    temperatureCoupledBase(patch(), ptf.KMethod(), ptf.kappaName()),
    neighbourFieldName_(ptf.neighbourFieldName_)
{}


thermalContactResistanceFvPatchScalarField::
thermalContactResistanceFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    neighbourFieldName_(dict.lookup("neighbourFieldName"))
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "thermalContactResistanceFvPatchScalarField::"
            "thermalContactResistanceFvPatchScalarField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<scalar, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if 
    (
        dict.found("R") && 
        !dict.found("d") && 
        !dict.found("K") &&
        !dict.found("h") 
    )
    {
        opMode_ = fixedThermalResistance;
        resistance_ = scalarField("R", dict, p.size());
    }
    else if 
    (
        !dict.found("R") && 
        dict.found("d") && 
        dict.found("K") &&
        !dict.found("h") 
    )
    {
        opMode_ = fixedThicknessAndKappa;
        thickness_ = scalarField("d", dict, p.size());
        kappa_ = scalarField("K", dict, p.size());
    }
    else if 
    (
        !dict.found("R") && 
        !dict.found("d") && 
        !dict.found("K") &&
        dict.found("h") 
    )
    {
        opMode_ = fixedThermalContactConductance;
        conductance_ = scalarField("h", dict, p.size());
    }
    else
    {
        FatalErrorIn
        (
            "thermalContactResistanceFvPatchScalarField::"
            "thermalContactResistanceFvPatchScalarField\n"
            "(\n"
            " const fvPatch& p,\n"
            " const DimensionedField<scalar, volMesh>& iF,\n"
            " const dictionary& dict\n"
            ")\n"
        ) << "\n patch type '" << p.type()
            << "' either R or d and K or h were not found '"
            << "\n for patch " << p.name()
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
            << exit(FatalError);
    }
    
    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


thermalContactResistanceFvPatchScalarField::
thermalContactResistanceFvPatchScalarField
(
    const thermalContactResistanceFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    temperatureCoupledBase(patch(), wtcsf.KMethod(), wtcsf.kappaName()),
    neighbourFieldName_(wtcsf.neighbourFieldName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void thermalContactResistanceFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const label samplePatchI = mpp.samplePolyPatch().index();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[samplePatchI];
        
    const scalar patchArea = gSum(nbrPatch.magSf());

    scalarField intFld = patchInternalField();


    // Calculate the temperature by harmonic averaging
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const thermalContactResistanceFvPatchScalarField& nbrField =
    refCast
    <
        const thermalContactResistanceFvPatchScalarField
    >
    (
        nbrPatch.lookupPatchField<volScalarField, scalar>
        (
            neighbourFieldName_
        )
    );

    // Swap to obtain full local values of neighbour internal field

    
    scalarField nbrIntFld(nbrField.patchInternalField());
    mpp.distribute(nbrIntFld);

    scalarField nbrFld =
        nbrPatch.lookupPatchField<volScalarField, scalar>
        (
            neighbourFieldName_
        );
    mpp.map().distribute(nbrFld);

    // Swap to obtain full local values of neighbour kappa*delta
    scalarField nbrKDelta(nbrField.kappa(nbrField)*nbrPatch.deltaCoeffs());
    mpp.distribute(nbrKDelta);

    scalarField myKDelta = kappa(*this)*patch().deltaCoeffs();
    
    scalarField KDeltaSolid(patch().size(), 0.0);
    
    
    if (opMode_ == fixedThermalResistance)
    {
        KDeltaSolid = 1/resistance_/patchArea;
    }
    else if (opMode_ == fixedThicknessAndKappa)
    {
        KDeltaSolid = kappa_/thickness_;
    }
    else
    {
        KDeltaSolid = conductance_;
    }
    
    
    
    const scalarField q
    (
        (intFld - nbrIntFld)/(1.0/myKDelta + 1.0/nbrKDelta + 1.0/KDeltaSolid)
    );

    this->refValue() = intFld - q/myKDelta;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar Q = gSum(kappa(*this)*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->dimensionedInternalField().name() << " <- "
            << nbrMesh.name() << ':'
            << nbrPatch.name() << ':'
            << this->dimensionedInternalField().name() << " :"
            << " heat transfer rate:" << Q
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }

    // Restore tag
    UPstream::msgType() = oldTag;
}


void thermalContactResistanceFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("neighbourFieldName")<< neighbourFieldName_
        << token::END_STATEMENT << nl;
    temperatureCoupledBase::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    thermalContactResistanceFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
