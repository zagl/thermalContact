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

Class
    Foam::compressible::thermalContactResistanceFvPatchScalarField

Description
    Mixed boundary condition for temperature, to be used for heat-transfer
    on back-to-back baffles, with thermal resistance. 

    Example usage:
        myInterfacePatchName
        {
            type        compressible::thermalContactResistance;
            neighbourFieldName  T;
            kappa               lookup;
            kappaName           kappa;
            d                   uniform 0.001;
            K                   uniform 10;
            value               uniform 300;
        }

    Needs to be on underlying mapped(Wall)FvPatch.

    Note: kappa : heat conduction at patch. Gets supplied how to lookup
        calculate kappa:
    - 'lookup' : lookup volScalarField (or volSymmTensorField) with name
    - 'fluidThermo' : use fluidThermo and compressible::RASmodel to calculate
        kappa
    - 'solidThermo' : use solidThermo kappa()
    - 'directionalSolidThermo' directionalKappa()

    Note: runs in parallel with arbitrary decomposition. Uses mapped
    functionality to calculate exchange.

SourceFiles
    thermalContactResistanceFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef thermalContactResistanceFvPatchScalarField_H
#define thermalContactResistanceFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
   Class thermalContactResistanceFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class thermalContactResistanceFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public temperatureCoupledBase
{

public:
        //- Operation mode enumeration
        enum operationMode
        {
            fixedThermalResistance,
            fixedThicknessAndKappa,
            fixedThermalContactConductance,
            unknown
        };
        static const NamedEnum<operationMode, 3> operationModeNames;

private:
    // Private data

        //- Name of field on the neighbour region
        const word neighbourFieldName_;
        
        //- Thermal resistance [K/W]
        mutable scalarField resistance_;
        
        //- Thickness [m]
        mutable scalarField thickness_;
        
        //- Kappa [W/m/K]
        mutable scalarField kappa_;
        
        //- Thermal contact conductance [W/m^2/K]
        mutable scalarField conductance_;
        
        //- Operation mode
        operationMode opMode_;


public:

    //- Runtime type information
    TypeName("compressible::thermalContactResistance");


    // Constructors

        //- Construct from patch and internal field
        thermalContactResistanceFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        thermalContactResistanceFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  thermalContactResistanceFvPatchScalarField onto a
        //  new patch
        thermalContactResistanceFvPatchScalarField
        (
            const thermalContactResistanceFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new thermalContactResistanceFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        thermalContactResistanceFvPatchScalarField
        (
            const thermalContactResistanceFvPatchScalarField&,
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
                new thermalContactResistanceFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
