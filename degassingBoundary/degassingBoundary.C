/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "degassingBoundary.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "fvmSup.H"

#include "volFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::degassingBoundary<Type>::degassingBoundary
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    patchNames_(dict.lookup("patchNames"))
{
    read(dict);

    forAll(patchNames_, patchi)
    {
        const fvPatch& patchField =
            mesh_.boundary()[patchNames_[patchi]];

        forAll(patchField.faceCells(), i)
        {
            cells_.append(patchField.faceCells()[i]);
            
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fv::degassingBoundary<Type>::addSup
(
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    const GeometricField<Type, fvPatchField, volMesh>& psi = eqn.psi();
    const scalar deltaT = mesh_.time().deltaT().value();

    typename GeometricField<Type, fvPatchField, volMesh>::Internal Su
    (
        IOobject
        (
            name_ + fieldNames_[fieldi] + "Su",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensioned<Type>
        (
            "zero",
            eqn.dimensions()/dimVolume,
            Zero
        ),
        false
    );


    const volScalarField& alpha = mesh_.objectRegistry::lookupObject<volScalarField>("alpha.air");
    const volScalarField&  rho  = mesh_.objectRegistry::lookupObject<volScalarField>("thermo:rho.air");

    if (psi.name() == "thermo:rho.air")
    {
        forAll(cells_, i)
        {
            Su[cells_[i]] = - alpha[cells_[i]]*psi[cells_[i]]/deltaT;
        }
    }else{
        forAll(cells_, i)
        {
            Su[cells_[i]] = - alpha[cells_[i]]*rho[cells_[i]]*psi[cells_[i]]/deltaT;
        }
    }


    //UIndirectList<Type>(Su, cells_) = alpha[cells_[i]]*rho[cells_[i]]/deltaT;

    eqn +=  Su;
}

template<class Type>
void Foam::fv::degassingBoundary<Type>::addSup
(
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    if (debug)
    {
        Info<< "degassingBoundary<" << pTraits<Type>::typeName
            << ">::addSup for source " << name_ << endl;
    }

    return this->addSup(eqn, fieldi);
}

// ************************************************************************* //
