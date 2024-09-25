/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Foam::fv::degassingBoundary

Description
    This option mimics a degassing boundary condition.

SourceFiles
    degassingBoundary.C

\*---------------------------------------------------------------------------*/

#ifndef degassingBoundary_H
#define degassingBoundary_H

#include "cellSet.H"
/*
#include "fvOption.H"
#include "surfaceFields.H"
#include "volFields.H"
*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

// Forward declaration of classes

template<class Type>
class SemiImplicitSource;


// Forward declaration of friend functions

template<class Type>
Ostream& operator<<
(
    Ostream&,
    const SemiImplicitSource<Type>&
);


/*---------------------------------------------------------------------------*\
                       Class degassingBoundary Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class degassingBoundary
:
    public option
{
protected:

    // Protected data

        //- Names of outlet boundary patches
        wordList patchNames_;

        //- Set of cells to apply source to
        labelList cells_;

    // Protected Member Functions

        //- Check if the phase is a degassing phase
        bool isDegassingPhase(const word& psiName) const;


public:

    //- Runtime type information
    TypeName("degassingBoundary");


    // Constructors

        //- Construct from components
        degassingBoundary
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        // Evaluation

            //- Add explicit contribution to equation
            virtual void addSup
            (
                fvMatrix<Type>& eqn,
                const label fieldi
            );

            //- Add explicit contribution to compressible equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<Type>& eqn,
                const label fieldi
            );

        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "degassingBoundary.C"
    #include "degassingBoundaryIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
