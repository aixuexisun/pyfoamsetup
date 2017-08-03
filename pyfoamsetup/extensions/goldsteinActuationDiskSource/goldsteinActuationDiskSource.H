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

Class
    Foam::fv::goldsetinActuationDiskSource

Description
    Actuation disk source based on Goldstein optimum

SourceFiles
    goldsteinActuationDiskSource.C
    goldsetinActuationDiskSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef goldsteinActuationDiskSource_H
#define goldsteinActuationDiskSource_H

#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class actuationDiskSource Declaration
\*---------------------------------------------------------------------------*/

class goldsteinActuationDiskSource
:
    public cellSetOption
{

protected:

    // Protected data

        //- Thrust coefficient
        scalar Thrust_;

        //- Torque coefficient
        scalar Torque_;

        //- Disk p1
        vector p1_;

        //- Disk p2
        vector p2_;

        // radius
        scalar Rp_;

        // hub radius
        scalar Rh_;

private:

    // Private Member Functions

        //- Check data
        void checkData() const;

        //- Add resistance to the UEqn
        template<class RhoFieldType>
        void addActuationDiskAxialInertialResistance
        (
            vectorField& Usource,
            const labelList& cells,
            const scalarField& V,
            const RhoFieldType& rho,
            const vectorField& U
        ) const;

        //- Disallow default bitwise copy construct
        goldsteinActuationDiskSource(const goldsteinActuationDiskSource&);

        //- Disallow default bitwise assignment
        void operator=(const goldsteinActuationDiskSource&);


public:

    //- Runtime type information
    TypeName("goldsteinActuationDiskSource");


    // Constructors

        //- Construct from components
        goldsteinActuationDiskSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~goldsteinActuationDiskSource()
    {}


    // Member Functions

        // Access

            //- Return Thrust
            scalar Thrust() const
            {
                return Thrust_;
            }

            //- Return Kq
            scalar Torque() const
            {
                return Torque_;
            }

            //- Return radiusPropeller
            scalar Rp() const
            {
                return Rp_;
            }

            //- Return radiusHub
            scalar radiusHub() const
            {
                return Rh_;
            }

            //- First point in actuator disk
            const vector& p1() const
            {
                return p1_;
            }

            //- Second point in actuator disk
            const vector& p2() const
            {
                return p2_;
            }

         // Add explicit and implicit contributions

            //- Source term to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "goldsteinActuationDiskSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
