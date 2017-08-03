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

#include "goldsteinActuationDiskSource.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "geometricOneField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(goldsteinActuationDiskSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        goldsteinActuationDiskSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::goldsteinActuationDiskSource::checkData() const
{
    if (Thrust_ <= VSMALL)
    {
        FatalErrorInFunction
           << "Thrust must be greater than zero"
           << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::goldsteinActuationDiskSource::goldsteinActuationDiskSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    p1_(coeffs_.lookup("p1")),
    p2_(coeffs_.lookup("p2")),
    Thrust_(readScalar(coeffs_.lookup("Thrust"))),
    Torque_(readScalar(coeffs_.lookup("Torque"))),
    Rp_(readScalar(coeffs_.lookup("Rp"))),
    Rh_(readScalar(coeffs_.lookup("Rh")))
{
    coeffs_.lookup("fieldNames") >> fieldNames_;
    applied_.setSize(fieldNames_.size(), false);

    Info<< "    - creating actuation disk zone: "
        << this->name() << endl;

    checkData();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::goldsteinActuationDiskSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    const scalarField& cellsV = mesh_.V();
    vectorField& Usource      = eqn.source();
    const vectorField& U      = eqn.psi();

    if (V() > VSMALL)
    {
        addActuationDiskAxialInertialResistance
        (
            Usource,
            cells_,
            cellsV,
            geometricOneField(),
            U
        );
    }
}


void Foam::fv::goldsteinActuationDiskSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    const scalarField& cellsV = mesh_.V();
    vectorField& Usource = eqn.source();
    const vectorField& U = eqn.psi();

    if (V() > VSMALL)
    {
        addActuationDiskAxialInertialResistance
        (
            Usource,
            cells_,
            cellsV,
            geometricOneField(),
            U
        );
    }
}


bool Foam::fv::goldsteinActuationDiskSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.readIfPresent("p1", p1_);
        coeffs_.readIfPresent("p2", p2_);
        coeffs_.readIfPresent("Thrust", Thrust_);
        coeffs_.readIfPresent("Torque", Torque_);
        coeffs_.readIfPresent("Rp", Rp_);
        coeffs_.readIfPresent("Rh", Rh_);

        checkData();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
