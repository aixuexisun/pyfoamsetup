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

\*---------------------------------------------------------------------------*/

#include "goldsteinActuationDiskSource.H"
#include "volFields.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::fv::goldsteinActuationDiskSource::addActuationDiskAxialInertialResistance
(
    vectorField& Usource,
    const labelList& cells,
    const scalarField& Vcells,
    const RhoFieldType& rho,
    const vectorField& U
) const
{
    vector diskDir = p2_ - p1_;
    scalar L = mag(diskDir);
    vector v = diskDir/L;

	scalar Delta = L;

	scalar factor = Delta*mathematical::pi*(3*Rh_ + 4*Rp_)*(Rp_ - Rh_);
    scalar Ax     = (105.0/8.0)*Thrust_/factor;
    scalar Atheta = (105.0/8.0)*Torque_/(Rp_*factor);

    scalar d      = 0.0;
    scalar r      = 0.0;
    scalar r_star = 0.0;
	scalar fbx    = 0.0;
	scalar fbt    = 0.0;

    scalar r_mark_hub = Rh_/Rp_;

    vector s(vector::zero);
    vector rVector(vector::zero);
    vector v_t(vector::zero);

    forAll(cells, i)
    {   
        s = mesh().cellCentres()[cells[i]] - p1_;
        d = s & v;

        if ((d > 0) && (d < L)) {
            rVector = s - d*v;

            r = mag(rVector);

            if ((r > Rh_) && (r < Rp_)) {
                r_star = (r - Rh_)/(Rp_-Rh_);

                // Add axial body force
				fbx = Ax*r_star*sqrt(1-r_star);
                Usource[cells[i]] += -fbx*Vcells[cells[i]]*v;

                // Add tangential body force
				fbt = Atheta*(r_star*sqrt(1 - r_star)/(r_star*(1-r_mark_hub) + r_mark_hub));
                v_t = v ^ rVector;
				v_t = v_t/mag(v_t);
                Usource[cells[i]] += -fbt*Vcells[cells[i]]*v_t;
            } 
        }

    }
}


// ************************************************************************* //
