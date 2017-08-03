/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "twoDimMotionPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "polyMesh.H"
#include "IFstream.H"
#include "Tuple2.H"
#include "interpolateSplineXY.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

twoDimMotionPointPatchVectorField::twoDimMotionPointPatchVectorField(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
):
    fixedValuePointPatchField<vector>(p, iF),
    origin_(vector::zero),
    p0_(p.localPoints())
{}


twoDimMotionPointPatchVectorField::twoDimMotionPointPatchVectorField (
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
):
    fixedValuePointPatchField<vector>(p, iF, dict),
    origin_(dict.lookup("origin")),
    timeDataFileName_(fileName(dict.lookup("timeDataFileName")).expand())
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }

    if (dict.found("p0"))
    {
        p0_ = vectorField("p0", dict , p.size());
    }
    else
    {
        p0_ = p.localPoints();
    }
}


twoDimMotionPointPatchVectorField::twoDimMotionPointPatchVectorField (
    const twoDimMotionPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
):
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    origin_(ptf.origin_),
    timeDataFileName_(ptf.timeDataFileName_),
    p0_(ptf.p0_, mapper)
{}


twoDimMotionPointPatchVectorField::twoDimMotionPointPatchVectorField (
    const twoDimMotionPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
):
    fixedValuePointPatchField<vector>(ptf, iF),
    origin_(ptf.origin_),
    timeDataFileName_(ptf.timeDataFileName_),
    p0_(ptf.p0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void twoDimMotionPointPatchVectorField::autoMap (const pointPatchFieldMapper& m) {
    fixedValuePointPatchField<vector>::autoMap(m);

    p0_.autoMap(m);
}

void twoDimMotionPointPatchVectorField::rmap (const pointPatchField<vector>& ptf, const labelList& addr) {
    const twoDimMotionPointPatchVectorField& aODptf = refCast<const twoDimMotionPointPatchVectorField>(ptf);

    fixedValuePointPatchField<vector>::rmap(aODptf, addr);

    p0_.rmap(aODptf.p0_, addr);
}

void twoDimMotionPointPatchVectorField::updateCoeffs() {
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = this->internalField().mesh()();

    // ----------- Read in the tabulated data ---------------------------------
    IFstream dataStream(timeDataFileName_);

    if (dataStream.good()) {
        List<Tuple2<scalar, vector> > timeValues
        (
            dataStream
        );

        times_.setSize(timeValues.size());
        x_.setSize(timeValues.size());
        y_.setSize(timeValues.size());
        omega_.setSize(timeValues.size());

        forAll(timeValues, i)
        {
            times_[i] = timeValues[i].first();
            x_[i]     = timeValues[i].second()[0];
            y_[i]     = timeValues[i].second()[1];
            omega_[i] = timeValues[i].second()[2];
        }
    }
    else {
        FatalErrorInFunction
            << "Cannot open time data file " << timeDataFileName_
            << exit(FatalError);
    }

    // --------------- Interpolate to the right time step ----------------------
    const scalar t = mesh.time().value();

    scalar x     = interpolateSplineXY(t, times_, x_);
    scalar y     = interpolateSplineXY(t, times_, y_);
    scalar omega = interpolateSplineXY(t, times_, omega_);

    // --------------- Calculate transformation --------------------------------
    vector axis  = vector(0.0, 0.0, 1.0);

    vector translation = vector(x, y, 0.0);

    vectorField p0Rel(p0_ - origin_);

    vectorField::operator=
    (
        p0Rel*(cos(omega) - 1)
      + (axis ^ p0Rel*sin(omega))
      + (axis & p0Rel)*(1 - cos(omega))*axis
      + (translation)
    );

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void twoDimMotionPointPatchVectorField::write(Ostream& os) const {
    pointPatchField<vector>::write(os);
    os.writeKeyword("origin")
        << origin_ << token::END_STATEMENT << nl;
    os.writeKeyword("timeDataFileName")
        << timeDataFileName_ << token::END_STATEMENT << nl;
    p0_.writeEntry("p0", os);
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField (pointPatchVectorField, twoDimMotionPointPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
