/*---------------------------------------------------------------------------*\
              Original work | Copyright (C) 2016-2017 DHI
              Modified work | Copyright (C) 2016-2017 OpenCFD Ltd.
              Modified work | Copyright (C) 2017-2018 Johan Roenby
-------------------------------------------------------------------------------

License
    This file is part of CLSAdvector which is an extension to OpenFOAM.

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

#include "CLSCutCell.H"
#include "scalarMatrices.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Debugging * * * * * * * * * * * * * //

#ifndef DebugPout
//Taken from OpenFOAM-plus/src/OpenFOAM/db/error/messageStream.H to make code
//compile with older OF versions.
#define DebugPout                                                              \
    if (debug) Pout
#endif

#ifndef InfoInFunction
#define InfoInFunction InfoIn(__func__)
#endif

#ifndef DebugInFunction
#define DebugInFunction                                                        \
    if (debug) InfoInFunction
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::CLSCutCell::debug = 0;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CLSCutCell::CLSCutCell(const fvMesh& mesh, const vectorField& nHat,scalarField& Psi)
:
    mesh_(mesh),
    cellI_(-1),
//    f_(f),
    nHat_(nHat),
    psi_(Psi),
    CLSValue_(0),
    CLSCutFace_(CLSCutFace(mesh_,nHat_)),
    CLSCutFaces_(10),
    CLSCutFacePoints_(10),
    CLSCutFaceCentres_(10),
    CLSCutFaceAreas_(10),
    CLSFaceEdges_(10),
    CLSFacePoints_(10),
    CLSFaceCentre_(vector::zero),
    CLSFaceArea_(vector::zero),
    subCellCentre_(vector::zero),
    subCellVolume_(-10),
    VOF_(-10),
    fullySubFaces_(10),
    cellStatus_(-1),
    subCellCentreAndVolumeCalculated_(false),
    CLSFaceCentreAndAreaCalculated_(false)
{
    clearStorage();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CLSCutCell::calcSubCellCentreAndVolume()
{
    if (cellStatus_ == 0)
    {
        subCellCentre_ = vector::zero;
        subCellVolume_ = 0.0;

        // Estimate the approximate cell centre as the average of face centres
        label nCellFaces(1 + CLSCutFaceCentres_.size() + fullySubFaces_.size());
        vector cEst = CLSFaceCentre_ + sum(CLSCutFaceCentres_);
        forAll(fullySubFaces_, facei)
        {
            cEst += mesh_.faceCentres()[fullySubFaces_[facei]];
        }
        cEst /= scalar(nCellFaces);


        // Contribution to subcell centre and volume from CLSface
        const scalar pyr3Vol0 =
            max(mag(CLSFaceArea_ & (CLSFaceCentre_ - cEst)), VSMALL);

        // Calculate face-pyramid centre
        const vector pc0 = 0.75*CLSFaceCentre_ + 0.25*cEst;

        // Accumulate volume-weighted face-pyramid centre
        subCellCentre_ += pyr3Vol0*pc0;

        // Accumulate face-pyramid volume
        subCellVolume_ += pyr3Vol0;

        // Contribution to subcell centre and volume from cut faces
        forAll(CLSCutFaceCentres_, facei)
        {
            // Calculate 3*face-pyramid volume
            scalar pyr3Vol =
                max
                (
                    mag
                    (
                        CLSCutFaceAreas_[facei]
                      & (CLSCutFaceCentres_[facei] - cEst)
                    ),
                    VSMALL
                );

            // Calculate face-pyramid centre
            vector pc = 0.75*CLSCutFaceCentres_[facei] + 0.25*cEst;

            // Accumulate volume-weighted face-pyramid centre
            subCellCentre_ += pyr3Vol*pc;

            // Accumulate face-pyramid volume
            subCellVolume_ += pyr3Vol;
        }

        // Contribution to subcell centre and volume from fully submerged faces
        forAll(fullySubFaces_, i)
        {
            const label facei = fullySubFaces_[i];
            const point& fCentre = mesh_.faceCentres()[facei];
            const vector& fArea = mesh_.faceAreas()[facei];

            // Calculate 3*face-pyramid volume
            scalar pyr3Vol = max(mag(fArea & (fCentre - cEst)), VSMALL);

            // Calculate face-pyramid centre
            vector pc = 0.75*fCentre + 0.25*cEst;

            // Accumulate volume-weighted face-pyramid centre
            subCellCentre_ += pyr3Vol*pc;

            // Accumulate face-pyramid volume
            subCellVolume_ += pyr3Vol;
        }

        subCellCentre_ /= subCellVolume_;
        subCellVolume_ /= scalar(3);
        VOF_ = subCellVolume_/mesh_.cellVolumes()[cellI_];

        subCellCentreAndVolumeCalculated_ = true;

        if (debug)
        {
            vector sumSf = CLSFaceArea_;
            scalar sumMagSf = mag(CLSFaceArea_);
            forAll(CLSCutFaceCentres_, facei)
            {
                sumSf += CLSCutFaceAreas_[facei];
                sumMagSf += mag(CLSCutFaceAreas_[facei]);
            }
            forAll(fullySubFaces_, facei)
            {
                sumSf += mesh_.faceAreas()[fullySubFaces_[facei]];
                sumMagSf += mag(CLSCutFaceAreas_[facei]);
            }
            if (mag(sumSf) > 1e-10)
            {
                Pout<< "Warninig: mag(sumSf)/magSumSf = "
                    << mag(sumSf)/sumMagSf << " for surface cell"
                    << cellI_ << endl;
            }
        }
    }
    else if (cellStatus_ == 1)
    {
        // Cell fully above CLSsurface
        subCellCentre_ = vector::zero;
        subCellVolume_ = 0;
        VOF_ = 0;
    }
    else if (cellStatus_ == -1)
    {
        // Cell fully below CLSsurface
        subCellCentre_ = mesh_.cellCentres()[cellI_];
        subCellVolume_ = mesh_.cellVolumes()[cellI_];
        VOF_ = 1;
    }
}


void Foam::CLSCutCell::calcCLSFaceCentreAndArea()
{
    // Initial guess of face centre from edge points
    point fCentre = vector::zero;
    label nEdgePoints = 0;
    forAll(CLSFaceEdges_, ei)
    {
        DynamicList<point>& edgePoints = CLSFaceEdges_[ei];
        forAll(edgePoints, pi)
        {
            fCentre += edgePoints[pi];
            nEdgePoints++;
        }
    }

    if (nEdgePoints > 0)
    {
        fCentre /= nEdgePoints;
    }
    else
    {
        DebugPout << "Warning: nEdgePoints = 0 for cell " << cellI_ << endl;
    }

    vector sumN = vector::zero;
    scalar sumA = 0;
    vector sumAc = vector::zero;

    forAll(CLSFaceEdges_, ei)
    {
        const DynamicList<point>& edgePoints = CLSFaceEdges_[ei];
        const label nPoints = edgePoints.size();
        for (label pi = 0; pi < nPoints-1; pi++)
        {
            const point& nextPoint = edgePoints[pi + 1];

            vector c = edgePoints[pi] + nextPoint + fCentre;
            vector n = (nextPoint - edgePoints[pi])^(fCentre - edgePoints[pi]);
            scalar a = mag(n);

            // Edges may have different orientation
            sumN += Foam::sign(n & sumN)*n;
            sumA += a;
            sumAc += a*c;
        }
    }

    // This is to deal with zero-area faces. Mark very small faces
    // to be detected in e.g., processorPolyPatch.
    if (sumA < ROOTVSMALL)
    {
        CLSFaceCentre_ = fCentre;
        CLSFaceArea_ = vector::zero;
    }
    else
    {
        CLSFaceCentre_ = sumAc/sumA/scalar(3);
        CLSFaceArea_ = 0.5*sumN;
    }


    // Check CLSFaceArea_ direction and change if not pointing out of subcell
    if ((CLSFaceArea_ & (CLSFaceCentre_ - subCellCentre())) < 0)
    {
        CLSFaceArea_ *= (-1);
    }

    CLSFaceCentreAndAreaCalculated_ = true;
}


void Foam::CLSCutCell::calcCLSFacePointsFromEdges()
{
    DebugPout
        << "Enter calcCLSFacePointsFromEdges() with CLSFaceArea_ = "
        << CLSFaceArea_ << " and CLSFaceCentre_ = " << CLSFaceCentre_
        << " and CLSFaceEdges_ = " << CLSFaceEdges_ << endl;

    // Defining local coordinates with zhat along CLSface normal and xhat from
    // CLSface centre to first point in CLSFaceEdges_
    const vector zhat = CLSFaceArea_/mag(CLSFaceArea_);
    vector xhat = CLSFaceEdges_[0][0] - CLSFaceCentre_;
    xhat = (xhat - (xhat & zhat)*zhat);
    xhat /= mag(xhat);
    vector yhat = zhat^xhat;
    yhat /= mag(yhat);

    DebugPout << "Calculated local coordinates" << endl;

    // Calculating CLSface point angles in local coordinates
    DynamicList<point> unsortedCLSFacePoints(3*CLSFaceEdges_.size());
    DynamicList<scalar> unsortedCLSFacePointAngles(3*CLSFaceEdges_.size());
    forAll(CLSFaceEdges_, ei)
    {
        const DynamicList<point>& edgePoints = CLSFaceEdges_[ei];
        forAll(edgePoints, pi)
        {
            const point& p = edgePoints[pi];
            unsortedCLSFacePoints.append(p);
            unsortedCLSFacePointAngles.append
            (
                Foam::atan2
                (
                    ((p - CLSFaceCentre_) & yhat),
                    ((p - CLSFaceCentre_) & xhat)
                )
            );
        }
    }

    DebugPout<< "Calculated CLSFace point angles" << endl;

    // Sorting CLSface points by angle and inserting into CLSFacePoints_
    labelList order(unsortedCLSFacePointAngles.size());
    Foam::sortedOrder(unsortedCLSFacePointAngles, order);
    CLSFacePoints_.append(unsortedCLSFacePoints[order[0]]);
    for (label pi = 1; pi < order.size(); pi++)
    {
        if
        (
            mag
            (
                unsortedCLSFacePointAngles[order[pi]]
               -unsortedCLSFacePointAngles[order[pi-1]]
            ) > 1e-8
        )
        {
            CLSFacePoints_.append(unsortedCLSFacePoints[order[pi]]);
        }
    }

    DebugPout<< "Sorted CLSface points by angle" << endl;
}


Foam::label Foam::CLSCutCell::calcSubCell
(
    const label celli,
    const scalar CLSValue
)
{
    // Populate CLSCutFaces_, CLSCutFacePoints_, fullySubFaces_, CLSFaceCentres_
    // and CLSFaceArea_.

    clearStorage();
    cellI_ = celli;
    CLSValue_ = CLSValue;
    const cell& c = mesh_.cells()[celli];
    forAll(c, fi)
    {
        const label facei = c[fi];

        const label faceStatus = CLSCutFace_.calcSubFace(celli,facei, CLSValue_);
        if (faceStatus == 0)
        {
            // Face is cut
            CLSCutFacePoints_.append(CLSCutFace_.subFacePoints());
            CLSCutFaceCentres_.append(CLSCutFace_.subFaceCentre());
            CLSCutFaceAreas_.append(CLSCutFace_.subFaceArea());
            CLSFaceEdges_.append(CLSCutFace_.surfacePoints());
        }
        else if (faceStatus == -1)
        {
            // Face fully below
            fullySubFaces_.append(facei);
        }
    }
    if (CLSCutFacePoints_.size())
    {
        // Cell cut at least at one face
        cellStatus_ = 0;
        calcCLSFaceCentreAndArea();
        // In the rare but occuring cases where a cell is only touched at a
        // point or a line the CLSFaceArea_ will have zero length and here the
        // cell should be treated as either completely empty or full.
        if (mag(CLSFaceArea_) < 10*SMALL)
        {
            if (fullySubFaces_.empty())
            {
                // Cell fully above CLSsurface
                cellStatus_ = 1;
            }
            else
            {
                // Cell fully below CLSsurface
                cellStatus_ = -1;
            }
        }
    }
    else if (fullySubFaces_.empty())
    {
        // Cell fully above CLSsurface
        cellStatus_ = 1;
    }
    else
    {
        // Cell fully below CLSsurface
        cellStatus_ = -1;
    }

    return cellStatus_;
}


const Foam::point& Foam::CLSCutCell::subCellCentre()
{
    if (!subCellCentreAndVolumeCalculated_)
    {
        calcSubCellCentreAndVolume();
    }

    return subCellCentre_;
}


Foam::scalar Foam::CLSCutCell::subCellVolume()
{
    if (!subCellCentreAndVolumeCalculated_)
    {
        calcSubCellCentreAndVolume();
    }

    return subCellVolume_;
}


const Foam::DynamicList<Foam::point>& Foam::CLSCutCell::CLSFacePoints()
{
    if (cellStatus_ == 0 && CLSFacePoints_.size() == 0)
    {
        calcCLSFacePointsFromEdges();
    }

    return CLSFacePoints_;
}


const Foam::point& Foam::CLSCutCell::CLSFaceCentre()
{
    if (!CLSFaceCentreAndAreaCalculated_)
    {
        calcCLSFaceCentreAndArea();
    }

    return CLSFaceCentre_;
}


const Foam::vector& Foam::CLSCutCell::CLSFaceArea()
{
    if (!CLSFaceCentreAndAreaCalculated_)
    {
        calcCLSFaceCentreAndArea();
    }

    return CLSFaceArea_;
}


Foam::scalar Foam::CLSCutCell::volumeOfFluid()
{
    if (!subCellCentreAndVolumeCalculated_)
    {
        calcSubCellCentreAndVolume();
    }

    return VOF_;
}


Foam::scalar Foam::CLSCutCell::CLSValue() const
{
    return CLSValue_;
}


void Foam::CLSCutCell::clearStorage()
{
    cellI_ = -1;
    CLSValue_ = 0;
    CLSCutFace_.clearStorage();
    CLSCutFaces_.clear();
    CLSCutFacePoints_.clear();
    CLSCutFaceCentres_.clear();
    CLSCutFaceAreas_.clear();
    CLSFaceEdges_.clear();
    CLSFacePoints_.clear();
    CLSFaceCentre_ = vector::zero;
    CLSFaceArea_ = vector::zero;
    subCellCentre_ = vector::zero;
    subCellVolume_ = -10;
    VOF_ = -10;
    fullySubFaces_.clear();
    cellStatus_ = -1;
    subCellCentreAndVolumeCalculated_ = false;
    CLSFaceCentreAndAreaCalculated_ = false;
}


Foam::label Foam::CLSCutCell::vofCutCell
(
    const label celli,
    const scalar alpha1,
    const scalar tol,
    const label maxIter
)
{
    DebugInFunction
        << "vofCutCell for cell " << celli << " with alpha1 = "
        << alpha1 << " ------" << endl;

    // Finding cell vertex extrema values
    const labelList& pLabels = mesh_.cellPoints(celli);
    scalarField fvert(pLabels.size());
    
    const vectorField& cellCentres = mesh_.cellCentres();
    const vector& CC = cellCentres[celli];
//    scalar& psi = psi_[celli];
    const pointField& points = mesh_.points();
    const scalar& aa = nHat_[celli].x();
    const scalar& bb = nHat_[celli].y();
    const scalar& cc = nHat_[celli].z();
    scalar d0 =  aa*CC.x() + bb*CC.y() + cc*CC.z();
// If fvert > 0 , the point is in phase2. vice versa.
    forAll(pLabels, pi)
    {
//        fvert[pi] = f_[pLabels[pi]];
        const vector& pp = points[pLabels[pi]];
        fvert[pi] = (aa*pp.x() + bb*pp.y() + cc*pp.z() - d0);
    }

    labelList order(fvert.size());
    sortedOrder(fvert, order);
    scalar f1 = fvert[order.first()];
    scalar f2 = fvert[order.last()];
//    DebugPout << "fvert = " << fvert << ", and order = " << order << endl;

    // Handling special case where method is handed an almost full/empty cell
    if (alpha1 < tol)
    {
//        psi = f1;
        return calcSubCell(celli, f1);
    }
    else if (1 - alpha1 < tol)
    {
//        psi = f2;
        return calcSubCell(celli, f2);
    }

    // Finding the two vertices inbetween which the CLSvalue giving alpha1 lies
    label L1 = 0;
    label L2 = fvert.size() - 1;
    scalar a1 = 1;
    scalar a2 = 0;
    scalar L3, f3, a3;

    while (L2 - L1 > 1)
    {
        L3 = round(0.5*(L1 + L2));
        f3 = fvert[order[L3]];
        calcSubCell(celli, f3);
        a3 = volumeOfFluid();
        if (a3 > alpha1)
        {
            L1 = L3; f1 = f3; a1 = a3;
        }
        else if (a3 < alpha1)
        {
            L2 = L3; f2 = f3; a2 = a3;
        }
    }
    if (mag(f1 - f2) < 10*SMALL)
    {
        DebugPout<< "Warning: mag(f1 - f2) < 10*SMALL" << endl;
//        psi = f1;
        return calcSubCell(celli, f1);
    }

    if (mag(a1 - a2) < tol)
    {
        DebugPout<< "Warning: mag(a1 - a2) < tol for cell " << celli << endl;
//        psi = 0.5*(f1+f2);
        return calcSubCell(celli, 0.5*(f1 + f2));
    }

    // Now we know that a(f) = alpha1 is to be found on the f interval
    // [f1, f2], i.e. alpha1 will be in the interval [a2,a1]
    DebugPout
        << "L1 = " << L1 << ", f1 = " << f1 << ", a1 = " << a1 << nl
        << "L2 = " << L2 << ", f2 = " << f2  << ", a2 = " << a2 << endl;


    // Finding coefficients in 3 deg polynomial alpha(f) from 4 solutions

    // Finding 2 additional points on 3 deg polynomial
    f3 = f1 + (f2 - f1)/scalar(3);
    calcSubCell(celli, f3);
    a3 = volumeOfFluid();

    scalar f4 = f1 + 2*(f2 - f1)/3;
    calcSubCell(celli, f4);
    scalar a4 = volumeOfFluid();

    // Building and solving Vandermonde matrix equation
    scalarField a(4), f(4), C(4);
    {
        a[0] = a1, a[1] = a3, a[2] = a4, a[3] = a2;
        f[0] = 0, f[1] = (f3-f1)/(f2-f1), f[2] = (f4-f1)/(f2-f1), f[3] = 1;
        scalarSquareMatrix M(4);
        forAll(f, i)
        {
            forAll(f, j)
            {
                M[i][j] = pow(f[i], 3 - j);
            }
        }

        // C holds the 4 polynomial coefficients
        C = a;
        LUsolve(M, C);
    }

    // Finding root with Newton method

    f3 = f[1]; a3 = a[1];
    label nIter = 0;
    scalar res = mag(a3 - alpha1);
    while (res > tol && nIter < 10*maxIter)
    {
        f3 -=
            (C[0]*pow3(f3) + C[1]*sqr(f3) + C[2]*f3 + C[3] - alpha1)
           /(3*C[0]*sqr(f3) + 2*C[1]*f3 + C[2]);
        a3 = C[0]*pow3(f3) + C[1]*sqr(f3) + C[2]*f3 + C[3];
        res = mag(a3 - alpha1);
        nIter++;
    }
    // Scaling back to original range
    f3 = f3*(f2 - f1) + f1;

    // Check result
    calcSubCell(celli, f3);
    
    const scalar VOF = volumeOfFluid();
    res = mag(VOF - alpha1);

    if (res > tol)
    {
        DebugPout
            << "Newton obtained f3 = " << f3 << " and a3 = " << a3
            << " with mag(a3-alpha1) = " << mag(a3-alpha1)
            << " but calcSubCell(celli,f3) gives VOF  = " << VOF << nl
            << "M(f)*C = a with " << nl
            << "f_scaled = " << f << nl
            << "f = " << f*(f2 - f1) + f1 << nl
            << "a = " << a << nl
            << "C = " << C << endl;
    }
    else
    {
        DebugPout<< "Newton did the job" << endl;
//        psi = f3;
        return cellStatus_;
    }

    // If tolerance not met use the secant method  with f3 as a hopefully very
    // good initial guess to crank res the last piece down below tol

    scalar x2 = f3;
    scalar g2 = VOF - alpha1;
    scalar x1 = max(1e-3*(f2 - f1), 100*SMALL);
    x1 = min(max(x1, f1), f2);
    calcSubCell(celli, x1);
    scalar g1 = volumeOfFluid() - alpha1;

    nIter = 0;
    scalar g0(0), x0(0);
    while (res > tol && nIter < maxIter && g1 != g2)
    {
        x0 = (x2*g1 - x1*g2)/(g1 - g2);
        calcSubCell(celli, x0);
        g0 = volumeOfFluid() - alpha1;
        res = mag(g0);
        x2 = x1; g2 = g1;
        x1 = x0; g1 = g0;
        nIter++;
    }

    if (debug)
    {
        if (res < tol)
        {
            Pout<< "Bisection finished the job in " << nIter << " iterations."
                << endl;
        }
        else
        {
            Pout<< "Warning: Bisection not converged " << endl;
            Pout<< "Leaving vofCutCell with f3 = " << f3 << " giving a3 = "
                << a3 << " so alpha1 - a3 = " << alpha1 - a3 << endl;
        }
    }
//    psi = x0;
    return cellStatus_;
}


void Foam::CLSCutCell::volumeOfFluid
(
    volScalarField& alpha1,
    const scalar f0
)
{
    // Setting internal field
    scalarField& alphaIn = alpha1;
    forAll(alphaIn, celli)
    {
        const label cellStatus = calcSubCell(celli, f0);
        if (cellStatus != 1)
        {
            // If cell not entirely above CLSsurface
            alphaIn[celli] = volumeOfFluid();
        }
    }

    // Setting boundary alpha1 values
    forAll(mesh_.boundary(), patchi)
    {
        if (mesh_.boundary()[patchi].size() > 0)
        {
            const label start = mesh_.boundary()[patchi].patch().start();
            scalarField& alphap = alpha1.boundaryFieldRef()[patchi];
            const scalarField& magSfp = mesh_.magSf().boundaryField()[patchi];
            const labelUList& Cb = mesh_.boundary()[patchi].patch().faceCells(); 

            forAll(alphap, patchFacei)
            {
                const label facei = patchFacei + start;
                const label celli = Cb[patchFacei];
                const label faceStatus = CLSCutFace_.calcSubFace(celli,facei, f0);

                if (faceStatus != 1)
                {
                    // Face not entirely above CLSsurface
                    alphap[patchFacei] =
                        mag(CLSCutFace_.subFaceArea())/magSfp[patchFacei];
                }
            }
        }
    }
}


// ************************************************************************* //
