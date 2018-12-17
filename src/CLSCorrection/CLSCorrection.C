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

#include "CLSCorrection.H"
#include "volFields.H"
#include "interpolationCellPoint.H"
#include "volPointInterpolation.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcGrad.H"
//#include "upwind.H"
#include "cellSet.H"
#include "meshTools.H"
#include "scalarMatrices.H"
#include "surfaceFields.H"
//#include "OBJstream.H"

// * * * * * * * * * * * * * * Debugging * * * * * * * * * * * * * //

#ifndef DebugInfo
//Taken from OpenFOAM-4.0/src/OpenFOAM/db/error/messageStream.H to make code
//compile with older OF versions.
#define DebugInfo                                                              \
    if (debug) Info
#endif

#ifndef InfoInFunction
#define InfoInFunction InfoIn(__func__)
#endif

#ifndef DebugInFunction
#define DebugInFunction                                                        \
    if (debug) InfoInFunction
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(CLSCorrection, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CLSCorrection::CLSCorrection
(
    volScalarField& alpha1,
    volScalarField& psi
)
:
    // General data
    mesh_(alpha1.mesh()),
    dict_(mesh_.solverDict(alpha1.name())),
    alpha1_(alpha1),
    alpha1In_(alpha1.ref()),

//    CorrectionTime_(0),

    psi_(psi),

//    epsilon_(epsilon),
    // Interpolation data
    ap_(mesh_.nPoints()),

    // Cell cutting data
//    surfCells_(label(0.2*mesh_.nCells())),
    isoCutCell_(mesh_, ap_)
{
//    CLSCutCell::debug = debug;
//    CLSCutFace::debug = debug;

    // Prepare lists used in parallel runs
/*    if (Pstream::parRun())
    {
        // Force calculation of required demand driven data (else parallel
        // communication may crash)
        mesh_.cellCentres();
        mesh_.cellVolumes();
        mesh_.faceCentres();
        mesh_.faceAreas();
        mesh_.magSf();
        mesh_.boundaryMesh().patchID();
        mesh_.cellPoints();
        mesh_.cellCells();
        mesh_.cells();

        // Get boundary mesh and resize the list for parallel comms
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        surfaceCellFacesOnProcPatches_.resize(patches.size());

        // Append all processor patch labels to the list
        forAll(patches, patchi)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchi])
             && patches[patchi].size() > 0
            )
            {
                procPatchLabels_.append(patchi);
            }
        }
    }
*/
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CLSCorrection::correct()
{
    ap_ = volPointInterpolation::New(mesh_).interpolate(psi_);
    const scalarField& V = mesh_.V();
    scalar psi2 = 0.05*max(psi_).value();
    scalar psi1 = 0.1*min(psi_).value();
    scalar psi3 = 0.05*min(psi_).value();
    scalar psi4 = 0;
    scalar  sumV = getSumVolume();
    scalar sumV1 = calculateVolume(psi1);
    scalar sumV2 = calculateVolume(psi2);
    scalar sumV3 = calculateVolume(psi3);
    scalar sumV4 = calculateVolume(psi4);

    scalarField a(4), f(4),C(4);
    {
        a[0] = sumV1, a[1] = sumV3, a[2] = sumV4, a[3] = sumV2;
        f[0] = 0, f[1] = (psi3 - psi1)/(psi2 - psi1), f[2] = (psi4-psi1)/(psi2 - psi1), f[3] = 1;
        scalarSquareMatrix M(4);
        forAll(f,i)
        {
            forAll(f,j)
            {
                M[i][j] = pow(f[i], 3 - j);
            }
        }
        C = a;
        LUsolve(M,C);
    }

    //Finding root with Newton Method
    psi3 = f[1]; sumV3 = a[1];
    label nIter = 0;
    scalar tol = 1e-7;
    scalar maxIter = 100;
    scalar res = mag(sumV3 - sumV);
    while (res > tol && nIter < 10*maxIter)
    {
        psi3 -=
            (C[0]*pow3(psi3) + C[1]*sqr(psi3) + C[2]*psi3 + C[3] - sumV)
           /(3*C[0]*sqr(psi3) + 2*C[1]*psi3 + C[2]);
        sumV3 = C[0]*pow3(psi3) + C[1]*sqr(psi3) + C[2]*psi3 + C[3];
        res = mag(sumV3 - sumV);
        nIter++;
    }

    // Scaling back to original range
    psi3 = psi3*(psi2 - psi1) + psi1;
//    Info <<"psi3 = "<<psi3<<tab<<"nIter = "<<nIter<<endl;
    scalar VOF = 0;
    // Check result
    forAll(alpha1_,celli)
    {
        if (alpha1_[celli]==1)
            VOF += V[celli];
        else if((alpha1_[celli]<1)&&(alpha1_[celli]>0))
        {
            isoCutCell_.calcSubCell(celli,psi3);
            scalar ss = (isoCutCell_.volumeOfFluid())*V[celli];
//            Info<<"volume = "<<vv<<tab<<"cut volume = "<<ss<<endl;
            VOF += ss;
        }
    }    
    scalar RES= mag(VOF - sumV);
//    Info << "RES = "<<RES<<tab<<"VOF = "<<VOF<<tab<<"sumV="<<sumV<<endl;
    if (RES > tol)
    {
        Info << "Newton failed"<<endl;
    }
    else
    {
        Info<< "Newton did the job" << endl;
    }

    psi_ = psi_ - psi3;
//    Info<<"psi1 = "<<psi1<<tab<<"psi2 = "<<psi2<<endl;
//    Info<<"sumV = "<<sumV<<tab<<"sumV0 = "<<sumV0<<tab<<"sumV1 = "<<sumV1<<tab<<"sumV2 = "<<sumV2<<endl;
}

Foam::scalar Foam::CLSCorrection::getSumVolume() const
{
    const scalarField& V = mesh_.V();
    scalar sumV = 0;
    forAll(alpha1_,celli)
    {
        if (alpha1_[celli]==1)
            sumV += V[celli];
        else if((alpha1_[celli]<1)&&(alpha1_[celli]>0))
        {
            sumV += alpha1_[celli]*V[celli];
        }
    } 
    return sumV;
}

Foam::scalar Foam::CLSCorrection::calculateVolume(const scalar& psi)
{
    const scalarField& V = mesh_.V();
    scalar sumV = 0;
    forAll(alpha1_,celli)
    {
        if (alpha1_[celli]==1)
            sumV += V[celli];
        else if((alpha1_[celli]<1)&&(alpha1_[celli]>0))
        {
            isoCutCell_.calcSubCell(celli,psi);
            scalar ss = (isoCutCell_.volumeOfFluid())*V[celli];
//            Info<<"volume = "<<vv<<tab<<"cut volume = "<<ss<<endl;
            sumV += ss;
        }
    }
    return sumV;
}

// ************************************************************************* //