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

Class
    Foam::CLSCutFace

Description
    Class for cutting a face, faceI, of an fvMesh, mesh_, at its intersection
    with an CLSsurface defined by the mesh point values f_ and the CLSvalue,
    CLSValue_.

    Original code supplied by Johan Roenby, DHI (2016)
    Modified by Linfan Zhang, THU (2018)

SourceFiles
    CLSCutFace.C

\*---------------------------------------------------------------------------*/

#ifndef CLSCutFace_H
#define CLSCutFace_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class CLSCutFaces Declaration
\*---------------------------------------------------------------------------*/

class CLSCutFace
{
    // Private data

        //- Mesh whose cells and faces to cut at their intersection with an
        //  CLSface
        const fvMesh& mesh_;

        //- CLSfunction values at mesh points. f_size() = mesh_.nPoints()
//        scalarField& f_;
        //-Normal vector
        const vectorField& nHat_;

        //- Point along first cut edge where CLSsurface cuts edge
        scalar firstEdgeCut_;

        //- Point along last cut edge where CLSsurface cuts edge
        scalar lastEdgeCut_;

        //- Index in mesh_.faces()[faceI_] of first fully submerged (f > f0)
        //  face point
        label firstFullySubmergedPoint_;

        //- Index in mesh_.faces()[faceI_] of last fully submerged (f > f0)
        //  face point
        label nFullySubmergedPoints_;

        //- Storage for centre of subface
        point subFaceCentre_;

        //- Storage for area vector of subface
        vector subFaceArea_;

        //- Storage for subFacePoints
        DynamicList<point> subFacePoints_;

        //- Storage for subFacePoints
        DynamicList<point> surfacePoints_;

        //- Boolean telling if subface centre and area have been calculated
        bool subFaceCentreAndAreaIsCalculated_;


    // Private Member Functions

        void calcSubFaceCentreAndArea();

        //- Calculate cut points along edges of face with values f[pLabels]
        //  Returns the face status, where:
        //  -1: face is fully below the CLSsurface
        //   0: face is cut, i.e. has values larger and smaller than CLSValue
        //  +1: face is fully above the CLSsurface
        label calcSubFace
        (
            const scalar CLSValue,
            const pointField& points,
            const scalarField& f,
            const labelList& pLabels
        );

        void subFacePoints(const pointField& points, const labelList& pLabels);

        void surfacePoints(const pointField& points, const labelList& pLabels);


public:

    // Constructors

        //- Construct from fvMesh and a scalarField
        //  Length of scalarField should equal number of mesh points
        CLSCutFace(const fvMesh& mesh, const vectorField& nHat);

    // Static data

        static int debug;


    // Member functions

        //- Calculate cut points along edges of faceI
        label calcSubFace
        (
            const label& celli,
            const label& faceI, 
            const scalar& CLSValue
        );

        //- Calculate cut points along edges of face with values f
        label calcSubFace
        (
            const pointField& points,
            const scalarField& f,
            const scalar& CLSValue
        );

        const point& subFaceCentre();

        const vector& subFaceArea();

        const DynamicList<point>& subFacePoints() const;

        const DynamicList<point>& surfacePoints() const;

        void clearStorage();

        //- Calculate time integrated area for a face
        scalar timeIntegratedArea
        (
            const pointField& fPts,
            const scalarField& pTimes,
            const scalar dt,
            const scalar magSf,
            const scalar Un0
        );

        void cutPoints
        (
            const pointField& pts,
            const scalarField& f,
            const scalar f0,
            DynamicList<point>& cutPoints
        );


        void quadAreaCoeffs
        (
            const DynamicList<point>& pf0,
            const DynamicList<point>& pf1,
            scalar& alpha,
            scalar& beta
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
