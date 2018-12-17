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
    Foam::CLSCutCell

Description
    Class for cutting a cell, celli, of an fvMesh, mesh_, at its intersection
    with an CLSsurface defined by the mesh point values f_ and the CLSvalue,
    CLSValue_.


    Reference:
        \verbatim
            Roenby, J., Bredmose, H. and Jasak, H. (2016).
            A computational method for sharp interface advection
            Royal Society Open Science, 3
            doi 10.1098/rsos.160405
        \endverbatim

    Original code supplied by Johan Roenby, DHI (2016)
    A little change for coupled with level set by Linfan Zhang,THU (2018)

SourceFiles
    CLSCutCell.C

\*---------------------------------------------------------------------------*/

#ifndef CLSCutCell_H
#define CLSCutCell_H

#include "fvMesh.H"
#include "volFieldsFwd.H"
#include "CLSCutFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class CLSCutCell Declaration
\*---------------------------------------------------------------------------*/

class CLSCutCell
{
    // Private data

        //- Mesh whose cells and faces to cut at their intersection with an
        //  CLSsurface.
        const fvMesh& mesh_;

        //- Cell to cut
        label cellI_;

        //- CLSfunction values at mesh points. f_size() = mesh_.nPoints().
//        scalarField& f_;
        //- Normal vector 
        const vectorField& nHat_;

        //- Initial psi
        scalarField& psi_;

        //- CLSvalue used to cut cell
        scalar CLSValue_;

        //- An CLSCutFace object to get access to its face cutting functionality
        CLSCutFace CLSCutFace_;

        //- List of face labels of CLSCutFaces
        DynamicList<label> CLSCutFaces_;

        //- List of point lists each defining an CLSCutFace
        DynamicList< DynamicList<point> > CLSCutFacePoints_;

        //- List of face centres for CLSCutFaces
        DynamicList<point> CLSCutFaceCentres_;

        //- List of face area vectors for CLSCutFaces
        DynamicList<vector> CLSCutFaceAreas_;

        //- Storage for subFace edges belonging to CLSFace
        DynamicList<DynamicList<point> > CLSFaceEdges_;

        //- Points constituting the cell-CLSsurface intersection (CLSface)
        DynamicList<point> CLSFacePoints_;

        //- Face centre of the CLSface
        point CLSFaceCentre_;

        //- Face normal of the CLSface by convention pointing from high to low
        //  values (i.e. opposite of the gradient vector).
        vector CLSFaceArea_;

        //- Cell centre of the subcell of celli which is "fully submerged", i.e.
        //  where the function value is higher than the CLSValue_
        point subCellCentre_;

        //- Volume of fully submerged subcell
        scalar subCellVolume_;

        //- Volume of Fluid for celli (subCellVolume_/mesh_.V()[celli])
        scalar VOF_;

        //- List of fully submerged faces
        DynamicList<label> fullySubFaces_;

        //- A cell status label taking one of the values:
        //
        //  - -1: cell is fully below the CLSsurface
        //  -  0: cell is cut
        //  - +1: cell is fully above the CLSsurface
        label cellStatus_;

        //- Boolean telling if subcell centre and volume have been calculated
        bool subCellCentreAndVolumeCalculated_;

        //- Boolean telling if CLSface centre and area have been calculated
        bool CLSFaceCentreAndAreaCalculated_;


    // Private Member Functions

            void calcSubCellCentreAndVolume();

            void calcCLSFaceCentreAndArea();

            void calcCLSFacePointsFromEdges();


public:

    // Constructors

        //- Construct from fvMesh and a scalarField
        //  Length of scalarField should equal number of mesh points
        CLSCutCell(const fvMesh&, const vectorField& ,scalarField&);

    // Static data

        static int debug;

    // Member functions

        label calcSubCell(const label celli, const scalar CLSValue);

        const point& subCellCentre();

        scalar subCellVolume();

        const DynamicList<point>& CLSFacePoints();

        const point& CLSFaceCentre();

        const vector& CLSFaceArea();

        scalar volumeOfFluid();

        scalar CLSValue() const;

        void clearStorage();

        label vofCutCell
        (
            const label celli,
            const scalar alpha1,
            const scalar tol,
            const label maxIter
        );

        void volumeOfFluid(volScalarField& alpha1, const scalar f0);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
