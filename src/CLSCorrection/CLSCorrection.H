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
    Foam::CLSCorrection

Description
    calculates the new VOF field after the calculation of level set Correction
    The fluid transport calculation is based on an idea of using the projection
    of vector from vertex to cell center on the normal vector which is defined by
    level set function.This is the coupled section to estimate the internal 
    distribution of fluid in cells.
    
    Original code supplied by Johan Roenby, DHI (2016)
    Modified by Linfan Zhang, THU (2018)

SourceFiles
    CLSCorrection.C
    CLSCorrectionTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef CLSCorrection_H
#define CLSCorrection_H
#include "fvMesh.H"
#include "volFieldsFwd.H"
#include "surfaceFields.H"
#include "className.H"
#include "isoCutCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class CLSCorrection
{

    // Private data types

        typedef DynamicList<label> DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point> DynamicPointList;

    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Dictionary for CLSCorrection controls
        const dictionary dict_;

        //- VOF field
        volScalarField& alpha1_;

        //- Often used reference to alpha1 internal field
        scalarField& alpha1In_;

        //- Time spent performing interface Correction
//        scalar CorrectionTime_;

        //Initial level set function: psi
        volScalarField& psi_;

        // Point interpolation data

            //- VOF field interpolated to mesh points
            scalarField ap_;


        // Switches and tolerances. Tolerances need to go into toleranceSwitches

            //- Tolerance for search of CLSFace giving specified VOF value
//            scalar CLSFaceTol_;

            //- Tolerance for marking of surface cells:
            //  Those with surfCellTol_ < alpha1 < 1 - surfCellTol_
//            scalar surfCellTol_;

        // Cell and face cutting

            //- List of surface cells
//            DynamicLabelList surfCells_;

            //- Cell cutting object
            isoCutCell isoCutCell_;

            //- Face cutting object
//            isoCutFace isoCutFace_;


        // Additional data for parallel runs

            //- List of processor patch labels
            DynamicLabelList procPatchLabels_;

            //- For each patch if it is a processor patch this is a list of the
            //  face labels on this patch that are downwind to a surface cell.
            //  For non-processor patches the list will be empty.
            List<DynamicLabelList> surfaceCellFacesOnProcPatches_;
    // Private Member Functions

        //- Disallow default bitwise copy construct
        CLSCorrection(const CLSCorrection&);

        //- Disallow default bitwise copy assignment
        void operator=(const CLSCorrection&);
/*
            //- Determine if a cell is a surface cell
            bool isASurfaceCell(const label celli) const
            {
                return
                (
                    surfCellTol_ < alpha1In_[celli]
                 && alpha1In_[celli] < 1 - surfCellTol_
                );
            }

            //- Clear out CLSFace data
            void clearCLSFaceData()
            {
                surfCells_.clear();
                bsFaces_.clear();
                bsx0_.clear();
                bsn0_.clear();
                bsUn0_.clear();
                bsf0_.clear();
                
                if (mesh_.topoChanging())
                {
                    // Introduced resizing to cope with changing meshes
                    checkBounding_.resize(mesh_.nCells());
                    cellIsBounded_.resize(mesh_.nCells());
                    ap_.resize(mesh_.nPoints());
                }
                checkBounding_ = false;
                cellIsBounded_ = false;

            }

        // Face value functions needed for random face access where the face
        // can be either internal or boundary face

            //- Return face value for a given Geometric surface field
            template<typename Type>
            Type faceValue
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>& f,
                const label facei
            ) const;

            //- Set face value for a given Geometric surface field
            template<typename Type>
            void setFaceValue
            (
                GeometricField<Type, fvsPatchField, surfaceMesh>& f,
                const label facei,
                const Type& value
            ) const;
*/
public:

    //- Runtime type information
    TypeName("CLSCorrection");

    //- Constructors

        //- Construct given alpha, phi and velocity field. Note: phi should be
        // divergence free up to a sufficient tolerance
        CLSCorrection
        (
            volScalarField& alpha1,
            volScalarField& psi
        );


    //- Destructor
    virtual ~CLSCorrection()
    {}


    // Member functions

        //- Advect the free surface. Updates alpha field, taking into account
        //  multiple calls within a single time step.
        void correct();

        // Access functions

            //- Return alpha field
            const volScalarField& alpha() const
            {
                return alpha1_;
            }

            //- Return the controls dictionary
            const dictionary& dict() const
            {
                return dict_;
            }
/*
            scalar CorrectionTime() const
            {
                return CorrectionTime_;
            }
*/
            scalar getSumVolume() const;
            scalar calculateVolume(const scalar& psi);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#   include "CLSCorrectionTemplates.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
