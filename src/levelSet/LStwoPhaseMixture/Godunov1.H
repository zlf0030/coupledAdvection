volScalarField psiNew = psi;
const faceList & ff =mesh.faces();
const pointField & pp = mesh.points();
const labelList& own = mesh.faceOwner();
const labelList& nei = mesh.faceNeighbour();

forAll(ff,facei)
{
    if(mesh.isInternalFace(ff[facei]))
    {    
        const cell& cellowner = mesh.cells()[own[facei]];
        const point& op = cellowner.centre(pp,ff);
        const cell& cellnei   = mesh.cells()[nei[facei]];
        const point& np = cellnei.centre(pp,ff);
        Dim = (op - np);
        nDim = Dim/mag(op - np);
        dc_ = (psi[own[facei]] - psi[nei[facei]])/mag(Dim);
        inpX = nDim&vector(1,0,0);
        inpY = nDim&vector(0,1,0);
        inpZ = nDim&vector(0,0,1);

        if (inpX>0)
        {
            owner.right() -= dc_*inpX;
            neigh.left() -= dc_*inpX;
        }
        else if (inpX<0)
        {
            owner.left() -= dc_*inpX;
            neigh.right() -= dc_*inpX;
        }
        
        if (inpY>0)
        {
            owner.north() -= dc_*inpY;
            neigh.sourth() -= dc_*inpY;
        }
        else if (inpY<0)
        {
            owner.sourth() -= dc_*inpY;
            neigh.north() -= dc_*inpY;
        }


        if (inpZ>0)
        {
            owner.top() -= dc_*inpZ;
            neigh.bottom() -= dc_*inpZ;
        }
        else if (inpZ<0)
        {
            owner.bottom() -= dc_*inpZ;
            neigh.top() -= dc_*inpZ;
        }
    }
    
    forAll(mesh.boundary(),patchi)
    {
        const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
        
        const fvsPatchField<scalar>& psisf = psi.boundaryField()[patchi];
        
        forAll(mesh.boundary()[patchi],facei)
        {
            
        }

        const point& fp = mesh.faceCentres()[ff[facei]];
        const cell& cellowner = mesh.cells()[own[ff[facei]]];
        const point& op = cellowner.centre(pp,ff);
        Dim = op -fp;
        nDim = Dim/mag(Dim);

        dc_ = (psi[own[facei]] - psi[nei[facei]])/mag(Dim);
    }

}

