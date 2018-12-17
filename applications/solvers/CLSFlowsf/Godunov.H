//B
//Reinitialization with Godunov scheme
//

psi.correctBoundaryConditions();

volScalarField psiNew = psi;

//const faceList & ff = mesh.faces();
//const pointField & pp = mesh.points();

const labelList& own = mesh.faceOwner();
const labelList& nei = mesh.faceNeighbour();

//Info << "Time0"<<endl;

//clock startTime;
//clock endTime;
//clock cellTime;
Info <<"NITER=" <<NITER<<endl;

const scalar& EPS = epsilon;

const volVectorField& CC= mesh.C();
//const surfaceVectorField& Cf=mesh.Cf();

for(int reinit = 0; reinit<NITER ; reinit++)
{
    volScalarField sign = psi/sqrt(pow(psi,2)+pow(EPS,2));
//    Info<<"start Time:"<<startTime.clockTimeIncrement()<<endl;
//    forAll(mesh.cells(),celli)
    forAll(band,celli)
    {
        if (band[celli] !=3)
        {
        const cell & cc = mesh.cells()[celli];
        const point& cellCenter = CC[celli];
        scalar dc_ = 0;
        vector Dim = vector(0,0,0);
        scalar l =0;
        scalar r =0;
        scalar s =0;
        scalar n =0;
        scalar b =0;
        scalar t =0;
        scalar G =0;
        scalar inipsi =0;
        
        forAll(cc,fi)
        {
            const label& facei = cc[fi];
            if(mesh.isInternalFace(facei))
            {
                if (own[facei]==celli)
                {
                    const point& np = CC[nei[facei]];
                    Dim = (cellCenter - np)/mag(cellCenter - np);
                    dc_ = (psi[celli]-psi[nei[facei]])/mag(Dim);
                }
                else if(nei[facei]==celli)
                {
                    const point& op = CC[own[facei]];
                    Dim = (cellCenter - op)/mag(cellCenter - op);
                    dc_ = (psi[celli]-psi[own[facei]])/mag(Dim);
                }
            }
/*            else
            {
                const point& fp = Cf[facei];
                Dim = (cellCenter - fp)/mag(cellCenter - fp);
                dc_ = scalar(0);
            }*/
            l += dc_*max(Dim&vector(-1,0,0),0);
            r -= dc_*max(Dim&vector(1,0,0),0);
            s += dc_*max(Dim&vector(0,-1,0),0);
            n -= dc_*max(Dim&vector(0,1,0),0);
            b += dc_*max(Dim&vector(0,0,-1),0);
            t -= dc_*max(Dim&vector(0,0,1),0);
        }
        
        if(reinit==0) inipsi = psi[celli];
        

        if(inipsi>0)
            G = Foam::sqrt(max(pow(max(l,0),2),pow(min(r,0),2))+max(pow(max(s,0),2),pow(min(n,0),2))+max(pow(max(b,0),2),pow(min(t,0),2))) - 1;
        else if(inipsi<0)
            G = Foam::sqrt(max(pow(min(l,0),2),pow(max(r,0),2))+max(pow(min(s,0),2),pow(max(n,0),2))+max(pow(min(b,0),2),pow(max(t,0),2))) - 1;
        else if(inipsi==0)
            G = 0;
        psi[celli] -= Dtau*sign[celli]*G;
//        cellTime=clock();
//        Info<<"cell Time:"<<cellTime.clockTimeIncrement()<<endl;
        }
    }
    psi.correctBoundaryConditions();
}
/*
scalar mmax = max(psi).value();
scalar mmin = min(psi).value();

forAll(psi,celli)
{
    if (psi[celli]>0)
        psi[celli] = psi[celli]*EPS/mmax;
    if (psi[celli]<0)
        psi[celli] = -psi[celli]*EPS/mmin;
}
*/

//    Info<<"end time:"<<endTime.clockTimeIncrement()<<endl;

//Info <<"Time0-1"<<endl;

