Info <<"Reinitialization Process" <<endl;
psi == (double(2.0)*alpha0 - double(1.0))*epsilon;
#include "makeBand.H"
dimensionedScalar one = dimensionedScalar("one",dimensionSet(0,1,0,0,0,0,0), 1.0);
scalar Dtau = CFL*DX;
scalar NITER = ENUM/CFL;

//Reinitialization with Godunov scheme

psi.correctBoundaryConditions();

volScalarField psiNew = psi;

const labelList& own = mesh.faceOwner();
const labelList& nei = mesh.faceNeighbour();

//Info << "Time0"<<endl;

//Info <<"NITER=" <<NITER<<endl;

const volVectorField& CC= mesh.C();

for(int reinit = 0; reinit<NITER ; reinit++)
{
    volScalarField sign = psi/sqrt(pow(psi,2)+pow(epsilon,2));
//    Info<<"start Time:"<<startTime.clockTimeIncrement()<<endl;
    forAll(band,celli)
    {
//        if (band[celli] !=3)
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
        }
//        else if (alpha1[celli] >= 0.9999)
//            psi[celli] = epsilon;
//        else if (alpha1[celli] <= 1e-6)
//            psi[celli] = -1*epsilon;
    }
    psi.correctBoundaryConditions();    
}

//Info <<"psi_max="<<max(psi).value()<<tab<<"psi_min="<<min(psi).value()<<endl;


// update Dirac function
//forAll(mesh.cells(),celli)
forAll (band,celli)
{
    if (band[celli] == 3)
//    if(mag(psi[celli]) > epsilon)
        {
            delta[celli] = double(0);
            if(psi[celli] < 0)
                H[celli] = double(0);
            else
                H[celli] = double(1);
        }
    else
    {
        H[celli] = double(1.0)/double(2.0)*(double(1.0)+psi[celli]/epsilon+Foam::sin(M_PI*psi[celli]/epsilon)/M_PI);
        delta[celli] = double(1.0)/(double(2.0)*epsilon)*(double(1.0)+Foam::cos(M_PI*psi[celli]/epsilon));
//      Info<<"celli = "<<celli<<tab<<"delta = "<<delta[celli]<<endl;
    }
};

const volScalarField limitedH
(
   "limitedH",
    min(max(H, scalar(0)), scalar(1))
);

H == limitedH;