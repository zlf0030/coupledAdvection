// solve Level-Set function as the re-initialization equation
Info<< "solve the reinitialization equation"<< nl << endl;
#include "makeBand.H"
dimensionedScalar one = dimensionedScalar("one",dimensionSet(0,1,0,0,0,0,0), 1.0);
scalar NITER = ENUM/CFL;
//    Info <<"EPS="<<EPS<<tab<<"NITER="<<NITER<<endl;

scalar Dtau = DX*CFL;
psi == psi0;
#include "Godunov.H"

/*
   psi == psi0;
    
   for (int corr=0; corr<NITER; corr++)
   {
      psi = psi + psi0/mag(psi0)*(double(1)-mag(fvc::grad(psi)*one))*Dtau;
      psi.correctBoundaryConditions();
   }


// update Dirac function
   forAll(mesh.cells(),celli)
   {
      if(mag(psi[celli]) > epsilon)
         delta[celli] = double(0);
      else
         delta[celli] = double(1.0)/(double(2.0)*epsilon)*(double(1.0)+Foam::cos(M_PI*psi[celli]/epsilon));
   };

// update Heaviside function

   forAll(mesh.cells(),celli)
   {
      if(psi[celli] < -epsilon)
         H[celli] = double(0);
      else if(epsilon < psi[celli])
         H [celli] = double(1);
      else
         H[celli] = double(1.0)/double(2.0)*(double(1.0)+psi[celli]/epsilon+Foam::sin(M_PI*psi[celli]/epsilon)/M_PI);
   };
//alpha1 = H;
//alpha2 = 1.0 - alpha1;
*/

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

rho == limitedH*rho1 + (1.0 - limitedH)*rho2;

//const_cast<volScalarField&>(mixture.nu()()) = limitedH*nu1 + (1.0 - limitedH)*nu2;

volScalarField& nuTemp = const_cast<volScalarField&>(mixture.nu()());
nuTemp == limitedH*nu1 + (1.0 - limitedH)*nu2;

H == limitedH;
