psi0 == psi;

volScalarField psi1(psi);
volScalarField psi2(psi);

dimensionedScalar dimChange
(
    dimensionedScalar("dimChange",dimTime, 1.0)
);

psi1= psi0 - dimChange*fvc::div(phi,psi0)*runTime.deltaTValue();
psi2= double(0.75)*psi0 +double(0.25)*psi1 - double(0.25)*dimChange*fvc::div(phi,psi1)*runTime.deltaTValue();
psi= double(1.0/3.0)*psi0+double(2.0/3.0)*psi2 - double(2.0/3.0)*dimChange*fvc::div(phi,psi2)*runTime.deltaTValue();
