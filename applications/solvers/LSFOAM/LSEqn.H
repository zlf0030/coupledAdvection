//Info <<"2nd max="<<max(d).value()<<tab<<"min="<<min(d).value()<<endl; 
d0 == d;

volScalarField d1(d);
volScalarField d2(d);

dimensionedScalar dimChange
(
    dimensionedScalar("dimChange",dimTime, 1.0)
);

Info <<"Before advection calculation  max="<<max(d).value()<<tab<<"min="<<min(d).value()<<endl; 

d1= d0 - dimChange*fvc::div(phi,d0)*runTime.deltaTValue();
d2= double(0.75)*d0 +double(0.25)*d1 - double(0.25)*dimChange*fvc::div(phi,d1)*runTime.deltaTValue();
d= double(1.0/3.0)*d0+double(2.0/3.0)*d2 - double(2.0/3.0)*dimChange*fvc::div(phi,d2)*runTime.deltaTValue();

Info <<"After advection calculation  max="<<max(d).value()<<tab<<"min="<<min(d).value()<<endl; 