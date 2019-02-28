// calculate the normal vector nHat
const volVectorField gradPsi = fvc::grad(psi);
forAll (band,celli)
{
    if (band[celli] > 2)
//    if(mag(psi[celli]) > epsilon)
    {
        nHat[celli] = vector(0,0,0);
    }
    else
    {
        nHat[celli] = gradPsi[celli]/(mag(gradPsi[celli])+scalar(1.0e-12));
    }
};
surfaceVectorField nVectfv(fvc::interpolate(nHat));
//surfaceVectorField nVecfv(gradPsif/(mag(gradPsif)+scalar(1.0e-12)/dimChange));
surfaceScalarField nVecf(nVectfv & mesh.Sf());

// calculate new curvature based on psi (LS function)
Curve == -fvc::div(nVecf);