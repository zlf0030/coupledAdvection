// calculate the normal vector nHat
const volVectorField gradPsi = fvc::grad(psi);
nHat = gradPsi/(mag(gradPsi)+deltaN_);
surfaceVectorField nVectfv(fvc::interpolate(nHat));
//surfaceVectorField nVecfv(gradPsif/(mag(gradPsif)+scalar(1.0e-12)/dimChange));
surfaceScalarField nVecf(nVectfv & mesh.Sf());

// calculate new curvature based on psi (LS function)
Curve == -fvc::div(nVecf);