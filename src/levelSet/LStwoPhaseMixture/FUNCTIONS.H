#ifndef FUNCTIONS_
#define FUNCTIONS_


namespace Foam
{
static const scalar Pi = constant::mathematical::pi;
volScalarField CALC_DIRAC(const volScalarField& D, const scalar EPS)
{
	const fvMesh& mesh = D.mesh();
	
	volScalarField DIRAC
	(
	  IOobject
	  (
	    "DIRAC",
	    D.instance(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	  ),
	  mesh,
	  dimensionedScalar("delta",dimensionSet(0,0,0,0,0,0,0),scalar(0)),
          "zeroGradient"
	);

	DIRAC = 0.5*(scalar(1.0)+Foam::cos(Pi*D/EPS))/EPS;

	forAll(mesh.cells(),cellI)
	{
		if (fabs(D[cellI])>EPS)
		{
			DIRAC[cellI] = 0.;
		}
	}

	DIRAC.correctBoundaryConditions();
	
    	forAll(mesh.boundary(), patchi)
    	{
	        if (!D.boundaryField()[patchi].coupled())
        	{
			    const fvPatchScalarField&  pD = D.boundaryField()[patchi];
                fvPatchScalarField&  pDIRAC = const_cast<Foam::fvPatchField<double>&>(DIRAC.boundaryField()[patchi]);

        		forAll(mesh.boundary()[patchi],facei)
        		{
				if (fabs(pD[facei])>EPS)
				{
					pDIRAC[facei]  = 0.;
				}
			}
		}
        }
	
	return DIRAC;
}

volScalarField CALC_HEAVISIDE(const volScalarField& D,const scalar EPS)
{
	const fvMesh& mesh = D.mesh();
	
	volScalarField H
	(
	    IOobject
	    (
	        "H",
	        D.instance(),
	        mesh,
	        IOobject::NO_READ,
	        IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("H",dimensionSet(0,0,0,0,0,0,0),scalar(0))
	);

	H = 0.5*(scalar(1)+(D/EPS)+(scalar(1)/Pi)*Foam::sin(Pi*D/EPS));

	forAll(mesh.cells(), cellI)
	{
		if (D[cellI] < - EPS)
		{
			H[cellI] = 0.;
		}
		if (D[cellI] > EPS)
		{
			H[cellI] = 1.0;
		}
	}

	H.correctBoundaryConditions();
	
    	forAll(mesh.boundary(), patchi)
    	{
	        if (!D.boundaryField()[patchi].coupled())
        	{
			    const fvPatchScalarField& pD = D.boundaryField()[patchi];
		        fvPatchScalarField& pH = const_cast<Foam::fvPatchField<double>&>(H.boundaryField()[patchi]);

        		forAll(mesh.boundary()[patchi],facei)
        		{
				if (pD[facei] < -EPS)
				{
					pH[facei]  = 0.;
				}
				if (pD[facei] > EPS)
				{
					pH[facei] = 1.0;
				}
			}
		}
	}
	
	return H;
}

volScalarField CALC_GRADIENT(const volScalarField& D,const scalar EPS)
{
    const fvMesh& mesh = D.mesh();
    dimensionedScalar one = dimensionedScalar("one", dimensionSet(0,1,0,0,0,0,0), 1.0);
    volScalarField gradd
    (
        IOobject
        (
            "G",
            D.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mag(one*fvc::grad(D))
    );

    forAll(D,celli)
    {
        if ((D[celli]>EPS)||(D[celli]<(-1.0*EPS)))
        {
            gradd[celli] = 1;
        }
    }

    gradd.correctBoundaryConditions();
    return gradd;
}


scalar CALC_Vol(const volScalarField& d, const scalar EPS)
{
  const fvMesh& mesh = d.mesh();
  
  const cellList& cells = mesh.cells();
  const faceList& faces = mesh.faces();
  const pointField& points = mesh.points();
  
  scalar V(0);
  
  forAll(cells, cellI)
  {
	if (fabs(d[cellI]) < EPS)
	{
	  V += cells[cellI].mag(points,faces);
	}
  }
  
  return V;
}

scalar CALC_ERR(const volScalarField& d, const volScalarField& NG, const scalar EPS)
{  
  const fvMesh& mesh = d.mesh();
  
  const cellList& cells = mesh.cells();
  const faceList& faces = mesh.faces();
  const pointField& points = mesh.points();
  
  const scalar tot_Vol = CALC_Vol(d, EPS);
  
  scalar magG(0);
  
  forAll(cells, cellI)
  {
	if (fabs(d[cellI]) < EPS)
	{
	  scalar V = cells[cellI].mag(points,faces);
	  magG += V*NG[cellI]/tot_Vol;
	}
  }
 
  return fabs(1-magG);
    
}
}
#endif
