/*---------------------------------------------------------------------------*\
|             isoAdvector | Copyright (C) 2016-2017 DHI                       |
-------------------------------------------------------------------------------

License
    This file is part of isoAdvector, which is an extension to OpenFOAM.

    IsoAdvector is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    IsoAdvector is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with IsoAdvector. If not, see <http://www.gnu.org/licenses/>.

Application
    calcAdvectErrors

Description
    Compares VOF field at a given time to an exact VOF solution for a plane,
    cylinder or sphere.

Author
    Johan Roenby, DHI, all rights reserved.

\*---------------------------------------------------------------------------*/

#include "timeSelector.H"
//#include "CLSAdvection.H"
//#include "CLSCorrection.H"
#include "fvCFD.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
//#include "isoCutCell.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    Foam::argList::noBanner();
    #include "setRootCase.H"

    // Reading user specified time(s)
    Foam::Time runTime(Foam::Time::controlDictName, args);
    instantList timeDirs = timeSelector::select0(runTime, args);
    // In case no times are specified we use latest time
/*    if (timeDirs.size() < 1)
    {
        instantList instList = runTime.times();
        timeDirs.append(instList.last());
    }*/
//    runTime.setTime(timeDirs[0], 0);
    instantList instList = runTime.times();
    runTime.setTime(instList[1],0);
//    Info<<"second="  << runTime.timeName()<<endl;
    // Reading mesh
    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"
//      Setting time to user specified time
//      runTime.setTime(timeDirs[0], 0);

    // Reading alpha1 field for user specified time
    immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

    volScalarField& alpha1(mixture.alpha1());
    volScalarField& alpha2(mixture.alpha2());

    volScalarField alpha0 = alpha1;
    const word& phase1Name=mixture.phase1Name();
    // Getting internal alpha1 field
    const scalarField& VOF_0 = alpha1.primitiveField();
    const dimensionedScalar& rho1 = mixture.rho1();
    const dimensionedScalar& rho2 = mixture.rho2();

    // Need to store rho for ddt(rho, U)
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        alpha1*rho1 + alpha2*rho2
    );
    rho.oldTime();

    scalar totalMass = 1;
    scalar totalMass0 = 1;
    const scalarField& V = mesh.V();
    totalMass0 = gSum(rho*V);
    const scalar V_0(sum(VOF_0*mesh.cellVolumes()));
    scalar V_Gamma=0;
    if (timeDirs.size() == 1)
    {
        // Setting time to second user specified time
        runTime.setTime(timeDirs[0],0);
        // Reading alpha1 field for user specified time
        volScalarField alphaT
        (
            IOobject
            (
                IOobject::groupName("alpha", phase1Name),
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        const scalarField& VOF_T = alphaT.primitiveField();
        // Calculating error measures
        const scalar V_T(sum(VOF_T*mesh.cellVolumes()));
    //    const scalar E1(sum(mag(VOF_calc - VOF_true)*mesh.cellVolumes()));
    //    const scalar E1rel(E1/V_calc);
        const scalar dV(mag(V_T - V_0));
        const scalar EV(dV/(V_0 + SMALL));

        rho == alphaT*rho1 + (1.0 - alphaT)*rho2;
//        const scalarField& V = mesh.V();
        totalMass = gSum(rho*V)/totalMass0;
        V_Gamma = 0;
        forAll(VOF_T,celli)
        {
            if (VOF_T[celli]>0.01 && VOF_T[celli]<0.99)
            {
                V_Gamma+=VOF_T[celli]*V[celli];
            }
        }
        const scalar ES(V_Gamma/V_T);
        const scalar aMin(min(VOF_T));
        const scalar aMaxMinus1(max(VOF_T) - 1);
        // Printing error measures
        Info<< "Advection errors for time: " << runTime.time().value() << endl;
        Info<< "dV       = " << dV << endl;
        Info<< "EV    = " << EV << endl;
        Info<< "aMin     = " << aMin << endl;
        Info<< "aMax - 1 = " << aMaxMinus1 << endl;
        Info<< "totalMass = "<<totalMass<<endl;
        Info<<"ES = "<<ES<<endl;
    }
    else
    {
        autoPtr<OFstream> errorFilePtr;
        forAll(timeDirs,i)
        {
            if (i!=0)
            {
                runTime.setTime(timeDirs[i],0);
    //            Info<<"second="  << runTime.timeName()<<endl;
                volScalarField alphaT
                (
                    IOobject
                    (
                        IOobject::groupName("alpha", phase1Name),
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh
                );
                const scalarField& VOF_T = alphaT.primitiveField();
                // Calculating error measures
                const scalar V_0(sum(VOF_0*mesh.cellVolumes()));
                const scalar V_T(sum(VOF_T*mesh.cellVolumes()));
            //    const scalar E1(sum(mag(VOF_calc - VOF_true)*mesh.cellVolumes()));
            //    const scalar E1rel(E1/V_calc);
                const scalar dV(mag(V_T - V_0));
                const scalar EV(dV/(V_0 + SMALL));
                rho == alphaT*rho1 + (1.0 - alphaT)*rho2;
                totalMass = gSum(rho*V)/totalMass0;
                V_Gamma = 0;
                forAll(VOF_T,celli)
                {
                    if (VOF_T[celli]>0.01 && VOF_T[celli]<0.99)
                    {
                        V_Gamma+=VOF_T[celli]*V[celli];
                    }
                }
                const scalar ES(V_Gamma/V_T);
                const scalar aMin(min(VOF_T));
                const scalar aMaxMinus1(max(VOF_T) - 1);
      // Create the residual file if not already created
                if (errorFilePtr.empty())
                {
                    Info<< "Creating error file." << endl;
                 // File update
                    if (Pstream::master())
                    {
                        fileName errors;
                        word name_ = "errors";
                        // Open new file at start up
                        errorFilePtr.reset(new OFstream("errors.dat"));
                    }
                }
                if (Pstream::master())
                {
                    errorFilePtr() << runTime.timeName()<<tab<<totalMass<<tab<<dV<<tab<<EV<<tab<<ES<<tab<<aMin<<tab<<aMaxMinus1<< endl;
                }
            }
        }
    }
    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //