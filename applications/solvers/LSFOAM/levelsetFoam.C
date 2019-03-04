/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    Made by Linfan Zhang, CNPE,2018

Application
    levelsetFoam

Description
    Two phase flow solver with levelset method capturing interface.

    
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "LSimmiscibleIncompressibleTwoPhaseMixture.H"
#include "subCycle.H"
#include "CrankNicolsonDdtScheme.H"
//#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"    // 后处理

    #include "setRootCase.H"    // 检查算例格式
    #include "createTime.H"     // 建立根目录
    #include "createMesh.H"     // 读取网格
    
    #include "createControl.H"  // 建立pimple算法控制
    #include "createTimeControls.H" // 读取控制时间参数
    #include "createFields.H"   // 建立量场
    #include "initContinuityErrs.H" // 建立累计误差
    #include "createFvOptions.H"    // 建立有限容积法的控制选项
    #include "correctPhi.H"         // 建立通量修正项

//    Info<<"1st:  " <<"max="<<max(d).value()<<tab<<"min="<<min(d).value()<<endl; 
//    turbulence->validate();

    if(!LTS)
    {
        #include "readTimeControls.H" //读取控制参数
        #include "CourantNo.H"  // 计算平均和最大库郎特数
        #include "setInitialDeltaT.H" // 设置变步长下的初始时间步长
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    autoPtr<OFstream> massFilePtr;

    scalar totalMass = 1;
    scalar totalMass0 = 1;
    const scalarField& V = mesh.V();
    totalMass0 = gSum(rho*V);

    Info<< "\nStarting time loop\n" << endl;
    
    while(runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"
        runTime++;

        Info<< "Time = "<< runTime.timeName() << nl <<endl;
        // --- Pressure-velocity PIMPLE corrector loop
//        while (pimple.loop())
        {
            #include "LSEqn.H"
//            Info<<"3rd:  " <<"max="<<max(d).value()<<tab<<"min="<<min(d).value()<<endl; 
            mixture.correct();
//            Info <<"really?" <<endl;
            rho == rho2+ (rho1-rho2)*H;
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }
/*
            if (pimple.turbCorr())
            {
                turbulence->correct();
            }*/
        }

        runTime.write();

        #include "writeMass.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl; 
    }

    Info<< "End\n" <<endl;

    return 0;
}
