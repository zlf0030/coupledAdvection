   
//    Info<<"reinitialization1: " <<"max="<<max(d_).value()<<tab<<"min="<<min(d_).value()<<tab<<"EPS="<<EPS<<endl;

   
//    Info<<"reinitialization2: " <<"max="<<max(d_).value()<<tab<<"min="<<min(d_).value()<<tab<<"Dtau="<<Dtau<<tab<<"NITER="<<NITER<<endl;
    volScalarField d0 = d_;
    volScalarField kk_=d_;
    for (int corr=0; corr <=NITER;corr++)
    {
        kk_ = d0/mag(d0)*(double(1)-mag(fvc::grad(d_)*one))*Dtau;
        Info<<"kk_max="<<max(kk_).value()<<tab<<"kk_min="<<min(kk_).value()<<endl;
/*        scalar mmax = max(d_).value();
        scalar mmin = min(d_).value();
        forAll(d_,celli)
        {
            if (d_[celli]>0)
                d_[celli]=d_[celli]*EPS/mmax;
            else
                d_[celli]=-1*d_[celli]*EPS/mmin;
        }*/
        d_=d_+kk_;
        d_.correctBoundaryConditions();
    }

/*        scalar mmax = max(d_).value();
        scalar mmin = min(d_).value();
        forAll(d_,celli)
        {
            if (d_[celli]>0)
                d_[celli]=d_[celli]*EPS/mmax;
            else
                d_[celli]=-1*d_[celli]*EPS/mmin;
        }
*/

//    Info<<"reinitialization3: " <<"max="<<max(d_).value()<<tab<<"min="<<min(d_).value()<<endl;

//    volScalarField Signd0 = d0/(sqrt(d0*d0+1e-08));
/*
    labelList d_L;
    forAll(d_,celli){
        if((d_[celli]<Thic)&&(d_[celli]>(-1.0*Thic))){
           d_L.append(celli);
        }
    }
*/
//    Info << "d_L="<<nl<<d_L<<endl;

/*
    for (int loop = 0; loop <= NITER; loop++)
    {
        // calculate distance function field
//        calculateGrad();
//        volVectorField Omiga = one*Signd0 * fvc::grad(d_)/G_;
        forAll(d_L,celli)
        {
//                Info<<"d_["<<celli<<"]=    "<<d_[celli]<<endl;
//                Info <<"G_["<<celli<<"]=    " <<G_[celli]<<endl;
//                Info <<"Sign["<<celli<<"]=    "<<Signd0[celli]<<endl;
                d_[celli] += Dtau*Signd0[celli]*(scalar(1)-G_[celli]);
//                Info<<"d2_["<<celli<<"]=    "<<d_[celli]<<endl;
//        d_ += Dtau*sign(d0)*(scalar(1)-one*mag(fvc::grad(d_)));
        }
    }
*/
//    scalar MAXD = max(d_).value();
//    Info<<"MAX_before="<<MAXD<<endl;
//    scalar MIND = min(d_).value();
//    Info<<"MIN_before="<<MIND<<endl;
/*
    forAll(d_,celli)
    {
        if (d_[celli]>Thic)
            {
//                Info << "d_["<<celli<<"]="<<d_[celli]<<endl;
                d_[celli]=Thic;
            }
        else if (d_[celli]<(-1.0*Thic))
            {
//                Info << "d_["<<celli<<"]="<<d_[celli]<<endl;
                d_[celli]=-1.0*Thic;
//                Info <<"modified d_="<<d_[celli]<<endl;
            }
    }

    forAll(mesh.boundary(),patchi)
    {
        fvPatchScalarField& pd_ = const_cast<Foam::fvPatchField<double>&>(d_.boundaryField()[patchi]);
        forAll(mesh.boundary()[patchi],facei)
        {
            forAll(pd_,facei)
            {
//                Info << "pd_["<<facei<<"]="<<d_[facei]<<endl;
                if (pd_[facei]>Thic)
                {
//                  Info << "pd_["<<facei<<"]="<<pd_[facei]<<endl;
                    pd_[facei]=Thic;
                }
                else if (pd_[facei]<(-1.0*Thic))
                {
//                  Info << "pd_["<<facei<<"]="<<pd_[facei]<<endl;
                    pd_[facei]=-1.0*Thic;
//                  Info <<"pmodified d_="<<pd_[facei]<<endl; 
                }
            }   
        }
    }
*/

//    Info <<"maxD_after=" <<max(d_).value()<<endl;
//    Info <<"minD_after=" <<min(d_).value()<<endl;

// volScalarField ss =mag(one*fvc::grad(d_)); 
// grad_d == ss;
//    volScalarField s = scalar(1) - mag(one*fvc::grad(d_));
//    Info<<"Dtau="<<Dtau<<"   min_s =" <<min(s).value()<<"   max_s="<<max(s).value()<<endl;

//    ERR = CALC_ERR(d_, mag(fvc::grad(d_)), EPS);
//    Info << "  Min(d) = " << min(d_).value() << "  Max(d) = " << max(d_).value()<< endl;
