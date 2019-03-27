scalarField& alpha1_ = alpha1.ref();
scalarField& f = psi.ref();
forAll(alpha1_,celli)
{
    if (alpha1_[celli]>0 || alpha1_[celli]<1)
    {
        scalar t = 0;
        scalar a = 0;
        scalar res = scalar(0.5)-alpha1_[celli];
        int maxIter = 100;
        int nIter = 0;
        scalar tol = 1e-8;
        while (res > tol && nIter < 10*maxIter)
        {
//            Info <<"t="<<t<<endl;
            t-=(pow3(t)-double(3.0)*t-double(2.0)+double(4.0)*alpha1_[celli])
                /(double(3.0)*sqr(t)-double(3.0));
            a = pow3(t)-double(3.0)*t-double(2.0)+double(4.0)*alpha1_[celli];
            res = mag(a);
            nIter++;
        }
        Info << "res = "<< res<<"t="<<t <<endl;
        f[celli] = 0.5*t*epsilon;
    }
    else
        f[celli]=(alpha1_[celli] - double(0.5))*epsilon;
}