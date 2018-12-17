    const labelListList pcs = mesh.pointCells();
    const labelListList cps = mesh.cellPoints();
    scalarField band (mesh.V().size(),2);
//	const pointField& pts = mesh.points();
//    const cell& Cells = mesh.cells();
    forAll(alpha1,celli)
    {
        if (alpha1[celli]>1e-6 && alpha1[celli]<0.999999)
        {
            band[celli] = 0;
//            Info <<"cellI="<<celli<<endl;
            forAll(cps[celli],pi)
            {
                const label& cp = cps[celli][pi];
                const labelList& pc=pcs[cp];
                forAll(pc, i)
                {
//                    Info <<"pcI="<<pc[i]<<endl;
                    if(pc[i]!=celli)
                    {
                        const label& l= pc[i];
                        if (alpha1[l]>1e-6 &&alpha1[l]<0.999999)
                        {
                            band[l] = 0;
                        }
                        else
                        {
                            band[l] = 1;
                        }
                    }
                }
            }
        }
    }


