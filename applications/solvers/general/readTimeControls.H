   adjustTimeStep =
        runTime.controlDict().lookupOrDefault("adjustTimeStep", false);
    
   maxCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);
    
   maxDeltaT =
        runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);