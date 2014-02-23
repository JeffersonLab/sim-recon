#ifndef CUDA_BARRIERFACTOR
#define CUDA_BARRIERFACTOR

#include "GPUManager/GPUCustomTypes.h"
#include "AMPTOOLS_AMPS/breakupMomentum.cuh"

// q     = breakup momentum
// spin  = angular momentum of the decay

static __device__ GDouble 
barrierFactor ( GDouble q, int spin ){

        GDouble barrier;
        GDouble z;

        z = ( (q*q) / (0.1973*0.1973) );

        switch (spin){

          case 0:
             barrier = 1.0;
             break;

          case 1:
             barrier = G_SQRT( (2.0*z) /
                              (z + 1.0) );
             break;
        
          case 2:
             barrier = G_SQRT( (13.0*z*z) /
                            ((z-3.0)*(z-3.0) + 9.0*z) );
             break;

          case 3:
             barrier = G_SQRT( (277.0*z*z*z) /
                            (z*(z-15.0)*(z-15.0) + 
                             9.0*(2.0*z-5.0)*(2.0*z-5.0)) );
             break;

          case 4:
             barrier = G_SQRT( (12746.0*z*z*z*z) /
                            ((z*z-45.0*z+105.0)*(z*z-45.0*z+105.0) +
                             25.0*z*(2.0*z-21.0)*(2.0*z-21.0)) );
             break;

          default:
             barrier = 0.0;
        }

        return barrier;
	
}

// mass0 = mass of parent
// spin  = angular momentum of the decay
// mass1 = mass of first daughter
// mass2 = mass of second daughter

static __device__ GDouble 
barrierFactor( GDouble mass0, int spin, GDouble mass1, GDouble mass2 ){
	
        GDouble q;

        q = breakupMomentum(mass0, mass1, mass2);

        return barrierFactor( q, spin );
}

#endif
