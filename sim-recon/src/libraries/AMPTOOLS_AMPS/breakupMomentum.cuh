#ifndef CUDA_BREAKUPMOMENTUM
#define CUDA_BREAKUPMOMENTUM

#include "GPUManager/GPUCustomTypes.h"

// mass0 = mass of parent
// mass1 = mass of first daughter
// mass2 = mass of second daughter

static __device__ GDouble 
breakupMomentum( GDouble mass0, GDouble mass1, GDouble mass2 ){
	
	// fabs -- correct?  consistent w/ previous E852 code
  return G_SQRT( G_FABS( mass0*mass0*mass0*mass0 + 
                         mass1*mass1*mass1*mass1 +
                         mass2*mass2*mass2*mass2 -
                         2.0*mass0*mass0*mass1*mass1 -
                         2.0*mass0*mass0*mass2*mass2 -
                         2.0*mass1*mass1*mass2*mass2  ) ) / (2.0 * mass0);
	
}

#endif
