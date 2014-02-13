
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

#include "GPUUtils/lorentzBoost.cuh"
#include "GPUUtils/threeVector.cuh"
#include "GPUUtils/wignerD.cuh"

__global__ void
GPUTwoPSAngles_kernel( GPU_AMP_PROTO, int j, int m, GDouble bigTheta, 
                       GDouble refFact ){

	int iEvent = GPU_THIS_EVENT;
  
  GDouble beam[4]   = GPU_P4(0);
  GDouble recoil[4] = GPU_P4(1);
	GDouble p1[4]     = GPU_P4(2);
	GDouble p2[4]     = GPU_P4(3);

  GDouble res[4];
  	
  for( int i = 0; i < 4; ++i ) res[i] = p1[i] + p2[i];

  boostToRest( beam   , res );
  boostToRest( recoil , res );
  boostToRest( p1     , res ); 

  GDouble z[3] = { beam[1], beam[2], beam[3] };
  makeUnit( z );
  
  GDouble y[3] = { recoil[1], recoil[2], recoil[3] };
  cross( y, z );
  makeUnit( y );
  
  // defines x and replaces it with the cross product
  // of y and z
  GDouble x[3] = { y[0], y[1], y[2] };
  cross( x, z );

  GDouble ang[3] = { dot( &(p1[1]), x ), 
                     dot( &(p1[1]), y ),
                     dot( &(p1[1]), z )  };

  GDouble cosTh  = cosTheta( ang );
  GDouble phiAng = phi( ang );    
 
  GDouble coef   = sqrt( ( 2. * j + 1 ) / ( 4 * 3.1416 ) );
       
  pcDevAmp[iEvent] = 
     ( coef * bigTheta * 
        ( wignerD( j, m, 0, cosTh, phiAng ) - 
          refFact * wignerD( j, -m, 0, cosTh, phiAng ) ) );

}

void
GPUTwoPSAngles_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
                     int j, int m, GDouble bigTheta, GDouble refFact )  
{  
  GPUTwoPSAngles_kernel<<< dimGrid, dimBlock >>>
     ( GPU_AMP_ARGS, j, m, bigTheta, refFact );
}
