
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

#include "AMPTOOLS_AMPS/breakupMomentum.cuh"
#include "AMPTOOLS_AMPS/barrierFactor.cuh"

__global__ void
GPUBreitWigner_kernel( GPU_AMP_PROTO, GDouble mass0, GDouble width0, 
                       GDouble orbitL, int daught1, int daught2 ){

	int iEvent = GPU_THIS_EVENT;

  // decode the list of daughter indices in each integer
  // be careful to handle the zero index correctly

  GDouble p1[4] = { 0, 0, 0, 0 };
  if( daught1 == 0 ){
  
    GDouble tmp[4] = GPU_P4(0);
    for( int i = 0; i < 4; ++i ) p1[i] += tmp[i];
  } 
  else if( daught1 > 0 ){
  
    while( daught1 > 0 ){
    
      int ind = daught1 % 10;
      GDouble tmp[4] = GPU_P4(ind);
      for( int i = 0; i < 4; ++i ) p1[i] += tmp[i];
      daught1 /= 10;
    }
  }
  
  
  GDouble p2[4] = { 0, 0, 0, 0 };
  if( daught2 == 0 ){
  
    GDouble tmp[4] = GPU_P4(0);
    for( int i = 0; i < 4; ++i ) p2[i] += tmp[i];
  } 
  else if( daught2 > 0 ){
  
    while( daught2 > 0 ){
    
      int ind = daught2 % 10;
      GDouble tmp[4] = GPU_P4(ind);
      for( int i = 0; i < 4; ++i ) p2[i] += tmp[i];
      daught2 /= 10;
    }
  }


  GDouble mass  = SQ( p1[0] + p2[0] );
  GDouble mass1 = SQ( p1[0] );
  GDouble mass2 = SQ( p2[0] );
   
  for( int i = 1; i <= 3; ++i ){
    
    mass  -= SQ( p1[i] + p2[i] );
    mass1 -= SQ( p1[i] );
    mass2 -= SQ( p2[i] );
  }
  
  mass  = G_SQRT( mass  );
  mass1 = G_SQRT( mass1 );
  mass2 = G_SQRT( mass2 );

  GDouble q  = fabs( breakupMomentum(  mass, mass1, mass2 ) );
  GDouble q0 = fabs( breakupMomentum( mass0, mass1, mass2 ) );

  GDouble F  = barrierFactor( q,  orbitL );
  GDouble F0 = barrierFactor( q0, orbitL );
  
  GDouble width = width0*(mass0/mass)*(q/q0)*((F*F)/(F0*F0));
//  GDouble width = width0;
 
  WCUComplex bwTop = { G_SQRT( mass0 * width0 / 3.1416 ), 0 };
  WCUComplex bwBot = { SQ( mass0 ) - SQ( mass ), -1.0 * mass0 * width };

  pcDevAmp[iEvent] = ( F * bwTop / bwBot );
}


void
GPUBreitWigner_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, 
                     GDouble mass, GDouble width, int orbitL,
                     int daught1, int daught2 )
{  

  GPUBreitWigner_kernel<<< dimGrid, dimBlock >>>
    ( GPU_AMP_ARGS, mass, width, orbitL, daught1, daught2 );
}
