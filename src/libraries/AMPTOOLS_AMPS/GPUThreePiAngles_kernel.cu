/*
 *  GPUThreePiAngles_kernel.cu
 *  GlueXTools
 *
 *  Created by Matthew Shepherd on 6/16/10.
 *  Copyright 2010 Home. All rights reserved.
 *
 */


#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

#include "AMPTOOLS_AMPS/breakupMomentum.cuh"

#include "GPUUtils/lorentzBoost.cuh"
#include "GPUUtils/threeVector.cuh"
#include "GPUUtils/wignerD.cuh"
#include "GPUUtils/clebsch.cuh"

__global__ void
GPUThreePiAngles_kernel( GPU_AMP_PROTO, int polBeam, GDouble polFrac, int jX, 
                         int parX, int iX, int lX, int jI, int iI, int iZ0, 
                         int iZ1, int iZ2 ){

	int iEvent = GPU_THIS_EVENT;
  
  GDouble beam[4]   = GPU_P4(0);
  GDouble recoil[4] = GPU_P4(1);
	GDouble p1[4]     = GPU_P4(2);
	GDouble p2[4]     = GPU_P4(3);
  GDouble p3[4]     = GPU_P4(4);

  GDouble alpha = phi( &(recoil[1]) );

  GDouble res[4];
  GDouble iso[4];
  	
  for( int i = 0; i < 4; ++i ){ 
    
    iso[i] = p1[i] + p2[i];
    res[i] = iso[i] + p3[i];
  }
  
  GDouble resMass = G_SQRT(res[0]*res[0]-res[1]*res[1]-res[2]*res[2]-res[3]*res[3]);
  GDouble isoMass = G_SQRT(iso[0]*iso[0]-iso[1]*iso[1]-iso[2]*iso[2]-iso[3]*iso[3]);
  GDouble p1Mass = G_SQRT(p1[0]*p1[0]-p1[1]*p1[1]-p1[2]*p1[2]-p1[3]*p1[3]);
  GDouble p2Mass = G_SQRT(p2[0]*p2[0]-p2[1]*p2[1]-p2[2]*p2[2]-p2[3]*p2[3]);
  GDouble p3Mass = G_SQRT(p3[0]*p3[0]-p3[1]*p3[1]-p3[2]*p3[2]-p3[3]*p3[3]);

  GDouble k = breakupMomentum( resMass, isoMass, p3Mass );
  GDouble q = breakupMomentum( isoMass, p1Mass, p2Mass );

  boostToRest( beam   , res );
  boostToRest( recoil , res );
  boostToRest( iso    , res ); 
  boostToRest( p1     , res );

  // now beam, recoil, iso, and p1 are at rest in the resonance frame

  // create the z axis in this frame
  GDouble zRes[3] = { -recoil[1], -recoil[2], -recoil[3] };
  makeUnit( zRes );
  
  // create the y axis from the cross product of the beam with z
  GDouble yRes[3] = { beam[1], beam[2], beam[3] };
  cross( yRes, zRes );
  makeUnit( yRes );
  
  // defines x and replaces it with the cross product
  // of y and z
  GDouble xRes[3] = { yRes[0], yRes[1], yRes[2] };
  cross( xRes, zRes );

  // rewrite the isobar direction in this coordinate system
  GDouble angRes[3] = { dot( &(iso[1]), xRes ), 
                        dot( &(iso[1]), yRes ),
                        dot( &(iso[1]), zRes )  };

  // and record the angles
  GDouble cosThRes  = cosTheta( angRes );
  GDouble phiAngRes = phi( angRes );
  
  boostToRest( p1 , iso );
  
  GDouble angIso[3] = { dot( &(p1[1]), xRes ),
                        dot( &(p1[1]), yRes ),
                        dot( &(p1[1]), zRes ) };
                        
  GDouble cosThIso = cosTheta( angIso );
  GDouble phiAngIso = phi( angIso );

  WCUComplex i = { 0, 1 };
  WCUComplex one = { 1, 0 };
  WCUComplex ans = { 0, 0 };
  
  // a prefactor the matrix elements that couple negative helicity
  // photons to the final state
   WCUComplex negResHelProd = ( polBeam == 0 ? 
    ( one * G_COS( 2 * alpha ) + i * G_SIN( 2 * alpha ) ) :
    ( one * G_COS( 2 * alpha ) + i * G_SIN( 2 * alpha ) ) * -1 );
  negResHelProd *= ( jX % 2 == 0 ? -parX : parX );
 
  // in general we also need a sum over resonance helicities here
  // however, we assume a production mechanism that only produces
  // resonance helicities +-1
 
  for( int mL = -lX; mL <= lX; ++mL ){
    
    WCUComplex term = { 0, 0 };
    
    for( int mI = -jI; mI <= jI; ++mI ){
      
        // CAREFUL!! ordering of arguments for GPU routine clebsch
        // is different from CPU routine clebschGordan
                                
        term += Y( jI, mI, cosThIso, phiAngIso ) *
        ( negResHelProd * clebsch( jI, mI, lX, mL, jX, -1 ) + 
          clebsch( jI, mI, lX, mL, jX,  1 ) );
    }
    
    term *= Y( lX, mL, cosThRes, phiAngRes );
    ans += term;
  }
  
  ans *= ( polBeam == 0 ? ( 1 + polFrac ) / 4 : ( 1 - polFrac ) / 4 );
  
  pcDevAmp[iEvent] = ans * 
      clebsch( 1, iZ0, 1, iZ1, iI, iZ0 + iZ1 ) *
      clebsch( iI, iZ0 + iZ1, 1, iZ2, iX, iZ0 + iZ1 + iZ2 ) *
      G_POW( k, lX ) * G_POW( q, jI );
}

void
GPUThreePiAngles_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
                       int polBeam, GDouble polFrac, int jX, int parX, int iX, 
                       int lX, int jI, int iI, int iZ0, int iZ1, int iZ2 )  
{  
  GPUThreePiAngles_kernel<<< dimGrid, dimBlock >>>
     ( GPU_AMP_ARGS, polBeam, polFrac, jX, parX, iX, lX, jI, iI, iZ0, iZ1, iZ2 );
}

