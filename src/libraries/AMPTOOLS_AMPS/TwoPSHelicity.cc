
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/TwoPSHelicity.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"

TwoPSHelicity::TwoPSHelicity( const vector< string >& args ) :
UserAmplitude< TwoPSHelicity >( args )
{
  assert( args.size() == 3 );
	
	m_j = atoi( args[0].c_str() );
	m_m = atoi( args[1].c_str() );
	m_e = atoi( args[2].c_str() );

  // make sure values are reasonable
  assert( abs( m_e ) == 1 );
	assert( m_m <= m_j );
  
  if( m_m == 0 ) m_bigTheta = 0.5;
	if( m_m  > 0 ) m_bigTheta = sqrt( 0.5 );
	if( m_m  < 0 ) m_bigTheta = 0;
  
	// the "reflectivity factor" is e*(-1)^(m)
	m_reflectivityFactor = ( m_m % 2 == 0 ? m_e : -m_e );
}


complex< GDouble >
TwoPSHelicity::calcAmplitude( GDouble** pKin ) const {
  
  TLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 
  TLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 
  TLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] ); 
  TLorentzVector p2     ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] ); 
  
  TLorentzVector resonance = p1 + p2;
  
  TLorentzRotation resRestBoost( -resonance.BoostVector() );
  
  TLorentzVector beam_res   = resRestBoost * beam;
  TLorentzVector recoil_res = resRestBoost * recoil;
  TLorentzVector p1_res = resRestBoost * p1;
  
  // helicity frame: z-axis is propagation of resonance X => opposite recoil proton in X rest frame
  TVector3 z = -1. * recoil_res.Vect().Unit();

  // y axis perpendicular to production plane
  TVector3 y = beam_res.Vect().Cross(z).Unit();

  TVector3 x = y.Cross(z);
  
  TVector3 angles( (p1_res.Vect()).Dot(x),
                   (p1_res.Vect()).Dot(y),
                   (p1_res.Vect()).Dot(z) );
  
  GDouble cosTheta = angles.CosTheta();
  GDouble phi = angles.Phi();
  
  GDouble coef = sqrt( ( 2. * m_j + 1 ) / ( 4 * 3.1416 ) );
  
  return complex< GDouble >( coef * m_bigTheta * 
                            ( wignerD( m_j, m_m, 0, cosTheta, phi ) - 
                             static_cast< GDouble>( m_reflectivityFactor ) * 
                             wignerD( m_j, -m_m, 0, cosTheta, phi ) ) );  
}

#ifdef GPU_ACCELERATION
void
TwoPSHelicity::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

  GPUTwoPSHelicity_exec( dimGrid,  dimBlock, GPU_AMP_ARGS,
                       m_j, m_m, m_bigTheta, 
                       static_cast< GDouble >( m_reflectivityFactor ) );
}
#endif //GPU_ACCELERATION
