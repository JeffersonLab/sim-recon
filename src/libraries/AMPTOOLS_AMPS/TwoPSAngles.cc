
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/TwoPSAngles.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/LorentzRotation.h"
#include "CLHEP/Vector/ThreeVector.h"

TwoPSAngles::TwoPSAngles( int j, int m, int e ) :
Amplitude(),
m_j( j ),
m_m( m ),
m_e( e )
{
  // make sure values are reasonable
  assert( abs( e ) == 1 );
	assert( m <= j );
  
  if( m_m == 0 ) m_bigTheta = 0.5;
	if( m_m  > 0 ) m_bigTheta = sqrt( 0.5 );
	if( m_m  < 0 ) m_bigTheta = 0;
  
	// the "reflectivity factor" is e*(-1)^(m)
	m_reflectivityFactor = ( m_m % 2 == 0 ? m_e : -m_e );
}


complex< GDouble >
TwoPSAngles::calcAmplitude( GDouble** pKin ) const {
  
  HepLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 
  HepLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 
  HepLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] ); 
  HepLorentzVector p2     ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] ); 
  
  HepLorentzVector resonance = p1 + p2;
  
  HepLorentzRotation resRestBoost( -resonance.boostVector() );
  
  HepLorentzVector beam_res   = resRestBoost * beam;
  HepLorentzVector recoil_res = resRestBoost * recoil;
  HepLorentzVector p1_res = resRestBoost * p1;
  
  Hep3Vector z = beam_res.vect().unit();
  Hep3Vector y = recoil_res.vect().cross(z).unit();
  Hep3Vector x = y.cross(z);
  
  Hep3Vector angles( (p1_res.vect()).dot(x),
                    (p1_res.vect()).dot(y),
                    (p1_res.vect()).dot(z) );
  
  GDouble cosTheta = angles.cosTheta();
  GDouble phi = angles.phi();
  
  GDouble coef = sqrt( ( 2. * m_j + 1 ) / ( 4 * 3.1416 ) );
  
  return complex< GDouble >( coef * m_bigTheta * 
                            ( wignerD( m_j, m_m, 0, cosTheta, phi ) - 
                             static_cast< GDouble>( m_reflectivityFactor ) * 
                             wignerD( m_j, -m_m, 0, cosTheta, phi ) ) );  
}

TwoPSAngles*
TwoPSAngles::newAmplitude( const vector< string >& args ) const {
  
	assert( args.size() == 3 );
	
	int j = atoi( args[0].c_str() );
	int m = atoi( args[1].c_str() );
	int e = atoi( args[2].c_str() );
	
	return new TwoPSAngles( j, m, e );
}

TwoPSAngles*
TwoPSAngles::clone() const {
  
	return ( isDefault() ? new TwoPSAngles() : 
          new TwoPSAngles( m_j, m_m, m_e ) );
}

#ifdef GPU_ACCELERATION
void
TwoPSAngles::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

  GPUTwoPSAngles_exec( dimGrid,  dimBlock, GPU_AMP_ARGS,
                       m_j, m_m, m_bigTheta, 
                       static_cast< GDouble >( m_reflectivityFactor ) );
}
#endif //GPU_ACCELERATION
