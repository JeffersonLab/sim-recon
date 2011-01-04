
#include <stdlib.h>

#include <cassert>
#include <iostream>
#include <string>
#include <sstream>

#include "IUAmpTools/AmpParameter.h"
#include "AMPTOOLS_AMPS/ThreePiAngles.h"
#include "AMPTOOLS_AMPS/clebschGordan.h"
#include "AMPTOOLS_AMPS/wignerD.h"
#include "AMPTOOLS_AMPS/breakupMomentum.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/LorentzRotation.h"
#include "CLHEP/Vector/ThreeVector.h"

ThreePiAngles::ThreePiAngles( int polBeam, const AmpParameter& polFrac, 
                              int jX, int parX, int iX, int lX, 
                              int jI, int iI, 
                              int iZ0, int iZ1, int iZ2 ) :
Amplitude(),
m_polBeam( polBeam ), // beam polarization component (X=0, Y=1)
m_polFrac( polFrac ), // fraction of polarization 0=0% 1=100%
m_jX( jX ),     // total J of produced resonance
m_parX( parX ), // parity of produced resonance
m_iX( iX ),     // total isospin of resonance
m_lX( lX ),     // l between bachelor and isobar
m_jI( jI ),     // total J of isobar
m_iI( iI ),     // total isospin of isobar
m_iZ0( iZ0 ),   // z component of isospin of final state particle 0
m_iZ1( iZ1 ),   // z component of isospin of final state particle 1
m_iZ2( iZ2 )    // z component of isospin of final state particle 2
{

  assert( ( polBeam == 0 ) || ( polBeam == 1 ) );
  assert( ( polFrac >= 0 ) && ( polFrac <= 1 ) );
  assert( jX >= 0  );
  assert( abs( (double)parX ) == 1 );
  assert( abs( (double)iX )   <= 1 );
  assert( lX <= jX );
  assert( jI >= 0  );
  assert( abs( (double)iI )  <= 1 );
  assert( abs( (double)iZ0 ) <= 1 );
  assert( abs( (double)iZ1 ) <= 1 );
  assert( abs( (double)iZ2 ) <= 1 );
    
  registerParameter( m_polFrac );
  
  // the first two elements are the beam and recoil
  m_iZ.push_back( 0 );
  m_iZ.push_back( 0 );
  m_iZ.push_back( iZ0 );
  m_iZ.push_back( iZ1 );
  m_iZ.push_back( iZ2 );
  
	setDefaultStatus( false );
}

complex< GDouble >
ThreePiAngles::calcAmplitude( GDouble** pKin ) const
{

  HepLorentzVector beam   ( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0] ); 
  HepLorentzVector recoil ( pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0] ); 
  HepLorentzVector p1     ( pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0] ); 
  HepLorentzVector p2     ( pKin[3][1], pKin[3][2], pKin[3][3], pKin[3][0] ); 
  HepLorentzVector p3     ( pKin[4][1], pKin[4][2], pKin[4][3], pKin[4][0] ); 
  
  HepLorentzVector isobar = p1 + p2;
  HepLorentzVector resonance = isobar + p3;
  
  // orientation of production plane in lab
  GDouble alpha = recoil.vect().phi();
  
  HepLorentzRotation resRestBoost( -resonance.boostVector() );
  
  HepLorentzVector beam_res   = resRestBoost * beam;
  HepLorentzVector recoil_res = resRestBoost * recoil;
  HepLorentzVector p3_res     = resRestBoost * p3;

  Hep3Vector zRes = -recoil_res.vect().unit();
  Hep3Vector yRes = beam_res.vect().cross(zRes).unit();
  Hep3Vector xRes = yRes.cross(zRes);
  
  Hep3Vector anglesRes( (p3_res.vect()).dot(xRes),
                        (p3_res.vect()).dot(yRes),
                        (p3_res.vect()).dot(zRes) );

  GDouble cosThetaRes = anglesRes.cosTheta();
  GDouble phiRes = anglesRes.phi();

  HepLorentzRotation isoRestBoost( -isobar.boostVector() );
  HepLorentzVector p1_iso = isoRestBoost * p1;
    
  Hep3Vector anglesIso( (p1_iso.vect()).dot(xRes),
                        (p1_iso.vect()).dot(yRes),
                        (p1_iso.vect()).dot(zRes) );

  GDouble cosThetaIso = anglesIso.cosTheta();
  GDouble phiIso = anglesIso.phi();
  
  GDouble k = breakupMomentum( resonance.m(), isobar.m(), p3.m() );
  GDouble q = breakupMomentum( isobar.m(), p1.m(), p2.m() );
  
  const vector< int >& perm = getCurrentPermutation();
  
  // get the z components of isospin (charges) for the pions
  int iZ0  = m_iZ[perm[2]];
  int iZ1  = m_iZ[perm[3]];
  int iZ2  = m_iZ[perm[4]];
            
  complex< GDouble > i( 0, 1 );
  complex< GDouble > ans( 0, 0 );
 
  // a prefactor the matrix elements that couple negative helicity
  // photons to the final state
  complex< GDouble > negResHelProd = ( m_polBeam == 0 ? 
     cos( 2 * alpha ) + i * sin( 2 * alpha ) : 
    -cos( 2 * alpha ) - i * sin( 2 * alpha ) );
  negResHelProd *= ( m_jX % 2 == 0 ? -m_parX : m_parX  );
 
  // in general we also need a sum over resonance helicities here
  // however, we assume a production mechanism that only produces
  // resonance helicities +-1
  
  for( int mL = -m_lX; mL <= m_lX; ++mL ){
    
    complex< GDouble > term( 0, 0 );
    
    for( int mI = -m_jI; mI <= m_jI; ++mI ){
              
      term += Y( m_jI, mI, cosThetaIso, phiIso ) *
       ( negResHelProd * clebschGordan( m_jI, m_lX, mI, mL, m_jX, -1 ) +
         clebschGordan( m_jI, m_lX, mI, mL, m_jX, 1 ) );
    }
    
    term *= Y( m_lX, mL, cosThetaRes, phiRes );
    ans += term;
  }

  ans *= ( m_polBeam == 0 ? ( 1 + m_polFrac ) / 4 : ( 1 - m_polFrac ) / 4 );

  ans *= clebschGordan( 1, 1, iZ0, iZ1, m_iI, iZ0 + iZ1 )    *
         clebschGordan( m_iI, 1, iZ0 + iZ1, iZ2, m_iX, iZ0 + iZ1 + iZ2 ) *
         pow( k, m_lX ) * pow( q, m_jI );
    
  return ans;
}

ThreePiAngles*
ThreePiAngles::newAmplitude( const vector< string >& args ) const {

	assert( args.size() == 11 );
	
	int polBeam = atoi( args[0].c_str() );
  AmpParameter polFrac( args[1] );
	int jX      = atoi( args[2].c_str() );
  int parX    = atoi( args[3].c_str() );
	int iX      = atoi( args[4].c_str() );
	int lX      = atoi( args[5].c_str() );
	int jI      = atoi( args[6].c_str() );
	int iI      = atoi( args[7].c_str() );
  int iZ0     = atoi( args[8].c_str() );
  int iZ1     = atoi( args[9].c_str() );
  int iZ2     = atoi( args[10].c_str() );
	
	return new ThreePiAngles( polBeam, polFrac, jX, parX, iX, lX, jI, iI, iZ0, iZ1, iZ2 );
}

ThreePiAngles*
ThreePiAngles::clone() const {

	return ( isDefault() ? new ThreePiAngles() : 
			 new ThreePiAngles( m_polBeam, m_polFrac, m_jX, m_parX, m_iX, m_lX, m_jI, 
                          m_iI, m_iZ0, m_iZ1, m_iZ2 ) );
}

#ifdef GPU_ACCELERATION

void
ThreePiAngles::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {

  const vector< int >& perm = getCurrentPermutation();
  
  // get the z components of isospin (charges) for the pions
  int iZ0  = m_iZ[perm[2]];
  int iZ1  = m_iZ[perm[3]];
  int iZ2  = m_iZ[perm[4]];

  GPUThreePiAngles_exec( dimGrid, dimBlock, GPU_AMP_ARGS,
                         m_polBeam, m_polFrac, m_jX, m_parX, m_iX, m_lX, 
                         m_jI, m_iI, iZ0, iZ1, iZ2 );
  
  
}

#endif
