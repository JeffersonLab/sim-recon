

#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>

#include "TLorentzVector.h"

#include "barrierFactor.h"
#include "breakupMomentum.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"

BreitWigner::BreitWigner( const vector< string >& args ) :
UserAmplitude< BreitWigner >( args )
{
  
  assert( args.size() == 5 );
	
	m_mass0 = AmpParameter( args[0] );
	m_width0 = AmpParameter( args[1] );
	m_orbitL = atoi( args[2].c_str() );
	m_daughters = pair< string, string >( args[3], args[4] );
  
  // need to register any free parameters so the framework knows about them
  registerParameter( m_mass0 );
  registerParameter( m_width0 );
  
  // make sure the input variables look reasonable
  assert( ( m_orbitL >= 0 ) && ( m_orbitL <= 4 ) );
}

complex< GDouble >
BreitWigner::calcAmplitude( GDouble** pKin ) const
{
  TLorentzVector P1, P2, Ptot, Ptemp;
  
  for( unsigned int i = 0; i < m_daughters.first.size(); ++i ){
    
    string num; num += m_daughters.first[i];
    int index = atoi(num.c_str());
    Ptemp.SetPxPyPzE( pKin[index][1], pKin[index][2],
                      pKin[index][3], pKin[index][0] );
    P1 += Ptemp;
    Ptot += Ptemp;
  }
  
  for( unsigned int i = 0; i < m_daughters.second.size(); ++i ){
    
    string num; num += m_daughters.second[i];
    int index = atoi(num.c_str());
    Ptemp.SetPxPyPzE( pKin[index][1], pKin[index][2],
                      pKin[index][3], pKin[index][0] );
    P2 += Ptemp;
    Ptot += Ptemp;
  }
  
  GDouble mass  = Ptot.M();
  GDouble mass1 = P1.M();
  GDouble mass2 = P2.M();
  
  // assert positive breakup momenta     
  GDouble q0 = fabs( breakupMomentum(m_mass0, mass1, mass2) );
  GDouble q  = fabs( breakupMomentum(mass,    mass1, mass2) );
  
  GDouble F0 = barrierFactor(q0, m_orbitL);
  GDouble F  = barrierFactor(q,  m_orbitL);
  
  GDouble width = m_width0*(m_mass0/mass)*(q/q0)*((F*F)/(F0*F0));
  //GDouble width = m_width0;
  
  // this first factor just gets normalization right for BW's that have
  // no additional s-dependence from orbital L
  complex<GDouble> bwtop( sqrt( m_mass0 * m_width0 / 3.1416 ), 0.0 );
  
  complex<GDouble> bwbottom( ( m_mass0*m_mass0 - mass*mass ) ,
                           -1.0 * ( m_mass0 * width ) );
  
  return( F * bwtop / bwbottom );
}

void
BreitWigner::updatePar( const AmpParameter& par ){
 
  // could do expensive calculations here on parameter updates
  
}

#ifdef GPU_ACCELERATION
void
BreitWigner::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {
  
  // use integers to endcode the string of daughters -- one index in each
  // decimal place
  
  int daught1 = atoi( m_daughters.first.c_str() );
  int daught2 = atoi( m_daughters.second.c_str() );
  
  GPUBreitWigner_exec( dimGrid,  dimBlock, GPU_AMP_ARGS, 
                       m_mass0, m_width0, m_orbitL, daught1, daught2 );

}
#endif //GPU_ACCELERATION

