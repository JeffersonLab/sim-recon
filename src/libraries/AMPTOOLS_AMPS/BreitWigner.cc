

#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>

#include "barrierFactor.h"
#include "breakupMomentum.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"
#include "CLHEP/Vector/LorentzVector.h"

BreitWigner::BreitWigner( const AmpParameter& mass0, 
                          const AmpParameter& width0, int orbitL,
                          pair<string,string> daughters) :
Amplitude(),
m_mass0( mass0 ),
m_width0( width0 ),
m_orbitL( orbitL ),
m_daughters( daughters )
{
  
  // this is not the default constructor
  setDefaultStatus( false );
  
  // need to register any free parameters so the framework knows about them
  registerParameter( m_mass0 );
  registerParameter( m_width0 );
  
  // make sure the input variables look reasonable
  assert((orbitL >= 0) && (orbitL <= 4));
}

complex< GDouble >
BreitWigner::calcAmplitude( GDouble** pKin ) const
{
  HepLorentzVector P1, P2, Ptot, Ptemp;
  
  for( unsigned int i = 0; i < m_daughters.first.size(); ++i ){
    
    string num; num += m_daughters.first[i];
    int index = atoi(num.c_str());
    Ptemp.set( pKin[index][1], pKin[index][2], 
               pKin[index][3], pKin[index][0] );                
    P1 += Ptemp;
    Ptot += Ptemp;
  }
  
  for( unsigned int i = 0; i < m_daughters.second.size(); ++i ){
    
    string num; num += m_daughters.second[i];
    int index = atoi(num.c_str());
    Ptemp.set( pKin[index][1], pKin[index][2], 
               pKin[index][3], pKin[index][0] );                
    P2 += Ptemp;
    Ptot += Ptemp;
  }
  
  GDouble mass  = Ptot.m();
  GDouble mass1 = P1.m();
  GDouble mass2 = P2.m();
  
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

BreitWigner*
BreitWigner::newAmplitude( const vector< string >& args ) const {
	
	assert( args.size() == 5 );
	
	AmpParameter mass( args[0] );
	AmpParameter width( args[1] );
	int orbitL = atoi( args[2].c_str() );
	pair< string, string > daughters( args[3], args[4] );
	
	return new BreitWigner( mass, width, orbitL, daughters );
}

BreitWigner*
BreitWigner::clone() const {
	
	return ( isDefault() ? new BreitWigner() : 
          new BreitWigner( m_mass0, m_width0, m_orbitL, m_daughters ) );
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

