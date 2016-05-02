

#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>

#include "TLorentzVector.h"

#include "barrierFactor.h"
#include "breakupMomentum.h"

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/BreitWigner3body.h"

BreitWigner3body::BreitWigner3body( const vector< string >& args ) :
UserAmplitude< BreitWigner3body >( args )
{
  
  assert( args.size() == 3 );
	
	m_mass0 = AmpParameter( args[0] );
	m_width0 = AmpParameter( args[1] );
	m_daughters = args[2];
  
  // need to register any free parameters so the framework knows about them
  registerParameter( m_mass0 );
  registerParameter( m_width0 );
  
}

complex< GDouble >
BreitWigner3body::calcAmplitude( GDouble** pKin ) const
{
  TLorentzVector Ptot, Ptemp;
  
  for( unsigned int i = 0; i < m_daughters.size(); ++i ){
    
    string num; num += m_daughters[i];
    int index = atoi(num.c_str());
    Ptemp.SetPxPyPzE( pKin[index][1], pKin[index][2],
                      pKin[index][3], pKin[index][0] );
    Ptot += Ptemp;
  }
  
  GDouble mass  = Ptot.M();
  
  GDouble width = m_width0;
  //GDouble width = m_width0;
  
  // this first factor just gets normalization right for BW's that have
  // no additional s-dependence from orbital L
  complex<GDouble> bwtop( sqrt( m_mass0 * m_width0 / 3.1416 ), 0.0 );
  
  complex<GDouble> bwbottom( ( m_mass0*m_mass0 - mass*mass ) ,
                           -1.0 * ( m_mass0 * width ) );
  
  return( bwtop / bwbottom );
}

void
BreitWigner3body::updatePar( const AmpParameter& par ){
 
  // could do expensive calculations here on parameter updates
  
}

