
#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>

#include "AMPTOOLS_AMPS/polCoef.h"

polCoef::polCoef( int polBeam, const AmpParameter& polFrac ) :
  Amplitude(),
  m_polBeam( polBeam ),  // beam polarization component (X=0, Y=1)
  m_polFrac( polFrac )  // fraction of polarization 0=0% 1=100%.
  
{
  assert( ( polBeam == 0 ) || ( polBeam == 1 ) );
}


complex< GDouble >
polCoef::calcAmplitude( GDouble** pKin ) const
{
  int pol=(m_polBeam==1 ? +1 : -1); // y and x-pol. respectively

  //(1+g) for x-pol, (1-g) for y-pol
  return complex<GDouble>((GDouble)sqrt((1.0-pol*m_polFrac)*0.5), 0);
}


//void polCoef::updatePar( const AmpParameter& par ){}


polCoef* polCoef::newAmplitude( const vector< string >& args ) const 
{
  int polBeam = atoi( args[0].c_str() );
  AmpParameter polFrac( args[1] );

  return new polCoef( polBeam, polFrac);
}

polCoef* polCoef::clone() const {
	
    return ( isDefault() ? new polCoef() :
	     new polCoef( m_polBeam, m_polFrac) );
}

