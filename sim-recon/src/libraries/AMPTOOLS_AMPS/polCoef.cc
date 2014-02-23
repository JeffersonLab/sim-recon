
#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>

#include "AMPTOOLS_AMPS/polCoef.h"

polCoef::polCoef( const vector< string >& args ) :
  UserAmplitude< polCoef >( args )
{
  
  m_polBeam = atoi( args[0].c_str() );
  m_polFrac = AmpParameter( args[1] );
  
  registerParameter( m_polFrac );
  
  assert( ( m_polBeam == 0 ) || ( m_polBeam == 1 ) );
}


complex< GDouble >
polCoef::calcAmplitude( GDouble** pKin ) const
{
  int pol=(m_polBeam==1 ? +1 : -1); // y and x-pol. respectively

  //(1+g) for x-pol, (1-g) for y-pol
  return complex<GDouble>((GDouble)sqrt((1.0-pol*m_polFrac)*0.5), 0);
}
