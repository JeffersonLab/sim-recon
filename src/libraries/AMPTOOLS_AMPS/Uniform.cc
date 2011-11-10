

#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>

#include "IUAmpTools/Kinematics.h"
#include "AMPTOOLS_AMPS/Uniform.h"

Uniform::Uniform( int arg) :
Amplitude()
{

}

complex< GDouble >
Uniform::calcAmplitude( GDouble** pKin ) const
{
  complex <GDouble> a(1,0);
  return a;
}

void
Uniform::updatePar( const AmpParameter& par ){
  
}

Uniform* Uniform::newAmplitude( const vector< string >& args ) const {
  int arg=0;
  return new Uniform(arg);
}

Uniform* Uniform::clone() const {
	
    int arg=0;
    return ( isDefault() ? new Uniform() :
	     new Uniform(arg) );
}

