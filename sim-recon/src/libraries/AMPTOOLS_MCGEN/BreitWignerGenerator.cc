
#include <math.h>
#include <cstdlib>
#include <cassert>

#include "AMPTOOLS_MCGEN/BreitWignerGenerator.h"

const double BreitWignerGenerator::kPi = 3.14159;

BreitWignerGenerator::BreitWignerGenerator() :
m_mass( 0 ),
m_width( 0 )
{}

BreitWignerGenerator::BreitWignerGenerator( double mass, double width ) :
m_mass( mass ),
m_width( width )
{}

pair< double, double >
BreitWignerGenerator::operator()() const
{
    
    // generate BW's with 100% efficiency by integrating
    // the normalized BW distribution from -infty to rho'
    // the value of this func at rho' ranges then from 0 -> 1
    // transform so rho = pi/2 - rho' * pi
    // rho ranges from -pi/2 -> pi/2
    // throw rho uniformly in this range and then
    // invert to get s
    // weight is the reciprocal of the normalized PDF and can be
    // used to reweight the events so they are flat in s
    // (i.e. flat in two-body phase space)
    
    double s = -1;
    double rho;
    
    assert( m_mass > 0 && m_width > 0 );
    
    // avoid potential funny business at extreme values of rho
    while( s < 0 ){
    
        rho = random( -kPi/2, kPi/2 );
        s = m_mass * m_mass + m_mass * m_width * tan( rho );
    }
    
    double weight = 1 / pdf( s );
    
    return pair< double, double >( sqrt(s), weight );
}

double
BreitWignerGenerator::pdf( double s ) const {
    
    return m_mass * m_width /
    (  kPi * ( ( s - m_mass * m_mass ) * 
               ( s - m_mass * m_mass ) +
              m_mass * m_mass * m_width * m_width ) );
}

double
BreitWignerGenerator::random( double low, double hi ) const {
	
	return( ( hi - low ) * drand48() + low );
}
