
#include <iostream>
#include <stdlib.h>

#include "AMPTOOLS_MCGEN/ProductionMechanism.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/LorentzRotation.h"

const double ProductionMechanism::kPi = 3.14159;

using namespace std;
using namespace CLHEP;

ProductionMechanism::ProductionMechanism( Recoil recoil, Type type, double slope ) :
m_type( type ),
m_lowMass( 0 ),
m_highMass( 0 ),
m_slope( slope ),
m_lastWeight( 1. )
{	
	switch( recoil ){
      
      // I'm sure the distinction between these doesn't matter!
      
		case kProton:  m_recMass = 0.9382; break;
		case kNeutron: m_recMass = 0.9395; break;
		default:       m_recMass = 0.9382; break;
	}
}

void
ProductionMechanism::setMassRange( double low, double high ){
	
	m_lowMass = low;
	m_highMass = high;
}

void 
ProductionMechanism::setGeneratorType( Type type ){
  
  m_type = type;
}


HepLorentzVector
ProductionMechanism::produceResonance( const HepLorentzVector& beam ){
	
	HepLorentzVector target( 0, 0, 0, 0.9382 );
	
	HepLorentzRotation lab2cmBoost( -( target + beam ).boostVector() );
	HepLorentzRotation cm2labBoost( ( target + beam ).boostVector() );
	
	double cmEnergy = ( lab2cmBoost * ( target + beam ) ).e();
	double beamMomCM = cmMomentum( cmEnergy, beam.m(), target.m() );

  double exptMax = exp(-1.)/m_slope;
  
  double t, tMax, resMass, resMomCM;

  do {
    resMass = generateMass();
    resMomCM  = cmMomentum( cmEnergy, resMass, m_recMass );
  
    tMax = 4. * beamMomCM * resMomCM;
    t = random( 0, tMax ); 
  } 
	while( random( 0., exptMax ) > t*exp(-m_slope*t) );
	
	Hep3Vector resonanceMomCM;
	resonanceMomCM.setRThetaPhi( resMomCM,
                              acos( 1. - 2.*t/tMax ),
                              random( -kPi, kPi ) );
	
	HepLorentzVector resonanceCM( resonanceMomCM, 
                               sqrt( resonanceMomCM.mag2() + 
                                    resMass * resMass ) );
	
	return cm2labBoost * resonanceCM;
}

void 
ProductionMechanism::addResonance( double mass, double width, double crossSec ){
  
  m_decGen.addChannel( m_bwGen.size(), crossSec );
  m_bwGen.push_back( BreitWignerGenerator( mass, width ) );
}

double
ProductionMechanism::generateMass(){
  
  if( m_type == kFlat ) return random( m_lowMass, m_highMass );
  
  double mass = 0;
  while( mass < m_lowMass || mass > m_highMass ){
    
    unsigned int channel = m_decGen();
    pair< double, double > bw = m_bwGen[channel]();
    
    mass = bw.first;
  }
  
  double prob = 0;
  for( unsigned int i = 0; i < m_bwGen.size(); ++i ){
    
    prob += m_bwGen[i].pdf( mass * mass ) * m_decGen.getProb( i );
  }
  
  // put in the factor of mass so resulting weights can be applied to 
  // obtain a distribution that is flat in mass instead of flat in s
  // (to get weights for reweighting flat in s remove the mass)
  
  m_lastWeight = 1 / ( prob * mass );
  
  return mass;
}

double
ProductionMechanism::cmMomentum( double M, double m1, double m2 ) const {
	
	// mini PDG Eq: 38.16
	
	double num1 = ( M * M - ( m1 + m2 ) * ( m1 + m2 ) );
	double num2 = ( M * M - ( m1 - m2 ) * ( m1 - m2 ) );
	
	return( sqrt( num1 * num2 ) / ( 2 * M ) );
}

double
ProductionMechanism::random( double low, double hi ) const {
	
	return( ( hi - low ) * drand48() + low );
}


