
#include <vector>
#include <stdlib.h>

#include "AMPTOOLS_MCGEN/ResonanceDecayFactory.h"

#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/LorentzRotation.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/RandBreitWigner.h"

using namespace CLHEP;

const double ResonanceDecayFactory::kPi = 3.14159;

ResonanceDecayFactory::ResonanceDecayFactory( double resMass, double isoMass, double isoWidth, double bachMass ) :
m_resMass( resMass ),
m_isoMass( isoMass ),
m_isoWidth( isoWidth ),
m_bachMass( bachMass ){}

vector< HepLorentzVector >
ResonanceDecayFactory::generateDecay() const {
    
    // initialize this high so we through a random number
    double c0Mass = 2 * m_resMass;
    
    // avoid threshold problems
    while( ( c0Mass + m_bachMass > 0.999 * m_resMass ) || ( c0Mass <= 0 ) ){
    
        c0Mass = RandBreitWigner::shoot( m_isoMass, m_isoWidth );
    }
    
    vector<HepLorentzVector> child( 2 );
	vector<Hep3Vector> childMom( 2 );

    childMom[0].setRThetaPhi( cmMomentum( m_resMass, c0Mass, m_bachMass ),
                              acos( random( -0.999999, 0.999999 ) ),
                              random( -kPi, kPi ) );
    childMom[1] = -childMom[0];
    
    child[0].setVect( childMom[0] );
    child[0].setE( sqrt( childMom[0].mag2() + c0Mass * c0Mass ) );
    
    child[1].setVect( childMom[1] );
    child[1].setE( sqrt( childMom[1].mag2() + m_bachMass * m_bachMass ) );
    
    return child;
}

double
ResonanceDecayFactory::cmMomentum( double M, double m1, double m2 ) const {
	
	// mini PDG Eq: 38.16
	
	double num1 = ( M * M - ( m1 + m2 ) * ( m1 + m2 ) );
	double num2 = ( M * M - ( m1 - m2 ) * ( m1 - m2 ) );
	
    if( ( num1 * num2 ) < 0 ){
        
        cout << "ERROR\t" << M << "\t" << m1 << "\t" << m2 << endl;
    }
    
	return( sqrt( num1 * num2 ) / ( 2 * M ) );
}

double
ResonanceDecayFactory::random( double low, double hi ) const {
	
	return( ( hi - low ) * drand48() + low );
}
