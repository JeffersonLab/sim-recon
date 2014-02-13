
#include <vector>
#include <stdlib.h>
#include <assert.h>

#include "AMPTOOLS_MCGEN/TwoBodyDecayFactory.h"

#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/LorentzRotation.h"
#include "CLHEP/Vector/ThreeVector.h"

using namespace CLHEP;

const double TwoBodyDecayFactory::kPi = 3.14159;

TwoBodyDecayFactory::TwoBodyDecayFactory( double parentMass, const vector<double>& childMass ) :
m_parentMass( parentMass ),
m_childMass( childMass )
{
	assert( childMass.size() == 2 );
}


vector<HepLorentzVector>
TwoBodyDecayFactory::generateDecay() const {
	
	vector<HepLorentzVector> child( 2 );
	vector<Hep3Vector> childMom( 2 );
	
    // let the X decay to isobar + 2 in the X CM
    // fill the isobar momentum vector and the bachelor momentum vector
    childMom[0].
        setRThetaPhi( cmMomentum( m_parentMass, m_childMass[0], m_childMass[1] ),
                      acos( random( -0.999999, 0.999999 ) ),
                      random( -kPi, kPi ) );
    childMom[1] = -childMom[0];
			
	// fill the final four-vectors
	for( int i = 0; i < 2; ++i ){
		
		child[i].setVect( childMom[i] );
		child[i].setE( sqrt( childMom[i].mag2() + 
							 m_childMass[i] * m_childMass[i] ) );
	}

	return child;
}

double
TwoBodyDecayFactory::cmMomentum( double M, double m1, double m2 ) const {
	
	// mini PDG Eq: 38.16
	
	double num1 = ( M * M - ( m1 + m2 ) * ( m1 + m2 ) );
	double num2 = ( M * M - ( m1 - m2 ) * ( m1 - m2 ) );
	
	return( sqrt( num1 * num2 ) / ( 2 * M ) );
}

double
TwoBodyDecayFactory::random( double low, double hi ) const {
	
	return( ( hi - low ) * drand48() + low );
}
