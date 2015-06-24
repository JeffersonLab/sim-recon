#if !defined(DALITZDECAYFACTORY)
#define DALITZDECAYFACTORY

#include <vector>

#include "TLorentzVector.h"

using namespace std;

class DalitzDecayFactory
{
	
public:
	
	DalitzDecayFactory( double parentMass, const vector<double>& childMass );
	
	vector<TLorentzVector> generateDecay() const;
	
private:
	
	static const double kPi;
	
	double cmMomentum( double M, double m1, double m2 ) const;
	double random( double low, double hi ) const;
	
	double m_parentMass;
	vector<double> m_childMass;
	
	double m_maxLorentzFactor;
};

#endif
