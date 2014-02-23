#if !defined(DALITZDECAYFFACTORY)
#define DALITZDECAYFACTORY

#include <vector>

#include "CLHEP/Vector/LorentzVector.h"

using namespace std;
using namespace CLHEP;

class DalitzDecayFactory
{
	
public:
	
	DalitzDecayFactory( double parentMass, const vector<double>& childMass );
	
	vector<HepLorentzVector> generateDecay() const;
	
private:
	
	static const double kPi;
	
	double cmMomentum( double M, double m1, double m2 ) const;
	double random( double low, double hi ) const;
	
	double m_parentMass;
	vector<double> m_childMass;
	
	double m_maxLorentzFactor;
};

#endif
