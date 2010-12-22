#if !defined(PRODUCTIONMECHANISM)
#define PRODUCTIONMECHANISM

#include <vector>

#include "CLHEP/Vector/LorentzVector.h"

#include "AMPTOOLS_MCGEN/BreitWignerGenerator.h"
#include "AMPTOOLS_MCGEN/DecayChannelGenerator.h"

using namespace CLHEP;
using namespace std;

class ProductionMechanism
{
	
public:
	
	enum Type { kResonant, kFlat };
	enum Recoil { kProton, kNeutron };
	
	ProductionMechanism( Recoil recoil, Type type, double slope = 5.0 );
	
	void setMassRange( double low, double high );
  void setGeneratorType( Type type );
	
	HepLorentzVector produceResonance( const HepLorentzVector& beam );
  
  // there may be a better way to do this, like pair< , >
  // but sometimes the user doesn't care about the weight
  double getLastGeneratedWeight() { return m_lastWeight; }
  
  void addResonance( double mass, double width, double crossSec );
	
private:
  
	static const double kPi;
  
  double generateMass();
  
	double cmMomentum( double M, double m1, double m2 ) const;
	double random( double low, double hi ) const;
  
	Type m_type;
  
	double m_lowMass;
	double m_highMass;
	double m_slope;
  
	double m_recMass;
  
  double m_lastWeight;
  
  vector< BreitWignerGenerator > m_bwGen;    
  DecayChannelGenerator m_decGen;
};

#endif
