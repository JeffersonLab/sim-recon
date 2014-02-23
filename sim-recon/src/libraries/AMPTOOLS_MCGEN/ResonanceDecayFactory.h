#if !defined(RESONANCEDECAYFACTORY)
#define RESONANCEDECAYFACTORY

#include <vector>

#include "CLHEP/Vector/LorentzVector.h"

using namespace std;
using namespace CLHEP;

class ResonanceDecayFactory
{
    
public:
    
    ResonanceDecayFactory( double resMass, double isoMass, double isoWidth, double bachMass );
    
    vector< HepLorentzVector > generateDecay() const;
   
private:
         
    static const double kPi;

    double cmMomentum( double M, double m1, double m2 ) const;
    double random( double low, double hi ) const;
        
    double m_resMass;
    double m_isoMass;
    double m_isoWidth;
    double m_bachMass;
    
};

#endif
