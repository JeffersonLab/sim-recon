#if !defined(RESONANCEDECAYFACTORY)
#define RESONANCEDECAYFACTORY

#include <vector>

#include "TLorentzVector.h"
#include "TRandom.h"

using namespace std;

class ResonanceDecayFactory
{
    
public:
    
    ResonanceDecayFactory( double resMass, double isoMass, double isoWidth, double bachMass );
    
    vector< TLorentzVector > generateDecay() const;
   
private:
         
    static const double kPi;

    double cmMomentum( double M, double m1, double m2 ) const;
    double random( double low, double hi ) const;
        
    double m_resMass;
    double m_isoMass;
    double m_isoWidth;
    double m_bachMass;
  
    mutable TRandom m_randGen;
};

#endif
