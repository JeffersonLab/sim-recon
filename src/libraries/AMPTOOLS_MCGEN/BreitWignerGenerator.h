#if !defined(BREITWIGNERGENERATOR)
#define BREITWIGNERGENERATOR

#include <utility>

using namespace std;

class BreitWignerGenerator
{
        
    public:
        
        BreitWignerGenerator();
        
        BreitWignerGenerator( double mass, double width );
        
        // output of the generation is a pair of doubles
        // the first is the mass and the second is the weight
        // to apply to this event to get back phase space
        pair< double, double > operator()() const;
    
        // returns the value of the PDF for some value of s
        double pdf( double s ) const;
        
    private:
    
        double random( double low, double hi ) const;
        
        static const double kPi;
    
        double m_mass;
        double m_width;
        bool m_flatMass;
};

#endif
