#if !(defined THREEPIPLOTGENERATORSCHILLING)
#define THREEPIPLOTGENERATORSCHILLING

#include <vector>
#include <string>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class FitResults;
class Kinematics;

class ThreePiPlotGeneratorSchilling : public PlotGenerator
{
    
public:
  
  // create an index for different histograms
  enum { k3PiMass = 0, kCosThetaPiMinus, kCosThetaPiPlus, kCosThetaPi0, kPhiPiPlus, kPhiPiMinus, kPhiPi0, kCosTheta, kPhi, kphi, kPsi, kt, kNumHists};
  
  ThreePiPlotGeneratorSchilling( const FitResults& results );
    
private:
        
  void projectEvent( Kinematics* kin );
  
};

#endif
