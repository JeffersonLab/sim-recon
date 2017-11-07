#if !(defined TWOPIPLOTGENERATOR)
#define TWOPIPLOTGENERATOR

#include <vector>
#include <string>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class FitResults;
class Kinematics;

class OmegaRadiativePlotGenerator : public PlotGenerator
{
    
public:
  
  // create an index for different histograms
  //enum { kOmegaMass = 0, kPi0CosTheta, kPi0Phi, kGammaCosTheta, kGammaPhi, kPhi, kphi, kPsi, kt, kNumHists};
  enum { kOmegaMass = 0, kCosThetaPi0, kCosThetaGamma, kPhiPi0, kPhiGamma, kCosTheta, kPhi, kphi, kPsi, kt, kNumHists};

  OmegaRadiativePlotGenerator( const FitResults& results );
    
private:
        
  void projectEvent( Kinematics* kin );
  
};

#endif
