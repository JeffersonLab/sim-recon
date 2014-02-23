#if !(defined THREEPIPLOTGENERATOR)
#define THREEPIPLOTGENERATOR

#include <vector>
#include <string>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class FitResults;
class Kinematics;

class ThreePiPlotGenerator : public PlotGenerator
{
    
public:
  
  // create an index for different histograms
  enum { k3PiMass = 0, kPiMPiP1Mass, kPiMPiP2Mass, kPiP1PiP2Mass,
         kAlpha, kCosThetaRes, kPhiRes, kNumHists };
  
  ThreePiPlotGenerator( const FitResults& results );
    
private:
        
  void projectEvent( Kinematics* kin );
  
};

#endif
