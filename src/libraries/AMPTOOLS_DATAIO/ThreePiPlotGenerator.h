#if !(defined THREEPIPLOTGENERATOR)
#define THREEPIPLOTGENERATOR

#include <vector>
#include <string>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class AmpToolsInterface;
class Kinematics;

class ThreePiPlotGenerator : public PlotGenerator
{
    
public:
  
  // create an index for different histograms
  enum { kHist1 = 0, kNumHists };
  
  ThreePiPlotGenerator( AmpToolsInterface& ati );
    
private:
        
  void projectEvent( Kinematics* kin );
  
};

#endif
