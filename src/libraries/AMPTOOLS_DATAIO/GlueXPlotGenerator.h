#if !(defined GLUEXPLOTGENERATOR)
#define GLUEXPLOTGENERATOR

#include <vector>
#include <string>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class AmpToolsInterface;
class Kinematics;

class GlueXPlotGenerator : public PlotGenerator
{
    
public:
  
  // create an index for different histograms
  enum { kHist1 = 0, kNumHists };
  
  GlueXPlotGenerator( AmpToolsInterface& ati );
    
private:
        
  void projectEvent( Kinematics* kin );
  
};

#endif
