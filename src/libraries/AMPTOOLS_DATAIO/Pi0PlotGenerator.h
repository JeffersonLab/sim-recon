#if !(defined PI0PLOTGENERATOR)
#define PI0PLOTGENERATOR

#include <vector>
#include <string>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class FitResults;
class Kinematics;

class Pi0PlotGenerator : public PlotGenerator
{
    
public:
  
  // create an index for different histograms
  enum { kCosTheta = 0, kPhi, kt, kCosTheta_phi, kt_phi, kNumHists};
  
  Pi0PlotGenerator( const FitResults& results );
  Pi0PlotGenerator( );
  
  void projectEvent( Kinematics* kin );
  
private:
  
  void createHistograms( );
  
};

#endif
