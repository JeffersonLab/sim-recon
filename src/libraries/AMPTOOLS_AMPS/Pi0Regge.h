#if !defined(PI0REGGE)
#define PI0REGGE

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include <string>
#include <complex>
#include <vector>

#ifdef GPU_ACCELERATION
void
GPUPi0Regge_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
                     int j, int m, GDouble bigTheta, GDouble refFact );

#endif // GPU_ACCELERATION

using std::complex;
using namespace std;

class Kinematics;

class Pi0Regge : public UserAmplitude< Pi0Regge >
{
    
public:
	
	Pi0Regge() : UserAmplitude< Pi0Regge >() { };
	Pi0Regge( const vector< string >& args );
	
	string name() const { return "Pi0Regge"; }
    
	complex< GDouble > calcAmplitude( GDouble** pKin ) const;
	
#ifdef GPU_ACCELERATION
  
  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;
  
	bool isGPUEnabled() const { return true; }
  
#endif // GPU_ACCELERATION
  
private:

};

#endif
