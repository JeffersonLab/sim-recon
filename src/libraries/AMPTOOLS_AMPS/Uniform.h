#if !defined(UNIFORM)
#define UNIFORM

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include <utility>
#include <string>
#include <complex>
#include <vector>


using std::complex;
using namespace std;

class Kinematics;

#ifdef GPU_ACCELERATION
void GPUUniform_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO);

#endif


class Uniform : public Amplitude
{  
public:
  
  Uniform() : Amplitude() { setDefaultStatus( true ); }
  Uniform(int arg);	
  ~Uniform(){}
  
  string name() const { return "Uniform"; }
  
  complex< GDouble > calcAmplitude( GDouble** pKin ) const;
  
  void updatePar( const AmpParameter& par );
  
  Uniform* newAmplitude( const vector< string >& args) const;
  
  Uniform* clone() const;
  
#ifdef GPU_ACCELERATION
  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const{
    GPUUniform_exec(dimGrid, dimBlock, GPU_AMP_ARGS);
  };
  
  bool isGPUEnabled() const { return true; }
#endif


};

#endif
