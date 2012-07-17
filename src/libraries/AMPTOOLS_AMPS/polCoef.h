#if !defined(POLCOEF)
#define POLCOEF

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
void GPUpolCoef_exec(dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
		     int polBeam, GDouble polFrac);

#endif


class polCoef : public Amplitude
{
  
public:
  
  polCoef() : Amplitude() { setDefaultStatus( true ); }
  polCoef( int polBeam, const AmpParameter& polFrac );	
  ~polCoef(){}
  
  string name() const { return "polCoef"; }
  
  complex< GDouble > calcAmplitude( GDouble** pKin=NULL) const;
  
  //void updatePar( const AmpParameter& par );
  
  polCoef* newAmplitude( const vector< string >& args) const;
  
  polCoef* clone() const;
  
private: 
  int m_polBeam;
  AmpParameter m_polFrac;

#ifdef GPU_ACCELERATION
  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const{
    GPUpolCoef_exec(dimGrid, dimBlock, GPU_AMP_ARGS, m_polBeam, m_polFrac);
  };

  bool isGPUEnabled() const { return true; }
#endif
  
};

#endif
