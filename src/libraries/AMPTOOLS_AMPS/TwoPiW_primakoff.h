#if !defined(TWOPIW_PRIMAKOFF)
#define TWOPIW_PRIMAKOFF

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "IUAmpTools/UserAmplitude.h"
#include "GPUManager/GPUCustomTypes.h"

#include <utility>
#include <string>
#include <complex>
#include <vector>

using std::complex;
using namespace std;

class Kinematics;

class TwoPiW_primakoff : public UserAmplitude< TwoPiW_primakoff >
{
  
public:
	
	TwoPiW_primakoff() : UserAmplitude< TwoPiW_primakoff >() {}
	TwoPiW_primakoff( const vector< string >& args );
	
	~TwoPiW_primakoff(){}
  
	string name() const { return "TwoPiW_primakoff"; }
  
  complex< GDouble > calcAmplitude( GDouble** pKin ) const;
	  
  void updatePar( const AmpParameter& par );
    
#ifdef GPU_ACCELERATION

  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;

	bool isGPUEnabled() const { return true; }

#endif // GPU_ACCELERATION
  
private:
	
  AmpParameter m_par1;    // for the moment assume W cross section has 2 parameters
  AmpParameter m_par2;
  
  pair< string, string > m_daughters;  
};

#endif
