#if !defined(BREITWIGNER3BODY)
#define BREITWIGNER3BODY

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "IUAmpTools/UserAmplitude.h"
#include "GPUManager/GPUCustomTypes.h"

#include <utility>
#include <string>
#include <complex>
#include <vector>

#ifdef GPU_ACCELERATION
void GPUBreitWigner3body_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
                          GDouble mass0, GDouble width0, int orbitL,
                          int daught1, int daught2 );

#endif // GPU_ACCELERATION

using std::complex;
using namespace std;

class Kinematics;

class BreitWigner3body : public UserAmplitude< BreitWigner3body >
{
  
public:
	
	BreitWigner3body() : UserAmplitude< BreitWigner3body >() {}
	BreitWigner3body( const vector< string >& args );
	
  ~BreitWigner3body(){}
  
	string name() const { return "BreitWigner3body"; }
  
  complex< GDouble > calcAmplitude( GDouble** pKin ) const;
	  
  void updatePar( const AmpParameter& par );
    
#ifdef GPU_ACCELERATION

	bool isGPUEnabled() const { return false; }

#endif // GPU_ACCELERATION
  
private:
	
  AmpParameter m_mass0;
  AmpParameter m_width0;
  
  string m_daughters;  
};

#endif
