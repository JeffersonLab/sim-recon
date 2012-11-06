#if !defined(TWOPSANGLES)
#define TWOPSANGLES

#include "IUAmpTools/Amplitude.h"
#include "GPUManager/GPUCustomTypes.h"

#include <string>
#include <complex>
#include <vector>

#ifdef GPU_ACCELERATION
void
GPUTwoPSAngles_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
                     int j, int m, GDouble bigTheta, GDouble refFact );

#endif // GPU_ACCELERATION

using std::complex;
using namespace std;

// A class for describing the angular portion of the decay in the
// reflectivity basis of R -> 1 2
// particles 1 and 2 are pseudoscalars
//
// j,m are the total and z projection of the spin of R
// e is reflectivity
// p is parity of R

class Kinematics;

class TwoPSAngles : public Amplitude
{
    
public:
	
	TwoPSAngles() : Amplitude() { setDefaultStatus( true ); }
	TwoPSAngles( int j, int m, int e ); 
	
	string name() const { return "TwoPSAngles"; }
    
	complex< GDouble > calcAmplitude( GDouble** pKin ) const;
	
	TwoPSAngles* newAmplitude( const vector< string >& args ) const;
	TwoPSAngles* clone() const;

#ifdef GPU_ACCELERATION
  
  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;
  
	bool isGPUEnabled() const { return true; }
  
#endif // GPU_ACCELERATION
  
private:
        
  int m_j;
	int m_m;
  int m_e;
	
	GDouble m_bigTheta;
	int m_reflectivityFactor;
};

#endif
