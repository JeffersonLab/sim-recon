#if !defined(THREEPIANGLES)
#define THREEPIANGLES

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/AmpParameter.h"

#include <string>
#include <complex>
#include <vector>

using std::complex;
using namespace std;

#ifdef GPU_ACCELERATION
void
GPUThreePiAngles_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
                       int polX, GDouble polFrac, int jX, int parX, int iX, 
                       int lX, int jI, int iI, int iZ0, int iZ1, int iZ2 );
#endif

class Kinematics;

class ThreePiAngles : public Amplitude
{

public:
	
	ThreePiAngles() : Amplitude() { setDefaultStatus( true ); }
  ThreePiAngles( int polX, const AmpParameter& polFrac, int jX, int parX, 
                 int iX, int lX, int jI, int iI, int iZ0, int iZ1, int iZ2 );
  
	string name() const { return "ThreePiAngles"; }
	bool containsFreeParameters() const { return false; }

	complex< GDouble > calcAmplitude( GDouble** pKin ) const;
	
	ThreePiAngles* newAmplitude( const vector< string >& args ) const;
	ThreePiAngles* clone() const;
	
#ifdef GPU_ACCELERATION
  
  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;
  
	bool isGPUEnabled() const { return true; }
  
#endif // GPU_ACCELERATION
  
private:
	
	int m_polBeam;
  AmpParameter m_polFrac;
	int m_jX;
  int m_parX;
	int m_iX;
	int m_lX;
	int m_jI;
	int m_iI;
  int m_iZ0;
  int m_iZ1;
  int m_iZ2;

  vector< int > m_iZ;
    
};

#endif
