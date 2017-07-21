#if !defined(TWOPIANGLES_PRIMAKOFF)
#define TWOPIANGLES_PRIMAKOFF

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include <string>
#include <complex>
#include <vector>

using std::complex;
using namespace std;

class Kinematics;

class TwoPiAngles_primakoff : public UserAmplitude< TwoPiAngles_primakoff >
{
    
public:
	
	TwoPiAngles_primakoff() : UserAmplitude< TwoPiAngles_primakoff >() { };
	TwoPiAngles_primakoff( const vector< string >& args );
	
	string name() const { return "TwoPiAngles_primakoff"; }
    
	complex< GDouble > calcAmplitude( GDouble** pKin ) const;
	
#ifdef GPU_ACCELERATION
	bool isGPUEnabled() const { return false; }
#endif // GPU_ACCELERATION
  
private:

  Double_t phipol;
  Int_t m_rho;
  Int_t PhaseFactor;
  AmpParameter polFrac;
  Int_t flat;

};

#endif
