#if !defined(PI0REGGE)
#define PI0REGGE

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include "TH1D.h"
#include <string>
#include <complex>
#include <vector>

#include "Pi0ReggeModel.h"

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
	
private:

	double polAngle;

	TH1D *totalFlux_vs_E;
	TH1D *polFlux_vs_E;
	TH1D *polFrac_vs_E;
};

#endif
