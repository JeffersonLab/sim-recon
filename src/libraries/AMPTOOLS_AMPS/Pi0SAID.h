#if !defined(PI0SAID)
#define PI0SAID

#include "TH2.h"

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

class Pi0SAID : public UserAmplitude< Pi0SAID >
{
    
public:
	
	Pi0SAID() : UserAmplitude< Pi0SAID >() { };
	Pi0SAID( const vector< string >& args );
	
	string name() const { return "Pi0SAID"; }
    
	complex< GDouble > calcAmplitude( GDouble** pKin ) const;
	
private:
	
	void FillDataTables();
	double DSG[31][41];
	double Sigma[31][41];

	TH2F *hCosTheta_Ebeam, *hSigma_Ebeam;
	GDouble Pgamma;
};

#endif
