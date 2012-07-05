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

class polCoef : public Amplitude
{
  
public:
  
  polCoef() : Amplitude() { setDefaultStatus( true ); }
  polCoef( int polBeam, const AmpParameter& polFrac );	
  ~polCoef(){}
  
  string name() const { return "polCoef"; }
  
  complex< GDouble > calcAmplitude( GDouble** pKin ) const;
  
  //void updatePar( const AmpParameter& par );
  
  polCoef* newAmplitude( const vector< string >& args) const;
  
  polCoef* clone() const;
  
private: 
  int m_polBeam;
  AmpParameter m_polFrac;

};

#endif
