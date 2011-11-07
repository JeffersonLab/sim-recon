#if !defined(B1PIANGAMP)
#define B1PIANGAMP

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/AmpParameter.h"

#include <string>
#include <complex>
#include <vector>

using std::complex;
using namespace std;


class Kinematics;

class b1piAngAmp : public Amplitude
{

public:
  
  b1piAngAmp() : Amplitude() { setDefaultStatus( true ); }
  b1piAngAmp(int polBeam, float polFrac, //const AmpParameter& polFrac,
	     int J_X, int Par_X, int L_X, int I_X, int epsilon_R, 
	     int Iz_b1, int Iz_pi,
	     float u_rho_1, float u_rho_3,
	     float u_omega_1, float u_omega_3,
	     float u_b1_0, float u_b1_2, bool orthocheck);
  
  string name() const { return "b1piAngAmp"; }
  bool containsFreeParameters() const { return false; }
  
  complex< GDouble > calcAmplitude( GDouble** pKin ) const;
  
  b1piAngAmp* newAmplitude( const vector< string >& args ) const;
  b1piAngAmp* clone() const;

  float u_rho(int J_rho) const;
  float u_omega(int L_omega) const;
  float u_b1(int L_b1) const;

private:
  
  int mpolBeam;
  //AmpParameter mpolFrac;
  float mpolFrac;
  int mJ_X, mPar_X, mL_X, mI_X, mepsilon_R, mIz_b1, mIz_pi;
  
  float m_u_rho_1, m_u_rho_3;
  float m_u_omega_1, m_u_omega_3;
  float m_u_b1_0, m_u_b1_2;
  float m_v_p, m_v_m;
  bool m_ORTHOCHECK;

};

#endif
