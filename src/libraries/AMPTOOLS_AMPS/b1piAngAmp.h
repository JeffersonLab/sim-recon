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
	     float u_b1_0, float u_b1_2, float G0_omega,float G0_b1,
	     bool orthocheck, bool fastCalc);
  
  string name() const { return "b1piAngAmp"; }
  bool containsFreeParameters() const { return false; }
  
  complex< GDouble > calcAmplitude( GDouble** pKin ) const;
  
  b1piAngAmp* newAmplitude( const vector< string >& args ) const;
  b1piAngAmp* clone() const;

  float u_rho(int J_rho) const;
  float u_omega(int L_omega) const;
  float u_b1(int L_b1) const;

  inline complex<GDouble> BreitWigner(GDouble m0, GDouble Gamma0, int L,
			       HepLorentzVector &P1,HepLorentzVector &P2) const;
  inline float CB(int j1, int j2, int m1, int m2, int J, int M) const;

  inline GDouble N(int J) const;

private:
  
  int mpolBeam;
  float mpolFrac;
  //AmpParameter mpolFrac;
  int mJ_X, mPar_X, mL_X, mI_X, mepsilon_R, mIz_b1, mIz_pi;
  
  float m_u_rho_1, m_u_rho_3;
  float m_u_omega_1, m_u_omega_3;
  float m_u_b1_0, m_u_b1_2, mG0_omega, mG0_b1;
  float m_v_p, m_v_m;
  bool m_ORTHOCHECK, m_fastCalc;
  bool m_disableBW_omega, m_disableBW_b1;

};

#endif
