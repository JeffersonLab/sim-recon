// b1piAmpCheck - a class for computing
// angles and invariant masses of b1pi events
//  by Igor Senderovich - 11/2011

#if !defined(B1PIAMPCHECK)
#define B1PIAMPCHECK

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include <utility>
#include <string>
#include <complex>
#include <vector>

#include <TLorentzRotation.h>
#include <TLorentzVector.h>
#include "CLHEP/Vector/LorentzVector.h"


using std::complex;
using namespace std;


/**
 * An object of the b1piAmpCheck class computes
 * kinamtical variables of a b1pi decay, such as
 * decay angles, and invariant masses of all resonances
 * in the decay tree
 */

class b1piAmpCheck 
{
  
public:
	
  b1piAmpCheck();
  b1piAmpCheck(Kinematics &evt);
  ~b1piAmpCheck(){}
    
  void SetEvent(Kinematics &evt);

  void ProcKin();
  
  //Getters
  double GetRhoPhi(){return m_rho_phi;}
  double GetRhoCosTheta(){return m_rho_cosTheta;}
  double GetOmegaPhi(){return m_omega_phi;}
  double GetOmegaCosTheta(){return m_omega_cosTheta;}
  double Getb1Phi(){return m_b1_phi;}
  double Getb1CosTheta(){return m_b1_cosTheta;}
  double GetXPhi(){return m_X_phi;}
  double GetXCosTheta(){return m_X_cosTheta;}  
  double GetMX(){return m_X.M();}
  double GetMb1(){return m_b1.M();}
  double GetMomega(){return m_omega.M();}
  double GetMrho(){return m_rho.M();}

  double GetAlpha(){return m_recoil.Phi();}

  const TLorentzVector& GetOmega(){return m_omega;}
  const TLorentzVector& GetXsPi(){return m_Xs_pi;}
  const TLorentzVector& Getb1sPi(){return m_b1s_pi;}

private:
  
  TLorentzVector& Hep2T(const HepLorentzVector &v1, TLorentzVector &v2);

  TLorentzVector& MoveToRF(TLorentzVector &parent,
			   TLorentzVector &daughter);
  
  TLorentzVector m_rhos_pip, m_rho,m_omega, m_b1, m_X, m_beam, m_recoil;
  TLorentzVector m_Xs_pi, m_b1s_pi, m_omegas_pi, m_rhos_pim;
  double m_rho_phi, m_rho_cosTheta;
  double m_omega_phi, m_omega_cosTheta;
  double m_b1_phi, m_b1_cosTheta;
  double m_X_phi, m_X_cosTheta;
  double m_alpha, m_X_M,m_b1_M,m_omega_M,m_rho_M;

};

#endif
