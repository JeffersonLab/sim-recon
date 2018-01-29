#ifndef TTCS_KINE_ifndef
#define TTCS_KINE_ifndef

#include <TVector3.h>
#include <TLorentzVector.h>

class TTCS_kine : public TObject
{
public:
  TTCS_kine(double m = 0.938, double E = 5.76);
  TTCS_kine( TLorentzVector, TLorentzVector,  TLorentzVector, 
	     double m = 0.938, double E = 5.76 ); // e-,e+,p, mp, Eg
  
  double GetPhi_cm() const;
  double GetTheta_cm() const;
  double Get_tM() const;
  double GetMinv() const;
  double GetMM2() const;
  double GetEg() const; // Energy of incpming photon
  double GetEq_prime() const; // Energy of the timelike photon
  double GetMis_mom() const;
  double GetPx_mis() const;
  double GetPy_mis() const;
  double GetPz_mis() const;
  double Get_L() const;
  double Get_L0() const;
  double GetQ2()const; // Q2 of the quasi-real photon
  void SetLemLepLp( TLorentzVector, TLorentzVector,  TLorentzVector);

  void Define_kinematic();
  //  ClassDef(TTCS_kine,1);

private:
  
  TLorentzVector Lp;    TLorentzVector Lp_cm;
  TLorentzVector Lp1;    TLorentzVector Lp1_cm;
  TLorentzVector Lem;   TLorentzVector Lem_cm;
  TLorentzVector Lep;   TLorentzVector Lep_cm;
  TLorentzVector Lg;   TLorentzVector Lg_cm;   
  TLorentzVector Lemep;   TLorentzVector Lemep_cm;
  TLorentzVector Lgemep;   TLorentzVector Lgemep_cm;
  TLorentzVector Lem_eep_cm;
  TLorentzVector Lep_eep_cm;
  TLorentzVector Lp1_eep_cm;
  TLorentzVector Lcm;
  TLorentzVector Lbeam;
  
  TLorentzVector L_mis;
  
  double Eq_prime;
  double mprot;
  double Egamma;
  double Eb;
  double Minv, tM;
  double MM2;
  double mis_mom;
  double px_mis;
  double py_mis;
  double pz_mis;
  double L, L0;
  double Q2; // Q2 of the quasi-real photon
  static constexpr double radian = 57.2957795130823229;
  static constexpr double PI = 3.14159265358979312;

  //  Lp.SetPxPyPzE(0., 0., 0., mprot);
  //  Lbeam.SetPxPyPzE(0, 0, Eb, Eb);
  
  TVector3 TV3_em, TV3_ep, TV3_p, TV3_p1;
  TVector3 TV3_emep_crs, TV3_pp1_crs; //cross products of (em X ep) and (p X p1)
  
  double phi_cm, theta_cm;

};

#endif
