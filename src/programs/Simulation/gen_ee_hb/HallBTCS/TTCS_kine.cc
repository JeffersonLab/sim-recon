#include "TTCS_kine.hh"
//#include "/home/rafayel/WORKDIR/e1-6/pass2_v1/tcs_analyse2/include/TTCS_kine.hh"

TTCS_kine::TTCS_kine(double m, double E)
{
  // Lem.SetPxPyPzE(0., 0., -1., 1.);
  // Lep.SetPxPyPzE(0., 0., -1., 1.);
  // Lp1.SetPxPyPzE(0., 0., -1., 1.);

  mprot = m;
  Eb = E;

  Lp.SetPxPyPzE(0., 0., 0., mprot);
  Lbeam.SetPxPyPzE(0, 0, Eb, Eb);
  
  //Define_kinematic();
}

TTCS_kine::TTCS_kine( TLorentzVector em, TLorentzVector ep, TLorentzVector p1, 
		      double m, double E )
{
  Lem = em;
  Lep = ep;
  Lp1 = p1;

  mprot = m;
  Eb = E;
  
  Lp.SetPxPyPzE(0., 0., 0., mprot);
  Lbeam.SetPxPyPzE(0, 0, Eb, Eb);

  Define_kinematic();
}

void TTCS_kine::SetLemLepLp(TLorentzVector em, TLorentzVector ep, TLorentzVector p1)
{
  Lem = em;
  Lep = ep;
  Lp1 = p1;

  Define_kinematic();
}

void TTCS_kine::Define_kinematic()
{
  L_mis = Lbeam + Lp - Lem - Lep - Lp1;
  Lg = Lp1 + Lem + Lep - Lp;
  Lemep = Lem + Lep;

  Minv = Lemep.M();
  Eq_prime = Lemep.E();
  tM = (Lp - Lp1).M2();
  Egamma = Lg.E();
  MM2 = L_mis.M2();
  mis_mom = L_mis.P();
  px_mis = L_mis.Px();
  py_mis = L_mis.Py();
  pz_mis = L_mis.Pz();
  
  Q2 = 2*Eb*(mis_mom - pz_mis);

  Lcm = Lg + Lp;
	  
  Lp_cm = Lp;
  Lp_cm.Boost( -Lcm.BoostVector() );
  Lp1_cm = Lp1;
  Lp1_cm.Boost( -Lcm.BoostVector() );
  Lem_cm = Lem;
  Lem_cm.Boost( -Lcm.BoostVector() );
  Lep_cm = Lep;
  Lep_cm.Boost( -Lcm.BoostVector() );
	  
  Lem_eep_cm = Lem;
  Lem_eep_cm.Boost( -Lemep.BoostVector() );
  Lep_eep_cm = Lep;
  Lep_eep_cm.Boost( -Lemep.BoostVector() );
  Lp1_eep_cm = Lp1;
  Lp1_eep_cm.Boost( -Lemep.BoostVector() );

  TV3_em = Lem_cm.Vect();
  TV3_ep = Lep_cm.Vect();
  TV3_p = Lp_cm.Vect();
  TV3_p1 = Lp1_cm.Vect();
  TV3_emep_crs = TV3_ep.Cross(TV3_em);
  TV3_pp1_crs = TV3_p.Cross(TV3_p1);

  if( TV3_em.Dot( TV3_pp1_crs ) > 0 )
    {
      phi_cm = TV3_emep_crs.Angle(TV3_pp1_crs)*radian;
    }
  else
    {
      phi_cm = (2*PI - TV3_emep_crs.Angle(TV3_pp1_crs))*radian;
    }

  theta_cm = (PI - Lem_eep_cm.Angle(Lp1_eep_cm.Vect()))*radian;
  
  double bb = 2*(Lem - Lep).Dot(Lp - Lp1);
  L = ((Minv*Minv - tM)*(Minv*Minv - tM) - bb*bb)/4.;
  L0 = Minv*Minv*Minv*Minv*sin(theta_cm/radian)*sin(theta_cm/radian);
}

double TTCS_kine::GetPhi_cm() const
{
  return phi_cm;
}

double TTCS_kine::GetTheta_cm() const
{
  return theta_cm;
}

double TTCS_kine::Get_tM() const
{
  return tM;
}

double TTCS_kine::GetMinv() const
{
  return Minv;
}

double TTCS_kine::GetMM2() const
{
  return MM2;
}

double TTCS_kine::GetEg() const
{
  return Egamma;
}

double TTCS_kine::GetEq_prime() const
{
  return Eq_prime;
}

double TTCS_kine::GetMis_mom() const
{
  return mis_mom;
}

double TTCS_kine::GetPx_mis() const
{
  return px_mis;
}

double TTCS_kine::GetPy_mis() const
{
  return py_mis;
}

double TTCS_kine::GetPz_mis() const
{
  return pz_mis;
}

double  TTCS_kine::GetQ2()const
{
  return Q2;
}

double TTCS_kine::Get_L() const
{
  return L;
}

double TTCS_kine::Get_L0() const
{
  return L0;
}
