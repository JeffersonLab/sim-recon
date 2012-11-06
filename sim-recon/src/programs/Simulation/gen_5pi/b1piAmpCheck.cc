#include <b1piAmpCheck.h>


b1piAmpCheck::b1piAmpCheck(){}

b1piAmpCheck::b1piAmpCheck(Kinematics &evt)
{
  SetEvent(evt);
}


TLorentzVector&
b1piAmpCheck::Hep2T(const HepLorentzVector &v1, TLorentzVector &v2)
{
  v2.SetXYZT(v1.px(),v1.py(),v1.pz(),v1.e());
  return v2;
}


/**
 * This function sets the 4 vectors of an event.
 * Seven vectors must be given in all:
 * beam photon, recoil proton and the 5 pions
 * starting from the bachelor and from there
 * in order down the tree. The calculation
 * are run automatically when a new event is loaded.
 */

void
b1piAmpCheck::SetEvent(Kinematics &evt) 
{
  assert(evt.particleList().size()==7);
  
  Hep2T(evt.particle(6), m_rhos_pip);
  Hep2T(evt.particle(5), m_rhos_pim);
  m_rho = m_rhos_pim + m_rhos_pip;

  Hep2T(evt.particle(4), m_omegas_pi);
  m_omega= m_rho + m_omegas_pi;

  Hep2T(evt.particle(3), m_b1s_pi);
  m_b1= m_omega + m_b1s_pi;

  Hep2T(evt.particle(2), m_Xs_pi);
  m_X= m_b1 + m_Xs_pi;

  Hep2T(evt.particle(0), m_beam);
  Hep2T(evt.particle(1), m_recoil);

  ProcKin();
}



TLorentzVector&
b1piAmpCheck::MoveToRF(TLorentzVector &parent,
		       TLorentzVector &daughter)
{
  daughter.RotateZ(-parent.Phi());
  daughter.RotateY(-parent.Theta());
  daughter.Boost(0,0,-parent.Beta());

  return daughter;
}

void
b1piAmpCheck::ProcKin()
{

  //Resonance RF, Godfried-Jackson frame
  TLorentzRotation XRFboost( -m_X.BoostVector() );

  TLorentzVector beam_XRF   = XRFboost * m_beam;
  TLorentzVector recoil_XRF = XRFboost * m_recoil;
  
  //Define coordinate system
  TVector3 zGJ = beam_XRF.Vect().Unit();
  TVector3 yGJ = zGJ.Cross(recoil_XRF.Vect()).Unit();
  TVector3 xGJ = yGJ.Cross(zGJ);
  
  TLorentzVector b1_XRF      = XRFboost * m_b1;
  TLorentzVector omega_XRF   = XRFboost * m_omega;
  TLorentzVector rho_XRF     = XRFboost * m_rho;
  TLorentzVector rhos_pip_XRF= XRFboost * m_rhos_pip;  


  TVector3 ang_b1( (b1_XRF.Vect()).Dot(xGJ),
		   (b1_XRF.Vect()).Dot(yGJ),
		   (b1_XRF.Vect()).Dot(zGJ) );
  
  m_X_phi=ang_b1.Phi();
  m_X_cosTheta=ang_b1.CosTheta();
  
  
  TLorentzVector omega_b1RF(MoveToRF(b1_XRF, omega_XRF));
  TLorentzVector rho_omegaRF(MoveToRF(omega_b1RF,
				      MoveToRF(b1_XRF, rho_XRF)));
  TLorentzVector rhos_pip_rhoRF(MoveToRF(rho_omegaRF,
					 MoveToRF(omega_b1RF,
						  MoveToRF(b1_XRF,rhos_pip_XRF))));
 
  m_b1_phi=omega_b1RF.Phi();
  m_b1_cosTheta=omega_b1RF.CosTheta();
  
  m_omega_phi=rho_omegaRF.Phi();
  m_omega_cosTheta=rho_omegaRF.CosTheta();
  
  m_rho_phi=rhos_pip_rhoRF.Phi();
  m_rho_cosTheta=rhos_pip_rhoRF.CosTheta();  
  
}
