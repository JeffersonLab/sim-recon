#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "AMPTOOLS_DATAIO/TwoPiPlotGenerator.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"

TwoPiPlotGenerator::TwoPiPlotGenerator( const FitResults& results ) :
PlotGenerator( results )
{
  // calls to bookHistogram go here
  
  // bookHistogram( k2PiMass, new Histogram1D( 200, 0., 2.0, "M2pi", "Invariant Mass of #pi^{+} #pi^{-}") );
  bookHistogram( k2PiMass, new Histogram1D( 200, 0.2, 0.8, "M2pi", "Invariant Mass of #pi^{+} #pi^{-}") );
  bookHistogram( kPiPCosTheta, new Histogram1D( 50, -1., 1., "cosTheta", "cos( #theta ) of Resonance Production") );

  bookHistogram( kPhiPiPlus,  new Histogram1D( 50, -1*PI, PI, "PhiPiPlus",  "#Phi_{#pi_{+}}" ) );
  bookHistogram( kPhiPiMinus, new Histogram1D( 50, -1*PI, PI, "PhiPiMinus", "#Phi_{#pi_{-}}" ) );
  bookHistogram( kPhi, new Histogram1D( 50, -1*PI, PI, "Phi", "#Phi" ) );
  bookHistogram( kphi, new Histogram1D( 50, -1*PI, PI, "phi", "#phi" ) );
  bookHistogram( kPsi, new Histogram1D( 50, -1*PI, PI, "psi", "#psi" ) );
  // bookHistogram( kt, new Histogram1D( 100, 0, 5, "t", "-t" ) );
  bookHistogram( kt, new Histogram1D( 100, 0, 0.05, "t", "-t" ) );
}

void
TwoPiPlotGenerator::projectEvent( Kinematics* kin ){
  
  TLorentzVector beam   = kin->particle( 0 );
  // TLorentzVector recoil = kin->particle( 1 );   // moved recoil to last particle.
  TLorentzVector p1 = kin->particle( 2 );
  TLorentzVector p2 = kin->particle( 3 );
  TLorentzVector target (0,0,0,208.);    // need to get target mass from configuration files
  TLorentzVector recoil = beam + target - p1 - p2;

  TLorentzVector resonance = p1 + p2; 
  TLorentzRotation resonanceBoost( -resonance.BoostVector() );

  TLorentzVector recoil_res = resonanceBoost * recoil;
  TLorentzVector p1_res = resonanceBoost * p1;

  // production plane is defined by the pi+ (neglect recoil)
  TVector3 zlab(0.,0.,1.0);     // z axis in lab
  TVector3 y = (p1.Vect().Cross(zlab)).Unit();    // perpendicular to decay plane. ensure that y is perpendicular to z

  double phipol = 0;     // should take this variable from the configuration file.
  TVector3 eps(cos(phipol), sin(phipol), 0.0); // beam polarization vector in lab
  TVector3 eps_perp = eps.Cross(zlab).Unit();         // perpendicular to plane defined by eps
  GDouble Phi_pip = atan2(y.Dot(eps),y.Dot(eps_perp));  // use this calculation to preserve sign of angle

  // choose helicity frame: z-axis opposite recoil target in rho rest frame. Note that for Primakoff recoil is never measured.
  y = (beam.Vect().Unit().Cross(-recoil.Vect().Unit())).Unit();   // redefine y normal to production plane
  
  // choose helicity frame: z-axis opposite recoil proton in rho rest frame
  TVector3 z = -1. * recoil_res.Vect().Unit();
  TVector3 x = y.Cross(z).Unit();
  TVector3 angles(   (p1_res.Vect()).Dot(x),
                     (p1_res.Vect()).Dot(y),
                     (p1_res.Vect()).Dot(z) );

  GDouble cosTheta = angles.CosTheta();
  
  GDouble phi = angles.Phi();
  
  GDouble Phi = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));

  // GDouble psi = phi - Phi;
  GDouble psi = Phi_pip;    // in the limit of forward scattering (Primakoff), Phi_pip is the angle between pip and the polarization
  if(psi < -1*PI) psi += 2*PI;
  if(psi > PI) psi -= 2*PI;

  // compute invariant t
  // GDouble t = - 2* recoil.M() * (recoil.E()-recoil.M());
  GDouble t = (beam - p1 - p2).M2();     // use measured particles to compute t
  
  // calls to fillHistogram go here
  cout << " mp1=" << p1.M() << " mp2=" << p2.M() << " mrecoil=" << recoil.M() << " m2pi=" << resonance.M() << endl;
  
  fillHistogram( k2PiMass, ( resonance ).M() );
  
  fillHistogram( kPiPCosTheta, cosTheta );

  fillHistogram( kPhiPiPlus,  p1.Phi() );
  fillHistogram( kPhiPiMinus, p2.Phi() );
  fillHistogram( kPhi, Phi );
  fillHistogram( kphi, phi );

  fillHistogram( kPsi, psi );
  fillHistogram( kt, -t );      // fill with -t to make positive
}
