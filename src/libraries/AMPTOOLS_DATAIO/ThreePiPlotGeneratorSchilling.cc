#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "AMPTOOLS_DATAIO/ThreePiPlotGeneratorSchilling.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"

ThreePiPlotGeneratorSchilling::ThreePiPlotGeneratorSchilling( const FitResults& results ) :
PlotGenerator( results )
{
  // calls to bookHistogram go here
  bookHistogram( k3PiMass, new Histogram1D( 75, 0.6, 0.9, "M3pi", "Invariant Mass of #pi^{+} #pi^{-} #pi^{0}") );
  bookHistogram( kCosThetaPiPlus, new Histogram1D( 50, -1., 1., "cosTheta", "cos( #theta ) of #pi^{+}") );
  bookHistogram( kCosThetaPiMinus, new Histogram1D( 50, -1., 1., "cosTheta", "cos( #theta ) of #pi^{-}") );
  bookHistogram( kCosThetaPi0, new Histogram1D( 50, -1., 1., "cosTheta", "cos( #theta ) of #pi^{0}") );
  bookHistogram( kPhiPiPlus,  new Histogram1D( 50, -1*PI, PI, "PhiPiPlus",  "#Phi_{#pi_{+}}" ) );
  bookHistogram( kPhiPiMinus, new Histogram1D( 50, -1*PI, PI, "PhiPiMinus", "#Phi_{#pi_{-}}" ) );
  bookHistogram( kPhiPi0, new Histogram1D( 50, -1*PI, PI, "PhiPi0", "#Phi_{#pi_{0}}" ) );
  bookHistogram( kCosTheta, new Histogram1D( 50, -1., 1., "CosTheta", "cos#theta" ) );
  bookHistogram( kPhi, new Histogram1D( 50, -1*PI, PI, "Phi", "#Phi" ) );
  bookHistogram( kphi, new Histogram1D( 50, -1*PI, PI, "phi", "#phi" ) );
  bookHistogram( kPsi, new Histogram1D( 50, -1*PI, PI, "psi", "#psi" ) );
  bookHistogram( kt, new Histogram1D( 100, 0, 1.0 , "t", "-t" ) );
}

void
ThreePiPlotGeneratorSchilling::projectEvent( Kinematics* kin ){
  
  TLorentzVector beam   = kin->particle( 0 );
  TLorentzVector recoil = kin->particle( 1 );
  TLorentzVector p1 = kin->particle( 2 );
  TLorentzVector p2 = kin->particle( 3 );
  TLorentzVector p3 = kin->particle( 4 );

  TLorentzVector resonance = p1 + p2 + p3; 
  TLorentzRotation resonanceBoost( -resonance.BoostVector() );

  TLorentzVector recoil_res = resonanceBoost * recoil;
  TLorentzVector p1_res = resonanceBoost * p1;
  TLorentzVector p2_res = resonanceBoost * p2;
  TLorentzVector p3_res = resonanceBoost * p3;

  // Three pi decay, use normal to decay plane 
  TVector3 norm = (p1_res.Vect().Cross(p2_res.Vect())).Unit();

  // normal to the production plane
  TVector3 y = (beam.Vect().Unit().Cross(-recoil.Vect().Unit())).Unit();
  
  // choose helicity frame: z-axis opposite recoil proton in rho rest frame
  TVector3 z = -1. * recoil_res.Vect().Unit();
  TVector3 x = y.Cross(z).Unit();
  TVector3 angles(   norm.Dot(x),
                     norm.Dot(y),
                     norm.Dot(z) );

  GDouble cosTheta = angles.CosTheta();
  
  GDouble phi = angles.Phi();
  
  TVector3 eps(1.0, 0.0, 0.0); // beam polarization vector
  GDouble Phi = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));

  GDouble psi = phi - Phi;
  if(psi < -1*PI) psi += 2*PI;
  if(psi > PI) psi -= 2*PI;

  // compute invariant t
  GDouble t = - 2* recoil.M() * (recoil.E()-recoil.M());

  // calls to fillHistogram go here
  
  fillHistogram( k3PiMass, ( resonance ).M() );
  fillHistogram( kCosThetaPiPlus, p1_res.CosTheta());
  fillHistogram( kCosThetaPiMinus, p2_res.CosTheta() );
  fillHistogram( kCosThetaPi0, p3_res.CosTheta() );
  fillHistogram( kPhiPiPlus,  p1.Phi() );
  fillHistogram( kPhiPiMinus, p2.Phi() );
  fillHistogram( kPhiPi0,     p3.Phi() );
  fillHistogram( kCosTheta,   cosTheta);
  fillHistogram( kPhi, Phi );
  fillHistogram( kphi, phi );
  fillHistogram( kPsi, psi );
  fillHistogram( kt, -t );      // fill with -t to make positive
}
