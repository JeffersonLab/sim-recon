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
  
  bookHistogram( k2PiMass, new Histogram1D( 200, 0.0, 2.0, "M2pi", "Invariant Mass of #pi^{+} #pi^{-}") );
  bookHistogram( kPiPCosTheta, new Histogram1D( 100, -1, 1, "cosTheta", "cos( #theta ) of Resonance Production") );

  bookHistogram( kPhi, new Histogram1D( 100, -1*PI, PI, "Phi", "#Phi" ) );
  bookHistogram( kphi, new Histogram1D( 100, -1*PI, PI, "phi","#phi" ) );

  bookHistogram( kPsi1, new Histogram1D( 100, -1*PI, PI, "psi1", "#psi 1") );
  bookHistogram( kPsi2, new Histogram1D( 100, -1*PI, PI, "psi2", "#psi 2") );
}

void
TwoPiPlotGenerator::projectEvent( Kinematics* kin ){
  
  TLorentzVector beam   = kin->particle( 0 );
  TLorentzVector recoil = kin->particle( 1 );
  TLorentzVector p1 = kin->particle( 2 );
  TLorentzVector p2 = kin->particle( 3 );

  TLorentzVector resonance = p1 + p2; 
  TLorentzRotation resonanceBoost( -resonance.BoostVector() );

  TLorentzVector beam_res = resonanceBoost * beam;
  TLorentzVector recoil_res = resonanceBoost * recoil;
  TLorentzVector p1_res = resonanceBoost * p1;

  TVector3 z = -recoil_res.Vect().Unit();
  TVector3 y = beam_res.Vect().Cross(z).Unit();
  TVector3 x = y.Cross(z).Unit();

  TVector3 angles(   (p1_res.Vect()).Dot(x),
                       (p1_res.Vect()).Dot(y),
                       (p1_res.Vect()).Dot(z) );

  GDouble cosTheta = angles.CosTheta();
  GDouble phi = angles.Phi();
  GDouble Phi = -1. * resonance.Vect().Phi();

  GDouble psi = p1_res.Phi();
  GDouble psi2 = phi - Phi;
  if(psi2 < -1*PI) psi2 += 2*PI;
  if(psi2 > PI) psi2 -= 2*PI;

  // calls to fillHistogram go here
  
  fillHistogram( k2PiMass, ( resonance ).M() );
  
  fillHistogram( kPiPCosTheta, cosTheta );
 
  fillHistogram( kphi, phi );
  fillHistogram( kPhi, Phi );

  fillHistogram( kPsi1, psi );
  fillHistogram( kPsi2, psi2 );
}
