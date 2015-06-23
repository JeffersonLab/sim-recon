
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "AMPTOOLS_DATAIO/ThreePiPlotGenerator.h"
#include "IUAmpTools/Histogram.h"
#include "IUAmpTools/Kinematics.h"

ThreePiPlotGenerator::ThreePiPlotGenerator( const FitResults& results ) :
PlotGenerator( results )
{
  // calls to bookHistogram go here
  
  bookHistogram( k3PiMass, new Histogram1D( 100, 0.7, 2, "3PiM", "Invariant Mass of #pi^{-} #pi^{+} #pi^{+}" ) );
  
  bookHistogram( kPiMPiP1Mass,  new Histogram1D( 100, 0.3, 1.8, "PiMPi1M", "Invariant Mass of #pi^{-} #pi^{+}_{1}" ) );
  bookHistogram( kPiMPiP2Mass,  new Histogram1D( 100, 0.3, 1.8, "PiMPi2M", "Invariant Mass of #pi^{-} #pi^{+}_{2}" ) );
  bookHistogram( kPiP1PiP2Mass, new Histogram1D( 100, 0.3, 1.8, "PiP1PiP2M", "Invariant Mass of #pi^{+}_{1} #pi^{+}_{2}" ) );
 
  bookHistogram( kAlpha, new Histogram1D( 100, -3.14, 3.14, "Alpha", "Laboratory Polar Angle of Production Plane" ) );
  bookHistogram( kCosThetaRes, new Histogram1D( 100, -1, 1, "CosThRes", "cos( #theta ) of Resonance Production" ) );
  bookHistogram( kPhiRes, new Histogram1D( 100, -3.14, 3.14, "phiRes", "#phi of Resonance Production" ) );
}

void
ThreePiPlotGenerator::projectEvent( Kinematics* kin ){
  
  TLorentzVector beam   = kin->particle( 0 );
  TLorentzVector recoil = kin->particle( 1 );
  TLorentzVector piP1 = kin->particle( 2 );
  TLorentzVector piM = kin->particle( 3 );
  TLorentzVector piP2 = kin->particle( 4 );
  
  TLorentzVector resonance = piM + piP1 + piP2;
  
  // orientation of production plane in lab
  GDouble alpha = recoil.Vect().Phi();
  
  TLorentzRotation resRestBoost( -resonance.BoostVector() );
  
  TLorentzVector beam_res   = resRestBoost * beam;
  TLorentzVector recoil_res = resRestBoost * recoil;
  TLorentzVector piP1_res     = resRestBoost * piP1;
  
  TVector3 zRes = -recoil_res.Vect().Unit();
  TVector3 yRes = beam_res.Vect().Cross(zRes).Unit();
  TVector3 xRes = yRes.Cross(zRes);
  
  TVector3 anglesRes( (piP1_res.Vect()).Dot(xRes),
                       (piP1_res.Vect()).Dot(yRes),
                       (piP1_res.Vect()).Dot(zRes) );
  
  GDouble cosThetaRes = anglesRes.CosTheta();
  GDouble phiRes = anglesRes.Phi();

  
  // calls to fillHistogram go here
  
  fillHistogram( k3PiMass, ( piM + piP1 + piP2 ).M() );
  
  fillHistogram( kPiMPiP1Mass,  ( piM + piP1 ).M() );
  fillHistogram( kPiMPiP2Mass,  ( piM + piP2 ).M() );
  fillHistogram( kPiP1PiP2Mass, ( piP1+ piP2 ).M() );

  fillHistogram( kAlpha, alpha );
  fillHistogram( kCosThetaRes, cosThetaRes );
  fillHistogram( kPhiRes, phiRes );
}
