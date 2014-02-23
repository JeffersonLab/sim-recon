
#include "AMPTOOLS_DATAIO/ThreePiPlotGenerator.h"
#include "IUAmpTools/Histogram.h"
#include "IUAmpTools/Kinematics.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/LorentzRotation.h"

ThreePiPlotGenerator::ThreePiPlotGenerator( const FitResults& results ) :
PlotGenerator( results )
{
  // calls to bookHistogram go here
  
  bookHistogram( k3PiMass, "Invariant Mass of #pi^{-} #pi^{+} #pi^{+}", Histogram( 100, 0.7, 2 ) );
  
  bookHistogram( kPiMPiP1Mass,  "Invariant Mass of #pi^{-} #pi^{+}_{1}", Histogram( 100, 0.3, 1.8 ) );
  bookHistogram( kPiMPiP2Mass,  "Invariant Mass of #pi^{-} #pi^{+}_{2}", Histogram( 100, 0.3, 1.8 ) );
  bookHistogram( kPiP1PiP2Mass, "Invariant Mass of #pi^{+}_{1} #pi^{+}_{2}", Histogram( 100, 0.3, 1.8 ) );
 
  bookHistogram( kAlpha, "Laboratory Polar Angle of Production Plane", Histogram( 100, -3.14, 3.14 ) );
  bookHistogram( kCosThetaRes, "cos( #theta ) of Resonance Production", Histogram( 100, -1, 1 ) );
  bookHistogram( kPhiRes, "#phi of Resonance Production", Histogram( 100, -3.14, 3.14 ) );
}

void
ThreePiPlotGenerator::projectEvent( Kinematics* kin ){
  
  HepLorentzVector beam   = kin->particle( 0 );
  HepLorentzVector recoil = kin->particle( 1 );
  HepLorentzVector piP1 = kin->particle( 2 );
  HepLorentzVector piM = kin->particle( 3 );
  HepLorentzVector piP2 = kin->particle( 4 );
  
  HepLorentzVector resonance = piM + piP1 + piP2;
  
  // orientation of production plane in lab
  GDouble alpha = recoil.vect().phi();
  
  HepLorentzRotation resRestBoost( -resonance.boostVector() );
  
  HepLorentzVector beam_res   = resRestBoost * beam;
  HepLorentzVector recoil_res = resRestBoost * recoil;
  HepLorentzVector piP1_res     = resRestBoost * piP1;
  
  Hep3Vector zRes = -recoil_res.vect().unit();
  Hep3Vector yRes = beam_res.vect().cross(zRes).unit();
  Hep3Vector xRes = yRes.cross(zRes);
  
  Hep3Vector anglesRes( (piP1_res.vect()).dot(xRes),
                       (piP1_res.vect()).dot(yRes),
                       (piP1_res.vect()).dot(zRes) );
  
  GDouble cosThetaRes = anglesRes.cosTheta();
  GDouble phiRes = anglesRes.phi();

  
  // calls to fillHistogram go here
  
  fillHistogram( k3PiMass, ( piM + piP1 + piP2 ).m() );
  
  fillHistogram( kPiMPiP1Mass,  ( piM + piP1 ).m() );
  fillHistogram( kPiMPiP2Mass,  ( piM + piP2 ).m() );
  fillHistogram( kPiP1PiP2Mass, ( piP1+ piP2 ).m() );

  fillHistogram( kAlpha, alpha );
  fillHistogram( kCosThetaRes, cosThetaRes );
  fillHistogram( kPhiRes, phiRes );
}
