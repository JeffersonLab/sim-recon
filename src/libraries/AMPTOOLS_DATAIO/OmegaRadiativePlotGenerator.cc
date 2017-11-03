#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "AMPTOOLS_DATAIO/OmegaRadiativePlotGenerator.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"

OmegaRadiativePlotGenerator::OmegaRadiativePlotGenerator( const FitResults& results ) :
PlotGenerator( results )
{
  // calls to bookHistogram go here
  
   bookHistogram( kOmegaMass, new Histogram1D( 200, 0.2, 0.8, "MOmega", "Invariant Mass of #pi^{0} #gamma") );
   bookHistogram( kCosThetaPi0, new Histogram1D( 50, -1., 1., "cosTheta", "cos( #theta ) of #pi^{0}") );
   bookHistogram( kCosThetaGamma, new Histogram1D( 50, -1., 1., "cosTheta", "cos( #theta ) of #Gamma") );
   bookHistogram( kPhiPi0,  new Histogram1D( 50, -1*PI, PI, "PhiPiPlus",  "#Phi_{#pi_{0}}" ) );
   bookHistogram( kPhiGamma, new Histogram1D( 50, -1*PI, PI, "PhiPiMinus", "#Phi_{#gamma}" ) );
   bookHistogram( kCosTheta, new Histogram1D( 50, -1., 1., "CosTheta", "cos#theta;cos#theta" ) );
   bookHistogram( kPhi, new Histogram1D( 50, -1*PI, PI, "Phi", "#Phi; #Phi[rad.]" ) );
   bookHistogram( kphi, new Histogram1D( 50, -1*PI, PI, "phi", "#phi; #phi[rad.]" ) );
   bookHistogram( kPsi, new Histogram1D( 50, -1*PI, PI, "psi", "#psi; #psi [rad.]" ) );
   bookHistogram( kt, new Histogram1D( 100, 0, 1.0 , "t", "-t" ) );
}

void
OmegaRadiativePlotGenerator::projectEvent( Kinematics* kin ){

   TLorentzVector beam   = kin->particle( 0 );
   TLorentzVector recoil = kin->particle( 1 );
   TLorentzVector p1 = kin->particle( 2 );
   TLorentzVector p2 = kin->particle( 3 );

   TLorentzVector resonance = p1 + p2; 
   TLorentzRotation resonanceBoost( -resonance.BoostVector() );

   TLorentzVector recoil_res = resonanceBoost * recoil;
   TLorentzVector p1_res = resonanceBoost * p1;
   TLorentzVector p2_res = resonanceBoost * p2;

   // normal to the production plane
   TVector3 y = (beam.Vect().Unit().Cross(-recoil.Vect().Unit())).Unit();

   // choose helicity frame: z-axis opposite recoil proton in rho rest frame
   TVector3 z = -1. * recoil_res.Vect().Unit();
   TVector3 x = y.Cross(z).Unit();
   TVector3 angles(   (p1_res.Vect()).Dot(x),
         (p1_res.Vect()).Dot(y),
         (p1_res.Vect()).Dot(z) );

   GDouble cosTheta = angles.CosTheta();

   GDouble phi = angles.Phi();

   //double polAngle = 0.0539258; // PARA Spring 2016
   double polAngle = 1.62927; // PERP Spring 2016
   //TVector3 eps(1.0, 0.0, 0.0); // beam polarization vector
   TVector3 eps(cos(polAngle), sin(polAngle), 0.0); // beam polarization vector
   GDouble Phi = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));

   GDouble psi = phi - Phi;
   if(psi < -1*PI) psi += 2*PI;
   if(psi > PI) psi -= 2*PI;

   // compute invariant t
   GDouble t = - 2* recoil.M() * (recoil.E()-recoil.M());

   // calls to fillHistogram go here

   fillHistogram( kOmegaMass, ( resonance ).M() );
   fillHistogram( kCosThetaPi0, p1_res.CosTheta());
   fillHistogram( kCosThetaGamma, p2_res.CosTheta() );
   fillHistogram( kPhiPi0,  p1.Phi() );
   fillHistogram( kPhiGamma, p2.Phi() );
   fillHistogram( kCosTheta,   cosTheta);
   fillHistogram( kPhi, Phi );
   fillHistogram( kphi, phi );
   fillHistogram( kPsi, psi );
   fillHistogram( kt, -t );      // fill with -t to make positive
}
