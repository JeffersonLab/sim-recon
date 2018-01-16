#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "AMPTOOLS_DATAIO/Pi0PlotGenerator.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"

/* Constructor to display FitResults */
Pi0PlotGenerator::Pi0PlotGenerator( const FitResults& results ) :
PlotGenerator( results )
{
	createHistograms();
}

/* Constructor for event generator (no FitResult) */
Pi0PlotGenerator::Pi0PlotGenerator( ) :
PlotGenerator( )
{
	createHistograms();
}

void
Pi0PlotGenerator::createHistograms( ) {
  // calls to bookHistogram go here 
  bookHistogram( kCosTheta, new Histogram1D( 50, -1., 1., "cosTheta", "cos(#theta)") );
  bookHistogram( kPhi, new Histogram1D( 50, -1*PI, PI, "Phi", "#Phi" ) );
  bookHistogram( kt, new Histogram1D( 100, 0.0, 5.0, "t", "-t" ) );
  bookHistogram( kCosTheta_phi, new Histogram2D( 180, -3.14, 3.14, 100, -1, 1, "cosTheta_phi", "cos#theta vs. #phi; #phi; cos#theta") );
  bookHistogram( kt_phi, new Histogram2D( 180, -3.14, 3.14, 100, 0.0, 5.0, "t_phi", "-t vs. #phi; #phi; -t") );
}

void
Pi0PlotGenerator::projectEvent( Kinematics* kin ){
  
  TLorentzVector beam   = kin->particle( 0 );
  TLorentzVector recoil = kin->particle( 1 );
  TLorentzVector p1 = kin->particle( 2 );

  TLorentzVector target  ( 0., 0., 0., 0.938);	
  
  TLorentzVector cm = recoil + p1;
  TLorentzRotation cmBoost( -cm.BoostVector() );
  TLorentzVector p1_cm = cmBoost * p1;
  
  GDouble t = (target - recoil).M2();
  GDouble cosTheta = p1_cm.CosTheta();
  GDouble phi = p1_cm.Phi();
  if(phi < -1*PI) phi += 2*PI;
  if(phi > PI) phi -= 2*PI;

  // calls to fillHistogram go here
  fillHistogram( kCosTheta, cosTheta );
  fillHistogram( kPhi, phi );
  fillHistogram( kt, -t );      // fill with -t to make positive
  fillHistogram( kCosTheta_phi, phi, cosTheta ); 
  fillHistogram( kt_phi, phi, -t ); 
}
