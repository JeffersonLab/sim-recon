
#include "AMPTOOLS_DATAIO/ThreePiPlotGenerator.h"
#include "IUAmpTools/Histogram.h"
#include "IUAmpTools/Kinematics.h"

ThreePiPlotGenerator::ThreePiPlotGenerator( AmpToolsInterface& ati ) :
PlotGenerator( ati )
{
  // calls to bookHistogram go here

  bookHistogram( kHist1, "Variable 1", Histogram( 100, 0, 10 ) );
  
}

void
ThreePiPlotGenerator::projectEvent( Kinematics* kin ){
  
  HepLorentzVector P1 = kin->particle( 0 );
  
  // calls to fillHistogram go here
  
  fillHistogram( kHist1, P1.e() );
}
