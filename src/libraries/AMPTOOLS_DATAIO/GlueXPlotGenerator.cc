
#include "AMPTOOLS_DATAIO/GlueXPlotGenerator.h"
#include "IUAmpTools/Histogram.h"
#include "IUAmpTools/Kinematics.h"

GlueXPlotGenerator::GlueXPlotGenerator( AmpToolsInterface& ati ) :
PlotGenerator( ati )
{
  // calls to bookHistogram go here

  bookHistogram( kHist1, "Variable 1", Histogram( 100, 0, 10 ) );
  
}

void
GlueXPlotGenerator::projectEvent( Kinematics* kin ){
  
  HepLorentzVector P1 = kin->particle( 0 );
  
  // calls to fillHistogram go here
  
  fillHistogram( kHist1, P1.e() );
}
