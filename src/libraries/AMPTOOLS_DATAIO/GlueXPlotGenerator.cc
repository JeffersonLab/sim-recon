
#include "AMPTOOLS_DATAIO/GlueXPlotGenerator.h"
#include "AMPTOOLS_AMPS/TwoPSAngles.h"
#include "AMPTOOLS_AMPS/ThreePiAngles.h"
#include "AMPTOOLS_AMPS/BreitWigner.h"

GlueXPlotGenerator::GlueXPlotGenerator( ConfigurationInfo* cfgInfo,
                                       const string& parFile ) :
PlotGenerator( cfgInfo, parFile ){
    
    initialize();
}


const vector< string >& 
GlueXPlotGenerator::availablePlots() const { 
    
    return m_histTitles; 
}

vector< Histogram > 
GlueXPlotGenerator::fillProjections( const string& fsName,
                                     PlotType type ){
    
    return vector< Histogram >();
}

void 
GlueXPlotGenerator::registerPhysics( AmplitudeManager* ampManager ){
    
  ampManager->registerAmplitudeFactor( BreitWigner() );
  ampManager->registerAmplitudeFactor( TwoPSAngles() );
  ampManager->registerAmplitudeFactor( ThreePiAngles() );
}

