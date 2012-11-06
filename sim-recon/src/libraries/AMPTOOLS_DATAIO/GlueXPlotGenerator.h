#if !(defined GLUEXPLOTGENERATOR)
#define GLUEXPLOTGENERATOR

#include <vector>
#include <string>

#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/PlotGenerator.h"
#include "IUAmpTools/Histogram.h"
#include "IUAmpTools/AmplitudeManager.h"

using namespace std;

class GlueXPlotGenerator : public PlotGenerator
{
    
public:
    
    GlueXPlotGenerator( ConfigurationInfo* cfgInfo,
                       const string& parFile );
    
    virtual ~GlueXPlotGenerator(){}
    
    const vector< string >& availablePlots() const;
    
private:
        
    vector< Histogram > fillProjections( const string& fsName,
                                         PlotType type );
    void registerPhysics( AmplitudeManager* ampManager );
    
    vector< string > m_histTitles;

};

#endif
