// $Id$
//
//    File: DApplication.cc
// Created: Mon Jul  3 21:46:01 EDT 2006
// Creator: davidl (on Darwin Harriet.local 8.6.0 powerpc)
//

#include "DApplication.h"
#include "HDDM/DEventSourceHDDMGenerator.h"
#include "TRACKING/DMagneticFieldMapGlueX.h"
#include "DFactoryGenerator.h"

//---------------------------------
// DApplication    (Constructor)
//---------------------------------
DApplication::DApplication(int narg, char* argv[]):JApplication(narg, argv)
{
	/// Add DEventSourceHDDMGenerator and
	/// DFactoryGenerator, which adds the default
	/// list of Hall-D factories
	AddEventSourceGenerator(new DEventSourceHDDMGenerator());
	AddFactoryGenerator(new DFactoryGenerator());
	
	if(const char *ptr = getenv("DANA_PLUGIN_PATH")){
		AddPluginPath(string(ptr));
	}
	if(const char *ptr = getenv("HALLD_MY")){
		AddPluginPath(string(ptr) + "/src/programs/Analysis/plugins");
	}
	if(const char *ptr = getenv("HALLD_HOME")){
		AddPluginPath(string(ptr) + "/src/programs/Analysis/plugins");
	}
	
	// Create magnetic field object for use by everyone
	bfield = new DMagneticFieldMapGlueX();
}

//---------------------------------
// ~DApplication    (Destructor)
//---------------------------------
DApplication::~DApplication()
{

}
