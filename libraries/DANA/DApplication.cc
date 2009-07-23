// $Id$
//
//    File: DApplication.cc
// Created: Mon Jul  3 21:46:01 EDT 2006
// Creator: davidl (on Darwin Harriet.local 8.6.0 powerpc)
//

#include <string>
using std::string;
#include <JANA/JVersion.h>

#include "DApplication.h"
#include "HDDM/DEventSourceHDDMGenerator.h"
#include "HDGEOMETRY/DMagneticFieldMapCalibDB.h"
#include "HDGEOMETRY/DMagneticFieldMapConst.h"
#include "HDGEOMETRY/DMagneticFieldMapSpoiled.h"
#include "HDGEOMETRY/DLorentzMapCalibDB.h"
//#include "HDGEOMETRY/DMaterialMapCalibDB.h"
#include "HDGEOMETRY/DRootGeom.h"
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
	
	// Add plugin paths to Hall-D specific binary directories
	const char *bms = getenv("BMS_OSNAME");
	string sbms(bms==NULL ? "":bms);
	
	if(const char *ptr = getenv("DANA_PLUGIN_PATH")){
		AddPluginPath(string(ptr));
	}
	if(const char *ptr = getenv("HALLD_MY")){
		AddPluginPath(string(ptr) + "/lib/" + sbms);
	}
	if(const char *ptr = getenv("HALLD_HOME")){
		AddPluginPath(string(ptr) + "/lib/" + sbms);
	}
	
	if(JVersion::minor<5)Init();
}

//---------------------------------
// Init
//---------------------------------
jerror_t DApplication::Init(void)
{
	this->JApplication::Init();

	// Create magnetic field object for use by everyone
	// Allow a trivial homogeneous map to be used if 
	// specified on the command line
	string bfield_type = "CalibDB";
	GetJParameterManager()->SetDefaultParameter("BFIELD_TYPE", bfield_type);
	if(bfield_type=="CalibDB"){
		bfield = new DMagneticFieldMapCalibDB(this);
		cout<<"Created Magnetic field map of type DMagneticFieldMapCalibDB."<<endl;
	}else if(bfield_type=="Const"){
		bfield = new DMagneticFieldMapConst(this);
		cout<<"Created Magnetic field map of type DMagneticFieldMapConst."<<endl;
	}else if(bfield_type=="Spoiled"){
		bfield = new DMagneticFieldMapSpoiled(this);
		cout<<"Created Magnetic field map of type DMagneticFieldMapSpoiled."<<endl;
	}else{
		_DBG_<<" Unknown DMagneticFieldMap subclass \"DMagneticFieldMap"<<bfield_type<<"\" !!"<<endl;
		exit(-1);
	}

	// Create Lorentz deflection object
	lorentz_def= new DLorentzMapCalibDB(this);
	
	// Create map of material properties
	//material = new DMaterialMapCalibDB(this);
	RootGeom= new DRootGeom(this);
	
	return NOERROR;
}

//---------------------------------
// ~DApplication    (Destructor)
//---------------------------------
DApplication::~DApplication()
{

}

//---------------------------------
// GetDGeometry
//---------------------------------
DGeometry* DApplication::GetDGeometry(unsigned int run_number)
{
	/// Get the DGeometry object for the specified run number.
	/// The DGeometry class is Hall-D specific. It uses the
	/// JGeometry class from JANA to access values in the HDDS
	/// XML files. However, it supplies some useful and more
	/// user friendly methods for getting at some of the values.
	///
	/// This will first look for the DGeometry object in a list
	/// kept internal to DApplication and return a pointer to the
	/// object if found there. If it is not found there, then
	/// a new DGeometry object will be created and added to the
	/// internal list before returning a pointer to it.
	///
	/// Note that since this method can change internal data
	/// members, a mutex is locked to ensure integrity. This
	/// means that it is <b>NOT</b> efficient to call this
	/// method for every event. The pointer should be obtained
	/// in a brun() method and kept in a local variable if
	/// needed outside of brun().
	
	// First, get the JGeometry object using our JApplication
	// base class. Then, use that to find the correct DGeometry
	// object if it exists.
	JGeometry *jgeom = GetJGeometry(run_number);
	if(!jgeom){
		_DBG_<<"ERROR: Unable get geometry for run "<<run_number<<"!"<<endl;
		_DBG_<<"Make sure you JANA_GEOMETRY_URL environment variable is set."<<endl;
		_DBG_<<"It should be set to something like:"<<endl;
		_DBG_<<endl;
		_DBG_<<"    xmlfile://${HALLD_HOME}/src/programs/Simulation/hdds/main_HDDS.xml"<<endl;
		_DBG_<<endl;
		_DBG_<<"Exiting now."<<endl;
		Quit();
		exit(-1);
		return NULL;
	}

	Lock();
	
	for(unsigned int i=0; i<geometries.size(); i++){
		if(geometries[i]->GetJGeometry() == jgeom){
			DGeometry *dgeom = geometries[i];
			Unlock();
			return dgeom;
		}
	}
	
	// Couldn't find a DGeometry object that uses this JGeometry object.
	// Create one and add it to the list.
	DGeometry *dgeom = new DGeometry(jgeom, this, run_number);
	geometries.push_back(dgeom);
	
	cout<<"Geometry created:"<<endl;
	cout<<"  Run requested:"<<jgeom->GetRunRequested()<<"  found:"<<jgeom->GetRunFound()<<endl;
	cout<<"  Run validity range: "<<jgeom->GetRunMin()<<"-"<<jgeom->GetRunMax()<<endl;
	cout<<"  URL=\""<<jgeom->GetURL()<<"\""<<"  context=\""<<jgeom->GetContext()<<"\""<<endl;
	
	Unlock();
	
	return dgeom;
}
