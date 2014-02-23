
#include <stdlib.h>
#include <dlfcn.h>
#include <unistd.h>

#include <iostream>
#include <vector>
#include <string>
#include <map>
using namespace std;

#include <JANA/JCalibrationFile.h>
#include <HDGEOMETRY/DMagneticFieldMapCalibDB.h>
#include <HDGEOMETRY/DMagneticFieldMapConst.h>
#include "HDGEOMETRY/DMagneticFieldMapSpoiled.h"
#include "HDGEOMETRY/DMagneticFieldMapParameterized.h"


//----------------------------------------------------------
// This file contains routines used to implement dynamic
// geometry linking in hdgeant. The routines generated 
// by the hdds-geant utility (part of the HDDS) package
// are what are linked here. Here is how it works:
// 
// The routines that were called from the hdgeant core
// have been replaced with a corresponding wrapper routine
// that is defined here. The wrapper simply calls the
// actual function using a function pointer. For example:
// 
// void md5geom_wrapper_(char *md5){
// 	
// 	(*md5geom_ptr)(md5);
// }
// 
// where the global variable md5geom_ptr is of type
// pointer to a function with format:
//     "void md5geom_(char *md5)"
// 
// Whereas previously one might simply call:
// 
//     md5geom_(md5);
// 
// now one calls:
// 
//     md5geom_wrapper(md5);
// 
// 
// The function pointers are all initialized to point to
// the statically linked routines so that they are what
// is used if the user does not specify the XML geometry
// be used. The pointers are initialized to point to the
// static routines with statements like this:
//
//     typeof(md5geom_) *md5geom_ptr = md5geom_;
//
// where the use of typeof() makes it so we don't have
// to enter the signature of the function twice.
// 
// If the user does specify that the XML source be used
// to dynamically regenerate the geometry, then the
// pointers are updated to point to routines found in
// the shared object. The shared object is generated at
// run time using the hdds-geant utility to generate a
// FORTRAN source file which is then compiled into a
// shared object and opened using the dl library.
// 
// This should work for most all modifications to the
// geometry. The one problem would be if a major change
// was made that required another routine be generated
// by hdds-geant. For example, a getsection_() routine
// is added. In that case, one would need to add the 
// declaration of a wrapper routine below and a function
// pointer for it. One would also need to add a call
// to GetRoutine() and then make sure all places in the
// code using the routine called the wrapper function.
//----------------------------------------------------------

extern "C" {
	
	void init_runtime_xml_(void);


	// Declare statically linked routines
	void hddsgeant3_(void);
	void md5geom_(char *md5);
	float guplsh_(int *medi0, int *medi1);
	void gufld_(float *r, float *B);
	void getoptical_(int *imat, float *E, float *refl, float *abs1, float *rind, float *plsh, float *eff);

	// Initialize routine pointers to use statically linked routines
	typeof(hddsgeant3_) *hddsgeant3_ptr = hddsgeant3_;
	typeof(md5geom_)    *md5geom_ptr    = md5geom_;
	typeof(guplsh_)     *guplsh_ptr     = guplsh_;
	typeof(gufld_)      *gufld_ptr      = gufld_;
	typeof(getoptical_) *getoptical_ptr = getoptical_;

	// Trivial wrapper routines use pointer to dispatch call
	void hddsgeant3_wrapper_(void){ (*hddsgeant3_ptr)(); }
	void md5geom_wrapper_(char *md5){ (*md5geom_ptr)(md5); }
	float guplsh_wrapper_(int *medi0, int *medi1){ return (*guplsh_ptr)(medi0, medi1); }
	void gufld_wrapper_(float *r, float *B){ (*gufld_ptr)(r, B); }
	void getoptical_wrapper_(int *imat, float *E, float *refl, float *abs1, float *rind, float *plsh, float *eff){ (*getoptical_ptr)(imat, E, refl, abs1, rind, plsh, eff);}

	// Below are several more routines which need to be implemented
	// using the same three steps as above. However, all of them
	// have the same format of returning an int and taking no
	// arguments. We compact things a bit using the following
	// macro which would expand to something like this:
	//
	//   int getcolumn_(void);
	//   int (*getcolumn_ptr)(void) = getcolumn_;
	//   int getcolumn_wrapper_(void) { return ( *getcolumn_ptr)(); }
	//
#define MakeDispatcherINT(N) \
	int N(void); \
	int (* N ## ptr)(void) = N; \
	int N ## wrapper_(void) { return ( *N ## ptr)(); }

	MakeDispatcherINT(getcolumn_);
	MakeDispatcherINT(getlayer_);
	MakeDispatcherINT(getmap_);
	MakeDispatcherINT(getmodule_);
	MakeDispatcherINT(getpackage_);
	MakeDispatcherINT(getplane_);
	MakeDispatcherINT(getring_);
	MakeDispatcherINT(getrow_);
	MakeDispatcherINT(getsector_);
}


void GetRoutine(void **ptr, const char *rname);

static void *dlgeom_handle=NULL;
string HDDS_XML = "$HDDS_HOME/main_HDDS.xml";


//------------------
// init_runtime_xml_
//------------------
void init_runtime_xml_(void)
{
	cout<<endl;
	cout<<"=============================================================="<<endl;
	cout<<"Enabling dynamic geometry rendering"<<endl;
	cout<<"- - - - - - - - - - - - - - - - - - -"<<endl;
	cout<<endl;
	cout<<"Please make sure the following environment variables are set:" <<endl;
	cout<<"   HDDS_HOME  "<< endl;
	cout<<"   BMS_OSNAME "<< endl;
	cout<<"   CERN       "<< endl;
	cout<<"   CERN_LEVEL "<< endl;

	// Generate FORTRAN code from XML
	cout<<endl;
	cout << "Generating FORTRAN from XML source ...." << endl;
	string cmd = "$HDDS_HOME/bin/$BMS_OSNAME/hdds-geant " + HDDS_XML + " > tmp.F";
	cout << cmd << endl;
	system(cmd.c_str());
	
	// Compile FORTRAN into shared object
	cout<<endl;
	cout << "Compiling FORTRAN into shared object ..." << endl;
	cmd = "gfortran -shared -fPIC -o tmp.so -I$CERN/$CERN_LEVEL/include tmp.F";
	cout << cmd << endl;
	system(cmd.c_str());

	// Attach shared object
	cout<<endl;
	cout << "Attaching shared object ..." << endl;
	string fname = "./tmp.so";
	void *handle = dlopen(fname.c_str(), RTLD_NOW | RTLD_GLOBAL);
	if(!handle){
		cerr<<"Unable to open \""<<fname<<"\"!"<<endl;
		cerr<<dlerror()<<endl;
		exit(-1);
	}
	
	// Copy handle to global variable
	dlgeom_handle = handle;
	
	// Get routines
	cout << "Linking routines ..." << endl;
	GetRoutine((void**)&hddsgeant3_ptr, "hddsgeant3_");
	GetRoutine((void**)&md5geom_ptr, "md5geom_");
	GetRoutine((void**)&guplsh_ptr, "guplsh_");
	GetRoutine((void**)&gufld_ptr, "gufld_");
	GetRoutine((void**)&getoptical_ptr, "getoptical_");

	GetRoutine((void**)&getcolumn_ptr, "getcolumn_");
	GetRoutine((void**)&getlayer_ptr, "getlayer_");
	GetRoutine((void**)&getmap_ptr, "getmap_");
	GetRoutine((void**)&getmodule_ptr, "getmodule_");
	GetRoutine((void**)&getpackage_ptr, "getpackage_");
	GetRoutine((void**)&getplane_ptr, "getplane_");
	GetRoutine((void**)&getring_ptr, "getring_");
	GetRoutine((void**)&getrow_ptr, "getrow_");
	GetRoutine((void**)&getsector_ptr, "getsector_");

	// Clean up (this won't delete the tmp.so file right away since
	// we still have it open).
	unlink("./tmp.F");
	unlink("./tmp.so");

	// do not close shared object since it may contain needed routines
	//dlclose(dlgeom_handle);

	cout<<endl;
	cout << "Geometry loaded successfully" << endl;
	cout<<"=============================================================="<<endl;
}

//------------------
// GetRoutine
//------------------
void GetRoutine(void **ptr, const char *rname)
{
	// This will extract the MD5 checksum of the 
	// geometry from the dynamically linked shared
	// object. It is called from the GetMD5Geom()
	// routine below. Use that routine to get the
	// checksum, not this one.

	// Create, compile and link shared object if needed.
//	if(!dlgeom_handle) init_runtime_xml();

	// Grab md5geom routine from shared object
	void *my_ptr = dlsym(dlgeom_handle, rname);
	char *err = dlerror();
	if(err != NULL){
		cout << "  ... " << rname << " not found. static version will be used" << endl;
	}else{
		*ptr = my_ptr;
		cout << "  ... linked " << rname << endl;
	}
}





