
#include <stdlib.h>
#include <dlfcn.h>
#include <unistd.h>

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoPcon.h>
#include <TGeoPgon.h>
#include <TGeoMatrix.h>

#include <DANA/DApplication.h>

//#include <hddsroot.h>

extern TGeoManager * hddsroot();
extern const char* md5geom(void);

void Usage(void);
void ParseCommandLineArguments(int narg, char *argv[]);

bool USE_XML = false;
bool SHOW_ANCESTORY=true;
bool SHOW_CCDB = true;
double X_LAB = 20.0;
double Y_LAB = 20.0;
double Z_LAB = 650.0;
int32_t RUN_NUMBER = 30000;

static void *dlgeom_handle=NULL;
string HDDS_XML = "$HDDS_HOME/main_HDDS.xml";

void init_runtime_xml(void);
void MakeSharedObjectFromXML(void);
TGeoManager* hddsroot_runtime(void);
const char* GetMD5Geom(void);
const char* md5geom_runtime(void);
const char* md5geom_xml(void);

//------------------
// main
//------------------
int main(int narg, char *argv[])
{
	
	ParseCommandLineArguments(narg, argv);
	
	double pos[3];
	pos[0] = X_LAB;
	pos[1] = Y_LAB;
	pos[2] = Z_LAB;

	if(!gGeoManager){
		new TGeoManager();
		cout<<"Created TGeoManager :"<<gGeoManager<<endl;
	}

	TGeoManager *DRGeom = NULL;
	if(USE_XML){
		DRGeom = hddsroot_runtime();
	}else{
		DRGeom = hddsroot();
	}
	if(!DRGeom){
		cerr<<"Can't get TGeoManager object!"<<endl;
		return -2;
	}

	TGeoNode *cnode = DRGeom->FindNode(X_LAB, Y_LAB, Z_LAB);
	if(!cnode){
		cerr<<"Can't get node object from TGeoManager!"<<endl;
		return -2;
	}

	TGeoVolume *vol = cnode->GetVolume();
	if(!vol){
		cerr<<"Can't get volume object from TGeoNode!"<<endl;
		return -2;
	}

  TGeoMaterial *mat = vol->GetMedium()->GetMaterial();
	if(!mat){
		cerr<<"Can't get material object from TGeoVolume!"<<endl;
		return -3;
	}
	
	DGeometry *geom = NULL;
	if(SHOW_CCDB){
		DApplication *dapp = new DApplication(narg, argv);
		geom = dapp->GetDGeometry(RUN_NUMBER);
	}
	
	cout<<endl;
	cout<<"Ignore above ROOT errors they are \"normal\" ;)" << endl;
	cout<<endl;

	double density=mat->GetDensity();
	double RadLen=mat->GetRadLen();
	double A=mat->GetA();
	double Z=mat->GetZ();
	
	cout<<endl;
	cout<<"   Location: (X, Y, Z) = ("<<X_LAB<<", "<<Y_LAB<<", "<<Z_LAB<<")"<<endl;
	cout<<"==============================================="<<endl;
	cout<<"     Volume: "<<vol->GetName()<<endl;
	cout<<"   material: "<<mat->GetName()<<endl;
	cout<<"    density: "<<density<<" g/cm^3"<<endl;
	cout<<"rad. length: "<<RadLen<<" cm"<<endl;
	cout<<"          A: "<<A<<endl;
	cout<<"          Z: "<<Z<<endl;
	
	if(SHOW_ANCESTORY){
	
		gGeoManager->GetCurrentNavigator()->SetCurrentPoint(pos);
		cout<<"  ancestory: ";
		
		for(int i=0; i<1000; i++){
			TGeoNode *node = gGeoManager->GetCurrentNavigator()->GetMother(i);
			if(!node)break;
			if(i>0) cout << " -> ";
			cout << node->GetVolume()->GetName();
		}
		cout<<endl;
	}
	
	// Optionally get info from material maps in CCDB and print
	if(SHOW_CCDB){
		DVector3 pos(X_LAB, Y_LAB, Z_LAB);
		const DMaterialMap *map = geom->FindDMaterialMap(pos);
		cout << endl;
		if(!map){
			cerr << "Unable to find material map in CCDB corresponding to this point" << endl;
		}else{
			const DMaterialMap::MaterialNode *node = map->FindNode(pos);
			
			// Some maps have numbers very close to, but not quite
			// zero fro rmin. For purposes of display, call anything 
			// within 1 um of r=0 to be actually at r=0
			double rmin = map->GetRmin();
			if( fabs(rmin) < 1.0E-6) rmin = 0.0;

			cout << "Material map info from CCDB for run " << RUN_NUMBER << endl;
			cout<<"==============================================="<<endl;
			cout << "    namepath: " << map->GetNamepath() << endl;
			cout << "R range (cm): " << rmin << " - " <<  map->GetRmax() << endl;
			cout << "Z range (cm): " << map->GetZmin() << " - " <<  map->GetZmax() << endl;
			cout << "     density: " << node->Density <<" g/cm^3"<< endl;
			cout << " rad. length: " << node->RadLen <<" cm"<< endl;
			cout << "           A: " << node->A << endl;
			cout << "           Z: " << node->Z << endl;
		}
	}
	
	cout << endl;

	return 0;
}


//------------------
// init_runtime_xml
//------------------
void init_runtime_xml(void)
{
	cout<<endl;
	cout<<"=============================================================="<<endl;
	cout<<"Enabling dynamic geometry rendering"<<endl;
	cout<<"- - - - - - - - - - - - - - - - - - -"<<endl;
	cout<<endl;
	
	// Check if a shared object file already exists and
	// if we can attach it. Do this by trying to get the
	// md5 checksum from it. If we succeed, get the md5
	// checksum for the XML so they can be compared
	string fname = "./tmp_hddsroot.so";
	string md5_shared = "";
	string md5_xml = "";
	void *handle = dlopen(fname.c_str(), RTLD_NOW | RTLD_GLOBAL);
	if(handle){
		const char* (*md5geom_ext)(void);
		*(void **) (&md5geom_ext) = dlsym(handle, "md5geom_ext");
		char *err = dlerror();
		if(err == NULL){
			md5_shared = (*md5geom_ext)();
			md5_xml = md5geom_xml();
		}
		dlclose(handle);
	}
	
	// Regenerate shared object if needed
	if(md5_shared.length()>0){
		cout<<"found existing shared object"<<endl;
		cout<<"   shared object checksum: "<<md5_shared<<endl;
		cout<<"             xml checksum: "<<md5_xml<<endl;
		if(md5_shared != md5_xml){
			cout<<"Checksums don't match. Shared object will be regenerated."<<endl;
			MakeSharedObjectFromXML();
		}
	}else{
		MakeSharedObjectFromXML();
	}

	// Attach shared object
	cout<<endl;
	cout << "Attaching shared object ..." << endl;
	handle = dlopen(fname.c_str(), RTLD_NOW | RTLD_GLOBAL);
	if(!handle){
		cerr<<"Unable to open \""<<fname<<"\"!"<<endl;
		cerr<<dlerror()<<endl;
		exit(-1);
	}
	
	// Copy handle to global variable
	dlgeom_handle = handle;
	
	// do not close shared object since it may contain needed routines
	//dlclose(dlgeom_handle);
}

//------------------
// MakeSharedObjectFromXML
//------------------
void MakeSharedObjectFromXML(void)
{
	cout<<"Please make sure root-config is in your PATH and that"<<endl;
	cout<<"the following environment variables are set:" <<endl;
	cout<<"   HDDS_HOME  "<< endl;
	cout<<"   BMS_OSNAME "<< endl;
	
	// Open temporary file and add some includes so it will compile
	ofstream ofs("tmp_hddsroot.cc");
	ofs << "#include <TSystem.h>"<<endl;
	ofs << "#include <TGeoManager.h>"<<endl;
	ofs << "#include <TGeoVolume.h>"<<endl;
	ofs << "#include <TGeoMaterial.h>"<<endl;
	ofs << "#include <TGeoMedium.h>"<<endl;
	ofs << "#include <TGeoPcon.h>"<<endl;
	ofs << "#include <TGeoPgon.h>"<<endl;
	ofs << "#include <TGeoMatrix.h>"<<endl;
	ofs << "extern \"C\" {"<<endl;
	ofs << "TGeoManager* hddsroot(void);"<<endl;
	ofs.close();

	// Generate C++ code from XML
	cout<<endl;
	cout << "Generating C++ from XML source ...." << endl;
	string cmd = "$HDDS_HOME/bin/$BMS_OSNAME/hdds-root_h " + HDDS_XML + " >> tmp_hddsroot.cc";
	cout << cmd << endl;
	int err = system(cmd.c_str());
	if(err!=0){ cerr << "Error running \""<< cmd << "\"" << endl; exit(-1); }
	
	// Close off file with externally accessible md5geom wrapper
	ofs.open("tmp_hddsroot.cc", ios_base::app);
	ofs << "const char* md5geom_ext(void){return md5geom();}"<<endl;
	ofs << "}"<<endl;	// close off extern "C"
	ofs.close();
	
	// Compile C++ into shared object
	cout<<endl;
	cout << "Compiling C++ into shared object ..." << endl;
	cmd = "c++ -shared -fPIC -o tmp_hddsroot.so `root-config --cflags --libs` -lGeom tmp_hddsroot.cc";
	cout << cmd << endl;
	err = system(cmd.c_str());
	if(err!=0){ cerr << "Error running \""<< cmd << "\"" << endl; exit(-1); }
	
	// Remove temporary C++ source file
	unlink("./tmp_hddsroot.cc");
}

//------------------
// hddsroot_runtime
//------------------
TGeoManager* hddsroot_runtime(void)
{
	// This should only be called if the user specifies the 
	// "-xml" command line switch. Otherwise, the built-in
	// geometry is used.
	
	// Create, compile and link shared object. That part is done
	// in a separate routine so md5geom_runtime() can use it too.
	if(!dlgeom_handle) init_runtime_xml();
	
	// Find hddsroot symbol inside shared object
	cout<<endl;
	cout << "Locating geometry ... " << endl;
	TGeoManager* (*my_hddsroot)(void);
	*(void **) (&my_hddsroot) = dlsym(dlgeom_handle, "hddsroot");
	char *err = dlerror();
	if(err != NULL){
		cerr << err << endl;
		exit(-1);
	}
	
	// Execute my_hddsroot
	cout<<endl;
	cout << "Loading geometry ... " << endl;
	TGeoManager *geo = (*my_hddsroot)();

	cout<<endl;
	cout << "Geometry loaded successfully" << endl;
	cout<<"=============================================================="<<endl;

	return geo;
}

//------------------
// md5geom_runtime
//------------------
const char* md5geom_runtime(void)
{
	// This will extract the MD5 checksum of the 
	// geometry from the dynamically linked shared
	// object. It is called from the GetMD5Geom()
	// routine below. Use that routine to get the
	// checksum, not this one.

	// Create, compile and link shared object if needed.
	if(!dlgeom_handle) init_runtime_xml();

	// Grab md5geom routine from shared object
	const char* (*md5geom_ext)(void);
	*(void **) (&md5geom_ext) = dlsym(dlgeom_handle, "md5geom_ext");
	char *err = dlerror();
	if(err != NULL){
		cerr << err << endl;
		exit(-1);
	}
	
	// Execute my_md5geom_
	return (*md5geom_ext)();
}

//------------------
// md5geom_xml
//------------------
const char* md5geom_xml(void)
{
	/// This will get the checksum from the XML directly by
	/// running the hdds-md5 program and parsing the output.
	static string md5_xml=""; 
	string fname = "tmp_hddsroot.md5";
	string cmd = "$HDDS_HOME/bin/$BMS_OSNAME/hdds-md5 " + HDDS_XML + " > " + fname;
	int err = system(cmd.c_str());
	if(err!=0){ cerr << "Error running \""<< cmd << "\"" << endl; exit(-1); }
	ifstream ifs(fname.c_str());
	if(ifs.is_open()){
		string str;
		while(ifs.good())ifs >> str;
		if(str.length()>=32)md5_xml = str.substr(str.length()-32);
		ifs.close();
	}
	unlink(fname.c_str());
	
	return md5_xml.c_str();
}


//------------------
// GetMD5Geom
//------------------
const char* GetMD5Geom(void)
{
	// Get the MD5 checksum of the geometry that will be
	// used for the simulation. This will retrieve the
	// geometry checksum from either what has been statically
	// linked in, or dynamically, whichever is being used.

	if(USE_XML){
		// Grab version from shared object
		return md5geom_runtime();
	}else{
		// Use compiled in version
		return md5geom();
	}
	
	return NULL;
}

//------------------
// ParseCommandLineArguments
//------------------
void ParseCommandLineArguments(int narg, char *argv[])
{
	vector<double> vals;
	bool print_xml_md5_checksum = false;

	for(int i=1; i<narg; i++){
		
		if(argv[i][0]=='-'){
			string arg(argv[i]);
			if(arg=="-h" || arg=="--help")Usage();

			if(arg=="-checksum" || arg=="--checksum")print_xml_md5_checksum = true;
			if(arg.find("-xml")==0){
				USE_XML = true;
				if(arg.find("=")!=string::npos){
					HDDS_XML = arg.substr(arg.find("=")+1);
				}
			}
			
		}else{
			vals.push_back(atof(argv[i]));
		}
	}

	// If user specified printing the checksum then 
	// do that and quit
	if(print_xml_md5_checksum){
		string checksum = GetMD5Geom();
		cout << "HDDS Geometry MD5 Checksum: " << checksum << endl;
		exit(0);
	}

	// Copy user-given coordinates to global variables
	if(vals.size() != 3)Usage();
	X_LAB = vals[0];
	Y_LAB = vals[1];
	Z_LAB = vals[2];
}

//------------------
// Usage
//------------------
void Usage(void)
{
	cout<<endl;
	cout<<"Usage:"<<endl;
	cout<<"    hd_geom_query [options] X Y Z"<<endl;
	cout<<endl;
	cout<<"Print the material properties for the specified point in lab"<<endl;
	cout<<" coordinates. Units of X,Y, and Z are cm."<<endl;
	cout<<endl;
	cout<<"By default, this uses the geometry in $HDDS/src/hddsroot.h"<<endl;
	cout<<"that was used to link this executable. The -xml switch may be used"<<endl;
	cout<<"to dynamically compile and link code generated from the XML at run"<<endl;
	cout<<"time. If an equals sign \"=\" follows the -xml switch then the"<<endl;
	cout<<"main_HDDS.xml file is taken from the remainder of that argument."<<endl;
	cout<<endl;
	cout<<"If the -xml switch is specified, then a file named \"tmp_hddsroot.so\""<<endl;
	cout<<"is searched for in the current directory. If found, it is opened and"<<endl;
	cout<<"the geometry checksum is read from it and compared to that of the XML"<<endl;
	cout<<"specified (which may be the default of $HDDS_HOME/main_HDDS.xml)."<<endl;
	cout<<"If the two match, then that shared object is used, bypassing the"<<endl;
	cout<<"(expensive) compilation phase. If the file is not present, is unreadable,"<<endl;
	cout<<"or the checksums don't match, then the shared object is automatically"<<endl;
	cout<<"(re)generated."<<endl;
	cout<<endl;
	cout<<" options:"<<endl;
	cout<<"    -h or --help          Print this usage statement"<<endl;
	cout<<"    -xml[=main_HDDS.xml]  Dynamically generate geometry"<<endl;
	cout<<"    -checksum             Print the MD5 checksum of the "<<endl;
	cout<<"                          geometry and exit"<<endl;
	cout<<endl;
	cout<<"If the -xml option is given and no file is specified,"<<endl;
	cout<<"then a value of: "<<HDDS_XML<<endl;
	cout<<"is used."<<endl;
	cout<<endl;

	exit(0);
}


