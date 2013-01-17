
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <string>
using namespace std;

#include <DANA/DApplication.h>


// These are defined in copytoplusplus.cc
extern string INFILE;
extern string OUTFILE;
extern bool POSTSMEAR;
extern string MCSMEAROPTS;
extern bool DELETEUNSMEARED;

// Defined in calibDB.cc
extern string HDDS_XML;

// Declare routines callable from FORTRAN
extern "C" int gxint_(void);
extern "C" void init_runtime_xml_(void); // defined in dl_routines.cc
extern "C" const char* GetMD5Geom(void); // defined in calibDB.cc

void Usage(void);

// Get access to FORTRAN common block with some control flags
#include "controlparams.h"

//------------------
// main
//------------------
int main(int narg, char *argv[])
{
	// This is needed so calibDB.cc can use it to get the
	// JCalibration object pointer. We want this to be done
	// in the same way as all other sim-recon software
	DApplication *dapp = new DApplication(narg, argv);
	dapp->Init();

	// Set some defaults. Note that most defaults related to the
	// simulation are set in uginit.F
	controlparams_.runtime_geom = 0;

	// Parse command line parameters
	bool print_xml_md5_checksum = false;
	for(int i=1; i<narg; i++){
		string arg = argv[i];
		string next = i<(narg+1) ? argv[i]:"";
		
		if(arg=="-h" || arg=="--help")Usage();
		if(arg=="-checksum" || arg=="--checksum")print_xml_md5_checksum = true;
		if(arg.find("-xml")==0){
			controlparams_.runtime_geom = 1;
			if(arg.find("=")!=string::npos){
				HDDS_XML = arg.substr(arg.find("=")+1);
			}
		}
	}
	
	// If specified to read in XML geometry, do necessary
	// initializations now.
	if(controlparams_.runtime_geom != 0){
		init_runtime_xml_();
	}
	
	// If user specified printing the checksum then 
	// do that and quit
	if(print_xml_md5_checksum){
	
		cout << "HDDS Geometry MD5 Checksum: " << GetMD5Geom() << endl;

		return 0;
	}
	
	// Run gxint_ proper
	int res = gxint_();

	return res;
}


//------------------
// Usage
//------------------
void Usage(void)
{
	cout<<endl;
	cout<<"Usage:"<<endl;
	cout<<"   hdgeant++ [options]"<<endl;
	cout<<endl;
	cout<<" Hall-D interactive Monte Carlo simulation based on GEANT3."<<endl;
	cout<<"Most configurable options are set using a control.in"<<endl;
	cout<<"file. A well-annotated example control.in can be"<<endl;
	cout<<"found in the HDGeant source code directory."<<endl;
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

