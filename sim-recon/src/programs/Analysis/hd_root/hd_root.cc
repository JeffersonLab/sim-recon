// Author: Edward Brash February 15, 2005
//
//
// hd_root.cc
//

#include <dlfcn.h>

#include <TFile.h>

#include "MyProcessor.h"
#include "DANA/DApplication.h"
using namespace std;

typedef void SetTFilePtrAddress_t(TFile **);
TFile* tfilePtr = NULL;
string OUTPUT_FILENAME = "hd_root.root";
string COMMAND_LINE_OUTPUT_FILENAME = "";
bool filename_from_command_line = false;

void ParseCommandLineArguments(int &narg, char *argv[]);
void DecideOutputFilename(void);
void Usage(void);


//-----------
// main
//-----------
int main(int narg, char *argv[])
{
	// Parse the command line
	ParseCommandLineArguments(narg, argv);

	// Instantiate our event processor
	MyProcessor myproc;

	// Instantiate an event loop object
	DApplication app(narg, argv);
	
	// Decide on the output filename
	DecideOutputFilename();
	
	// Run though all events, calling our event processor's methods
	app.monitor_heartbeat = 0;
	app.Run(&myproc);
	
	return 0;
}


//-----------
// ParseCommandLineArguments
//-----------
void ParseCommandLineArguments(int &narg, char *argv[])
{
	if(narg==1)Usage();
	
	for(int i=1;i<narg;i++){
		if(argv[i][0] != '-')continue;
		switch(argv[i][1]){
			case 'h':
				Usage();
				break;
			case 'D':
				toprint.push_back(&argv[i][2]);
				break;
			case 'A':
				ACTIVATE_ALL = 1;
			case 'o':
				if(i>=narg-1){
					cerr<<"\"-o\" requires a filename!"<<endl;
					exit(-1);
				}
				COMMAND_LINE_OUTPUT_FILENAME = argv[i+1];
				filename_from_command_line = true;
				
				// Remove the "-o fname" arguments from file list so
				// JANA won't think the "fname" is an input file.
				for(int j=i; j<(narg-2); j++)argv[j] = argv[j+2];
				narg -= 2;
				break;
		}
	}
}

//-----------
// DecideOutputFilename
//-----------
void DecideOutputFilename(void)
{
	/// Decide on the output filename to use based on the command line
	/// input and configuration parameter input. The command line takes
	/// precedence. This also makes sure to copy the filename that is
	/// being used into the configuration parameter.

	// Set the default output filename (avoids later warnings from JANA)
	gPARMS->SetDefaultParameter("OUTPUT_FILENAME", OUTPUT_FILENAME,"Output filename used by hd_root");
	
	// If the user specified an output filename on the command line,
	// use it to overwrite the config. parameter/default one
	if(filename_from_command_line){
		OUTPUT_FILENAME = COMMAND_LINE_OUTPUT_FILENAME;
	
		// Set the actual output filename in config. param.
		gPARMS->SetParameter("OUTPUT_FILENAME", OUTPUT_FILENAME);
	}

	jout<<"OUTPUT_FILENAME: "<<OUTPUT_FILENAME<<endl;

}

//-----------
// Usage
//-----------
void Usage(void)
{
	cout<<"Usage:"<<endl;
	cout<<"       hd_root [options] source1 source2 ..."<<endl;
	cout<<endl;
	cout<<"Process events from a Hall-D data source (e.g. a file)"<<endl;
	cout<<"This will create a ROOT file that plugins or debug histos"<<endl;
	cout<<"can write into."<<endl;
	cout<<endl;
	cout<<"Options:"<<endl;
	cout<<endl;
	cout<<"   -h        Print this message"<<endl;
	cout<<"   -Dname    Activate factory for data of type \"name\" (can be used multiple times)"<<endl;
	cout<<"   -A        Activate factories (overrides and -DXXX options)"<<endl;
	cout<<"   -o fname  Set output filename (default is \"hd_root.root\")"<<endl;
	cout<<endl;

	exit(0);
}
