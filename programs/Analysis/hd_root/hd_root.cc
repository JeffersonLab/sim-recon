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

void ParseCommandLineArguments(int &narg, char *argv[]);
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
				OUTPUT_FILENAME = argv[i+1];
_DBG_<<"OUTPUT_FILENAME="<<OUTPUT_FILENAME<<endl;
				break;
		}
	}
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
