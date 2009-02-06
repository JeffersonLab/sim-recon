
// This is s simple example program that shows how to use a DTrackResolution
// class to apply detector derived resolutions to charged tracks.
// There are detailed comments that should step you through the program
// with enough information to understand it and allow you to modify it to
// your own needs.


#include <iostream>
using namespace std;

#include <particleType.h>
#include <DANA/DApplication.h>

#include "DTrackingResolutionGEANT.h"
#include "DEventProcessor_HDParSim.h"

string OUTPUT_FILENAME = "hdparsim.root";

void Usage(DApplication &app);
void ParseCommandLineArguments(int &narg, char *argv[], DApplication &app);

//------------
// main
//------------
int main(int narg, char *argv[])
{
	// Instantiate an event loop object
	DApplication app(narg, argv);

	ParseCommandLineArguments(narg, argv, app);
	if(narg<=1)Usage(app);

	// Run though all events, calling our event processor's methods
	app.Run(new DEventProcessor_HDParSim(OUTPUT_FILENAME.c_str()), 1);

	return 0;
}

//-----------
// ParseCommandLineArguments
//-----------
void ParseCommandLineArguments(int &narg, char *argv[], DApplication &app)
{
	for(int i=1;i<narg;i++){
		if(argv[i][0] != '-')continue;
		switch(argv[i][1]){
			case 'h':
				Usage(app);
				break;
			case 'o':
				if(i>=narg-1){
					cerr<<"\"-o\" requires a filename!"<<endl;
					exit(-1);
				}
				OUTPUT_FILENAME = argv[i+1];
				break;
		}
	}
}

//-----------
// Usage
//-----------
void Usage(DApplication &app)
{
	cout<<"Usage:"<<endl;
	cout<<"       hdparsim [options] source1 source2 ..."<<endl;
	cout<<endl;
	cout<<"Parmetric simulation of Hall-D GlueX detector."<<endl;
	cout<<"Read generated events from an HDDM file and apply resolutions"<<endl;
	cout<<"and acceptances to "<<endl;
	cout<<"can write into."<<endl;
	cout<<endl;
	cout<<"Options:"<<endl;
	cout<<endl;
	cout<<"   -h        Print this message"<<endl;
	cout<<"   -o fname  Set output filename (default is \"hdparsim.root\")"<<endl;
	cout<<endl;
	
	cout<<"JANA options:"<<endl;
	app.Usage();
	cout<<endl;

	exit(0);
}
