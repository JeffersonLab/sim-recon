// Author: David Lawrence  June 25, 2004
//
// modified June 6, 2006
//
// mctrk_ana.cc
//

#include <string>
using namespace std;

#include <dlfcn.h>

#include <TFile.h>

#include "DANA/DApplication.h"

#include "MyProcessor.h"

const char* OUTPUTFILE = "mctrk_ana.root";

//-----------
// main
//-----------
int main(int narg, char *argv[])
{
	// Parse command line
	for(int i=1; i<narg; i++){
		if(string("-o") == argv[i]){
			i++;
			if(i>=narg){
				cout<<"You must supply an output filename when using \"-o\"!"<<endl;
				exit(-1);
			}
			if(argv[i][0]=='-'){
				cout<<"You must supply an output filename when using \"-o\"!"<<endl;
				exit(-1);
			}
			OUTPUTFILE = argv[i];
		}
	}

	// Instantiate an event loop object
	DApplication app(narg, argv);
	//app.SetShowTicker(0);

	// Instantiate our event processor
	//MyProcessor myproc;
	//app.AddProcessor(&myproc);
	
	// Add plugins
	app.AddPlugin("trackeff_hists");
	app.AddPlugin("acceptance_hists");
	app.AddPlugin("mcthrown_hists");

	// Run though all events, calling our event processor's methods
	app.monitor_heartbeat = false;
	app.Run(NULL,1);
	
	return 0;
}

