// Author: David Lawrence  June 25, 2004
//
// modified June 6, 2006
//
// mctrk_ana.cc
//

#include <dlfcn.h>

#include <TFile.h>

#include "DANA/DApplication.h"

#include "MyProcessor.h"

//-----------
// main
//-----------
int main(int narg, char *argv[])
{	
	// Instantiate an event loop object
	DApplication app(narg, argv);

	// Instantiate our event processor
	MyProcessor myproc;
	app.AddProcessor(&myproc);
	
	// Always use the track_hists.so plugin
	app.AddPlugin("track_hists");

	// Run though all events, calling our event processor's methods
	app.monitor_heartbeat = false;
	app.Run(NULL,1);
	
	return 0;
}

