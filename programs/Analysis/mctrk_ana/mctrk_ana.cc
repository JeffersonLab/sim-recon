// Author: David Lawrence  June 25, 2004
//
// modified June 6, 2006
//
// mctrk_ana.cc
//

#include <dlfcn.h>

#include <TFile.h>

#include "MyProcessor.h"
#include "DApplication.h"

typedef void SetTFilePtrAddress_t(TFile **);
TFile* tfilePtr = NULL;

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

	// Get the list of shared objects (if any) and try setting their
	// TFile pointer address globals so that they use the same
	// ROOT file for output.
	vector<void*> sohandles = app.GetSharedObjectHandles();
	for(unsigned int i=0; i<sohandles.size(); i++){
		SetTFilePtrAddress_t *SetTFilePtrAddress = (SetTFilePtrAddress_t*)dlsym(sohandles[i], "SetTFilePtrAddress");
		if(SetTFilePtrAddress)(*SetTFilePtrAddress)(&tfilePtr);
	}

	// Run though all events, calling our event processor's methods
	app.monitor_heartbeat = false;
	app.Run(NULL,1);
	
	return 0;
}

