// Author: Edward Brash February 15, 2005
//
//
// hd_root.cc
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
	// Instantiate our event processor
	MyProcessor myproc;

	// Instantiate an event loop object
	DApplication app(narg, argv);
	
	// Get the list of shared objects (if any) and try setting their
	// TFile pointer address globals so that they use the same
	// ROOT file for output.
	vector<void*> sohandles = app.GetSharedObjectHandles();
	for(unsigned int i=0; i<sohandles.size(); i++){
		SetTFilePtrAddress_t *SetTFilePtrAddress = (SetTFilePtrAddress_t*)dlsym(sohandles[i], "SetTFilePtrAddress");
		if(SetTFilePtrAddress){
			(*SetTFilePtrAddress)(&tfilePtr);
		}
	}
	tfilePtr = NULL;

	// Run though all events, calling our event processor's methods
	app.Run(NULL);
	
	return 0;
}

