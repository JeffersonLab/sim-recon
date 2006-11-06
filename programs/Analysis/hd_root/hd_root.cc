// Author: Edward Brash February 15, 2005
//
//
// hd_root.cc
//

#include <dlfcn.h>

#include <TFile.h>

#include "MyProcessor.h"
#include "DANA/DApplication.h"

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

	// Run though all events, calling our event processor's methods
	app.Run(&myproc);
	
	return 0;
}

