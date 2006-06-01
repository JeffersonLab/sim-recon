// Author: Edward Brash February 15, 2005
//
//
// hd_root.cc
//

#include "MyProcessor.h"
#include "DApplication.h"

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

