#include "FillHitsProc.h"
#include "DANA/DApplication.h"

//-----------
// main
//-----------
int main(int narg, char *argv[])
{
	// Instantiate our event processor
	FillHitsProc myproc;
	
	// Instantiate an event loop object
	DApplication app(narg, argv);
	
	// Run though all events, calling our event processor's methods
	app.Run(&myproc, 1);
	
	return 0;
}
