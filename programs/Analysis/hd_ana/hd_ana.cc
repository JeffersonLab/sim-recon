// Author: David Lawrence  June 25, 2004
//
//
// hd_ana.cc
//

#include "DANA/DApplication.h"

//-----------
// main
//-----------
int main(int narg, char *argv[])
{
	// Instantiate an event loop object
	DApplication app(narg, argv);

	// Run though all events, calling our event processor's methods
	app.Run(NULL, 1);
	
	return 0;
}

