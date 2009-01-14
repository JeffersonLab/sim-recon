#include "DANA/DApplication.h"
#include "MyProcessor.h"

//-----------
// main
//-----------
int main(int narg, char *argv[])
{
	// Instantiate an event loop object
	DApplication app(narg, argv);
	app.monitor_heartbeat = false;
	app.SetShowTicker(0);

	// Run though all events, calling our event processor's methods
        app.AddProcessor(new MyProcessor());
	app.Run(NULL, 1);
	
	return 0;
}

// end of source file
