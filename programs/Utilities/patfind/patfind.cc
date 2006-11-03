// Author: David Lawrence  June 25, 2004
//
//
// hd_ana.cc
//

#include <TApplication.h>

#include "MyMainFrame.h"
#include "MyProcessor.h"
#include "JANA/JEventLoop.h"


extern DApplication *japp;
JEventLoop *eventloop = NULL;
MyProcessor *myproc = NULL;
MyMainFrame *mmf = NULL;

//-----------
// main
//-----------
int main(int narg, char *argv[])
{
	// Instantiate a DApplication object this has to be done BEFORE
	// creating the TApplication object since that modifes the argument list.
	japp = new DApplication(narg, argv);
	
	// Open Window
	TApplication app("PatFind", &narg, argv);
	mmf = new MyMainFrame(gClient->GetRoot(), 1200, 1200);
	
	// Create a MyProcessor
	myproc = new MyProcessor();
	japp->AddProcessor(myproc);
	eventloop = new JEventLoop(japp);
	japp->Init();

	// Hand control to ROOT event loop
	app.Run();
	
	// Close out event loop
	japp->Fini();
	
	// clean up
	delete myproc;
	delete mmf;
	delete japp;
	
	return 0;
}

