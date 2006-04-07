// Author: David Lawrence  June 25, 2004
//
//
// hd_ana.cc
//

#include <TApplication.h>

#include "MyMainFrame.h"
#include "MyProcessor.h"
#include "DEventLoop.h"


DApplication *dapp = NULL;
DEventLoop *eventloop = NULL;
MyProcessor *myproc = NULL;
MyMainFrame *mmf = NULL;

//-----------
// main
//-----------
int main(int narg, char *argv[])
{
	// Instantiate a DApplication object this has to be done BEFORE
	// creating the TApplication object since that modifes the argument list.
	dapp = new DApplication(narg, argv);
	
	// Open Window
	TApplication app("PatFind", &narg, argv);
	mmf = new MyMainFrame(gClient->GetRoot(), 1200, 1200);
	
	// Create a MyProcessor
	myproc = new MyProcessor();
	dapp->AddProcessor(myproc);
	eventloop = new DEventLoop(dapp);
	dapp->Init();

	// Hand control to ROOT event loop
	app.Run();
	
	// Close out event loop
	dapp->Fini();
	
	// clean up
	delete myproc;
	delete mmf;
	delete dapp;
	
	return 0;
}

