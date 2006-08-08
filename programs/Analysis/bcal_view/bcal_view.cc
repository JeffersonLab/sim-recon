

#include <unistd.h>

#include <TApplication.h>
#include <TCanvas.h>
#include <TGeometry.h>


#include <unistd.h>
#include <pthread.h>
#include <TGApplication.h>

#include "bcview.h"
#include "bcv_mainframe.h"
#include "MyProcessor.h"

bool GO = false; // true=continuously display events false=wait for user

TCanvas *maincanvas=NULL;
extern JApplication *japp;
JEventLoop *eventloop =NULL;
MyProcessor *myproc = NULL;
bcv_mainframe *bcvmf=NULL;

//-------------------
// main
//-------------------
int main(int narg, char *argv[])
{
	// Instantiate a JApplication object this has to be done BEFORE
	// creating the TApplication object since that modifes the argument list.
	japp = new DApplication(narg, argv);
	
	// Open Window
	TApplication app("BCAL_Module_Viewer", &narg, argv);
	bcvmf = new bcv_mainframe(gClient->GetRoot(), 600, 400);
	
	// This is done AFTER creating the TApplication object so when the
	// init routine is called, the window will be mapped and it can
	// draw the detectors.
	myproc = new MyProcessor();
	japp->AddProcessor(myproc);
	eventloop = new JEventLoop(japp);
	japp->Init();

	// Hand control to ROOT event loop
	app.Run();
	
	// Clean-up the app (call erun and fini methods, delete sources)
	japp->Fini();
		
	delete japp;

	return 0;
}


