
#include <unistd.h>
#include <pthread.h>
#include <TGApplication.h>

#include "hdview.h"
#include "hdv_mainframe.h"
#include "MyProcessor.h"

int GO = 0; // 1=continuously display events 0=wait for user

TCanvas *maincanvas=NULL;
extern JApplication *japp;
JEventLoop *eventloop =NULL;
MyProcessor *myproc = NULL;
hdv_mainframe *hdvmf=NULL;

//-------------------
// main
//-------------------
int main(int narg, char *argv[])
{
	// Instantiate a DApplication object this has to be done BEFORE
	// creating the TApplication object since that modifes the argument list.
	japp = new DApplication(narg, argv);
	
	// Open Window
	TApplication app("HDView", &narg, argv);
	hdvmf = new hdv_mainframe(gClient->GetRoot(), 1000, 600);
	
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

