
#include <TGApplication.h>

#include "hdview.h"
#include "hdv_mainframe.h"
#include "MyProcessor.h"

TCanvas *maincanvas=NULL;
DEventLoop *eventloop=NULL;
MyProcessor *myproc = NULL;
hdv_mainframe *hdvmf=NULL;

//-------------------
// main
//-------------------
int main(int narg, char *argv[])
{
	// Instantiate an event loop object and initialize it
	eventloop = new DEventLoop(narg, argv);
	myproc = new MyProcessor();
	eventloop->AddProcessor(myproc);
	eventloop->Init();

	// Open Window
	TApplication app("HDView", &narg, argv);
	hdvmf = new hdv_mainframe(gClient->GetRoot(), 1000, 600);
	
	// Hand control to ROOT event loop
	app.Run();
	
	// Close out event loop
	eventloop->Fini();

	return 0;
}

//-------------------
// hdv_getevent
//-------------------
derror_t hdv_getevent(void)
{
	// Read in next event. 
	derror_t err;
	err = eventloop->OneEvent();
	if(err!=NOERROR)return err;
		

	return NOERROR;
}
