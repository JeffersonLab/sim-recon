


#include "hdview.h"
#include "hdv_mainframe.h"
#include "MyProcessor.h"

TCanvas *maincanvas=NULL;
DEventLoop *eventloop=NULL;
MyProcessor *myproc = NULL;

int main(int narg, char *argv[])
{
	// Instantiate an event loop object and initialize it
	eventloop = new DEventLoop(narg, argv);
	myproc = new MyProcessor();
	eventloop->AddProcessor(myproc);
	eventloop->Init();

	// Open Window
	TApplication app("HDView", &narg, argv);
	hdv_mainframe hdvmf(gClient->GetRoot(), 600, 500);
	
	// Hand control to ROOT event loop
	app.Run();
	
	// Close out event loop
	eventloop->Fini();

	return 0;
}
