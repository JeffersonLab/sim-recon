
#include <unistd.h>
#include <pthread.h>
#include <TGApplication.h>

#include "hdview.h"
#include "hdv_mainframe.h"
#include "MyProcessor.h"

void* EventUpdateThread(void *arg);
int DONE = 0;
int GO = 1;

TCanvas *maincanvas=NULL;
DEventLoop *eventloop=NULL;
MyProcessor *myproc = NULL;
hdv_mainframe *hdvmf=NULL;

//-------------------
// main
//-------------------
int main(int narg, char *argv[])
{
	// Instantiate an event loop object this has to be done BEFORE
	// creating the TApplication object since that modifes the argument list.
	eventloop = new DEventLoop(narg, argv);

	// Open Window
	TApplication app("HDView", &narg, argv);
	hdvmf = new hdv_mainframe(gClient->GetRoot(), 1000, 600);
	
	// This is done AFTER creating the TApplication object so when the
	// init routine is called, the window will be mapped and it can
	// draw the detectors.
	myproc = new MyProcessor();
	eventloop->AddProcessor(myproc);
	eventloop->Init();
	
	// create a thread to auto-advance events
	pthread_t thr;
	pthread_create(&thr, NULL,EventUpdateThread, NULL);

	// Hand control to ROOT event loop
	app.Run();
	
	DONE=1;
	void *retval;
	pthread_join(thr, &retval);
	
	// Close out event loop
	eventloop->Fini();

	return 0;
}

//-------------------
// EventUpdateThread
//-------------------
void* EventUpdateThread(void *arg)
{
	int Niterations = 0;
	while(!DONE){
		sleep(1);
		
		if(GO)hdv_getevent();
	}
	
	pthread_exit(NULL);
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
