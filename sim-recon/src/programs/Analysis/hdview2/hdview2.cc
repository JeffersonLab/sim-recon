
#include <unistd.h>
#include <pthread.h>
#include <TGApplication.h>

#include "hdview2.h"
#include "MyProcessor.h"

int GO = 0; // 1=continuously display events 0=wait for user
bool PRINT_FACTORY_LIST = false;

TCanvas *maincanvas=NULL;
extern JApplication *japp;
JEventLoop *eventloop =NULL;
MyProcessor *myproc = NULL;

void PrintFactoryList(JApplication *japp);
void ParseCommandLineArguments(int &narg, char *argv[], JApplication *japp);
void Usage(JApplication *japp);

//-------------------
// main
//-------------------
int main(int narg, char *argv[])
{
	// Parse the command line arguments
	ParseCommandLineArguments(narg, argv, japp);
	
	// Instantiate a DApplication object. This has to be done BEFORE
	// creating the TApplication object since that modifies the argument list.
	japp = new DApplication(narg, argv);
	
	// Create a ROOT TApplication object
	TApplication app("HDView", &narg, argv);
	
	// This is done AFTER creating the TApplication object so when the
	// init routine is called, the window will be mapped and it can
	// draw the detectors.
	myproc = new MyProcessor();
	japp->AddProcessor(myproc);

	// Call the JApplication's Init() routine to attach any plugins
	japp->Init();

	// Create the JEventLoop object explicitly.
	eventloop = new JEventLoop(japp);
	eventloop->SetAutoFree(0); // prevent auto-freeing of event after OneEvent is called



	// We need to re-call myproc->init here (it was already called from the japp->Init()
	// call above). This is because the previous call was done before the JEventLoop
	// object existed and so the hdv_mainframe object wasn't created etc... We have
	// to call japp->Init() before instantiated the JEventLoop so that the plugins
	// are attached and the JEventLoop can capture the full list of factories and
	// processors.
	myproc->init();

	// If the PRINT_FACTORY_LIST flag was set, then print the factory list
	if(PRINT_FACTORY_LIST)PrintFactoryList(japp);

	// Process the first event
	eventloop->OneEvent();

	// Hand control to the ROOT "event" loop
	japp->GetJParameterManager()->PrintParameters();
	app.SetReturnFromRun(true);
	app.Run();
	
	// Clean-up the app (call erun and fini methods, delete sources)
	//japp->Fini(); // This now actually done in hdv_mainframe::DoQuit()
		
	if(japp)delete japp;

	return 0;
}

//-----------
// PrintFactoryList
//-----------
void PrintFactoryList(JApplication *japp)
{
	// When we get here, the Run() method hasn't been
	// called so the JEventLoop objects haven't
	// been created yet and cansequently the factory objects
	// don't yet exist. Since we want the "list factories"
	// option to work even without an input file, we need
	// to first make the factories before we can list them.
	// To do this we only need to instantiate a JEventLoop object
	// passing it our "app" pointer. The JEventLoop will automatically
	// register itself with the DApplication and the factories
	// will be made, even ones from plugins passed on the command
	// line.
	vector<JEventLoop*> loops = japp->GetJEventLoops();
	if(loops.size()<1){
		_DBG_<<"No JEventLoops in japp!"<<endl;
		return;
	}
	JEventLoop *loop = loops[0];
	
	// Print header
	cout<<endl;
	cout<<"  Factory List"<<endl;
	cout<<"-------------------------"<<endl;
	
	// Get list of factories from the JEventLoop and loop over them
	// Printing out the data types and tags.
	vector<JFactory_base*> factories = loop->GetFactories();
	vector<JFactory_base*>::iterator iter = factories.begin();
	for(; iter!=factories.end(); iter++){
		cout<<" "<<(*iter)->GetDataClassName();
		if(strlen((*iter)->Tag()) !=0){
			cout<<" : "<<(*iter)->Tag();
		}
		cout<<endl;
	}
	cout<<endl;
	cout<<" "<<factories.size()<<" factories registered"<<endl;
	cout<<endl;
	
}


//-----------
// ParseCommandLineArguments
//-----------
void ParseCommandLineArguments(int &narg, char *argv[], JApplication *japp)
{
	if(narg==1)Usage(japp);

	for(int i=1;i<narg;i++){
		if(argv[i][0] != '-')continue;
		switch(argv[i][1]){
			case 'h':
				Usage(japp);
				break;
			case 'L':
				PRINT_FACTORY_LIST = true;
				//PrintFactoryList(japp);
				break;
		}
	}
}

//-----------
// Usage
//-----------
void Usage(JApplication *japp)
{
	cout<<"Usage:"<<endl;
	cout<<"       hdview2 [options] source1 source2 ..."<<endl;
	cout<<endl;
	cout<<"Print the contents of a Hall-D data source (e.g. a file)"<<endl;
	cout<<"to the screen."<<endl;
	cout<<endl;
	cout<<"Options:"<<endl;
	cout<<endl;
	cout<<"   -h        Print this message"<<endl;
	cout<<"   -Dname    Print the data of type \"name\" (can be used multiple times)"<<endl;
	cout<<"   -A        Print ALL data types (overrides and -DXXX options)"<<endl;
	cout<<"   -L        List available factories and exit"<<endl;
	cout<<"   -p        Don't pause for keystroke between events (def. is to pause)"<<endl;
	cout<<"   -s        Skip events which don't have any of the specified data types"<<endl;
	cout<<endl;
	cout<<"JANA options:"<<endl;
	cout<<endl;
	japp->Usage();
	cout<<endl;

	exit(0);
}
