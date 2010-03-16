// Author: David Lawrence  June 25, 2004
//
//
// hd_ana.cc
//

#include <iostream>
using namespace std;

#include <termios.h>

#include "MyProcessor.h"
#include "DANA/DApplication.h"
#include "HDDM/DEventSourceHDDMGenerator.h"

void PrintFactoryList(DApplication *app);
void ParseCommandLineArguments(int &narg, char *argv[]);
void Usage(void);

bool LIST_FACTORIES = false;

//-----------
// main
//-----------
int main(int narg, char *argv[])
{
	// Parse the command line
	ParseCommandLineArguments(narg, argv);

	// Instantiate our event processor
	MyProcessor myproc;

	// Instantiate an event loop object
	DApplication *app = new DApplication(narg, argv);
	
	// Set tag prefix for JANA streams to empty
	jout.SetTag("");
	
	// If LIST_FACTORIES is set, print all factories and exit
	if(LIST_FACTORIES){
		PrintFactoryList(app);
		return 0;
	}

	// This monkeyshines is needed to get getchar() to return single
	// characters without waiting for the user to hit return
	struct termios t;
	tcgetattr(fileno(stdin), &t);
	t.c_lflag &= (~ICANON);
	//t.c_cc[VMIN] = 1;
	tcsetattr(fileno(stdin), TCSANOW, &t);


	// Run though all events, calling our event processor's methods
	app->SetShowTicker(0);
	app->monitor_heartbeat = false;
	app->Run(&myproc);
	
	delete app;

	return 0;
}

//-----------
// PrintFactoryList
//-----------
void PrintFactoryList(DApplication *app)
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
	app->Init();
	JEventLoop *loop = new JEventLoop(app);
	
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
void ParseCommandLineArguments(int &narg, char *argv[])
{
	if(narg==1)Usage();

	for(int i=1;i<narg;i++){
		if(argv[i][0] != '-')continue;
		switch(argv[i][1]){
			case 'h':
				Usage();
				break;
			case 'D':
				toprint.push_back(&argv[i][2]);
				break;
			case 'p':
				PAUSE_BETWEEN_EVENTS = 0;
				break;
			case 's':
				SKIP_BORING_EVENTS = 1;
				break;
			case 'A':
				PRINT_ALL = 1;
				break;
			case 'L':
				LIST_FACTORIES = 1;
				break;
			case 'a':
				LIST_ASSOCIATED_OBJECTS = true;
				break;
			case 'f':
				PRINT_SUMMARY_HEADER = false;
				break;
		}
	}
}

//-----------
// Usage
//-----------
void Usage(void)
{
	DApplication dapp(0,NULL);

	cout<<"Usage:"<<endl;
	cout<<"       hd_dump [options] source1 source2 ..."<<endl;
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
	cout<<"   -a        List types and number of associated objects"<<endl;
	cout<<"   -f        Skip printing summary header lisiting all factories"<<endl;
	cout<<"             (This disables auto-activating every single factory)"<<endl;
	cout<<endl;

	dapp.Usage();

	exit(0);
}


