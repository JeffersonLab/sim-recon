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

//int COUNT;

bool LIST_FACTORIES = false;

int VERBOSE;
float MASS;
float SMEARWEIGHT;
char *OUTNAME;

//-----------
// main
//-----------
int main(int narg, char *argv[])
{
  VERBOSE = 0;
  SMEARWEIGHT = 0.01;
  // Parse the command line
  ParseCommandLineArguments(narg, argv);

  // Instantiate our event processor
  MyProcessor myproc;

  // Instantiate an event loop object
  DApplication *app = new DApplication(narg, argv);

  // If LIST_FACTORIES is set, print all factories and exit
  if(LIST_FACTORIES){
    PrintFactoryList(app);
    return 0;
  }

  COUNT = 0;
  // Run though all events, calling our event processor's methods
  app->SetShowTicker(0);
  app->monitor_heartbeat = false;
  app->Run(&myproc);

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
    cout<<" "<<(*iter)->dataClassName();
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
      case 'v':
        VERBOSE++;
        cerr << "Verbosity: " << VERBOSE << endl;
        break;
      case 'M':
        MASS = atof(&argv[i][2]);
        cerr << "Mass to which to fit: " << MASS << endl;
        break;
      case 'm':
        MAX_EVENTS = atof(&argv[i][2]);
        cerr << "Max events: " << MAX_EVENTS << endl;
        break;
      case 'p':
        PAUSE_BETWEEN_EVENTS = 0;
        break;
      case 's':
        SKIP_BORING_EVENTS = 1;
        break;
      case 'A':
        PRINT_ALL = 1;
      case 'W':
        SMEARWEIGHT = atof(&argv[i][2]);
        cerr << "Inverse of the error matrix weight: " << SMEARWEIGHT << endl;
        break;
        break;
      case 'L':
        LIST_FACTORIES = 1;
        break;
      case 'o':
        OUTNAME = &argv[i][2];
        break;
    }
  }
}

//-----------
// Usage
//-----------
void Usage(void)
{
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
  cout<<endl;

  exit(0);
}


