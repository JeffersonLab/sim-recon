// Author: David Lawrence  June 25, 2004
//
//
// hd_ana.cc
//

#include <iostream>
using namespace std;

#include <termios.h>

#include "MyProcessor.h"
#include "DApplication.h"

void ParseCommandLineArguments(int &narg, char *argv[]);
void Usage(void);


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
	DApplication app(narg, argv);

	// This monkey shines is needed to get getchar() to return single
	// characters without waiting for the user to hit return
	struct termios t;
	tcgetattr(fileno(stdin), &t);
	t.c_lflag &= (~ICANON);
	//t.c_cc[VMIN] = 1;
	tcsetattr(fileno(stdin), TCSANOW, &t);


	// Run though all events, calling our event processor's methods
	app.SetShowTicker(0);
	app.Run(&myproc);
	
	return 0;
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
	cout<<"   -p        Don't pause for keystroke between events (def. is to pause)"<<endl;
	cout<<"   -s        Skip events which don't have any of the specified data types"<<endl;
	cout<<endl;

	exit(0);
}


