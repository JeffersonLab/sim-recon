// Author: David Lawrence  June 25, 2004
//
//
// hd_ana.cc
//

#include "DANA/DApplication.h"
using namespace std;

void Usage(JApplication &app);


//-----------
// main
//-----------
int main(int narg, char *argv[])
{
	// Instantiate an event loop object
	DApplication app(narg, argv);

	if(narg<=1)Usage(app);

	// Run though all events, calling our event processor's methods
	app.Run(NULL, 1);
	
	return 0;
}

//-----------
// Usage
//-----------
void Usage(JApplication &app)
{
	cout<<endl;
	cout<<"Usage:"<<endl;
	cout<<"    hd_ana [options] source1 source2 source3 ..."<<endl;
	cout<<endl;
	app.Usage();
	cout<<endl;
	
	exit(0);
}

