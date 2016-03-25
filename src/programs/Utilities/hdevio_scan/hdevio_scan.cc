
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <thread>
using namespace std;


#include <DAQ/HDEVIO.h>

class WorkerThread;


vector<string> filenames;



//----------------
// Usage
//----------------
void Usage(string mess="")
{


	if(mess != "") cout << endl << mess << endl << endl;
	
	exit(0);
}

//----------------
// ParseCommandLineArguments
//----------------
void ParseCommandLineArguments(int narg, char *argv[])
{

	if(narg<2) Usage("You must supply a filename!");

	for(int i=1; i<narg; i++){
	
		filenames.push_back(argv[i]);
	}
}

//----------------
// main
//----------------
int main(int narg, char *argv[])
{

	ParseCommandLineArguments(narg, argv);
	
	// Loop over input files
	for(uint32_t i=0; i<filenames.size(); i++){
		string &filename = filenames[i];
		cout << "Processing file " << (i+1) << "/" << filenames.size() << " : " << filename << endl;

		HDEVIO *hdevio = new HDEVIO(filename);
		if(!hdevio->is_open){
			cout << hdevio->err_mess.str() << endl;
			continue;
		}
		
		time_t start_time = time(NULL);
		hdevio->PrintFileSummary();
		time_t end_time = time(NULL);
		
		delete hdevio;

		cout << (end_time - start_time) << " sec " << endl;
	}

	return 0;
}


