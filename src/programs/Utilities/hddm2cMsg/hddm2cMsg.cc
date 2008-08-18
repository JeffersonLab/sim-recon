// $Id$
//
// Created June 22, 2005  David Lawrence

#include <iostream>
#include <iomanip>
using namespace std;

#include <signal.h>
#include <time.h>

#include <cMsg.hxx>

#include "hddm_s.h"

void ParseCommandLineArguments(int narg, char* argv[]);
void Usage(void);
void ctrlCHandle(int x);

char *INFILENAME = NULL;
char *UDL = "cMsg:cMsg//localhost:3456";
int QUIT = 0;

//-----------
// main
//-----------
int main(int narg,char* argv[])
{
	// Set up to catch SIGINTs for graceful exits
	signal(SIGINT,ctrlCHandle);

	ParseCommandLineArguments(narg, argv);
	
	cout<<" input file: "<<INFILENAME<<endl;
	cout<<" output UDL: "<<UDL<<endl;
	
	// Open cMsg output
	cMsg cMsgSys(UDL, "hddm2cMsg", "Feed events from HDDM file to cMsg");
	try{
		cout<<"Connecting to cMsg system ..."<<endl;
		cMsgSys.connect();
	}catch(cMsgException e){
		cout<<e.toString()<<endl;
	}
	
	cMsgSys.disconnect();
	exit(0);

	// Open Input file
	s_iostream_t *fin = open_s_HDDM(INFILENAME);
	if(!fin){
		cout<<" Error opening input file \""<<INFILENAME<<"\"!"<<endl;
		exit(-1);
	}
		
	// Loop over events in input file
	s_HDDM_t *hddm_s;
	int NEvents = 0;
	time_t last_time = time(NULL);
	while((hddm_s = read_s_HDDM(fin))){
		NEvents++;
		time_t now = time(NULL);
		if(now != last_time){
			cout<<"  "<<NEvents<<" events processed      \r";cout.flush();
			last_time = now;
		}
		
		
		if(QUIT)break;
	}
	
	// close input and output files
	close_s_HDDM(fin);

	cout<<" "<<NEvents<<" events read"<<endl;

	return 0;
}

//-----------
// ParseCommandLineArguments
//-----------
void ParseCommandLineArguments(int narg, char* argv[])
{

	for(int i=1; i<narg; i++){
		char *ptr = argv[i];
		
		if(ptr[0] == '-'){
			switch(ptr[1]){
				case 'h': Usage();													break;
			}
		}else{
			INFILENAME = argv[i];
		}
	}

	if(!INFILENAME){
		cout<<endl<<"You must enter a filename!"<<endl<<endl;
		Usage();
	}
}


//-----------
// Usage
//-----------
void Usage(void)
{
	cout<<endl<<"Usage:"<<endl;
	cout<<"     hddm2cMsg [options] file.hddm"<<endl;
	cout<<endl;
	cout<<" Read the given, Geant-produced HDDM file as input and send"<<endl;
	cout<<"the events to the specified cMsg UDL."<<endl;
	cout<<endl;
	cout<<"  options:"<<endl;
	cout<<"    -h       Print this usage statement."<<endl;
	cout<<endl;
	cout<<" Example:"<<endl;
	cout<<endl;
	cout<<"     hddm2cMsg hdgeant.hddm"<<endl;
	cout<<endl;
	cout<<endl;

	exit(0);
}

//-----------------------------------------------------------------
// ctrlCHandle
//-----------------------------------------------------------------
void ctrlCHandle(int x)
{
	QUIT++;
	cerr<<endl<<"SIGINT received ("<<QUIT<<")....."<<endl;
}
