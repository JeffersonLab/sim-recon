// $Id: mcsmear.cc 2388 2007-01-10 16:46:03Z davidl $
//
// Created August 24, 2007  David Lawrence

#include <iostream>
#include <iomanip>
using namespace std;

#include <signal.h>
#include <time.h>

#include "HDDM/hddm_s.h"

bool Filter(s_HDDM_t *hddm_s);
void ParseCommandLineArguments(int narg, char* argv[]);
void Usage(void);
void ctrlCHandle(int x);

char *INFILENAME = NULL;
char *OUTFILENAME = NULL;
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
	cout<<" output file: "<<OUTFILENAME<<endl;
	
	// Open Input file
	s_iostream_t *fin = open_s_HDDM(INFILENAME);
	if(!fin){
		cout<<" Error opening input file \""<<INFILENAME<<"\"!"<<endl;
		exit(-1);
	}
	
	// Output file
	s_iostream_t *fout = init_s_HDDM(OUTFILENAME);
	if(!fout){
		cout<<" Error opening output file \""<<OUTFILENAME<<"\"!"<<endl;
		exit(-1);
	}
	
	// Loop over events in input file
	s_HDDM_t *hddm_s;
	int NEvents_read = 0;
	int NEvents_written = 0;
	time_t last_time = time(NULL);
	while((hddm_s = read_s_HDDM(fin))){
		NEvents_read++;
		time_t now = time(NULL);
		if(now != last_time){
			cout<<" "<<NEvents_read<<" events read -- "<<NEvents_written<<" events written      \r";cout.flush();
			last_time = now;
		}
		
		// Write or don't depending on return value of Filter()
		if(Filter(hddm_s)){
			flush_s_HDDM(hddm_s, fout);
			NEvents_written++;
		}
		
		if(QUIT)break;
	}
	
	// close input and output files
	close_s_HDDM(fin);
	close_s_HDDM(fout);

	cout<<endl<<"FINAL:"<<endl;
	cout<<" "<<NEvents_read<<" events read -- "<<NEvents_written<<" events written"<<endl;
	cout<<"Output file has "<<100.0*(double)NEvents_written/(double)NEvents_read<<"% of the events that were in the input file."<<endl;

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
	
	// Generate output filename based on input filename
	char *ptr, *path_stripped;
	path_stripped = ptr = strdup(INFILENAME);
	while((ptr = strstr(ptr, "/")))path_stripped = ++ptr;
	ptr = strstr(path_stripped, ".hddm");
	if(ptr)*ptr=0;
	char str[256];
	sprintf(str, "%s_filtered.hddm", path_stripped);
	OUTFILENAME = strdup(str);
}


//-----------
// Usage
//-----------
void Usage(void)
{
	cout<<endl<<"Usage:"<<endl;
	cout<<"     filtergen [options] file.hddm"<<endl;
	cout<<endl;
	cout<<" Read filter events from the given HDDM file while copying it"<<endl;
	cout<<"to another HDDM file. This was written to provide an easy way"<<endl;
	cout<<"to filter events from a file produced by pythia that one does"<<endl;
	cout<<"not want to waste time tracking through the whole detector."<<endl;
	cout<<endl;
	cout<<"At this point, the filtering conditions are hardwired and"<<endl;
	cout<<"can only be changed by editing the source. In the future, we"<<endl;
	cout<<"will want command line options so the filter can be tuned"<<endl;
	cout<<"without recompiling."<<endl;
	cout<<endl;
	cout<<"  options:"<<endl;
	cout<<"    -h       Print this usage statement."<<endl;
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
