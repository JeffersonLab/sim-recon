// $Id$
//
// Created Dec 22, 2007  David Lawrence

#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

#include <signal.h>
#include <time.h>

#include "HDDM/hddm_s.h"

void Smear(s_HDDM_t *hddm_s);
void ParseCommandLineArguments(int narg, char* argv[]);
void Usage(void);
void ctrlCHandle(int x);

vector<char*> INFILENAMES;
char *OUTFILENAME = NULL;
int QUIT = 0;


#define _DBG_ cout<<__FILE__<<":"<<__LINE__<<" "
#define _DBG__ cout<<__FILE__<<":"<<__LINE__<<endl


//-----------
// main
//-----------
int main(int narg,char* argv[])
{
	// Set up to catch SIGINTs for graceful exits
	signal(SIGINT,ctrlCHandle);

	ParseCommandLineArguments(narg, argv);
	
	
	// Output file
	cout<<" output file: "<<OUTFILENAME<<endl;
	s_iostream_t *fout = init_s_HDDM(OUTFILENAME);
	if(!fout){
		cout<<" Error opening output file \""<<OUTFILENAME<<"\"!"<<endl;
		exit(-1);
	}

	// Loop over input files
	int NEvents = 0;
	int NEvents_read = 0;
	time_t last_time = time(NULL);
	for(unsigned int i=0; i<INFILENAMES.size(); i++){
		cout<<" input file: "<<INFILENAMES[i]<<endl;
		s_iostream_t *fin = open_s_HDDM(INFILENAMES[i]);
		if(!fin){
			cout<<" Error opening input file \""<<INFILENAMES[i]<<"\"!"<<endl;
			exit(-1);
		}
			
		// Loop over all events in input
		while(true){
			s_HDDM_t *hddm_s = read_s_HDDM(fin);
			if(!hddm_s)break;
			NEvents_read++;
			
			// Write this output event to file and free its memory
			flush_s_HDDM(hddm_s, fout);
			NEvents++;
		
			// Update ticker
			time_t now = time(NULL);
			if(now != last_time){
				cout<<"  "<<NEvents_read<<" events read     ("<<NEvents<<" event written) \r";cout.flush();
				last_time = now;
			}

			if(QUIT)break;
		}

		// Close input file
		close_s_HDDM(fin);
	}
		
	// Close output file
	close_s_HDDM(fout);
	
	cout<<" "<<NEvents<<" events read"<<endl;

	return 0;
}

//-----------
// ParseCommandLineArguments
//-----------
void ParseCommandLineArguments(int narg, char* argv[])
{
	INFILENAMES.clear();

	for(int i=1; i<narg; i++){
		char *ptr = argv[i];
		
		if(ptr[0] == '-'){
			switch(ptr[1]){
				case 'h': Usage();						break;
				case 'o': OUTFILENAME=&ptr[2];		break;
			}
		}else{
			INFILENAMES.push_back(argv[i]);
		}
	}

	if(INFILENAMES.size()==0){
		cout<<endl<<"You must enter a filename!"<<endl<<endl;
		Usage();
	}
	
	if(OUTFILENAME==NULL){
		char *default_outfile = (char*)malloc(20);
		OUTFILENAME = strcpy(default_outfile,"merged_files.hddm");
	}
}


//-----------
// Usage
//-----------
void Usage(void)
{
	cout<<endl<<"Usage:"<<endl;
	cout<<"     hddm_merge_files [-oOutputfile] file1.hddm file2.hddm ..."<<endl;
	cout<<endl;
	cout<<"options:"<<endl;
	cout<<"    -oOutputfile  Set output filename (def. merged_files.hddm)"<<endl;
	cout<<endl;
	cout<<" This will merge 1 or more HDDM files into a single HDDM file."<<endl;
	cout<<" "<<endl;
	cout<<" "<<endl;
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
