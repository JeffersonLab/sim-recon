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
unsigned int EVENTS_TO_SKIP = 0;
unsigned int EVENTS_TO_KEEP = 1;
unsigned int SPECIFIC_EVENT_TO_KEEP = 0;
unsigned int SPECIFIC_OFFSET_TO_KEEP = 0;

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
	
	cout<<"Skipping "<<EVENTS_TO_SKIP<<endl;
	cout<<"Keeping  "<<EVENTS_TO_KEEP<<endl;
	
	// Output file
	cout<<" output file: "<<OUTFILENAME<<endl;
	s_iostream_t *fout = init_s_HDDM(OUTFILENAME);
	if(!fout){
		cout<<" Error opening output file \""<<OUTFILENAME<<"\"!"<<endl;
		exit(-1);
	}

	// Loop over input files
	unsigned int NEvents = 0;
	unsigned int NEvents_read = 0;
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
			
			bool write_this_event = false;
			
			// Loop over physics events within this event and see if one
			// has the run number of interest
			if(hddm_s->physicsEvents!=HDDM_NULL){
				for(unsigned int i=0; i<hddm_s->physicsEvents->mult; i++){
					int eventNo = hddm_s->physicsEvents->in[i].eventNo;
					if((unsigned int)eventNo == SPECIFIC_EVENT_TO_KEEP)write_this_event = true;
				}
			}
			
			// Check if we're in the range of offsets to write out
			if(NEvents_read>EVENTS_TO_SKIP)write_this_event = true;
			
			// Write this output event to file and free its memory
			if(write_this_event){
				flush_s_HDDM(hddm_s, fout);
				NEvents++;
			}else{
				flush_s_HDDM(hddm_s, NULL);
			}
		
			// Update ticker
			time_t now = time(NULL);
			if(now != last_time){
				cout<<"  "<<NEvents_read<<" events read     ("<<NEvents<<" event written) \r";cout.flush();
				last_time = now;
			}
			
			// Quit as soon as we wrote all of the events we're going to
			if(NEvents_read>=(EVENTS_TO_SKIP+EVENTS_TO_KEEP))break;

			if(QUIT)break;
		}

		// Close input file
		close_s_HDDM(fin);
	}
		
	// Close output file
	close_s_HDDM(fout);
	
	cout<<endl;
	cout<<" "<<NEvents_read<<" events read, "<<NEvents<<" events written"<<endl;

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
				case 'h': Usage();										break;
				case 'o': OUTFILENAME=&ptr[2];						break;
				case 's': EVENTS_TO_SKIP=atoi(&ptr[2]);			break;
				case 'k': EVENTS_TO_KEEP=atoi(&ptr[2]);			break;
				case 'e': SPECIFIC_OFFSET_TO_KEEP=atoi(&ptr[2]);	break;
				case 'E': SPECIFIC_EVENT_TO_KEEP=atoi(&ptr[2]);	break;
			}
		}else{
			INFILENAMES.push_back(argv[i]);
		}
	}

	if(INFILENAMES.size()==0){
		cout<<endl<<"You must enter a filename!"<<endl<<endl;
		Usage();
	}
	
	if(SPECIFIC_OFFSET_TO_KEEP>0){
		EVENTS_TO_KEEP=1;
		EVENTS_TO_SKIP=SPECIFIC_OFFSET_TO_KEEP-1;
	}

	if(SPECIFIC_EVENT_TO_KEEP>0){
		EVENTS_TO_KEEP=1;
		EVENTS_TO_SKIP=1000000000; // Large number of events to read it while looking for the specified event
	}
	
	if(OUTFILENAME==NULL){
		if(SPECIFIC_OFFSET_TO_KEEP>0){
			OUTFILENAME = new char[256];
			sprintf(OUTFILENAME,"evt%d.hddm", SPECIFIC_OFFSET_TO_KEEP);
		}else{
			OUTFILENAME = (char*)"culled.hddm";
		}
	}
}


//-----------
// Usage
//-----------
void Usage(void)
{
	cout<<endl<<"Usage:"<<endl;
	cout<<"     hddm_cull_events [-oOutputfile] [-sNeventsToSkip] [-kNeventsToKeep] file1.hddm file2.hddm ..."<<endl;
	cout<<endl;
	cout<<"options:"<<endl;
	cout<<"    -oOutputfile     Set output filename (def. culled.hddm)"<<endl;
	cout<<"    -sNeventsToSkip  Set number of events to skip (def. 0)"<<endl;
	cout<<"    -kNeventsToKeep  Set number of events to keep (def. 1)"<<endl;
	cout<<"    -eSingleEvent    Keep only the single, specified event (file pos.)"<<endl;
	cout<<"    -ESingleEvent    Keep only the single, specified event (event number)"<<endl;
	cout<<endl;
	cout<<" This will copy a continguous set of events from the combined event streams"<<endl;
	cout<<" into a seperate output file. The primary use for this would be to copy"<<endl;
	cout<<" a single, problematic event into a seperate file for easier debugging."<<endl;
	cout<<endl;
	cout<<" If the -eNNN option is used then only a single event is extracted"<<endl;
	cout<<" (the NNN-th event) and written to a file with the name evtNNN.hddm."<<endl;
	cout<<" Note that the event is counted from the begining of the file, starting"<<endl;
	cout<<" with \"1\". This does NOT look at the event number stored in the event itself."<<endl;
	cout<<" "<<endl;
	cout<<" If the -ENNN option is used then only a single event is extracted"<<endl;
	cout<<" (the specified event number) and written to a file with the name evtNNN.hddm."<<endl;
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
