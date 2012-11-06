// $Id: mcsmear.cc 2388 2007-01-10 16:46:03Z davidl $
//
// Created June 22, 2005  David Lawrence

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

vector<int> Nevents_to_merge;
vector<bool> loop_source;
vector<char*> INFILENAMES;
char *OUTFILENAME = NULL;
int QUIT = 0;
int MAXEVENTS=-1;


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
	
	// Dummy check
	if(Nevents_to_merge.size() != INFILENAMES.size()){
		_DBG_<<"Size of Nevents_to_merge and INFILENAMES vectors not the same!"<<endl;
		_DBG_<<"This indicates a bug in the program. Exiting ..."<<endl;
		exit(-1);
	}
	if(Nevents_to_merge.size() != loop_source.size()){
		_DBG_<<"Size of Nevents_to_merge and loop_source vectors not the same!"<<endl;
		_DBG_<<"This indicates a bug in the program. Exiting ..."<<endl;
		exit(-1);
	}
	
	// Open Input file(s)
	vector<s_iostream_t*> instreams;
	for(unsigned int i=0; i<INFILENAMES.size(); i++){
		cout<<" input file: "<<INFILENAMES[i]<<"  ("<<Nevents_to_merge[i]<<" events)"<<endl;
		s_iostream_t *fin = open_s_HDDM(INFILENAMES[i]);
		if(!fin){
			cout<<" Error opening input file \""<<INFILENAMES[i]<<"\"!"<<endl;
			exit(-1);
		}
		
		instreams.push_back(fin);
	}
	
	// Output file
	cout<<" output file: "<<OUTFILENAME<<endl;
	s_iostream_t *fout = init_s_HDDM(OUTFILENAME);
	if(!fout){
		cout<<" Error opening output file \""<<OUTFILENAME<<"\"!"<<endl;
		exit(-1);
	}
	
	// Loop over output events until one of the inputs does not
	// have enough events to make an output
	int NEvents = 0;
	int NEvents_read = 0;
	time_t last_time = time(NULL);
	while(true){
		// First, read in all events we want to combine
		vector<s_HDDM_t*> events;
		
		// Loop over inputs
		bool done=false;
		events.clear();
		for(unsigned int i=0; i<instreams.size(); i++){
			// Loop over events to merge for this input
			for(unsigned int j=0; j<(unsigned int)Nevents_to_merge[i]; j++){
				s_HDDM_t *hddm_s = read_s_HDDM(instreams[i]);
				if(!hddm_s){
					// Looks like this source is out of events. Check if we need
					// to re-open this input to keep reading events from it.
					if(loop_source[i]){
						close_s_HDDM(instreams[i]);
						instreams[i] = open_s_HDDM(INFILENAMES[i]);
						cout<<endl<<"Reopened \""<<INFILENAMES[i]<<"\" ..."<<endl;
						hddm_s = read_s_HDDM(instreams[i]);
					}
				}
				if(!hddm_s){
					done = true;
					break;
				}
				events.push_back(hddm_s);
				NEvents_read++;
			}
			if(done)break;
		}
		if(done)break;
		
		// Make a new output event
		s_HDDM_t* output_hddm_s = make_s_HDDM();
		
		// An input file's event may have more than 1 PhysicsEvent already.
		// Loop through and find the total number of Physics events we
		// need to allocate for.
		int Nphysics_events = 0;
		for(unsigned int i=0; i<events.size(); i++)Nphysics_events += events[i]->physicsEvents->mult;
		
		// Add enough PhysicsEvents to hold all of our inputs
		output_hddm_s->physicsEvents = make_s_PhysicsEvents(Nphysics_events);
		
		// Copy PhysicsEvents for inputs, transferring ownership to
		// output event and freeing memory from the "heads" of the
		// input events.
		unsigned int &mult = output_hddm_s->physicsEvents->mult;
		mult = 0;
		for(unsigned int i=0; i<events.size(); i++){
			s_HDDM_t* input_hddm_s = events[i];
			
			// Loop over PhysicsEvents inside this input event
			for(unsigned int j=0; j<input_hddm_s->physicsEvents->mult; j++){
				output_hddm_s->physicsEvents->in[mult++] = input_hddm_s->physicsEvents->in[j];
				
				// Set the 2 pointers in PhysicsEvent to HDDM_NULL so they aren't freed
				// when the input event is freed.
				input_hddm_s->physicsEvents->in[j].reactions = (s_Reactions_t*)HDDM_NULL;
				input_hddm_s->physicsEvents->in[j].hitView = (s_HitView_t*)HDDM_NULL;
			}
			
			// Free this input event
			flush_s_HDDM(input_hddm_s, NULL);
		}
		
		// Write this output event to file and free its memory
		flush_s_HDDM(output_hddm_s, fout);
		NEvents++;
		
		// Update ticker
		time_t now = time(NULL);
		if(now != last_time){
			cout<<"  "<<NEvents_read<<" events read     ("<<NEvents<<" event written) \r";cout.flush();
			last_time = now;
		}

		if(QUIT)break;
		if(MAXEVENTS>0){
			if(NEvents>=MAXEVENTS)break;
		}
	}
	
	// Close all inputs
	for(unsigned int i=0; i<instreams.size(); i++)close_s_HDDM(instreams[i]);
	close_s_HDDM(fout);
	
	cout<<" "<<NEvents<<" events read"<<endl;

	return 0;
}

//-----------
// ParseCommandLineArguments
//-----------
void ParseCommandLineArguments(int narg, char* argv[])
{
	int Nmerge = 1;
	bool loop = false;

	INFILENAMES.clear();
	loop_source.clear();
	Nevents_to_merge.clear();

	for(int i=1; i<narg; i++){
		char *ptr = argv[i];
		
		if(ptr[0] == '-'){
			switch(ptr[1]){
				case 'h': Usage();						break;
				case 'N': Nmerge=atoi(&ptr[2]);		break;
				case 'M': MAXEVENTS=atoi(&ptr[2]);	break;
				case 'o': OUTFILENAME=&ptr[2];		break;
				case 'l': loop=true;						break;
				case 's': loop=false;					break;
			}
		}else{
			INFILENAMES.push_back(argv[i]);
			loop_source.push_back(loop);
			Nevents_to_merge.push_back(Nmerge);
		}
	}

	if(INFILENAMES.size()==0){
		cout<<endl<<"You must enter a filename!"<<endl<<endl;
		Usage();
	}
	
	if(OUTFILENAME==NULL){
		char *default_outfile = (char*)malloc(20);
		OUTFILENAME = strcpy(default_outfile,"merged.hddm");
	}
}


//-----------
// Usage
//-----------
void Usage(void)
{
	cout<<endl<<"Usage:"<<endl;
	cout<<"     hddm_merge_events [-Mmaxevents] [-oOutputfile] [-l|s -Nnum] file1.hddm [-l|s -Nnum] file2.hddm ..."<<endl;
	cout<<endl;
	cout<<"options:"<<endl;
	cout<<"    -oOutputfile  Set output filename (def. merged.hddm)"<<endl;
	cout<<"    -Mmaxevents   Do not write out more than maxevents events"<<endl;
	cout<<"    -Nnum         Merge together num events from each of the following input files"<<endl;
	cout<<"    -l            Continually loop over events from the following inputs"<<endl;
	cout<<"    -s            Stop when all events from the following inputs have been read"<<endl;
	cout<<endl;
	cout<<" This will take events from 1 or more HDDM files and merge them"<<endl;
	cout<<" into events written to an output HDDM file. This works by simply"<<endl;
	cout<<" copying the s_PhysicsEvent objects into a single event that gets"<<endl;
	cout<<" written to the output file."<<endl;
	cout<<" "<<endl;
	cout<<" Multiple input files can be specified and for each, the number of"<<endl;
	cout<<" events to merge. This allows one to do things such as merge 3 events"<<endl;
	cout<<" from a file of background events with a single event from a file of"<<endl;
	cout<<" signal events to get a sample of events with 3 times the backaground"<<endl;
	cout<<" rate overlayed on the signal."<<endl;
	cout<<" "<<endl;
	cout<<" By specifying the -l flag, certain input sources can be looped over"<<endl;
	cout<<" to recycle thier events. This is useful if one has a sample of"<<endl;
	cout<<" background events that is smaller than the sample of signal events."<<endl;
	cout<<" See the examples below."<<endl;
	cout<<" "<<endl;
	cout<<" NOTE: Events are not merged such that double hits are"<<endl;
	cout<<" merged. Hits are only copied so it is possible for"<<endl;
	cout<<" the output event to have two hits on the same wire"<<endl;
	cout<<" at the same time."<<endl;
	cout<<" "<<endl;
	cout<<" Example 1:"<<endl;
	cout<<"     hddm_merge_events -N2 file.hddm"<<endl;
	cout<<" "<<endl;
	cout<<"  This will combine every 2 events from file.hddm into a single in the"<<endl;
	cout<<" output. Since no output file name is specified, the file \"merged.hddm\""<<endl;
	cout<<" will be used."<<endl;
	cout<<" "<<endl;
	cout<<" "<<endl;
	cout<<" Example 2:"<<endl;
	cout<<"     hddm_merge_events signal.hddm -N10 em_bkgnd.hddm -N1 hadronic_bkgnd.hddm"<<endl;
	cout<<" "<<endl;
	cout<<" This will combine 1 event from signal.hddm, 10 events from em_bkgnd.hddm,"<<endl;
	cout<<" and 1 event from hadronic_bkgnd.hdm into a single output event."<<endl;
	cout<<" "<<endl;
	cout<<" "<<endl;
	cout<<" Example 3:"<<endl;
	cout<<"     hddm_merge_events -N3 file1.hddm file2.hddm file3.hddm -N1 file4.hddm"<<endl;
	cout<<" "<<endl;
	cout<<" This will combine 3 events from each of file1.hddm, file2.hddm, and"<<endl;
	cout<<" file3.hddm with 1 event from file4.hddm."<<endl;
	cout<<" "<<endl;
	cout<<" "<<endl;
	cout<<" Example 4:"<<endl;
	cout<<"     hddm_merge_events signal.hddm -l -N10 background.hddm"<<endl;
	cout<<" "<<endl;
	cout<<" This will combine 1 signal event with 10 background events to create"<<endl;
	cout<<" and output event. The events in the background file will be looped"<<endl;
	cout<<" over continuosly as needed until all of the signal events have been"<<endl;
	cout<<" processed."<<endl;
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
