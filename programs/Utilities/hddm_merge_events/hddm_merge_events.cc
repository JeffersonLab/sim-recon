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

char *INFILENAME = NULL;
char *OUTFILENAME = NULL;
int QUIT = 0;
int Nevents_to_merge=2;

s_HDDM_t *merged_event=NULL;

void AddEvent(s_HDDM_t *new_event);
void AddCDCHits(s_HDDM_t *new_event);
void AddFDCHits(s_HDDM_t *new_event);
void AddThrowns(s_HDDM_t *new_event);


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
	int NEvents = 0;
	time_t last_time = time(NULL);
	while((hddm_s = read_s_HDDM(fin))){
		NEvents++;
		time_t now = time(NULL);
		if(now != last_time){
			cout<<"  "<<NEvents<<" events processed      \r";cout.flush();
			last_time = now;
		}
		
		// Accumulate events
		if(NEvents%Nevents_to_merge == 1){
			merged_event = hddm_s;
		}else{
			AddEvent(hddm_s);
			flush_s_HDDM(hddm_s, NULL);
		}
		
		// Every N events, write them out
		if(NEvents%Nevents_to_merge == 0){
			if(merged_event)flush_s_HDDM(merged_event, fout);
			merged_event=NULL;
		}		
		
		if(QUIT)break;
	}
	
	// close input and output files
	close_s_HDDM(fin);
	close_s_HDDM(fout);

	cout<<" "<<NEvents<<" events read"<<endl;

	return 0;
}

//-----------
// AddEvent
//-----------
void AddEvent(s_HDDM_t *new_event)
{
	AddCDCHits(new_event);
	AddFDCHits(new_event);
	AddThrowns(new_event);
}

//-----------
// AddCDCHits
//-----------
void AddCDCHits(s_HDDM_t *new_event)
{
	// Make sure structures exist in merged_event
	if(merged_event->physicsEvents == HDDM_NULL)merged_event->physicsEvents = make_s_PhysicsEvents(1);
	s_HitView_t* &old_hits = merged_event->physicsEvents->in[0].hitView;
	s_HitView_t* &new_hits = new_event->physicsEvents->in[0].hitView;
	if(old_hits == HDDM_NULL)old_hits = make_s_HitView();
	if(old_hits->centralDC==HDDM_NULL)old_hits->centralDC = make_s_CentralDC();
	if(old_hits->centralDC->cdcStraws==HDDM_NULL)old_hits->centralDC->cdcStraws = make_s_CdcStraws(0);
	if(old_hits->centralDC->cdcTruthPoints==HDDM_NULL)old_hits->centralDC->cdcTruthPoints = make_s_CdcTruthPoints(0);

	if(new_hits->centralDC==HDDM_NULL)return;
	if(new_hits->centralDC->cdcStraws==HDDM_NULL)return;
	if(new_hits->centralDC->cdcTruthPoints==HDDM_NULL)return;
	
	s_CdcStraws_t* &old_cdcStraws = old_hits->centralDC->cdcStraws;
	s_CdcStraws_t* &new_cdcStraws = new_hits->centralDC->cdcStraws;
	s_CdcTruthPoints_t* &old_cdcTruthPoints = old_hits->centralDC->cdcTruthPoints;
	s_CdcTruthPoints_t* &new_cdcTruthPoints = new_hits->centralDC->cdcTruthPoints;

	int Nold = old_cdcStraws->mult;
	int Nnew = new_cdcStraws->mult;
	if(Nnew>0){
		s_CdcStraws_t *cdcStraws = make_s_CdcStraws(Nold+Nnew);

		for(int i=0; i<Nold; i++){
			cdcStraws->in[cdcStraws->mult++] = old_cdcStraws->in[i];
		}
		for(int i=0; i<Nnew; i++){
			cdcStraws->in[cdcStraws->mult++] = new_cdcStraws->in[i];
		}
		
		// Here we need to replace the new_fdcChambers structure
		// so when new_event is freed, the underlying structures
		// are not.
		free(new_cdcStraws);
		new_cdcStraws = (s_CdcStraws_t*)HDDM_NULL;
		
		free(old_cdcStraws);
		old_cdcStraws = cdcStraws;
	}

	Nold = old_cdcTruthPoints->mult;
	Nnew = new_cdcTruthPoints->mult;
	if(Nnew>0){
		s_CdcTruthPoints_t *cdcTruthPoints = make_s_CdcTruthPoints(Nold+Nnew);

		for(int i=0; i<Nold; i++){
			cdcTruthPoints->in[cdcTruthPoints->mult++] = old_cdcTruthPoints->in[i];
		}
		for(int i=0; i<Nnew; i++){
			cdcTruthPoints->in[cdcTruthPoints->mult++] = new_cdcTruthPoints->in[i];
		}
		
		// Here we need to replace the new_fdcChambers structure
		// so when new_event is freed, the underlying structures
		// are not.
		free(new_cdcTruthPoints);
		new_cdcTruthPoints = (s_CdcTruthPoints_t*)HDDM_NULL;
		
		free(old_cdcTruthPoints);
		old_cdcTruthPoints = cdcTruthPoints;
	}

}

//-----------
// AddFDCHits
//-----------
void AddFDCHits(s_HDDM_t *new_event)
{
	// Make sure structures exist in merged_event
	if(merged_event->physicsEvents == HDDM_NULL)merged_event->physicsEvents = make_s_PhysicsEvents(1);
	s_HitView_t* &old_hits = merged_event->physicsEvents->in[0].hitView;
	s_HitView_t* &new_hits = new_event->physicsEvents->in[0].hitView;
	if(old_hits == HDDM_NULL)old_hits = make_s_HitView();
	if(old_hits->forwardDC==HDDM_NULL)old_hits->forwardDC = make_s_ForwardDC();
	if(old_hits->forwardDC->fdcChambers==HDDM_NULL)old_hits->forwardDC->fdcChambers = make_s_FdcChambers(0);
	
	if(new_hits->forwardDC==HDDM_NULL)return;
	if(new_hits->forwardDC->fdcChambers==HDDM_NULL)return;
	
	s_FdcChambers_t* &old_fdcChambers = old_hits->forwardDC->fdcChambers;
	s_FdcChambers_t* &new_fdcChambers = new_hits->forwardDC->fdcChambers;
	int Nold = old_fdcChambers->mult;
	int Nnew = new_fdcChambers->mult;
	if(Nnew>0){
		s_FdcChambers_t *fdcChambers = make_s_FdcChambers(Nold+Nnew);

		for(int i=0; i<Nold; i++){
			fdcChambers->in[fdcChambers->mult++] = old_fdcChambers->in[i];
		}
		for(int i=0; i<Nnew; i++){
			fdcChambers->in[fdcChambers->mult++] = new_fdcChambers->in[i];
		}
		
		// Here we need to replace the new_fdcChambers structure
		// so when new_event is freed, the underlying structures
		// are not.
		free(new_fdcChambers);
		new_fdcChambers = (s_FdcChambers_t*)HDDM_NULL;
		
		free(old_fdcChambers);
		old_fdcChambers = fdcChambers;
	}
}


//-----------
// AddThrowns
//-----------
void AddThrowns(s_HDDM_t *new_event)
{
	// Make sure structures exist in merged_event
	if(merged_event->physicsEvents == HDDM_NULL)merged_event->physicsEvents = make_s_PhysicsEvents(1);
	s_Reactions_t* &old_reactions = merged_event->physicsEvents->in[0].reactions;
	s_Reactions_t* &new_reactions = new_event->physicsEvents->in[0].reactions;
	if(old_reactions == HDDM_NULL)old_reactions = make_s_Reactions(0);
	
	if(new_reactions==HDDM_NULL)return;
	
	int Nold = old_reactions->mult;
	int Nnew = new_reactions->mult;
	if(Nnew>0){
		s_Reactions_t *reactions = make_s_Reactions(Nold+Nnew);

		for(int i=0; i<Nold; i++){
			reactions->in[reactions->mult++] = old_reactions->in[i];
		}
		for(int i=0; i<Nnew; i++){
			reactions->in[reactions->mult++] = new_reactions->in[i];
		}
		
		// Here we need to replace the new_reactions structure
		// so when new_event is freed, the underlying structures
		// are not.
		free(new_reactions);
		new_reactions = (s_Reactions_t*)HDDM_NULL;
		
		free(old_reactions);
		old_reactions = reactions;
	}
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
				case 'N': Nevents_to_merge=atoi(&ptr[2]);						break;
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
	sprintf(str, "%s_merged.hddm", path_stripped);
	OUTFILENAME = strdup(str);
}


//-----------
// Usage
//-----------
void Usage(void)
{
	cout<<endl<<"Usage:"<<endl;
	cout<<"     hddm_merge_events [-Nnum] file.hddm"<<endl;
	cout<<endl;
	cout<<"options:"<<endl;
	cout<<"    -Nnum    Merge together num events"<<endl;
	cout<<endl;
	cout<<" This will merge hits from multiple events into one"<<endl;
	cout<<" event and write the merged events to an output file."<<endl;
	cout<<" At this time, only FDC, CDC, and thrown particle "<<endl;
	cout<<" information is copied."<<endl;
	cout<<" NOTE: Events are not merged such that double hits are"<<endl;
	cout<<" merged. Hits are only copied so it is possible for"<<endl;
	cout<<" the output event to have two hits on the same wire"<<endl;
	cout<<" at the same time."<<endl;
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
