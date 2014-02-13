// $Id$
//
// Created Oct 10, 2013  David Lawrence

#include "hddm_cull_events.h"

#include <HDDM/hddm_s.h>



//-----------
// Process_s  --  HDDM simulation format
//-----------
void Process_s(unsigned int &NEvents, unsigned int &NEvents_read)
{
	// Output file
	cout<<" output file: "<<OUTFILENAME<<endl;
	s_iostream_t *fout = init_s_HDDM(OUTFILENAME);
	if(!fout){
		cout<<" Error opening output file \""<<OUTFILENAME<<"\"!"<<endl;
		exit(-1);
	}

	// Loop over input files
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
			// has the event number of interest
			if(EVENT_TO_KEEP_MODE && hddm_s->physicsEvents!=HDDM_NULL){
				for(unsigned int i=0; i<hddm_s->physicsEvents->mult; i++){
					int eventNo = hddm_s->physicsEvents->in[i].eventNo;
					if((unsigned int)eventNo == SPECIFIC_EVENT_TO_KEEP){
						write_this_event = true;
						QUIT = true;
					}
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
}
