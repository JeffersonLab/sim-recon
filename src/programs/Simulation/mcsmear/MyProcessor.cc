// Author: David Lawrence  Sat Jan 29 09:37:37 EST 2011
//
//
// MyProcessor.cc
//

#include <iostream>
using namespace std;

#include "MyProcessor.h"

#include <JANA/JEvent.h>

#include <HDDM/DEventSourceHDDM.h>
#include <TRACKING/DMCThrown.h>

extern void Smear(s_HDDM_t *hddm_s);
extern char *OUTFILENAME;

//------------------------------------------------------------------
// init   -Open output file 
//------------------------------------------------------------------
jerror_t MyProcessor::init(void)
{
	// open HDDM file
	file = init_s_HDDM(OUTFILENAME);
	if(!file){
		cout<<" Error opening output file \""<<OUTFILENAME<<"\"!"<<endl;
		exit(-1);
	}

	Nevents_written = 0;

	pthread_mutex_init(&output_file_mutex, NULL);	

	return NOERROR;
}

//------------------------------------------------------------------
// evnt   -Fill histograms here
//------------------------------------------------------------------
jerror_t MyProcessor::evnt(JEventLoop *loop, int eventnumber)
{
	// This is a little complicated. We need to get a hold of the s_HDDM_t
	// structure pointer for this event so we can pass it to flush_s_HDDM()
	// along with our ouput stream pointer. The flush routine frees up the
	// memory in the s_HDDM_t structure. When the framework tries "flush"ing
	// a second time, we get a seg fault. To prevent the framework from
	// flushing, we have to clear the free_on_flush flag (by default set
	// to true). This means we need to get the DEventSource pointer and
	// downcast to a DEventSourceHDDM structure. It's a little strange setting
	// this for every event, but we have no way of knowing when the event
	// source changes and this at least guarantees it for all event sources.
	JEvent& event = loop->GetJEvent();
	JEventSource *source = event.GetJEventSource();
	DEventSourceHDDM *hddm_source = dynamic_cast<DEventSourceHDDM*>(source);
	if(!hddm_source){
		cerr<<" This program MUST be used with an HDDM file as input!"<<endl;
		exit(-1);
	}
	s_HDDM_t *hddm = (s_HDDM_t*)event.GetRef();
	if(!hddm)return NOERROR;
	hddm_source->flush_on_free = false;
	
	// Smear values and add noise hits
	Smear(hddm);
	
	// Write event to output file
	pthread_mutex_lock(&output_file_mutex);
	flush_s_HDDM(hddm, file);
	Nevents_written++;
	pthread_mutex_unlock(&output_file_mutex);

	return NOERROR;
}

//------------------------------------------------------------------
// fini   -Close output file here
//------------------------------------------------------------------
jerror_t MyProcessor::fini(void)
{
	if(file){
		close_s_HDDM(file);
		cout<<endl<<"Closed HDDM file"<<endl;
	}
	cout<<" "<<Nevents_written<<" event written to "<<OUTFILENAME<<endl;
	
	return NOERROR;
}

