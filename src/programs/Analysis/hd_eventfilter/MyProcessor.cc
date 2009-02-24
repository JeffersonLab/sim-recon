// Author: David Lawrence  June 25, 2004
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

//------------------------------------------------------------------
// init   -Open output file here (e.g. a ROOT file)
//------------------------------------------------------------------
jerror_t MyProcessor::init(void)
{
	// open HDDM file
	filename = "filtered.hddm";
	file = init_s_HDDM((char*)filename.c_str());
	Nevents_written = 0;

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
	// downcast to a DEventSourceHDDM structure. It a little strange setting
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
	
	
	// Initialize write_out flag. We set this to true to write the event
	// out and false to ignore it. Initialize it to false so that we can
	// pick the conditions needed to keep it.
	bool write_out=false;

	// Here we do whatever calculations are needed to determine if we keep
	// the event. Since this is basically a skeleton meant as an example
	// we do a trivial check on the momentum of the thrown particles.
	// In practice, one could request objects that require full reconstruction
	// as well so that filters could be built on those quantities as well.

	//---------------------- Filter Code Start ----------------------
	// Get data
	vector<const DMCThrown*> mcthrowns;
	loop->Get(mcthrowns);

	// Loop over thrown tracks
	for(unsigned int i=0;i<mcthrowns.size();i++){
		const DMCThrown *mcthrown = mcthrowns[i];
		
		// keep tracks with at least 1 thrown particle greater than 1GeV/c
		if(mcthrown->momentum().Mag()>1.0)write_out = true;

	}
	//----------------------- Filter Code End -----------------------
	
	// If write_out flag is set, write this event to our output file
	if(write_out){
		flush_s_HDDM(hddm, file);
		Nevents_written++;
	}
	
	// If we write the event out, then tell source not to free it
	hddm_source->flush_on_free = !write_out;

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
	cout<<" "<<Nevents_written<<" event written to "<<filename<<endl;
	
	return NOERROR;
}

