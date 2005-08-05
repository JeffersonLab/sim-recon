// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.cc
//

#include <iostream>
using namespace std;

#include "MyProcessor.h"
#include "DEvent.h"
#include "hddm_s.h"

#include "DMCTrackEfficiency.h"
#include "DEventSourceHDDM.h"
#include "DMCTrackEfficiency.h"
#include "DMCThrown.h"
#include "DMCReconstructed.h"

//------------------------------------------------------------------
// init   -Open output file here (e.g. a ROOT file)
//------------------------------------------------------------------
derror_t MyProcessor::init(void)
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
derror_t MyProcessor::evnt(DEventLoop *loop, int eventnumber)
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
	DEvent& event = loop->GetDEvent();
	DEventSource *source = event.GetDEventSource();
	DEventSourceHDDM *hddm_source = dynamic_cast<DEventSourceHDDM*>(source);
	if(!hddm_source){
		cerr<<" This program MUST be used with an HDDM file as input!"<<endl;
		exit(-1);
	}
	s_HDDM_t *hddm = (s_HDDM_t*)event.GetRef();
	if(!hddm)return NOERROR;
	
	// Get data
	vector<const DMCTrackEfficiency*> mctrackefficiencies;
	loop->Get(mctrackefficiencies);

	// Loop over thrown tracks
	bool write_out=false;
	for(unsigned int i=0;i<mctrackefficiencies.size();i++){
		const DMCTrackEfficiency *trkeff = mctrackefficiencies[i];
		
		if(trkeff->fittable){
			if(trkeff->Nhits_found){
				float fraction_from_thrown = (float)trkeff->Nhits_thrown_and_found/(float)trkeff->Nhits_found;
				if(fraction_from_thrown <0.70){
					write_out=true;
				}
			}else{
				write_out=true;
			}
		}
	}
	
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
derror_t MyProcessor::fini(void)
{
	if(file){
		close_s_HDDM(file);
		cout<<endl<<"Closed HDDM file"<<endl;
	}
	cout<<" "<<Nevents_written<<" event written to "<<filename<<endl;
	
	return NOERROR;
}

