// $Id$
//
//    File: JEventProcessor_DAQTree.cc
// Created: Tue Oct 22 14:55:40 EDT 2013
// Creator: dalton (on Linux gluon45.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//

#include <iostream>
#include <algorithm>
#include "JEventProcessor_DAQTree.h"
#include <JANA/JFactory.h>
using namespace jana;

#include <stdint.h>
#include <DAQ/Df125WindowRawData.h>
#include <DAQ/DF1TDCHit.h>
#include <DAQ/DF1TDCTriggerTime.h>
#include <DAQ/Df250WindowRawData.h>
#include <DAQ/Df250PulseRawData.h>
#include <DAQ/Df250PulseIntegral.h>
#include <DAQ/Df250PulseTime.h>
#include <DAQ/Df250PulsePedestal.h>
#include <DAQ/Df250TriggerTime.h>

/// Define a comparison operator for sorting objets of all the DAQ types.
bool Df125WindowRawData_cmp(const Df125WindowRawData *a,const Df125WindowRawData *b){
	// sort by crate, then by slot, then by channel, then by trigger number
	if (a->rocid   != b->rocid)   return a->rocid < b->rocid;
	if (a->slot    != b->slot )   return a->slot < b->slot;
	if (a->channel != b->channel) return a->channel < b->channel;
	return a->itrigger < b->itrigger;
}
bool DF1TDCHit_cmp(const DF1TDCHit *a,const DF1TDCHit *b){
	// sort by crate, then by slot, then by channel, then by trigger number
	if (a->rocid   != b->rocid)   return a->rocid < b->rocid;
	if (a->slot    != b->slot )   return a->slot < b->slot;
	if (a->channel != b->channel) return a->channel < b->channel;
	return a->itrigger < b->itrigger;
}
bool DF1TDCTriggerTime_cmp(const DF1TDCTriggerTime *a,const DF1TDCTriggerTime *b){
	// sort by crate, then by slot, then by channel, then by trigger number
	if (a->rocid   != b->rocid)   return a->rocid < b->rocid;
	if (a->slot    != b->slot )   return a->slot < b->slot;
	return a->itrigger < b->itrigger;
}
bool Df250WindowRawData_cmp(const Df250WindowRawData *a,const Df250WindowRawData *b){
	// sort by crate, then by slot, then by channel, then by trigger number
	if (a->rocid   != b->rocid)   return a->rocid < b->rocid;
	if (a->slot    != b->slot )   return a->slot < b->slot;
	if (a->channel != b->channel) return a->channel < b->channel;
	return a->itrigger < b->itrigger;
}
bool Df250PulseRawData_cmp(const Df250PulseRawData *a,const Df250PulseRawData *b){
	// sort by crate, then by slot, then by channel, then by trigger number
	if (a->rocid   != b->rocid)   return a->rocid < b->rocid;
	if (a->slot    != b->slot )   return a->slot < b->slot;
	if (a->channel != b->channel) return a->channel < b->channel;
	return a->itrigger < b->itrigger;
}
bool Df250PulseIntegral_cmp(const Df250PulseIntegral *a,const Df250PulseIntegral *b){
	// sort by crate, then by slot, then by channel, then by trigger number
	if (a->rocid   != b->rocid)   return a->rocid < b->rocid;
	if (a->slot    != b->slot )   return a->slot < b->slot;
	if (a->channel != b->channel) return a->channel < b->channel;
	return a->itrigger < b->itrigger;
}
bool Df250PulseTime_cmp(const Df250PulseTime *a,const Df250PulseTime *b){
	// sort by crate, then by slot, then by channel, then by trigger number
	if (a->rocid   != b->rocid)   return a->rocid < b->rocid;
	if (a->slot    != b->slot )   return a->slot < b->slot;
	if (a->channel != b->channel) return a->channel < b->channel;
	return a->itrigger < b->itrigger;
}
bool Df250PulsePedestal_cmp(const Df250PulsePedestal *a,const Df250PulsePedestal *b){
	// sort by crate, then by slot, then by channel, then by trigger number
	if (a->rocid   != b->rocid)   return a->rocid < b->rocid;
	if (a->slot    != b->slot )   return a->slot < b->slot;
	if (a->channel != b->channel) return a->channel < b->channel;
	return a->itrigger < b->itrigger;
}
bool Df250TriggerTime_cmp(const Df250TriggerTime *a,const Df250TriggerTime *b){
	// sort by crate, then by slot, then by channel, then by trigger number
	if (a->rocid   != b->rocid)   return a->rocid < b->rocid;
	if (a->slot    != b->slot )   return a->slot < b->slot;
	return a->itrigger < b->itrigger;
}


// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
extern "C"{
	void InitPlugin(JApplication *app){
		InitJANAPlugin(app);
		app->AddProcessor(new JEventProcessor_DAQTree());
	}
} // "C"


//------------------
// JEventProcessor_DAQTree (Constructor)
//------------------
JEventProcessor_DAQTree::JEventProcessor_DAQTree()
{

}

//------------------
// ~JEventProcessor_DAQTree (Destructor)
//------------------
JEventProcessor_DAQTree::~JEventProcessor_DAQTree()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_DAQTree::init(void)
{
	// This is called once at program startup. If you are creating
	// and filling historgrams in this plugin, you should lock the
	// ROOT mutex like this:
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();
	//

	printf("JEventProcessor_DAQTree::init()\n");

	/// Initialize the flags
	f125WRDtree_exists = 0;
	F1TDCHtree_exists = 0;
	F1TDCTTtree_exists = 0;
	f250WRDtree_exists = 0;
	f250PRDtree_exists = 0;
	f250PItree_exists = 0;
	f250PTtree_exists = 0;
	f250PPtree_exists = 0;
	f250TTtree_exists = 0;

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_DAQTree::brun(JEventLoop *eventLoop, int runnumber)
{
	// This is called whenever the run number changes
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_DAQTree::evnt(JEventLoop *loop, int eventnumber)
{
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// loop->Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.
	// Here's an example:
	//

	/// Trees are filled with data
	japp->RootWriteLock();


	/// Df125WindowRawData
	const uint32_t numDf125WRDpedsamps = 10;
	const Int_t Df125WRDminpeakheight = 100;
	/// Get a vector of Df125WindowRawData objects for this event (1 object for each crate/slot/channel above threshold)
	vector<const Df125WindowRawData*> f125WindowRawData_vec;
	loop->Get(f125WindowRawData_vec);
	sort(f125WindowRawData_vec.begin(), f125WindowRawData_vec.end(), Df125WindowRawData_cmp);
	unsigned int num_f125WRD = f125WindowRawData_vec.size();
	/// Create tree if doesn't exist and objects found.
	if (!f125WRDtree_exists && num_f125WRD>0) {
		printf("DAQTree >>eventnum %i, found %4i Df125WindowRawData objects\n",eventnumber,num_f125WRD);
		printf("DAQTree >>Creating tree Df125WindowRawData_tree\n");
		Df125WindowRawData_tree = new TTree("Df125WindowRawData",
											"tree of flash 125 raw window data (waveform samples) for each channel and event");
		Df125WindowRawData_tree->Branch("channelnum",&channelnum,"channelnum/i");
		Df125WindowRawData_tree->Branch("eventnum",&eventnum,"eventnum/i");
		Df125WindowRawData_tree->Branch("rocid",&rocid,"rocid/i");
		Df125WindowRawData_tree->Branch("slot",&slot,"slot/i");
		Df125WindowRawData_tree->Branch("channel",&channel,"channel/i");
		Df125WindowRawData_tree->Branch("itrigger",&itrigger,"itrigger/i");
		Df125WindowRawData_tree->Branch("waveform",&waveform);
		Df125WindowRawData_tree->Branch("nsamples",&nsamples,"nsamples/i");
		Df125WindowRawData_tree->Branch("w_integral",&w_integral,"w_integral/i");
		Df125WindowRawData_tree->Branch("w_min",&w_min,"w_min/i");
		Df125WindowRawData_tree->Branch("w_max",&w_max,"w_max/i");
		Df125WindowRawData_tree->Branch("w_samp1",&w_samp1,"w_samp1/i");
		Df125WindowRawData_tree->Branch("w_ped",&w_ped,"w_ped/i");
		Df125WindowRawData_tree->Branch("w_time",&w_time,"w_time/f");
		Df125WindowRawData_tree->Branch("invalid_samples",&invalid_samples,"invalid_samples/b");
		Df125WindowRawData_tree->Branch("overflow",&overflow,"overflow/b");
		f125WRDtree_exists = 1;
	}
	eventnum = eventnumber;
	// Loop over all objects in this event
	for(unsigned int c_chan=0; c_chan<num_f125WRD; c_chan++){
		waveform.clear();
		channelnum = c_chan;
		const Df125WindowRawData *f125WindowRawData = f125WindowRawData_vec[c_chan];
		rocid = f125WindowRawData->rocid;
		slot = f125WindowRawData->slot;
		channel = f125WindowRawData->channel;
		itrigger = f125WindowRawData->itrigger;
		invalid_samples = f125WindowRawData->invalid_samples;
		overflow = f125WindowRawData->overflow;

		// Get a vector of the samples for this channel
		const vector<uint16_t> &samplesvector = f125WindowRawData->samples;
		nsamples=samplesvector.size();
		/// loop over the samples to calculate integral, min, max
		for (uint16_t c_samp=0; c_samp<nsamples; c_samp++) {
			if (samplesvector[c_samp]<4096) {
				waveform.push_back(samplesvector[c_samp]); // push the sample into the waveform vector
			} else {
				waveform.push_back(4096); // push the sample into the waveform vector
			}
			if (c_samp==0) {  // use first sample for initialization
				w_integral = samplesvector[0]; 
				w_min = samplesvector[0];
				w_max = samplesvector[0];
				w_samp1 = samplesvector[0];
				w_ped = samplesvector[0]; 
			} else {
				if (c_samp<numDf125WRDpedsamps) {
					w_ped += samplesvector[c_samp];
				}
				w_integral += samplesvector[c_samp];
				if (w_min > samplesvector[c_samp]) w_min = samplesvector[c_samp];
				if (w_max < samplesvector[c_samp]) w_max = samplesvector[c_samp];
			}		
		}
		/// find the time to cross half peak height
		Int_t lastbelowsamp=0, peakheight = w_max-w_min;
		Float_t threshold = w_min + peakheight/2.0;
		Float_t  firstaboveheight=0, lastbelowheight=0;
		w_time=0;
		if (peakheight > Df125WRDminpeakheight) { 
			for (uint16_t c_samp=0; c_samp<nsamples; c_samp++) {
				if (samplesvector[c_samp]>threshold) { 
					firstaboveheight = samplesvector[c_samp];
					lastbelowsamp = c_samp-1;
					lastbelowheight = samplesvector[c_samp-1];
					break;
				}
			}
			w_time = lastbelowsamp + (threshold-lastbelowheight)/(firstaboveheight-lastbelowheight);
		}
		Df125WindowRawData_tree->Fill();
	}


	/// DF1TDCHit
	/// Get a vector of DF1TDCHit objects for this event (1 object for each crate/slot/channel above threshold)
	vector<const DF1TDCHit*> F1TDCHit_vec;
	loop->Get(F1TDCHit_vec);
	sort(F1TDCHit_vec.begin(), F1TDCHit_vec.end(), DF1TDCHit_cmp);
	unsigned int num_F1TDCH = F1TDCHit_vec.size();
	/// Create tree if doesn't exist and objects found.
	if (!F1TDCHtree_exists && num_F1TDCH>0) {
		printf("DAQTree >>eventnum %i, found %4i DF1TDCHit objects\n",eventnumber,num_F1TDCH);
		printf("DAQTree >>Creating tree F1TDCHit_tree\n");
		DF1TDCHit_tree = new TTree("DF1TDCHit",	"tree of F1 TDC hit times for each channel and event");
		DF1TDCHit_tree->Branch("channelnum",&channelnum,"channelnum/i");
		DF1TDCHit_tree->Branch("eventnum",&eventnum,"eventnum/i");
		DF1TDCHit_tree->Branch("rocid",&rocid,"rocid/i");
		DF1TDCHit_tree->Branch("slot",&slot,"slot/i");
		DF1TDCHit_tree->Branch("channel",&channel,"channel/i");
		DF1TDCHit_tree->Branch("itrigger",&itrigger,"itrigger/i");
		DF1TDCHit_tree->Branch("trig_time",&trig_time,"trig_time/i");
		DF1TDCHit_tree->Branch("time",&time,"time/I");
		F1TDCHtree_exists = 1;
	}
	eventnum = eventnumber;
	/// Loop over all DF1TDCHit objects in this event
	for(unsigned int c_chan=0; c_chan<num_F1TDCH; c_chan++){
		channelnum = c_chan;
		const DF1TDCHit *F1TDCHit = F1TDCHit_vec[c_chan];
		rocid = F1TDCHit->rocid;
		slot = F1TDCHit->slot;
		channel = F1TDCHit->channel;
		itrigger = F1TDCHit->itrigger;
		trig_time= F1TDCHit->trig_time;
		time = F1TDCHit->time;
		DF1TDCHit_tree->Fill();
	}


	/// DF1TDCTriggerTime
	/// Get a vector of DF1TDCTriggerTime objects for this event (1 object for each crate/slot)
	vector<const DF1TDCTriggerTime*> F1TDCTriggerTime_vec;
	loop->Get(F1TDCTriggerTime_vec);
	sort(F1TDCTriggerTime_vec.begin(), F1TDCTriggerTime_vec.end(), DF1TDCTriggerTime_cmp);
	unsigned int num_F1TDCTT = F1TDCTriggerTime_vec.size();
	/// Create tree if doesn't exist and objects found.
	if (!F1TDCTTtree_exists && num_F1TDCTT>0) {
		printf("DAQTree >>eventnum %i, found %4i DF1TDCTriggerTime objects\n",eventnumber,num_F1TDCTT);
		printf("DAQTree >>Creating tree DF1TDCTriggerTime_tree\n");
		DF1TDCTriggerTime_tree = new TTree("DF1TDCTriggerTime", "tree of F1 TDC trigger times for each slot and event");
		DF1TDCTriggerTime_tree->Branch("eventnum",&eventnum,"eventnum/i");
		DF1TDCTriggerTime_tree->Branch("rocid",&rocid,"rocid/i");
		DF1TDCTriggerTime_tree->Branch("slot",&slot,"slot/i");
		DF1TDCTriggerTime_tree->Branch("itrigger",&itrigger,"itrigger/i");
		DF1TDCTriggerTime_tree->Branch("time",&time,"time/I");
		F1TDCTTtree_exists = 1;
	}
	eventnum = eventnumber;
	/// Loop over all DF1TDCTriggerTime objects in this event
	for(unsigned int c_chan=0; c_chan<num_F1TDCTT; c_chan++){
		channelnum = c_chan;
		const DF1TDCTriggerTime *F1TDCTriggerTime = F1TDCTriggerTime_vec[c_chan];
		rocid = F1TDCTriggerTime->rocid;
		slot = F1TDCTriggerTime->slot;
		itrigger = F1TDCTriggerTime->itrigger;
		time = F1TDCTriggerTime->time;
		DF1TDCTriggerTime_tree->Fill();
	}



	/// Df250WindowRawData
	const uint32_t numDf250WRDpedsamps= 10;
	const Int_t Df250WRDminpeakheight = 11;
	/// Get a vector of Df250WindowRawData objects for this event (1 object for each crate/slot/channel above threshold)
	vector<const Df250WindowRawData*> f250WindowRawData_vec;
	loop->Get(f250WindowRawData_vec);
	sort(f250WindowRawData_vec.begin(), f250WindowRawData_vec.end(), Df250WindowRawData_cmp);
	unsigned int num_f250WRD = f250WindowRawData_vec.size();
	/// Create tree if doesn't exist and objects found.
	if (!f250WRDtree_exists && num_f250WRD>0) {
		printf("DAQTree >>eventnum %i, found %4i Df250WindowRawData objects\n",eventnumber,num_f250WRD);
		printf("DAQTree >>Creating tree Df250WindowRawData_tree\n");
		Df250WindowRawData_tree = new TTree("Df250WindowRawData",
											"tree of flash 250 raw window data (waveform samples) for each channel and event");
		Df250WindowRawData_tree->Branch("channelnum",&channelnum,"channelnum/i");
		Df250WindowRawData_tree->Branch("eventnum",&eventnum,"eventnum/i");
		Df250WindowRawData_tree->Branch("rocid",&rocid,"rocid/i");
		Df250WindowRawData_tree->Branch("slot",&slot,"slot/i");
		Df250WindowRawData_tree->Branch("channel",&channel,"channel/i");
		Df250WindowRawData_tree->Branch("itrigger",&itrigger,"itrigger/i");
		Df250WindowRawData_tree->Branch("waveform",&waveform);
		Df250WindowRawData_tree->Branch("nsamples",&nsamples,"nsamples/i");
		Df250WindowRawData_tree->Branch("w_integral",&w_integral,"w_integral/i");
		Df250WindowRawData_tree->Branch("w_min",&w_min,"w_min/i");
		Df250WindowRawData_tree->Branch("w_max",&w_max,"w_max/i");
		Df250WindowRawData_tree->Branch("w_samp1",&w_samp1,"w_samp1/i");
		Df250WindowRawData_tree->Branch("w_ped",&w_ped,"w_ped/i");
		Df250WindowRawData_tree->Branch("w_time",&w_time,"w_time/f");
		Df250WindowRawData_tree->Branch("invalid_samples",&invalid_samples,"invalid_samples/b");
		Df250WindowRawData_tree->Branch("overflow",&overflow,"overflow/b");
		f250WRDtree_exists = 1;
	}
	eventnum = eventnumber;
	// Loop over all  objects in this event
	for(unsigned int c_chan=0; c_chan<num_f250WRD; c_chan++){
		waveform.clear();
		channelnum = c_chan;
		const Df250WindowRawData *f250WindowRawData = f250WindowRawData_vec[c_chan];
		rocid = f250WindowRawData->rocid;
		slot = f250WindowRawData->slot;
		channel = f250WindowRawData->channel;
		itrigger = f250WindowRawData->itrigger;
		invalid_samples = f250WindowRawData->invalid_samples;
		overflow = f250WindowRawData->overflow;
		// Get a vector of the samples for this channel
		const vector<uint16_t> &samplesvector = f250WindowRawData->samples;
		nsamples=samplesvector.size();
		/// loop over the samples to calculate integral, min, max
		for (uint16_t c_samp=0; c_samp<nsamples; c_samp++) {
			if (samplesvector[c_samp]<4096) {
				waveform.push_back(samplesvector[c_samp]); // push the sample into the waveform vector
			} else {
				waveform.push_back(4096); // push the sample into the waveform vector
			}
			if (c_samp==0) {  // use first sample for initialization
				w_integral = samplesvector[0]; 
				w_min = samplesvector[0];
				w_max = samplesvector[0];
				w_samp1 = samplesvector[0];
				w_ped = samplesvector[0]; 
			} else {
				if (c_samp<numDf250WRDpedsamps) {
					w_ped += samplesvector[c_samp];
				}
				w_integral += samplesvector[c_samp];
				if (w_min > samplesvector[c_samp]) w_min = samplesvector[c_samp];
				if (w_max < samplesvector[c_samp]) w_max = samplesvector[c_samp];
			}		
		}
		/// find the time to cross half peak height
		Int_t lastbelowsamp=0, peakheight = w_max-w_min;
		Float_t threshold = w_min + peakheight/2.0;
		Float_t  firstaboveheight=0, lastbelowheight=0;
		w_time=0;
		if (peakheight > Df250WRDminpeakheight) { 
			for (uint16_t c_samp=0; c_samp<nsamples; c_samp++) {
				if (samplesvector[c_samp]>threshold) { 
					firstaboveheight = samplesvector[c_samp];
					lastbelowsamp = c_samp-1;
					lastbelowheight = samplesvector[c_samp-1];
					break;
				}
			}
			w_time = lastbelowsamp + (threshold-lastbelowheight)/(firstaboveheight-lastbelowheight);
		}
		Df250WindowRawData_tree->Fill();
	}


	/// Df250PulseRawData
	const uint32_t numDf250PRDpedsamps= 4;
	const Int_t Df250PRDminpeakheight = 11;
	/// Get a vector of Df250PulseRawData objects for this event (1 object for each crate/slot/channel above threshold)
	vector<const Df250PulseRawData*> f250PulseRawData_vec;
	loop->Get(f250PulseRawData_vec);
	sort(f250PulseRawData_vec.begin(), f250PulseRawData_vec.end(), Df250PulseRawData_cmp);
	unsigned int num_f250PRD = f250PulseRawData_vec.size();
	/// Create tree if doesn't exist and objects found.
	if (!f250PRDtree_exists && num_f250PRD>0) {
		printf("DAQTree >>eventnum %i, found %4i Df250PulseRawData objects\n",eventnumber,num_f250PRD);
		printf("DAQTree >>Creating tree Df250PulseRawData_tree\n");
		Df250PulseRawData_tree = new TTree("Df250PulseRawData",
										   "tree of flash 250 pulse raw data for each channel and event");
		Df250PulseRawData_tree->Branch("channelnum",&channelnum,"channelnum/i");
		Df250PulseRawData_tree->Branch("eventnum",&eventnum,"eventnum/i");
		Df250PulseRawData_tree->Branch("rocid",&rocid,"rocid/i");
		Df250PulseRawData_tree->Branch("slot",&slot,"slot/i");
		Df250PulseRawData_tree->Branch("channel",&channel,"channel/i");
		Df250PulseRawData_tree->Branch("itrigger",&itrigger,"itrigger/i");
		Df250PulseRawData_tree->Branch("pulse_number",&pulse_number,"pulse_number/i");
		Df250PulseRawData_tree->Branch("first_sample_number",&first_sample_number,"first_sample_number/i");
		Df250PulseRawData_tree->Branch("waveform",&waveform);
		Df250PulseRawData_tree->Branch("nsamples",&nsamples,"nsamples/i");
		Df250PulseRawData_tree->Branch("w_integral",&w_integral,"w_integral/i");
		Df250PulseRawData_tree->Branch("w_min",&w_min,"w_min/i");
		Df250PulseRawData_tree->Branch("w_max",&w_max,"w_max/i");
		Df250PulseRawData_tree->Branch("w_samp1",&w_samp1,"w_samp1/i");
		Df250PulseRawData_tree->Branch("w_ped",&w_ped,"w_ped/i");
		Df250PulseRawData_tree->Branch("w_time",&w_time,"w_time/f");
		Df250PulseRawData_tree->Branch("invalid_samples",&invalid_samples,"invalid_samples/b");
		Df250PulseRawData_tree->Branch("overflow",&overflow,"overflow/b");
		f250PRDtree_exists = 1;
	}
	eventnum = eventnumber;
	// Loop over all  objects in this event
	for(unsigned int c_chan=0; c_chan<num_f250PRD; c_chan++){
		waveform.clear();
		channelnum = c_chan;
		const Df250PulseRawData *f250PulseRawData = f250PulseRawData_vec[c_chan];
		rocid = f250PulseRawData->rocid;
		slot = f250PulseRawData->slot;
		channel = f250PulseRawData->channel;
		itrigger = f250PulseRawData->itrigger;
		pulse_number  = f250PulseRawData->pulse_number;
		first_sample_number = f250PulseRawData->first_sample_number;
		invalid_samples = f250PulseRawData->invalid_samples;
		overflow = f250PulseRawData->overflow;
		// Get a vector of the samples for this channel
		const vector<uint16_t> &samplesvector = f250PulseRawData->samples;
		nsamples=samplesvector.size();
		/// loop over the samples to calculate integral, min, max
		for (uint16_t c_samp=0; c_samp<nsamples; c_samp++) {
			/// Deal with overflow by setting sample to max val
			if (samplesvector[c_samp]<4096) {
				waveform.push_back(samplesvector[c_samp]); // push the sample into the waveform vector
			} else {
				waveform.push_back(4096); // push the sample into the waveform vector
			}
			if (c_samp==0) {  // use first sample for initialization
				w_integral = samplesvector[0]; 
				w_min = samplesvector[0];
				w_max = samplesvector[0];
				w_samp1 = samplesvector[0]; 
				w_ped = samplesvector[0]; 
			} else {
				if (c_samp<numDf250PRDpedsamps) {
					w_ped += samplesvector[c_samp];
				}
				w_integral += samplesvector[c_samp];
				if (w_min > samplesvector[c_samp]) w_min = samplesvector[c_samp];
				if (w_max < samplesvector[c_samp]) w_max = samplesvector[c_samp];
			}		
		}
		/// find the time to cross half peak height
		Int_t lastbelowsamp=0, peakheight = w_max-w_min;
		Float_t threshold = w_min + peakheight/2.0;
		Float_t  firstaboveheight=0, lastbelowheight=0;
		w_time=0;
		if (peakheight > Df250PRDminpeakheight) { 
			for (uint16_t c_samp=0; c_samp<nsamples; c_samp++) {
				if (samplesvector[c_samp]>threshold) { 
					firstaboveheight = samplesvector[c_samp];
					lastbelowsamp = c_samp-1;
					lastbelowheight = samplesvector[c_samp-1];
					break;
				}
			}
			w_time = lastbelowsamp + (threshold-lastbelowheight)/(firstaboveheight-lastbelowheight);
		}
		Df250PulseRawData_tree->Fill();
	}

	/// Df250PulseIntegral
	/// Get a vector of Df250PulseIntegral objects for this event (1 object for each crate/slot/channel above threshold)
	vector<const Df250PulseIntegral*> f250PulseIntegral_vec;
	loop->Get(f250PulseIntegral_vec);
	sort(f250PulseIntegral_vec.begin(), f250PulseIntegral_vec.end(), Df250PulseIntegral_cmp);
	unsigned int num_f250PI = f250PulseIntegral_vec.size();
	/// Create tree if doesn't exist and objects found.
	if (!f250PItree_exists && num_f250PI>0) {
		printf("DAQTree >>eventnum %i, found %4i Df250PulseIntegral objects\n",eventnumber,num_f250PI);
		printf("DAQTree >>Creating tree Df250PulseIntegral_tree\n");
		Df250PulseIntegral_tree = new TTree("Df250PulseIntegral",
											"tree of flash 250 pulse integral for each channel and event");
		Df250PulseIntegral_tree->Branch("channelnum",&channelnum,"channelnum/i");
		Df250PulseIntegral_tree->Branch("eventnum",&eventnum,"eventnum/i");
		Df250PulseIntegral_tree->Branch("rocid",&rocid,"rocid/i");
		Df250PulseIntegral_tree->Branch("slot",&slot,"slot/i");
		Df250PulseIntegral_tree->Branch("channel",&channel,"channel/i");
		Df250PulseIntegral_tree->Branch("itrigger",&itrigger,"itrigger/i");
		Df250PulseIntegral_tree->Branch("pulse_number",&pulse_number,"pulse_number/i");
		Df250PulseIntegral_tree->Branch("quality_factor",&quality_factor,"quality_factor/i");
		Df250PulseIntegral_tree->Branch("integral",&integral,"integral/I");
		Df250PulseIntegral_tree->Branch("pedestal",&pedestal,"pedestal/I");
		f250PItree_exists = 1;
	}
	eventnum = eventnumber;
	/// Loop over all Df250PulseIntegral objects in this event
	for(unsigned int c_chan=0; c_chan<num_f250PI; c_chan++){
		channelnum = c_chan;
		const Df250PulseIntegral *f250PulseIntegral = f250PulseIntegral_vec[c_chan];
		rocid = f250PulseIntegral->rocid;
		slot = f250PulseIntegral->slot;
		channel = f250PulseIntegral->channel;
		itrigger = f250PulseIntegral->itrigger;
		pulse_number = f250PulseIntegral->pulse_number;
		quality_factor= f250PulseIntegral->quality_factor;
		integral = f250PulseIntegral->integral;
		pedestal = f250PulseIntegral->pedestal;
		Df250PulseIntegral_tree->Fill();
	}


	/// Df250PulseTime
	/// Get a vector of Df250PulseTime objects for this event (1 object for each crate/slot/channel above threshold)
	vector<const Df250PulseTime*> f250PulseTime_vec;
	loop->Get(f250PulseTime_vec);
	sort(f250PulseTime_vec.begin(), f250PulseTime_vec.end(), Df250PulseTime_cmp);
	unsigned int num_f250PT = f250PulseTime_vec.size();
	/// Create tree if doesn't exist and objects found.
	if (!f250PTtree_exists && num_f250PT>0) {
		printf("DAQTree >>eventnum %i, found %4i Df250PulseTime objects\n",eventnumber,num_f250PT);
		printf("DAQTree >>Creating tree f250PulseTime_tree\n");
		Df250PulseTime_tree = new TTree("Df250PulseTime",
										"tree of flash 250 pulse times for each channel and event");
		Df250PulseTime_tree->Branch("channelnum",&channelnum,"channelnum/i");
		Df250PulseTime_tree->Branch("eventnum",&eventnum,"eventnum/i");
		Df250PulseTime_tree->Branch("rocid",&rocid,"rocid/i");
		Df250PulseTime_tree->Branch("slot",&slot,"slot/i");
		Df250PulseTime_tree->Branch("channel",&channel,"channel/i");
		Df250PulseTime_tree->Branch("itrigger",&itrigger,"itrigger/i");
		Df250PulseTime_tree->Branch("pulse_number",&pulse_number,"pulse_number/i");
		Df250PulseTime_tree->Branch("quality_factor",&quality_factor,"quality_factor/i");
		Df250PulseTime_tree->Branch("time",&time,"time/I");
		f250PTtree_exists = 1;
	}
	eventnum = eventnumber;
	/// Loop over all Df250PulseTime objects in this event
	for(unsigned int c_chan=0; c_chan<num_f250PT; c_chan++){
		channelnum = c_chan;
		const Df250PulseTime *f250PulseTime = f250PulseTime_vec[c_chan];
		rocid = f250PulseTime->rocid;
		slot = f250PulseTime->slot;
		channel = f250PulseTime->channel;
		itrigger = f250PulseTime->itrigger;
		pulse_number = f250PulseTime->pulse_number;
		quality_factor= f250PulseTime->quality_factor;
		time = f250PulseTime->time;
		Df250PulseTime_tree->Fill();
	}

	/// Df250PulsePedestal
	/// Get a vector of Df250PulsePedestal objects for this event (1 object for each crate/slot/channel above threshold)
	vector<const Df250PulsePedestal*> f250PulsePedestal_vec;
	loop->Get(f250PulsePedestal_vec);
	sort(f250PulsePedestal_vec.begin(), f250PulsePedestal_vec.end(), Df250PulsePedestal_cmp);
	unsigned int num_f250PP = f250PulsePedestal_vec.size();
	/// Create tree if doesn't exist and objects found.
	if (!f250PPtree_exists && num_f250PP>0) {
		printf("DAQTree >>eventnum %i, found %4i Df250PulsePedestal objects\n",eventnumber,num_f250PP);
		printf("DAQTree >>Creating tree f250PulsePedestal_tree\n");
		Df250PulsePedestal_tree = new TTree("Df250PulsePedestal",
										"tree of flash 250 pulse times for each channel and event");
		Df250PulsePedestal_tree->Branch("channelnum",&channelnum,"channelnum/i");
		Df250PulsePedestal_tree->Branch("eventnum",&eventnum,"eventnum/i");
		Df250PulsePedestal_tree->Branch("rocid",&rocid,"rocid/i");
		Df250PulsePedestal_tree->Branch("slot",&slot,"slot/i");
		Df250PulsePedestal_tree->Branch("channel",&channel,"channel/i");
		Df250PulsePedestal_tree->Branch("itrigger",&itrigger,"itrigger/i");
		Df250PulsePedestal_tree->Branch("pulse_number",&pulse_number,"pulse_number/i");
		Df250PulsePedestal_tree->Branch("pedestal",&pedestal,"pedestal/i");
		Df250PulsePedestal_tree->Branch("pulse_peak",&pulse_peak,"pulse_peak/I");
		f250PPtree_exists = 1;
	}
	eventnum = eventnumber;
	/// Loop over all Df250PulsePedestal objects in this event
	for(unsigned int c_chan=0; c_chan<num_f250PP; c_chan++){
		channelnum = c_chan;
		const Df250PulsePedestal *f250PulsePedestal = f250PulsePedestal_vec[c_chan];
		rocid = f250PulsePedestal->rocid;
		slot = f250PulsePedestal->slot;
		channel = f250PulsePedestal->channel;
		itrigger = f250PulsePedestal->itrigger;
		pulse_number = f250PulsePedestal->pulse_number;
		pedestal= f250PulsePedestal->pedestal;
		pulse_peak = f250PulsePedestal->pulse_peak;
		Df250PulsePedestal_tree->Fill();
	}


	/// Df250TriggerTime
	/// Get a vector of Df250TriggerTime objects for this event (1 object for each crate/slot)
	vector<const Df250TriggerTime*> f250TriggerTime_vec;
	loop->Get(f250TriggerTime_vec);
	sort(f250TriggerTime_vec.begin(), f250TriggerTime_vec.end(), Df250TriggerTime_cmp);
	unsigned int num_f250TT = f250TriggerTime_vec.size();
	/// Create tree if doesn't exist and objects found.
	if (!f250TTtree_exists && num_f250TT>0) {
		printf("DAQTree >>eventnum %i, found %4i Df250TriggerTime objects\n",eventnumber,num_f250TT);
		printf("DAQTree >>Creating tree Df250TriggerTime_tree\n");
		Df250TriggerTime_tree = new TTree("Df250TriggerTime",
										  "tree of flash 250 trigger times for each slot and event");
		Df250TriggerTime_tree->Branch("eventnum",&eventnum,"eventnum/i");
		Df250TriggerTime_tree->Branch("rocid",&rocid,"rocid/i");
		Df250TriggerTime_tree->Branch("slot",&slot,"slot/i");
		Df250TriggerTime_tree->Branch("itrigger",&itrigger,"itrigger/i");
		Df250TriggerTime_tree->Branch("time",&time,"time/I");
		f250TTtree_exists = 1;
	}
	eventnum = eventnumber;
	/// Loop over all Df250TriggerTime objects in this event
	for(unsigned int c_chan=0; c_chan<num_f250TT; c_chan++){
		channelnum = c_chan;
		const Df250TriggerTime *f250TriggerTime = f250TriggerTime_vec[c_chan];
		rocid = f250TriggerTime->rocid;
		slot = f250TriggerTime->slot;
		itrigger = f250TriggerTime->itrigger;
		time = f250TriggerTime->time;
		Df250TriggerTime_tree->Fill();
	}

	japp->RootUnLock();
	
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_DAQTree::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_DAQTree::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}



/* emacs
 * Local Variables:
 * mode:C++
 * mode:font-lock
 * c-file-style: "stroustrup"
 * tab-width: 4
 * End:
 */
