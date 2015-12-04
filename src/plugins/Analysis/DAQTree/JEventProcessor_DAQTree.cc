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
#include <DAQ/Df125PulseRawData.h>
#include <DAQ/Df125PulseIntegral.h>
#include <DAQ/Df125PulseTime.h>
#include <DAQ/Df125PulsePedestal.h>
#include <DAQ/Df125TriggerTime.h>
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
bool Df125PulseRawData_cmp(const Df125PulseRawData *a,const Df125PulseRawData *b){
	// sort by crate, then by slot, then by channel, then by trigger number
	if (a->rocid   != b->rocid)   return a->rocid < b->rocid;
	if (a->slot    != b->slot )   return a->slot < b->slot;
	if (a->channel != b->channel) return a->channel < b->channel;
	return a->itrigger < b->itrigger;
}
bool Df125PulseIntegral_cmp(const Df125PulseIntegral *a,const Df125PulseIntegral *b){
	// sort by crate, then by slot, then by channel, then by trigger number
	if (a->rocid   != b->rocid)   return a->rocid < b->rocid;
	if (a->slot    != b->slot )   return a->slot < b->slot;
	if (a->channel != b->channel) return a->channel < b->channel;
	return a->itrigger < b->itrigger;
}
bool Df125PulseTime_cmp(const Df125PulseTime *a,const Df125PulseTime *b){
	// sort by crate, then by slot, then by channel, then by trigger number
	if (a->rocid   != b->rocid)   return a->rocid < b->rocid;
	if (a->slot    != b->slot )   return a->slot < b->slot;
	if (a->channel != b->channel) return a->channel < b->channel;
	return a->itrigger < b->itrigger;
}
bool Df125PulsePedestal_cmp(const Df125PulsePedestal *a,const Df125PulsePedestal *b){
	// sort by crate, then by slot, then by channel, then by trigger number
	if (a->rocid   != b->rocid)   return a->rocid < b->rocid;
	if (a->slot    != b->slot )   return a->slot < b->slot;
	if (a->channel != b->channel) return a->channel < b->channel;
	return a->itrigger < b->itrigger;
}
bool Df125TriggerTime_cmp(const Df125TriggerTime *a,const Df125TriggerTime *b){
	// sort by crate, then by slot, then by channel, then by trigger number
	if (a->rocid   != b->rocid)   return a->rocid < b->rocid;
	if (a->slot    != b->slot )   return a->slot < b->slot;
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
	f125PRDtree_exists = 0;
	f125PItree_exists = 0;
	f125PTtree_exists = 0;
	f125PPtree_exists = 0;
	f125TTtree_exists = 0;
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
jerror_t JEventProcessor_DAQTree::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_DAQTree::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// loop->Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.
	// Here's an example:
	//

	vector<const Df125WindowRawData*> f125WindowRawData_vec;
	loop->Get(f125WindowRawData_vec);
	sort(f125WindowRawData_vec.begin(), f125WindowRawData_vec.end(), Df125WindowRawData_cmp);

	vector<const Df125PulseRawData*> f125PulseRawData_vec;
	loop->Get(f125PulseRawData_vec);
	sort(f125PulseRawData_vec.begin(), f125PulseRawData_vec.end(), Df125PulseRawData_cmp);

	vector<const Df125PulseIntegral*> f125PulseIntegral_vec;
	loop->Get(f125PulseIntegral_vec);
	sort(f125PulseIntegral_vec.begin(), f125PulseIntegral_vec.end(), Df125PulseIntegral_cmp);

	vector<const Df125PulseTime*> f125PulseTime_vec;
	loop->Get(f125PulseTime_vec);
	sort(f125PulseTime_vec.begin(), f125PulseTime_vec.end(), Df125PulseTime_cmp);

	vector<const Df125PulsePedestal*> f125PulsePedestal_vec;
	loop->Get(f125PulsePedestal_vec);
	sort(f125PulsePedestal_vec.begin(), f125PulsePedestal_vec.end(), Df125PulsePedestal_cmp);

	vector<const Df125TriggerTime*> f125TriggerTime_vec;
	loop->Get(f125TriggerTime_vec);
	sort(f125TriggerTime_vec.begin(), f125TriggerTime_vec.end(), Df125TriggerTime_cmp);

	vector<const DF1TDCHit*> F1TDCHit_vec;
	loop->Get(F1TDCHit_vec);
	sort(F1TDCHit_vec.begin(), F1TDCHit_vec.end(), DF1TDCHit_cmp);

	vector<const DF1TDCTriggerTime*> F1TDCTriggerTime_vec;
	loop->Get(F1TDCTriggerTime_vec);
	sort(F1TDCTriggerTime_vec.begin(), F1TDCTriggerTime_vec.end(), DF1TDCTriggerTime_cmp);

	vector<const Df250WindowRawData*> f250WindowRawData_vec;
	loop->Get(f250WindowRawData_vec);
	sort(f250WindowRawData_vec.begin(), f250WindowRawData_vec.end(), Df250WindowRawData_cmp);

	vector<const Df250PulseRawData*> f250PulseRawData_vec;
	loop->Get(f250PulseRawData_vec);
	sort(f250PulseRawData_vec.begin(), f250PulseRawData_vec.end(), Df250PulseRawData_cmp);

	vector<const Df250PulseIntegral*> f250PulseIntegral_vec;
	loop->Get(f250PulseIntegral_vec);
	sort(f250PulseIntegral_vec.begin(), f250PulseIntegral_vec.end(), Df250PulseIntegral_cmp);

	vector<const Df250PulseTime*> f250PulseTime_vec;
	loop->Get(f250PulseTime_vec);
	sort(f250PulseTime_vec.begin(), f250PulseTime_vec.end(), Df250PulseTime_cmp);

	vector<const Df250PulsePedestal*> f250PulsePedestal_vec;
	loop->Get(f250PulsePedestal_vec);
	sort(f250PulsePedestal_vec.begin(), f250PulsePedestal_vec.end(), Df250PulsePedestal_cmp);

	vector<const Df250TriggerTime*> f250TriggerTime_vec;
	loop->Get(f250TriggerTime_vec);
	sort(f250TriggerTime_vec.begin(), f250TriggerTime_vec.end(), Df250TriggerTime_cmp);

	/// Trees are filled with data
	japp->RootWriteLock();


	/// Df125WindowRawData
	const uint32_t numDf125WRDpedsamps = 10;
	const Int_t Df125WRDminpeakheight = 100;
	/// Get a vector of Df125WindowRawData objects for this event (1 object for each crate/slot/channel above threshold)
	unsigned int num_f125WRD = f125WindowRawData_vec.size();
	/// Create tree if doesn't exist and objects found.
	if (!f125WRDtree_exists && num_f125WRD>0) {
		printf("DAQTree >>eventnum %llu, found %4i Df125WindowRawData objects\n",(long long unsigned int)eventnumber,num_f125WRD);
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


	/// Df125PulseRawData
	const uint32_t numDf125PRDpedsamps= 4;
	const Int_t Df125PRDminpeakheight = 11;
	/// Get a vector of Df125PulseRawData objects for this event (1 object for each crate/slot/channel above threshold)
	unsigned int num_f125PRD = f125PulseRawData_vec.size();
	/// Create tree if doesn't exist and objects found.
	if (!f125PRDtree_exists && num_f125PRD>0) {
		printf("DAQTree >>eventnum %llu, found %4i Df125PulseRawData objects\n",(long long unsigned int)eventnumber,num_f125PRD);
		printf("DAQTree >>Creating tree Df125PulseRawData_tree\n");
		Df125PulseRawData_tree = new TTree("Df125PulseRawData",
										   "tree of flash 125 pulse raw data for each channel and event");
		Df125PulseRawData_tree->Branch("channelnum",&channelnum,"channelnum/i");
		Df125PulseRawData_tree->Branch("eventnum",&eventnum,"eventnum/i");
		Df125PulseRawData_tree->Branch("rocid",&rocid,"rocid/i");
		Df125PulseRawData_tree->Branch("slot",&slot,"slot/i");
		Df125PulseRawData_tree->Branch("channel",&channel,"channel/i");
		Df125PulseRawData_tree->Branch("itrigger",&itrigger,"itrigger/i");
		Df125PulseRawData_tree->Branch("pulse_number",&pulse_number,"pulse_number/i");
		Df125PulseRawData_tree->Branch("first_sample_number",&first_sample_number,"first_sample_number/i");
		Df125PulseRawData_tree->Branch("waveform",&waveform);
		Df125PulseRawData_tree->Branch("nsamples",&nsamples,"nsamples/i");
		Df125PulseRawData_tree->Branch("w_integral",&w_integral,"w_integral/i");
		Df125PulseRawData_tree->Branch("w_min",&w_min,"w_min/i");
		Df125PulseRawData_tree->Branch("w_max",&w_max,"w_max/i");
		Df125PulseRawData_tree->Branch("w_samp1",&w_samp1,"w_samp1/i");
		Df125PulseRawData_tree->Branch("w_ped",&w_ped,"w_ped/i");
		Df125PulseRawData_tree->Branch("w_time",&w_time,"w_time/f");
		Df125PulseRawData_tree->Branch("invalid_samples",&invalid_samples,"invalid_samples/b");
		Df125PulseRawData_tree->Branch("overflow",&overflow,"overflow/b");
		f125PRDtree_exists = 1;
	}
	eventnum = eventnumber;
	// Loop over all  objects in this event
	for(unsigned int c_chan=0; c_chan<num_f125PRD; c_chan++){
		waveform.clear();
		channelnum = c_chan;
		const Df125PulseRawData *f125PulseRawData = f125PulseRawData_vec[c_chan];
		rocid = f125PulseRawData->rocid;
		slot = f125PulseRawData->slot;
		channel = f125PulseRawData->channel;
		itrigger = f125PulseRawData->itrigger;
		pulse_number  = f125PulseRawData->pulse_number;
		first_sample_number = f125PulseRawData->first_sample_number;
		invalid_samples = f125PulseRawData->invalid_samples;
		overflow = f125PulseRawData->overflow;
		// Get a vector of the samples for this channel
		const vector<uint16_t> &samplesvector = f125PulseRawData->samples;
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
				if (c_samp<numDf125PRDpedsamps) {
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
		if (peakheight > Df125PRDminpeakheight) { 
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
		Df125PulseRawData_tree->Fill();
	}

	/// Df125PulseIntegral
	/// Get a vector of Df125PulseIntegral objects for this event (1 object for each crate/slot/channel above threshold)
	unsigned int num_f125PI = f125PulseIntegral_vec.size();
	/// Create tree if doesn't exist and objects found.
	if (!f125PItree_exists && num_f125PI>0) {
		printf("DAQTree >>eventnum %llu, found %4i Df125PulseIntegral objects\n",(long long unsigned int)eventnumber,num_f125PI);
		printf("DAQTree >>Creating tree Df125PulseIntegral_tree\n");
		Df125PulseIntegral_tree = new TTree("Df125PulseIntegral",
											"tree of flash 125 pulse integral for each channel and event");
		Df125PulseIntegral_tree->Branch("channelnum",&channelnum,"channelnum/i");
		Df125PulseIntegral_tree->Branch("eventnum",&eventnum,"eventnum/i");
		Df125PulseIntegral_tree->Branch("rocid",&rocid,"rocid/i");
		Df125PulseIntegral_tree->Branch("slot",&slot,"slot/i");
		Df125PulseIntegral_tree->Branch("channel",&channel,"channel/i");
		Df125PulseIntegral_tree->Branch("itrigger",&itrigger,"itrigger/i");
		Df125PulseIntegral_tree->Branch("pulse_number",&pulse_number,"pulse_number/i");
		Df125PulseIntegral_tree->Branch("quality_factor",&quality_factor,"quality_factor/i");
		Df125PulseIntegral_tree->Branch("integral",&integral,"integral/i");
		Df125PulseIntegral_tree->Branch("pedestal",&pedestal,"pedestal/i");
		f125PItree_exists = 1;
	}
	eventnum = eventnumber;
	/// Loop over all Df125PulseIntegral objects in this event
	for(unsigned int c_chan=0; c_chan<num_f125PI; c_chan++){
		channelnum = c_chan;
		const Df125PulseIntegral *f125PulseIntegral = f125PulseIntegral_vec[c_chan];
		rocid = f125PulseIntegral->rocid;
		slot = f125PulseIntegral->slot;
		channel = f125PulseIntegral->channel;
		itrigger = f125PulseIntegral->itrigger;
		pulse_number = f125PulseIntegral->pulse_number;
		quality_factor= f125PulseIntegral->quality_factor;
		integral = f125PulseIntegral->integral;
		pedestal = f125PulseIntegral->pedestal;
		Df125PulseIntegral_tree->Fill();
	}


	/// Df125PulseTime
	/// Get a vector of Df125PulseTime objects for this event (1 object for each crate/slot/channel above threshold)
	unsigned int num_f125PT = f125PulseTime_vec.size();
	/// Create tree if doesn't exist and objects found.
	if (!f125PTtree_exists && num_f125PT>0) {
		printf("DAQTree >>eventnum %llu, found %4i Df125PulseTime objects\n",(long long unsigned int)eventnumber,num_f125PT);
		printf("DAQTree >>Creating tree f125PulseTime_tree\n");
		Df125PulseTime_tree = new TTree("Df125PulseTime",
										"tree of flash 125 pulse times for each channel and event");
		Df125PulseTime_tree->Branch("channelnum",&channelnum,"channelnum/i");
		Df125PulseTime_tree->Branch("eventnum",&eventnum,"eventnum/i");
		Df125PulseTime_tree->Branch("rocid",&rocid,"rocid/i");
		Df125PulseTime_tree->Branch("slot",&slot,"slot/i");
		Df125PulseTime_tree->Branch("channel",&channel,"channel/i");
		Df125PulseTime_tree->Branch("itrigger",&itrigger,"itrigger/i");
		Df125PulseTime_tree->Branch("pulse_number",&pulse_number,"pulse_number/i");
		Df125PulseTime_tree->Branch("quality_factor",&quality_factor,"quality_factor/i");
		Df125PulseTime_tree->Branch("time",&time,"time/i");
		f125PTtree_exists = 1;
	}
	eventnum = eventnumber;
	/// Loop over all Df125PulseTime objects in this event
	for(unsigned int c_chan=0; c_chan<num_f125PT; c_chan++){
		channelnum = c_chan;
		const Df125PulseTime *f125PulseTime = f125PulseTime_vec[c_chan];
		rocid = f125PulseTime->rocid;
		slot = f125PulseTime->slot;
		channel = f125PulseTime->channel;
		itrigger = f125PulseTime->itrigger;
		pulse_number = f125PulseTime->pulse_number;
		quality_factor= f125PulseTime->quality_factor;
		time = f125PulseTime->time;
		Df125PulseTime_tree->Fill();
	}

	/// Df125PulsePedestal
	/// Get a vector of Df125PulsePedestal objects for this event (1 object for each crate/slot/channel above threshold)
	unsigned int num_f125PP = f125PulsePedestal_vec.size();
	/// Create tree if doesn't exist and objects found.
	if (!f125PPtree_exists && num_f125PP>0) {
		printf("DAQTree >>eventnum %llu, found %4i Df125PulsePedestal objects\n",(long long unsigned int)eventnumber,num_f125PP);
		printf("DAQTree >>Creating tree f125PulsePedestal_tree\n");
		Df125PulsePedestal_tree = new TTree("Df125PulsePedestal",
										"tree of flash 125 pulse times for each channel and event");
		Df125PulsePedestal_tree->Branch("channelnum",&channelnum,"channelnum/i");
		Df125PulsePedestal_tree->Branch("eventnum",&eventnum,"eventnum/i");
		Df125PulsePedestal_tree->Branch("rocid",&rocid,"rocid/i");
		Df125PulsePedestal_tree->Branch("slot",&slot,"slot/i");
		Df125PulsePedestal_tree->Branch("channel",&channel,"channel/i");
		Df125PulsePedestal_tree->Branch("itrigger",&itrigger,"itrigger/i");
		Df125PulsePedestal_tree->Branch("pulse_number",&pulse_number,"pulse_number/i");
		Df125PulsePedestal_tree->Branch("pedestal",&pedestal,"pedestal/i");
		Df125PulsePedestal_tree->Branch("pulse_peak",&pulse_peak,"pulse_peak/i");
		f125PPtree_exists = 1;
	}
	eventnum = eventnumber;
	/// Loop over all Df125PulsePedestal objects in this event
	for(unsigned int c_chan=0; c_chan<num_f125PP; c_chan++){
		channelnum = c_chan;
		const Df125PulsePedestal *f125PulsePedestal = f125PulsePedestal_vec[c_chan];
		rocid = f125PulsePedestal->rocid;
		slot = f125PulsePedestal->slot;
		channel = f125PulsePedestal->channel;
		itrigger = f125PulsePedestal->itrigger;
		pulse_number = f125PulsePedestal->pulse_number;
		pedestal= f125PulsePedestal->pedestal;
		pulse_peak = f125PulsePedestal->pulse_peak;
		Df125PulsePedestal_tree->Fill();
	}


	/// Df125TriggerTime
	/// Get a vector of Df125TriggerTime objects for this event (1 object for each crate/slot)
	unsigned int num_f125TT = f125TriggerTime_vec.size();
	/// Create tree if doesn't exist and objects found.
	if (!f125TTtree_exists && num_f125TT>0) {
		printf("DAQTree >>eventnum %llu, found %4i Df125TriggerTime objects\n",(long long unsigned int)eventnumber,num_f125TT);
		printf("DAQTree >>Creating tree Df125TriggerTime_tree\n");
		Df125TriggerTime_tree = new TTree("Df125TriggerTime",
										  "tree of flash 125 trigger times for each slot and event");
		Df125TriggerTime_tree->Branch("eventnum",&eventnum,"eventnum/i");
		Df125TriggerTime_tree->Branch("rocid",&rocid,"rocid/i");
		Df125TriggerTime_tree->Branch("slot",&slot,"slot/i");
		Df125TriggerTime_tree->Branch("itrigger",&itrigger,"itrigger/i");
		Df125TriggerTime_tree->Branch("time",&time,"time/l");
		f125TTtree_exists = 1;
	}
	eventnum = eventnumber;
	/// Loop over all Df125TriggerTime objects in this event
	for(unsigned int c_chan=0; c_chan<num_f125TT; c_chan++){
		channelnum = c_chan;
		const Df125TriggerTime *f125TriggerTime = f125TriggerTime_vec[c_chan];
		rocid = f125TriggerTime->rocid;
		slot = f125TriggerTime->slot;
		itrigger = f125TriggerTime->itrigger;
		time = f125TriggerTime->time;
		Df125TriggerTime_tree->Fill();
	}


	/// DF1TDCHit
	/// Get a vector of DF1TDCHit objects for this event (1 object for each crate/slot/channel above threshold)
	unsigned int num_F1TDCH = F1TDCHit_vec.size();
	/// Create tree if doesn't exist and objects found.
	if (!F1TDCHtree_exists && num_F1TDCH>0) {
		printf("DAQTree >>eventnum %llu, found %4i DF1TDCHit objects\n",(long long unsigned int)eventnumber,num_F1TDCH);
		printf("DAQTree >>Creating tree F1TDCHit_tree\n");
		DF1TDCHit_tree = new TTree("DF1TDCHit",	"tree of F1 TDC hit times for each channel and event");
		DF1TDCHit_tree->Branch("channelnum",&channelnum,"channelnum/i");
		DF1TDCHit_tree->Branch("eventnum",&eventnum,"eventnum/i");
		DF1TDCHit_tree->Branch("rocid",&rocid,"rocid/i");
		DF1TDCHit_tree->Branch("slot",&slot,"slot/i");
		DF1TDCHit_tree->Branch("channel",&channel,"channel/i");
		DF1TDCHit_tree->Branch("itrigger",&itrigger,"itrigger/i");
		DF1TDCHit_tree->Branch("trig_time",&trig_time,"trig_time/i");
		DF1TDCHit_tree->Branch("time",&time,"time/i");
		DF1TDCHit_tree->Branch("data_word",&data_word,"data_word/i");
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
		data_word = F1TDCHit->data_word;
		DF1TDCHit_tree->Fill();
	}


	/// DF1TDCTriggerTime
	/// Get a vector of DF1TDCTriggerTime objects for this event (1 object for each crate/slot)
	unsigned int num_F1TDCTT = F1TDCTriggerTime_vec.size();
	/// Create tree if doesn't exist and objects found.
	if (!F1TDCTTtree_exists && num_F1TDCTT>0) {
		printf("DAQTree >>eventnum %llu, found %4i DF1TDCTriggerTime objects\n",(long long unsigned int)eventnumber,num_F1TDCTT);
		printf("DAQTree >>Creating tree DF1TDCTriggerTime_tree\n");
		DF1TDCTriggerTime_tree = new TTree("DF1TDCTriggerTime", "tree of F1 TDC trigger times for each slot and event");
		DF1TDCTriggerTime_tree->Branch("eventnum",&eventnum,"eventnum/i");
		DF1TDCTriggerTime_tree->Branch("rocid",&rocid,"rocid/i");
		DF1TDCTriggerTime_tree->Branch("slot",&slot,"slot/i");
		DF1TDCTriggerTime_tree->Branch("itrigger",&itrigger,"itrigger/i");
		DF1TDCTriggerTime_tree->Branch("time",&time,"time/l");
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
	unsigned int num_f250WRD = f250WindowRawData_vec.size();
	/// Create tree if doesn't exist and objects found.
	if (!f250WRDtree_exists && num_f250WRD>0) {
		printf("DAQTree >>eventnum %llu, found %4i Df250WindowRawData objects\n",(long long unsigned int)eventnumber,num_f250WRD);
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
	unsigned int num_f250PRD = f250PulseRawData_vec.size();
	/// Create tree if doesn't exist and objects found.
	if (!f250PRDtree_exists && num_f250PRD>0) {
		printf("DAQTree >>eventnum %llu, found %4i Df250PulseRawData objects\n",(long long unsigned int)eventnumber,num_f250PRD);
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
	unsigned int num_f250PI = f250PulseIntegral_vec.size();
	/// Create tree if doesn't exist and objects found.
	if (!f250PItree_exists && num_f250PI>0) {
		printf("DAQTree >>eventnum %llu, found %4i Df250PulseIntegral objects\n",(long long unsigned int)eventnumber,num_f250PI);
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
		Df250PulseIntegral_tree->Branch("integral",&integral,"integral/i");
		Df250PulseIntegral_tree->Branch("pedestal",&pedestal,"pedestal/i");
		Df250PulseIntegral_tree->Branch("nsamples_integral",&nsamples_integral,"nsamples_integral/i");
		Df250PulseIntegral_tree->Branch("nsamples_pedestal",&nsamples_pedestal,"nsamples_pedestal/i");
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
		nsamples_integral = f250PulseIntegral->nsamples_integral;
		nsamples_pedestal = f250PulseIntegral->nsamples_pedestal;
		Df250PulseIntegral_tree->Fill();
	}


	/// Df250PulseTime
	/// Get a vector of Df250PulseTime objects for this event (1 object for each crate/slot/channel above threshold)
	unsigned int num_f250PT = f250PulseTime_vec.size();
	/// Create tree if doesn't exist and objects found.
	if (!f250PTtree_exists && num_f250PT>0) {
		printf("DAQTree >>eventnum %llu, found %4i Df250PulseTime objects\n",(long long unsigned int)eventnumber,num_f250PT);
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
		Df250PulseTime_tree->Branch("time",&time,"time/i");
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
	unsigned int num_f250PP = f250PulsePedestal_vec.size();
	/// Create tree if doesn't exist and objects found.
	if (!f250PPtree_exists && num_f250PP>0) {
		printf("DAQTree >>eventnum %llu, found %4i Df250PulsePedestal objects\n",(long long unsigned int)eventnumber,num_f250PP);
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
		Df250PulsePedestal_tree->Branch("pulse_peak",&pulse_peak,"pulse_peak/i");
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
	unsigned int num_f250TT = f250TriggerTime_vec.size();
	/// Create tree if doesn't exist and objects found.
	if (!f250TTtree_exists && num_f250TT>0) {
		printf("DAQTree >>eventnum %llu, found %4i Df250TriggerTime objects\n",(long long unsigned int)eventnumber,num_f250TT);
		printf("DAQTree >>Creating tree Df250TriggerTime_tree\n");
		Df250TriggerTime_tree = new TTree("Df250TriggerTime",
										  "tree of flash 250 trigger times for each slot and event");
		Df250TriggerTime_tree->Branch("eventnum",&eventnum,"eventnum/i");
		Df250TriggerTime_tree->Branch("rocid",&rocid,"rocid/i");
		Df250TriggerTime_tree->Branch("slot",&slot,"slot/i");
		Df250TriggerTime_tree->Branch("itrigger",&itrigger,"itrigger/i");
		Df250TriggerTime_tree->Branch("time",&time,"time/l");
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
