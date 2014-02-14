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
#include <DAQ/Df250WindowRawData.h>
#include <DAQ/Df250PulseIntegral.h>


bool Df250WindowRawData_cmp(const Df250WindowRawData *a,const Df250WindowRawData *b){
	// sort by crate, then by slot, then by channel
	if(a->rocid != b->rocid)return a->rocid < b->rocid;
	if(a->slot  != b->slot )return a->slot < b->slot;
	return a->channel < b->channel;
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

	printf("Welcome to JEventProcessor_DAQTree\n");

	japp->RootWriteLock();

	/// Trees are created
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

	/// Trees are created
	Df250PulseIntegral_tree = new TTree("Df250PulseIntegral",
										"tree of flash 250 pulse integral for each channel and event");
	Df250PulseIntegral_tree->Branch("channelnum",&channelnum,"channelnum/i");
	Df250PulseIntegral_tree->Branch("eventnum",&eventnum,"eventnum/i");

	japp->RootUnLock();
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
	eventnum = eventnumber;

	// Get a vector of Df250WindowRawData objects for this event (1 object for each crate/slot/channel)
	vector<const Df250WindowRawData*> f250WindowRawData_vec;
	loop->Get(f250WindowRawData_vec);
	sort(f250WindowRawData_vec.begin(), f250WindowRawData_vec.end(), Df250WindowRawData_cmp);
	uint32_t Nchannels = f250WindowRawData_vec.size();

	// Loop over all channels in this event
	for(unsigned int c_chan=0; c_chan<Nchannels; c_chan++){
		waveform.clear();
		channelnum = c_chan;
		const Df250WindowRawData *f250WindowRawData = f250WindowRawData_vec[c_chan];
		rocid = f250WindowRawData->rocid;
		slot = f250WindowRawData->slot;
		channel = f250WindowRawData->channel;
		itrigger = f250WindowRawData->itrigger;

		// Get a vector of the samples for this channel
		const vector<uint16_t> &samplesvector = f250WindowRawData->samples;
		nsamples=samplesvector.size();
		// loop over the samples to calculate integral, min, max
		for (uint16_t c_samp=0; c_samp<nsamples; c_samp++) {
			waveform.push_back(samplesvector[c_samp]); // push the sample into the waveform vector
			if (c_samp==0) {  // use first sample for initialization
				w_integral = samplesvector[0]; 
				w_min = samplesvector[0];
				w_max = samplesvector[0];
				w_samp1 = samplesvector[0];
			} else {
				w_integral += samplesvector[c_samp];
				if (w_min > samplesvector[c_samp]) w_min = samplesvector[c_samp];
				if (w_max < samplesvector[c_samp]) w_max = samplesvector[c_samp];
			}
		}
		Df250WindowRawData_tree->Fill();
	}

	// Get a vector of Df250WindowRawData objects for this event (1 object for each crate/slot/channel)
	vector<const Df250PulseIntegral*> f250PulseIntegral_vec;
	loop->Get(f250PulseIntegral_vec);
	Nchannels = f250PulseIntegral_vec.size();
	// Loop over all channels in this event
	for(unsigned int c_chan=0; c_chan<Nchannels; c_chan++){
		channelnum = c_chan;
		Df250PulseIntegral_tree->Fill();
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
