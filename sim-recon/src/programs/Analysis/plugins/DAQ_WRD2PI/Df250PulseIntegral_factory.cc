// $Id$
//
//    File: Df250PulseIntegral_factory.cc
// Created: Thu Feb 13 12:49:12 EST 2014
// Creator: dalton (on Linux gluon104.jlab.org 2.6.32-358.23.2.el6.x86_64 x86_64)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include <JANA/JApplication.h>
using namespace jana;
#include <DAQ/Df250WindowRawData.h>

#include "Df250PulseIntegral_factory.h"
#include "JFactoryGenerator_Df250PulseIntegral.h"


// Routine used to create our JEventProcessor

extern "C"{
	void InitPlugin(JApplication *app){
		InitJANAPlugin(app);
		app->AddFactoryGenerator(new JFactoryGenerator_Df250PulseIntegral());
	}
} // "C"

//------------------
// init
//------------------
jerror_t Df250PulseIntegral_factory::init(void)
{
	printf("Welcome to Df250PulseIntegral_factory\n");
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t Df250PulseIntegral_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t Df250PulseIntegral_factory::evnt(JEventLoop *loop, int eventnumber)
{

	// Code to generate factory data goes here. Add it like:
	//
	// Df250PulseIntegral *myDf250PulseIntegral = new Df250PulseIntegral;
	// myDf250PulseIntegral->x = x;
	// myDf250PulseIntegral->y = y;
	// ...
	// _data.push_back(myDf250PulseIntegral);
	//
	// Note that the objects you create here will be deleted later
	// by the system and the _data vector will be cleared automatically.

	mypulse_number = 0;
	myquality_factor = 0;

	// Get a vector of objects for this event (1 object for each crate/slot/channel)
	vector<const Df250WindowRawData*> f250WindowRawData_vec;
	loop->Get(f250WindowRawData_vec);
//	sort(f250WindowRawData_vec.begin(), f250WindowRawData_vec.end(), Df250WindowRawData_cmp);
	uint32_t Nchannels = f250WindowRawData_vec.size();
	printf("channles: %i,  ", Nchannels);

	// Loop over all channels in this event
	for(unsigned int c_chan=0; c_chan<Nchannels; c_chan++){
		// get Df250WindowRawData object
		const Df250WindowRawData *f250WindowRawData = f250WindowRawData_vec[c_chan];
		// create new Df250PulseIntegral object
		Df250PulseIntegral *myDf250PulseIntegral = new Df250PulseIntegral;
		myDf250PulseIntegral->rocid =f250WindowRawData->rocid;
		myDf250PulseIntegral->slot = f250WindowRawData->slot;
		myDf250PulseIntegral->channel = f250WindowRawData->channel;
		myDf250PulseIntegral->itrigger = f250WindowRawData->itrigger;

		// Get a vector of the samples for this channel
		const vector<uint16_t> &samplesvector = f250WindowRawData->samples;
		pedestal = 0;
		nsamples=samplesvector.size();
		// loop over the first X samples to calculate pedestal
		for (uint16_t c_samp=0; c_samp<ped_samples; c_samp++) {
			pedestal += samplesvector[c_samp];
		}
		pedestal /= ped_samples;
		// loop over the remaining samples to calculate integral
		for (uint16_t c_samp=ped_samples; c_samp<nsamples; c_samp++) {
			myintegral += samplesvector[c_samp];
		}
		myintegral = myintegral - ((pedestal * (nsamples - ped_samples))/ped_samples);
		myDf250PulseIntegral->pulse_number = mypulse_number;
		myDf250PulseIntegral->quality_factor = myquality_factor;
		myDf250PulseIntegral->integral = myintegral;
		_data.push_back(myDf250PulseIntegral);
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t Df250PulseIntegral_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t Df250PulseIntegral_factory::fini(void)
{
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
