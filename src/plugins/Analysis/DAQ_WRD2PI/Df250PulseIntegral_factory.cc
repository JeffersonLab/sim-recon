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
	printf("Df250PulseIntegral_factory::init()\n");
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t Df250PulseIntegral_factory::brun(jana::JEventLoop *eventLoop, int32_t runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t Df250PulseIntegral_factory::evnt(JEventLoop *loop, uint64_t eventnumber)
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

	uint32_t pulse_number = 0;
	uint32_t quality_factor = 0;

	// Get a vector of objects for this event (1 object for each crate/slot/channel)
	vector<const Df250WindowRawData*> f250WindowRawData_vec;
	loop->Get(f250WindowRawData_vec);
	uint32_t Nchannels = f250WindowRawData_vec.size();

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
		uint32_t nsamples=samplesvector.size();
		uint32_t pedestalsum = 0;
		uint32_t signalsum = 0;

		// loop over the first X samples to calculate pedestal
		for (uint16_t c_samp=0; c_samp<ped_samples; c_samp++) {
			pedestalsum += samplesvector[c_samp];
		}
		// loop over the remaining samples to calculate integral
		for (uint16_t c_samp=ped_samples; c_samp<nsamples; c_samp++) {
			signalsum += samplesvector[c_samp];
		}
		//uint32_t pedestaleffect = ((pedestalsum * isamples)/ped_samples);
		// myintegral = signalsum - pedestaleffect;
		// if (myintegral > 10000) jout << myintegral << "  " <<  signalsum << "  " <<  pedestaleffect << "  " <<  pedestalsum << "  " << nsamples  << "  " << ped_samples  << "  " << isamples << " \n";
		myDf250PulseIntegral->pulse_number = pulse_number;
		myDf250PulseIntegral->quality_factor = quality_factor;
		myDf250PulseIntegral->integral = signalsum;
		myDf250PulseIntegral->pedestal = pedestalsum;
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
