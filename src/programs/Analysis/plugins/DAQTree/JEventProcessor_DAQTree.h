// $Id$
//
//    File: JEventProcessor_DAQTree.h
// Created: Tue Oct 22 14:55:40 EDT 2013
// Creator: dalton (on Linux gluon45.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_DAQTree_
#define _JEventProcessor_DAQTree_

#include <JANA/JEventProcessor.h>
#include "TH1D.h"
#include "TTree.h"

#include <stdint.h>

// Doxygen documentation
/** 
 This plugin is designed to make root trees of low level data from the DAQ.
 Each low level data type will have its own tree in the file.
 Each event and each channel is a new entry in the tree.

 Trees will be named after the low level data type used to fill them.
 Currently the data types that are supported are: /n
   Df250WindowRawData /n
 The data types that will ultimately be supported are: /n
   Df125PulseIntegral /n
   Df125PulseTime /n
   Df125TriggerTime /n
   DF1TDCHit /n
   DF1TDCTriggerTime /n
   Df250PulseIntegral /n
   Df250PulseRawData /n
   Df250PulseTime /n
   Df250StreamingRawData /n
   Df250TriggerTime /n
   Df250WindowRawData /n
   Df250WindowSum /n

 Wiki documentation can be found here: https://halldweb1.jlab.org/wiki/index.php/DAQTree_plugin
*/

class JEventProcessor_DAQTree:public jana::JEventProcessor{
	public:
		JEventProcessor_DAQTree();
		~JEventProcessor_DAQTree();
		const char* className(void){return "JEventProcessor_DAQTree";}
		TTree *Df250WindowRawData_tree;
		uint32_t channelnum;         ///< Arbitrary global channel number (sorted by crate, slot, channel)
		uint32_t eventnum;	         ///< Event number	
		uint32_t rocid;              ///< (from DDAQAddress) Crate number
		uint32_t slot;               ///< (from DDAQAddress) Slot number in crate
		uint32_t channel;            ///< (from DDAQAddress) Channel number in slot
		uint32_t itrigger;           ///< (from DDAQAddress) Trigger number for cases when this hit was read in a multi-event block (from DDAQAddress)
		vector<uint32_t> waveform;   ///< STL vector of samples in the waveform for the event
		uint32_t nsamples;           ///< Number of samples in the waveform
		uint32_t w_integral;         ///< Sum of all samples in the waveform
		uint32_t w_min;              ///< Minimum sample in the waveform
		uint32_t w_max;              ///< Maximum sample in the waveform
		uint32_t w_samp1;            ///< First sample in the waveform  (for simple analysis in case the STL vector is difficult to access)
        TTree *Df250PulseIntegral_tree;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_DAQTree_
