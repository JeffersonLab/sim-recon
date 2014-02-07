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


//-----------------------------------------
/// This plugin is designed to make a tree of low level data from the DAQ.
/// Each event and each channel is a new entry in the tree.
//-----------------------------------------

class JEventProcessor_DAQTree:public jana::JEventProcessor{
	public:
		JEventProcessor_DAQTree();
		~JEventProcessor_DAQTree();
		const char* className(void){return "JEventProcessor_DAQTree";}
		TH1D *histogram;
		TTree *sampletree;
		uint32_t channelnum;         /// Arbitrary global channel number
		uint32_t eventnum;	     /// Event number	
		uint32_t rocid;              /// Crate number
		uint32_t slot;               /// Slot number in crate
		uint32_t channel;            /// Channel number in slot
		uint32_t itrigger; 
///trigger number for cases when this hit was read in a multi-event block (DDAQAddress)
		vector<uint32_t> waveform;   /// STL vector of samples in the waveform for the event
		uint32_t nsamples;           /// Number of samples in the waveform
		uint32_t w_integral;         /// Sum of all samples in the waveform
		uint32_t w_min;              /// Minimum sample in the waveform
		uint32_t w_max;              /// Maximum sample in the waveform
		uint32_t w_samp1;            /// First sample in the waveform


	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_DAQTree_
