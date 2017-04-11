// $Id$
//
//    File: DL3Trigger_factory.h
// Created: Wed Jul 31 14:34:24 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DL3Trigger_factory_
#define _DL3Trigger_factory_

#include <mutex>
using std::mutex;

#include <JANA/JFactory.h>
#include "DL3Trigger.h"

#ifdef HAVE_TMVA
#include <TMVA/Reader.h>
#endif

class DL3Trigger_factory:public jana::JFactory<DL3Trigger>{
	public:
		DL3Trigger_factory(){};
		~DL3Trigger_factory(){};

		double FRACTION_TO_KEEP;
		bool DO_WIRE_BASED_TRACKING;
		bool DO_BCAL_CLUSTER;
		uint32_t L1_TRIG_MASK;
		uint32_t L1_FP_TRIG_MASK;
		string MVA_WEIGHTS;
		double MVA_CUT;
		
		
#ifdef HAVE_TMVA
		TMVA::Reader *mvareader;
#endif
		Float_t Nstart_counter;
		Float_t Ntof;
		Float_t Nbcal_points;
		Float_t Nbcal_clusters;
		Float_t Ebcal_points;
		Float_t Ebcal_clusters;
		Float_t Nfcal_clusters;
		Float_t Efcal_clusters;
		Float_t Ntrack_candidates;
		Float_t Ptot_candidates;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DL3Trigger_factory_

