// $Id$
//
//    File: DEventProcessor_p4pi_hists.h
// Created: Mon Aug 29 16:20:58 EDT 2016
// Creator: aaustreg (on Linux halld01.jlab.org 2.6.32-642.3.1.el6.x86_64 x86_64)
//

#ifndef _DEventProcessor_p4pi_hists_
#define _DEventProcessor_p4pi_hists_

#include <JANA/JEventProcessor.h>
#include <JANA/JApplication.h>

#include <ANALYSIS/DEventWriterROOT.h>
#include <HDDM/DEventWriterREST.h>
#include <ANALYSIS/DHistogramActions.h>

#include "DFactoryGenerator_p4pi_hists.h"

using namespace jana;
using namespace std;

class DEventProcessor_p4pi_hists : public jana::JEventProcessor
{
	public:
		const char* className(void){return "DEventProcessor_p4pi_hists";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop* locEventLoop, int32_t locRunNumber);	///< Called every time a new run number is detected.
		jerror_t evnt(jana::JEventLoop* locEventLoop, uint64_t locEventNumber);	///< Called every event.
		jerror_t erun(void);						///< Called every time run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		//For non-custom reaction-independent histograms, it is recommended that you simultaneously run the monitoring_hists plugin instead of defining them here

		//For making EventStore skims
		int Get_FileNumber(JEventLoop* locEventLoop) const;
		ofstream dEventStoreSkimStream;
};

#endif // _DEventProcessor_p4pi_hists_

