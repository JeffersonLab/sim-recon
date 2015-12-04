// $Id$
//
//    File: DEventProcessor_run_summary.h
// Created: Tue Nov 18 15:44:17 EST 2014
// Creator: sdobbs (on Linux ifarm1102 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#ifndef _DEventProcessor_run_summary_
#define _DEventProcessor_run_summary_

#include <JANA/JEventProcessor.h>
#include <JANA/JApplication.h>

#include <ANALYSIS/DEventWriterROOT.h>
#include <HDDM/DEventWriterREST.h>
#include <ANALYSIS/DHistogramActions.h>

#include <TTree.h>

#include "DEPICSstore.h"

///#include "DFactoryGenerator_run_summary.h"

using namespace jana;
using namespace std;

class DEventProcessor_run_summary : public jana::JEventProcessor
{
	public:
		const char* className(void){return "DEventProcessor_run_summary";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop* locEventLoop, int locRunNumber);	///< Called every time a new run number is detected.
		jerror_t evnt(jana::JEventLoop* locEventLoop, uint64_t locEventNumber);	///< Called every event.
		jerror_t erun(void);						///< Called every time run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		int current_run_number;
		TTree *conditions_tree;
		DEPICSstore *epics_info;
};

#endif // _DEventProcessor_run_summary_

