// $Id$
//
//    File: DEventProcessor_b1pi_hists.h
// Created: Mon Apr  3 11:38:03 EDT 2006
// Creator: pmatt (on Darwin swire-b241.jlab.org 8.4.0 powerpc)
//

#ifndef _DEventProcessor_b1pi_hists_
#define _DEventProcessor_b1pi_hists_

#include "TFile.h"
#include "TROOT.h"

#include "JANA/JEventProcessor.h"

#include "DANA/DApplication.h"
#include "ANALYSIS/DAnalysisResults.h"
#include "ANALYSIS/DEventWriterROOT.h"
#include "HDDM/DEventWriterREST.h"

#include "DFactoryGenerator_b1pi_hists.h"

using namespace jana;

class DEventProcessor_b1pi_hists : public JEventProcessor
{
	public:
		DEventProcessor_b1pi_hists(){};
		~DEventProcessor_b1pi_hists(){};
		const char* className(void){return "DEventProcessor_b1pi_hists";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};


#endif // _DEventProcessor_b1pi_hists_

