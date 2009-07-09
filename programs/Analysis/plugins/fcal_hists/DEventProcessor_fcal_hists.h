// $Id$
//
//    File: DEventProcessor_bcal_hists.h
// Created: Mon Apr  3 11:38:03 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.4.0 powerpc)
//

#ifndef _DEventProcessor_fcal_hists_
#define _DEventProcessor_fcal_hists_

#include <JANA/JEventProcessor.h>
using namespace jana;

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

class DEventProcessor_fcal_hists:public JEventProcessor{
	public:
		DEventProcessor_fcal_hists(){};
		~DEventProcessor_fcal_hists(){};
		const char* className(void){return "DEventProcessor_fcal_hists";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
		
		TH2D *dE_over_E_vs_E;
};

#endif // _DEventProcessor_fcal_hists_

