// $Id$
//
//    File: DEventProcessor_track_hists.h
// Created: Sat Jun  3 13:58:03 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.4.0 powerpc)
//

#ifndef _DEventProcessor_track_hists_
#define _DEventProcessor_track_hists_

#include "DEventProcessor.h"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

class DEventProcessor_track_hists:public DEventProcessor{
	public:
		DEventProcessor_track_hists(){};
		~DEventProcessor_track_hists(){};
		const char* className(void){return "DEventProcessor_track_hists";}

	private:
		derror_t init(void);						///< Called once at program start.
		derror_t brun(DEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		derror_t evnt(DEventLoop *eventLoop, int eventnumber);	///< Called every event.
		derror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		derror_t fini(void);						///< Called after last event of last event source has been processed.

		TFile *ROOTfile;
		
		TH1F *track_p, *track_theta, *track_phi;
		TH1F *trackCandidate_p;
};

#endif // _DEventProcessor_track_hists_

