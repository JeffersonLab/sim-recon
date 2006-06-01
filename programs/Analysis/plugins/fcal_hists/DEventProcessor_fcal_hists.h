// $Id$
//
//    File: DEventProcessor_bcal_hists.h
// Created: Mon Apr  3 11:38:03 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.4.0 powerpc)
//

#ifndef _DEventProcessor_fcal_hists_
#define _DEventProcessor_fcal_hists_

#include "DEventProcessor.h"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

class DEventProcessor_fcal_hists:public DEventProcessor{
	public:
		DEventProcessor_fcal_hists(){};
		~DEventProcessor_fcal_hists(){};
		const char* className(void){return "DEventProcessor_fcal_hists";}

	private:
		derror_t init(void);						///< Called once at program start.
		derror_t brun(DEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		derror_t evnt(DEventLoop *eventLoop, int eventnumber);	///< Called every event.
		derror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		derror_t fini(void);						///< Called after last event of last event source has been processed.

		TFile *ROOTfile;
		
		TH1F* two_gamma_mass, *two_gamma_mass_corr;
		TH1F* pi0_zdiff;
		TH2F* xy_shower;
		TH1F* z_shower;
		TH1F* E_shower;
		TH2F* E_over_Erec_vs_E;
		TH2F* E_over_Erec_vs_R;
		TH2F* E_over_Ereccorr_vs_z;
};

#endif // _DEventProcessor_fcal_hists_

