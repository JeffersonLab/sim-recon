// $Id$
//
//    File: DEventProcessor_bcal_hists.h
// Created: Mon Apr  3 11:38:03 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.4.0 powerpc)
//

#ifndef _DEventProcessor_bcal_hists_
#define _DEventProcessor_bcal_hists_

#include <JANA/JEventProcessor.h>
using namespace jana;

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

class DEventProcessor_bcal_hists:public JEventProcessor{
	public:
		DEventProcessor_bcal_hists(){};
		~DEventProcessor_bcal_hists(){};
		const char* className(void){return "DEventProcessor_bcal_hists";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
		
		TH1F* two_gamma_mass, *two_gamma_mass_corr, *two_gamma_mass_cut;
		TH1F* bcal_fcal_two_gamma_mass, *bcal_fcal_two_gamma_mass_cut;
		TH2F* xy_shower;
		TH1F* z_shower;
		TH1F* E_shower;
		TH2F* Erec_over_Ethrown_vs_z;
		TH2F* Ecorr_over_Erec_vs_z;
		TH2F* Ereconstructed_vs_Ethrown;
		TH1F* Etot_truth, *Etot_hits;
		TH2F* Etruth_over_Ethrown_vs_z;
		TH2F *Edeposited_over_Ethrown_vs_z;
};

#endif // _DEventProcessor_bcal_hists_

