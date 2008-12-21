// $Id: DEventProcessor_acceptance_hists.h 1816 2006-06-06 14:38:18Z davidl $
//
//    File: DEventProcessor_acceptance_hists.h
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DEventProcessor_acceptance_hists_
#define _DEventProcessor_acceptance_hists_

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include <JANA/JFactory.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>
using namespace jana;

class DEventProcessor_acceptance_hists:public JEventProcessor{

	public:
		DEventProcessor_acceptance_hists();
		~DEventProcessor_acceptance_hists();
		
		TH2F *CDC, *FDC, *CDC_FDC;
		TH2F *BCAL, *FCAL, *TOF;
		TH2F *thrown_charged, *thrown_photon;
		
		TH1D *FDC_anode_hits_per_event;
		TH1D *FDC_anode_hits_per_layer;
		TH1D *FDC_anode_hits_per_wire;
		
		TH1D *CDC_nhits_vs_pthrown;
		TH1D *FDC_nhits_vs_pthrown;
		TH1D *pthrown;

	private:
		jerror_t init(void);	///< Invoked via DEventProcessor virtual method
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
		jerror_t erun(void);					///< Invoked via DEventProcessor virtual method
		jerror_t fini(void);					///< Invoked via DEventProcessor virtual method

};

#endif // _DEventProcessor_acceptance_hists_

