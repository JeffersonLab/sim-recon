// $Id: DEventProcessor_mcthrown_hists.h 1816 2006-06-06 14:38:18Z davidl $
//
//    File: DEventProcessor_mcthrown_hists.h
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DEventProcessor_mcthrown_hists_
#define _DEventProcessor_mcthrown_hists_

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include "DFactory.h"
#include "DEventProcessor.h"
#include "DEventLoop.h"

class DEventProcessor_mcthrown_hists:public DEventProcessor{

	public:
		DEventProcessor_mcthrown_hists();
		~DEventProcessor_mcthrown_hists();
		
		TFile* ROOTfile;
		TH1F *pmom, *theta, *phi, *energy;
		TH1F *Nparticles_per_event, *particle_type;
		TH3F *vertex;

	private:
		derror_t init(void);	///< Invoked via DEventProcessor virtual method
		derror_t evnt(DEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
		derror_t erun(void);					///< Invoked via DEventProcessor virtual method
		derror_t fini(void);					///< Invoked via DEventProcessor virtual method

};

#endif // _DEventProcessor_mcthrown_hists_

