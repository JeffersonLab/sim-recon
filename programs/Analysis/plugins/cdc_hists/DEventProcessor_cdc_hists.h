// $Id: $
//
//    File: DEventProcessor_cdc_hists.h
//

#ifndef _DEventProcessor_cdc_hists_
#define _DEventProcessor_cdc_hists_

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include "DFactory.h"
#include "DEventProcessor.h"
#include "DEventLoop.h"

class DEventProcessor_cdc_hists:public DEventProcessor{

	public:
		DEventProcessor_cdc_hists();
		~DEventProcessor_cdc_hists();
		
		TFile* ROOTfile;
		TH1F *dE;

	private:
		derror_t init(void);	///< Invoked via DEventProcessor virtual method
		derror_t evnt(DEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
		derror_t erun(void);					///< Invoked via DEventProcessor virtual method
		derror_t fini(void);					///< Invoked via DEventProcessor virtual method

};

#endif // _DEventProcessor_cdc_hists_

