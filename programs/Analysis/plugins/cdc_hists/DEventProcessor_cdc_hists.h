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

#include <JANA/JFactory.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>

class DEventProcessor_cdc_hists:public JEventProcessor{

	public:
		DEventProcessor_cdc_hists();
		~DEventProcessor_cdc_hists();
		
		TH1F *dE;

	private:
		jerror_t init(void);	///< Invoked via DEventProcessor virtual method
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
		jerror_t erun(void);					///< Invoked via DEventProcessor virtual method
		jerror_t fini(void);					///< Invoked via DEventProcessor virtual method

};

#endif // _DEventProcessor_cdc_hists_

