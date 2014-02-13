// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.h
//
/// Example program for a Hall-D analyzer which uses DANA
///

#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>


class MyProcessor:public JEventProcessor
{
	public:
		jerror_t init(void);					///< Called once at program start.
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t fini(void);					///< Called after last event of last event source has been processed.
		
		TH2F *R_vs_theta, *R_over_sintheta_vs_theta;
		
		TTree *fit_parms;
		float val[100]; // holds values used to fill fit_parms
		
		string CANDIDATE_TAG;
};
