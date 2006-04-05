// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.h
//
/// Example program for a Hall-D analyzer which uses DANA
///

#include "DEventProcessor.h"
#include "DEventLoop.h"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>


class MyProcessor:public DEventProcessor
{
	public:
		derror_t init(void);					///< Called once at program start.
		derror_t evnt(DEventLoop *eventLoop, int eventnumber);	///< Called every event.
		derror_t fini(void);					///< Called after last event of last event source has been processed.

		TFile *ROOTfile;
		TH1F *stats, *frac, *h4_dist, *h4_dist_primary;
		TH2F *delta_p;
		TH1F *delta_p_over_p;
		TH1F *hits_per_thrown_track; 
		
};
