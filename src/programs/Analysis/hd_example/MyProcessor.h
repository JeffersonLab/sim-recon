// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.h
//
/// Example program for a Hall-D analyzer which uses DANA
///

#ifndef _MyProcessor_
#define _MyProcessor_

#include "JANA/JEventProcessor.h"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>


using namespace jana;

class MyProcessor:public JEventProcessor
{
	public:
		jerror_t init(void);										///< Called once at program start.
		jerror_t brun(JEventLoop *eventLoop, int runnumber){return NOERROR;}	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);						///< Called every event.
		jerror_t erun(void){return NOERROR;}				///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);										///< Called after last event of last event source has been processed.

		TFile *ROOTfile;
		TH2 *cdc_ring_vs_straw;
		TH1 *fcal_y_vs_x, *fcalhitE;
};

#endif
