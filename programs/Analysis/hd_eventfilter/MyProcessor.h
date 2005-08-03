// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.h
//
/// Example program for a Hall-D analyzer which uses DANA
///

#include "DEventProcessor.h"
#include "DEventLoop.h"
#include "hddm_s.h"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>


class MyProcessor:public DEventProcessor
{
	public:
		derror_t init(void);										///< Called once at program start.
		derror_t brun(DEventLoop *eventLoop, int runnumber){return NOERROR;}	///< Called everytime a new run number is detected.
		derror_t evnt(DEventLoop *eventLoop, int eventnumber);						///< Called every event.
		derror_t erun(void){return NOERROR;}				///< Called everytime run number changes, provided brun has been called.
		derror_t fini(void);										///< Called after last event of last event source has been processed.

		string filename;
		s_iostream_t *file;
		unsigned long Nevents_written;
};
