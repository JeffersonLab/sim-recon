// $Id$
// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.h
//
/// hd_dump print event info to screen
///

#include "DEventProcessor.h"
#include "DEventLoop.h"
#include "DFactory.h"

extern int PAUSE_BETWEEN_EVENTS;
extern int SKIP_BORING_EVENTS;
extern int PRINT_ALL;

extern vector<string> toprint;

class MyProcessor:public DEventProcessor
{
	public:
		derror_t init(void){return NOERROR;};				///< Called once at program start.
		derror_t brun(DEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		derror_t evnt(DEventLoop *eventLoop, int eventnumber);						///< Called every event.
		derror_t erun(void){return NOERROR;};				///< Called everytime run number changes, provided brun has been called.
		derror_t fini(void){return NOERROR;};				///< Called after last event of last event source has been processed.

};
