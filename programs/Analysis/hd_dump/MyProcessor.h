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

extern int PAUSE_BETWEEN_EVENTS;
extern int SKIP_BORING_EVENTS;

extern int Ntoprint;
extern char *toprint[1024];


class MyProcessor:public DEventProcessor
{
	public:
		derror_t init(void){};				///< Called once at program start.
		derror_t brun(int runnumber){};	///< Called everytime a new run number is detected.
		derror_t evnt(int eventnumber);	///< Called every event.
		derror_t erun(void){};				///< Called everytime run number changes, provided brun has been called.
		derror_t fini(void){};				///< Called after last event of last event source has been processed.

};
