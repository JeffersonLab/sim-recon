// Author: David Lawrence  June 25, 2004
//
//
// DTOF.h
//
/// Example program for a Hall-D analyzer which uses DANA
///

#include "DEventProcessor.h"

class DTOF:public DEventProcessor
{
	public:
		derror_t init(void);					///< Called once at program start.
		derror_t brun(int runnumber);		///< Called everytime a new run number is detected.
		derror_t evnt(int eventnumber);	///< Called every event.
		derror_t erun(void);					///< Called everytime run number changes, provided brun has been called.
		derror_t fini(void);					///< Called after last event of last event source has been processed.
};
