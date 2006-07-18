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
#include <TTree.h>

//#include "CDCHitCopy.h"
//#include "FCALHitCopy.h"

class MyProcessor:public JEventProcessor
{
	public:
		jerror_t init(void);					///< Called once at program start.
		jerror_t brun(JEventLoop *eventLoop, int runnumber){return NOERROR;};	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void){return NOERROR;};				///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);					///< Called after last event of last event source has been processed.

		TFile *ROOTfile;
		TTree *ROOTtree;
		//CDCHitCopy *ROOTcdchit;
		//FCALHitCopy *ROOTfcalhit;
};
