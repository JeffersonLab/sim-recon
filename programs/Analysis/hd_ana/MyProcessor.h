// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.h
//
/// Example program for a Hall-D analyzer which uses DANA
///

#include "DEventProcessor.h"

#include <TFile.h>
#include <TH1.h>


class MyProcessor:public DEventProcessor
{
	public:
		derror_t init(void);					///< Called once at program start.
		derror_t brun(int runnumber);		///< Called everytime a new run number is detected.
		derror_t evnt(int eventnumber);	///< Called every event.
		derror_t erun(void);					///< Called everytime run number changes, provided brun has been called.
		derror_t fini(void);					///< Called after last event of last event source has been processed.

		TFile *ROOTfile;
		TH1 *Photon_Energy;
		TH1 *FCAL_Energy;
};
