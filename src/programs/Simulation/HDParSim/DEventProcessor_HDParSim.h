// $Id$
//
//    File: DEventProcessor_HDParSim.h
// Created: Tue Feb  3 08:39:28 EST 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 i386)
//

#ifndef _DEventProcessor_HDParSim_
#define _DEventProcessor_HDParSim_

#include <pthread.h>

#include <TFile.h>
#include <TH2D.h>

#include <JANA/JEventProcessor.h>
using namespace jana;

#include "DTrackingResolutionGEANT.h"


class DEventProcessor_HDParSim:public JEventProcessor{
	public:
		DEventProcessor_HDParSim(const char *fname);
		~DEventProcessor_HDParSim(){};
		const char* className(void){return "DEventProcessor_HDParSim";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		
		const char *rootFileName;
		TFile *rootFile;
		TH2D *dp_over_p_vs_p;
		
		pthread_mutex_t root_mutex;
		DTrackingResolution *res;
};

#endif // _DEventProcessor_HDParSim_

