// $Id$
//
//    File: DEventProcessor_eloss.h
// Created: Mon Mar 22 15:58:47 EDT 2010
// Creator: davidl (on Darwin harriet.jlab.org 9.8.0 i386)
//

#ifndef _DEventProcessor_eloss_
#define _DEventProcessor_eloss_

#include <TTree.h>

#include <JANA/JEventProcessor.h>

class DEventProcessor_eloss:public jana::JEventProcessor{
	public:
		DEventProcessor_eloss();
		~DEventProcessor_eloss();
		const char* className(void){return "DEventProcessor_eloss";}
		
		typedef struct{
			int event;
			float s;
			float x,y,z;
			float P;
			float dP;
			int mech;
		}event_t;
		
		TTree *geant;
		TTree *dana;
		event_t geant_event;
		event_t dana_event;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DEventProcessor_eloss_

