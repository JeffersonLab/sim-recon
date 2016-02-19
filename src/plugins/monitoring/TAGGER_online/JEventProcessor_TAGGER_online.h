// $Id$
//
//    File: JEventProcessor_TAGGER_online.h
// Created: Thu Feb 18 07:45:18 EST 2016
// Creator: jrsteven (on Linux gluon110.jlab.org 2.6.32-358.23.2.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_TAGGER_online_
#define _JEventProcessor_TAGGER_online_

#include <JANA/JEventProcessor.h>

#include <TAGGER/DTAGMHit.h>
#include <START_COUNTER/DSCHit.h>
#include <PID/DBeamPhoton.h>

#include "TH1.h"
#include "TH2.h"
#include "TDirectoryFile.h"

class JEventProcessor_TAGGER_online:public jana::JEventProcessor{
	public:
		JEventProcessor_TAGGER_online();
		~JEventProcessor_TAGGER_online();
		const char* className(void){return "JEventProcessor_TAGGER_online";}

	private:
		TH2D *dTaggerEnergy_DeltaTSC;

		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_TAGGER_online_

