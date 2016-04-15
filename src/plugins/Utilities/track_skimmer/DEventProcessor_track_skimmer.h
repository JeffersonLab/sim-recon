// $Id$
//
//    File: DEventProcessor_track_skimmer.h
// Created: Tue Jan 13 11:08:16 EST 2015
// Creator: Paul (on Darwin Pauls-MacBook-Pro.local 14.0.0 i386)
//

#ifndef _DEventProcessor_track_skimmer_
#define _DEventProcessor_track_skimmer_

#include <map>
#include <fstream>

#include <JANA/JEventProcessor.h>
#include <JANA/JApplication.h>

#include <ANALYSIS/DEventWriterROOT.h>
#include <HDDM/DEventWriterREST.h>
#include <ANALYSIS/DHistogramActions.h>

#include "DFactoryGenerator_track_skimmer.h"

using namespace jana;
using namespace std;

class DEventProcessor_track_skimmer : public jana::JEventProcessor
{
	public:
		const char* className(void){return "DEventProcessor_track_skimmer";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop* locEventLoop, int locRunNumber);	///< Called every time a new run number is detected.
		jerror_t evnt(jana::JEventLoop* locEventLoop, uint64_t locEventNumber);	///< Called every event.
		jerror_t erun(void);						///< Called every time run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		int Get_FileNumber(JEventLoop* locEventLoop) const;

		map<string, ofstream*> dIDXAStreamMap;
};

#endif // _DEventProcessor_track_skimmer_

