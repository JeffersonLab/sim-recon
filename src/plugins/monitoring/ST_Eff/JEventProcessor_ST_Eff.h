// $Id$
//
// File: JEventProcessor_ST_Eff.h
//

#ifndef _JEventProcessor_ST_Eff_
#define _JEventProcessor_ST_Eff_

#include <JANA/JEventProcessor.h>
#include <JANA/JApplication.h>

#include "TH1I.h"
#include "TH2I.h"

#include "START_COUNTER/DSCHit.h"
#include "TRIGGER/DTrigger.h"
#include "TRACKING/DTrackTimeBased.h"

#include "PID/DChargedTrack.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DParticleID.h"
#include "PID/DDetectorMatches.h"
#include "ANALYSIS/DCutActions.h"
#include "ANALYSIS/DTreeInterface.h"

#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <set>
#include <thread>

using namespace jana;
using namespace std;

class JEventProcessor_ST_Eff : public jana::JEventProcessor
{
	public:
		JEventProcessor_ST_Eff(){};
		~JEventProcessor_ST_Eff(){};
		const char* className(void){return "JEventProcessor_ST_Eff";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop* locEventLoop, int locRunNumber);	///< Called every time a new run number is detected.
		jerror_t evnt(jana::JEventLoop* locEventLoop, uint64_t locEventNumber);	///< Called every event.
		jerror_t erun(void);						///< Called every time run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		//TRACK REQUIREMENTS
		double dMinTrackingFOM, dMinPIDFOM;
		unsigned int dMinNumTrackHits;
		int dMinHitRingsPerCDCSuperlayer, dMinHitPlanesPerFDCPackage;
		DCutAction_TrackHitPattern* dCutAction_TrackHitPattern;

		//HISTOGRAMS
		double locHistMaxDeltaPhi;
		TH2I* dHist_HitFound;
		TH2I* dHist_HitTotal;

		//TREE
		DTreeInterface* dTreeInterface;
		//thread_local: Each thread has its own object: no lock needed
			//important: manages it's own data internally: don't want to call new/delete every event!
		static thread_local DTreeFillData dTreeFillData;
};

#endif // _JEventProcessor_ST_Eff_

