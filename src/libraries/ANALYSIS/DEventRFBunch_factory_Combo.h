// $Id$
//
//    File: DEventRFBunch_factory_Combo.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt
//

#ifndef _DEventRFBunch_factory_Combo_
#define _DEventRFBunch_factory_Combo_

#include <iostream>
#include <deque>

#include "TFile.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TDirectoryFile.h"
#include "TROOT.h"

#include "JANA/JFactory.h"
#include "particleType.h"

#include "TRACKING/DTrackTimeBased.h"
#include <BCAL/DBCALShower.h>
#include <FCAL/DFCALShower.h>

#include "PID/DEventRFBunch.h"
#include "PID/DParticleID.h"
#include "PID/DVertex.h"
#include "PID/DDetectorMatches.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DChargedTrack.h"
#include "PID/DNeutralShower.h"
#include "PID/DDetectorMatches.h"
#include "PID/DEventRFBunch_factory.h"
#include "RF/DRFTime_factory.h"

#include "ANALYSIS/DMCThrownMatching.h"
#include "ANALYSIS/DParticleComboBlueprint.h"

using namespace jana;
using namespace std;

class DEventRFBunch_factory_Combo:public jana::JFactory<DEventRFBunch>
{
	public:
		DEventRFBunch_factory_Combo(){use_factory = 1;}; //prevents JANA from searching the input file for these objects
		~DEventRFBunch_factory_Combo(){};
		const char* Tag(void){return "Combo";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		DEventRFBunch_factory* dEventRFBunchFactory;
		DRFTime_factory* dRFTimeFactory;

		const DParticleID* dParticleID;

		string dShowerSelectionTag;
		string dTrackSelectionTag;

		double dBeamBunchPeriod;
		double dTargetCenterZ;
		double dMinThrownMatchFOM;

		bool Get_StartTime(JEventLoop* locEventLoop, const DTrackTimeBased* locTrackTimeBased, double& locStartTime);
		double Calc_StartTime(const DNeutralShower* locNeutralShower, const DVertex* locVertex);
		int Find_BestRFBunchShift(double locRFHitTime, const vector<double>& locTimes, int& locBestNumVotes);

		bool Is_AllTruePID(const DMCThrownMatching* locMCThrownMatching, const DParticleComboBlueprint* locParticleComboBlueprint);

		map<const DReaction*, TH1I*> dHistMap_RFParticleDeltaT;

		map<const DReaction*, TH1I*> dHistMap_DeltaRFTime_TruePID; //given that the PIDs are all correct, diff between selected & true RF times
		map<const DReaction*, TH1I*> dHistMap_DeltaRFTime; //diff between selected & true RF times (all combos)
};

#endif // _DEventRFBunch_factory_Combo_

