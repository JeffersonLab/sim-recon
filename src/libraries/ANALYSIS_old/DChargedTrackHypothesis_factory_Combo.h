// $Id$
//
//    File: DChargedTrackHypothesis_factory_Combo.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DChargedTrackHypothesis_factory_Combo_
#define _DChargedTrackHypothesis_factory_Combo_

#include "JANA/JFactory.h"
#include "TRACKING/DTrackTimeBased.h"
#include "PID/DChargedTrack.h"
#include "PID/DDetectorMatches.h"
#include "PID/DEventRFBunch.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DChargedTrackHypothesis_factory.h"
#include "ANALYSIS/DParticleComboBlueprint.h"
#include "ANALYSIS/DReaction.h"

using namespace jana;
using namespace std;

class DChargedTrackHypothesis_factory_Combo : public jana::JFactory<DChargedTrackHypothesis>
{
	public:
		DChargedTrackHypothesis_factory_Combo(){};
		~DChargedTrackHypothesis_factory_Combo(){};
		const char* Tag(void){return "Combo";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		void Build_RFPIDMap(map<Particle_t, vector<const DReaction*> >& locPIDMap, map<const DReaction*, set<const DEventRFBunch*> >& locRFBunchReactionMap, map<Particle_t, set<const DEventRFBunch*> >& locRFBunchPIDMap);
		void Create_TrackHypos(JEventLoop* locEventLoop, const DChargedTrack* locChargedTrack, map<Particle_t, vector<const DReaction*> >& locPIDMap, map<Particle_t, set<const DEventRFBunch*> >& locRFBunchPIDMap);
		DChargedTrackHypothesis* Create_TrackHypo(JEventLoop* locEventLoop, const DEventRFBunch* locEventRFBunch, const DChargedTrack* locChargedTrack, Particle_t locPID);

		DChargedTrackHypothesis_factory* dChargedTrackHypothesisFactory;
		const DDetectorMatches* dDetectorMatches;
		vector<const DReaction*> dReactions;
		map<pair<const DChargedTrack*, Particle_t>, const DTrackTimeBased*> dTimeBasedSourceMap;

		map<Particle_t, vector<const DReaction*> > dPositivelyChargedPIDs; //vector: reactions for which they are needed
		map<Particle_t, vector<const DReaction*> > dNegativelyChargedPIDs; //vector: reactions for which they are needed

		string dTrackSelectionTag;
};

#endif // _DChargedTrackHypothesis_factory_Combo_
