// $Id$
//
//    File: DReaction_factory_ReactionFilter.h
// Created: Mon Nov 21 17:54:40 EST 2016
// Creator: pmatt (on Darwin Pauls-MacBook-Pro-2.local 13.4.0 i386)
//

#ifndef _DReaction_factory_ReactionFilter_
#define _DReaction_factory_ReactionFilter_

#include <iostream>
#include <iomanip>

#include <JANA/JFactory.h>
#include <ANALYSIS/DReaction.h>
#include <ANALYSIS/DHistogramActions.h>
#include <ANALYSIS/DCutActions.h>

#include "DCustomAction_dEdxCut.h"
#include "FSInfo.h"

using namespace std;
using namespace jana;

class DReaction_factory_ReactionFilter : public jana::JFactory<DReaction>
{
	public:
		DReaction_factory_ReactionFilter()
		{
			// This is so that the created DReaction objects persist throughout the life of the program instead of being cleared each event. 
			SetFactoryFlag(PERSISTANT);
		}
		const char* Tag(void){return "ReactionFilter";}

	private:
		jerror_t init(void);
		jerror_t brun(JEventLoop* locEventLoop, int32_t locRunNumber);
		jerror_t evnt(JEventLoop* locEventLoop, uint64_t locEventNumber);
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		// Create Reaction Steps
		void Create_FirstStep(DReaction* locReaction, FSInfo* locFSInfo);
		void Create_DecaySteps(DReaction* locReaction, FSInfo* locFSInfo);
		void Create_DecayStep(DReaction* locReaction, FSInfo* locFSInfo, Particle_t locPID);

		// Actions & cuts
		void Define_LooseCuts(void);
		void Add_PreComboCuts(DReaction* locReaction, FSInfo* locFSInfo);
		void Add_PIDActions(DReaction* locReaction);
		void Add_MassHistograms(DReaction* locReaction, FSInfo* locFSInfo, bool locUseKinFitResultsFlag, string locBaseUniqueName = "");

		// User-input channels
		deque<FSInfo*> dFSInfos;

		// Keep track of DReactionSteps that have been created for decaying particles
		// These can be re-used between DReactions, allowing the analysis library to save memory for combo steps
		map<pair<Particle_t, bool>, set<DReactionStep*> > dDecayStepMap_All; //bool: true for mass constrained by kinfit, false if not

		// Cuts
		map<Particle_t, pair<double, double> > dMissingMassCuts; //Unknown = none missing //if negative, uses missing mass squared instead
		map<Particle_t, pair<double, double> > dInvariantMassCuts;
		map<Particle_t, map<DetectorSystem_t, double> > dPIDTimingCuts;

		double dBeamBunchPeriod;
		deque<DReactionStep*> dReactionStepPool; //to prevent memory leaks
};

#endif // _DReaction_factory_ReactionFilter_

