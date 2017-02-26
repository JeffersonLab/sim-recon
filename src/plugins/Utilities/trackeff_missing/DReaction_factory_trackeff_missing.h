// $Id$
//
//    File: DReaction_factory_trackeff_missing.h
// Created: Wed Feb 25 08:58:19 EST 2015
// Creator: pmatt (on Linux pmattdesktop.jlab.org 2.6.32-504.8.1.el6.x86_64 x86_64)
//

#ifndef _DReaction_factory_trackeff_missing_
#define _DReaction_factory_trackeff_missing_

#include <iostream>
#include <iomanip>

#include <JANA/JFactory.h>
#include <ANALYSIS/DReaction.h>
#include <ANALYSIS/DHistogramActions.h>
#include <ANALYSIS/DCutActions.h>

using namespace std;
using namespace jana;

class DReaction_factory_trackeff_missing : public jana::JFactory<DReaction>
{
	public:
		DReaction_factory_trackeff_missing()
		{
			// This is so that the created DReaction objects persist throughout the life of the program instead of being cleared each event. 
			SetFactoryFlag(PERSISTANT);
		}
		const char* Tag(void){return "trackeff_missing";}

	private:
		jerror_t init(void);
		jerror_t brun(JEventLoop* locEventLoop, int32_t locRunNumber);
		jerror_t evnt(JEventLoop* locEventLoop, uint64_t locEventNumber);
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		// Actions & cuts
		void Define_LooseCuts(void);
		void Add_PIDActions(DReaction* locReaction);
		bool Add_MassCuts(DReaction* locReaction, bool locKinFitFlag = false); //returns false if no cuts placed
		void Add_MassHistograms(DReaction* locReaction, bool locKinFitFlag, string locBaseUniqueName = "");

		// Utilities
		set<Particle_t> Get_InvariantMassPIDs(DReaction* locReaction, bool locKinFitFlag = false);
		map<Particle_t, pair<int, deque<Particle_t> > > Get_MissingMassPIDs(DReaction* locReaction, bool locKinFitFlag = false);

		double dBeamBunchPeriod;
		deque<DReactionStep*> dReactionStepPool; //to prevent memory leaks

		// Cuts
		map<Particle_t, pair<double, double> > dMissingMassCuts; //Unknown = none missing //if negative, uses missing mass squared instead
		map<Particle_t, pair<double, double> > dInvariantMassCuts;
		map<Particle_t, pair<double, double> > dMissingMassCuts_KinFit; //Unknown = none missing //if negative, uses missing mass squared instead
		map<Particle_t, pair<double, double> > dInvariantMassCuts_KinFit;
		map<Particle_t, map<DetectorSystem_t, double> > dPIDTimingCuts;
};

#endif // _DReaction_factory_trackeff_missing_
