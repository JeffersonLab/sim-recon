// $Id$
//
//    File: DCustomAction_HistMass_X_2000.h
// Created: Mon Dec  2 12:18:54 EST 2013
// Creator: pmatt (on Darwin pmattLaptop 10.8.0 i386)
//

#ifndef _DCustomAction_HistMass_X_2000_
#define _DCustomAction_HistMass_X_2000_

#include <string>
#include <iostream>

#include "JANA/JEventLoop.h"
#include "JANA/JApplication.h"

#include "ANALYSIS/DAnalysisAction.h"
#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DAnalysisUtilities.h"

using namespace std;
using namespace jana;

class DCustomAction_HistMass_X_2000 : public DAnalysisAction
{
	public:

		DCustomAction_HistMass_X_2000(const DReaction* locReaction, bool locUseKinFitResultsFlag, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Custom_HistMass_X_2000", locUseKinFitResultsFlag, locActionUniqueString) {}

		void Initialize(JEventLoop* locEventLoop);

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		// Optional: Useful utility functions.
		const DAnalysisUtilities* dAnalysisUtilities;

		//Store any histograms as member variables here
		TH1D* dMassHist;

		//Used for determining when the group of particles used for the invariant mass is identical to a previous combo (to prevent double counting)
		deque<set<const DKinematicData*> > dPastParticles;
};

#endif // _DCustomAction_HistMass_X_2000_

