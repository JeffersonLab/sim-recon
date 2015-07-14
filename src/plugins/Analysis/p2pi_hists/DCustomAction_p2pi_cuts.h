// $Id$
//
//    File: DCustomAction_p2pi_cuts.h
// Created: Wed Jan 21 16:53:41 EST 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _DCustomAction_p2pi_cuts_
#define _DCustomAction_p2pi_cuts_

#include <string>
#include <iostream>

#include "TH1.h"

#include "JANA/JEventLoop.h"
#include "JANA/JApplication.h"

#include "TAGGER/DTAGHGeometry.h"
#include "ANALYSIS/DAnalysisAction.h"
#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DAnalysisUtilities.h"

using namespace std;
using namespace jana;

class DCustomAction_p2pi_cuts : public DAnalysisAction
{
	public:

		DCustomAction_p2pi_cuts(const DReaction* locReaction, bool locUseKinFitResultsFlag, string locActionUniqueString = "") : 
	        DAnalysisAction(locReaction, "Custom_p2pi_cuts", locUseKinFitResultsFlag, locActionUniqueString) {}

		void Initialize(JEventLoop* locEventLoop);

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		// Optional: Useful utility functions.
		const DAnalysisUtilities* dAnalysisUtilities;

		Particle_t dMissingPID;
};

#endif // _DCustomAction_p2pi_cuts_

