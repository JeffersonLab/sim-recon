// $Id$
//
//    File: DCustomAction_p2pi_taggerCoincidence.h
// Created: Thu Jan 22 08:06:18 EST 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _DCustomAction_p2pi_taggerCoincidence_
#define _DCustomAction_p2pi_taggerCoincidence_

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

class DCustomAction_p2pi_taggerCoincidence : public DAnalysisAction
{
	public:

                DCustomAction_p2pi_taggerCoincidence(const DReaction* locReaction, bool locUseKinFitResultsFlag, double locDeltaTMax, string locActionUniqueString = "") : 
	        DAnalysisAction(locReaction, "Custom_p2pi_taggerCoincidence", locUseKinFitResultsFlag, locActionUniqueString), dDeltaTMax(locDeltaTMax){}

		void Initialize(JEventLoop* locEventLoop);

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		double dDeltaTMax;

		// Optional: Useful utility functions.
		// const DAnalysisUtilities* dAnalysisUtilities;

		// need PID algos for SC matching
		const DParticleID* dParticleID;

		//Store any histograms as member variables here
		TH2I *dMatch_E_DeltaT_All, *dMatch_E_DeltaT_SC;
};

#endif // _DCustomAction_p2pi_taggerCoincidence_

