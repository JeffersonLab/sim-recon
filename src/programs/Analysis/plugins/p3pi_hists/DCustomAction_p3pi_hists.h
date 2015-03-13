// $Id$
//
//    File: DCustomAction_p3pi_hists.h
// Created: Wed Jan 21 16:53:41 EST 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _DCustomAction_p3pi_hists_
#define _DCustomAction_p3pi_hists_

#include <string>
#include <iostream>

#include "TH1.h"

#include "JANA/JEventLoop.h"
#include "JANA/JApplication.h"

#include "ANALYSIS/DAnalysisAction.h"
#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DAnalysisUtilities.h"

using namespace std;
using namespace jana;

class DCustomAction_p3pi_hists : public DAnalysisAction
{
	public:

		DCustomAction_p3pi_hists(const DReaction* locReaction, bool locUseKinFitResultsFlag, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Custom_p3pi_hists", locUseKinFitResultsFlag, locActionUniqueString) {}

		void Initialize(JEventLoop* locEventLoop);

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		// Optional: Useful utility functions.
		const DAnalysisUtilities* dAnalysisUtilities;

		Particle_t dMissingPID;

		//Store any histograms as member variables here
		TH1I *dEgamma;
		TH2I *dMM_M3pi;

		TH2I *dMM2_M3pi, *dProton_dEdx_P, *dProton_P_Theta, *dDeltaE_M3pi;
};

#endif // _DCustomAction_p3pi_hists_

