// $Id$
//
//    File: DCustomAction_p2pi_hists.h
// Created: Wed Jan 21 16:53:41 EST 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _DCustomAction_p2pi_hists_
#define _DCustomAction_p2pi_hists_

#include <string>
#include <iostream>

#include "TH1.h"
#include "TLorentzRotation.h"

#include "JANA/JEventLoop.h"
#include "JANA/JApplication.h"

#include "ANALYSIS/DAnalysisAction.h"
#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DAnalysisUtilities.h"

using namespace std;
using namespace jana;

class DCustomAction_p2pi_hists : public DAnalysisAction
{
	public:

		DCustomAction_p2pi_hists(const DReaction* locReaction, bool locUseKinFitResultsFlag, string locActionUniqueString = "") : 
	        DAnalysisAction(locReaction, "Custom_p2pi_hists", locUseKinFitResultsFlag, locActionUniqueString) {}

		void Initialize(JEventLoop* locEventLoop);

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		// Optional: Useful utility functions.
		const DAnalysisUtilities* dAnalysisUtilities;

		// need PID algos for SC matching
                const DParticleID* dParticleID;

		Particle_t dMissingPID;

		//Store any histograms as member variables here
		TH1I *dEgamma;
		TH2I *dMM_M2pi, *dMM_M2pi_noEle, *dt_M2pi_noEle;
		TH2I *dMM2_M2pi, *dProton_dEdx_P, *dProton_P_Theta, *dDeltaE_M2pi, *dDeltaE_M2pi_ProtonTag;
		TH2I *dDalitz_p2pi, *dMppiplus_M2pi, *dMppiminus_M2pi, *dEgamma_M2pi;
		TH2I *dPiPlusPsi_t;
		TH2I *dPiPlusPsi_Egamma, *dProtonPhi_Egamma;
		TH2I *dBaryonM_CosTheta_Egamma1, *dBaryonM_CosTheta_Egamma2, *dBaryonM_CosTheta_Egamma3; 
};

#endif // _DCustomAction_p2pi_hists_

