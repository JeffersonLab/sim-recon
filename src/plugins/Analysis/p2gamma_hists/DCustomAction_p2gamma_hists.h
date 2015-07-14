// $Id$
//
//    File: DCustomAction_p2gamma_hists.h
// Created: Wed Jan 21 16:53:41 EST 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _DCustomAction_p2gamma_hists_
#define _DCustomAction_p2gamma_hists_

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

class DCustomAction_p2gamma_hists : public DAnalysisAction
{
	public:

		DCustomAction_p2gamma_hists(const DReaction* locReaction, bool locUseKinFitResultsFlag, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Custom_p2gamma_hists", locUseKinFitResultsFlag, locActionUniqueString) {}

		void Initialize(JEventLoop* locEventLoop);

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		// Optional: Useful utility functions.
//		const DAnalysisUtilities* dAnalysisUtilities;
		
		set<set<pair<const JObject*, Particle_t> > > dPreviousSourceObjects;

		//Store any histograms as member variables here
		TH1I *dEgamma;

		TH2I *dMM2_M2g, *dProton_dEdx_P, *dProton_P_Theta, *dProtonPhi_Egamma, *dProtonPhi_Theta, *dProtonPhi_t;
		TH2I *dPi0Phi_Egamma, *dPi0Phi_Theta, *dDeltaE_M2g, *dPi0EgammaCorr;
		TH2I *dMM2_M2g_ProtonTag, *dDeltaE_M2g_ProtonTag, *dMM2_DeltaE_ProtonTag;
		TH2I *dMM2_M2g_CoplanarTag, *dDeltaE_M2g_CoplanarTag, *dMM2_DeltaE_CoplanarTag;
		TH2I *dDeltaPhi_M2g, *dPhi2g_PhiP;
		TH2I *dEgamma_M2g_ProtonTag;
		TH2I *dMinEgamma_M2g;
};

#endif // _DCustomAction_p2gamma_hists_

