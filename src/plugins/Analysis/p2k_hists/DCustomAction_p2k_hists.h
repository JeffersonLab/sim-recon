// $Id$
//
//    File: DCustomAction_p2k_hists.h
// Created: Wed Jan 21 16:53:41 EST 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _DCustomAction_p2k_hists_
#define _DCustomAction_p2k_hists_

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

class DCustomAction_p2k_hists : public DAnalysisAction
{
	public:

		DCustomAction_p2k_hists(const DReaction* locReaction, bool locUseKinFitResultsFlag, string locActionUniqueString = "") : 
	        DAnalysisAction(locReaction, "Custom_p2k_hists", locUseKinFitResultsFlag, locActionUniqueString) {}

		void Initialize(JEventLoop* locEventLoop);

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		// Parameters for event selection to fill histograms
		int endpoint_energy_bins;
		double cohmin_energy, cohedge_energy, endpoint_energy;
		double dEdxCut, minMM2Cut, maxMM2Cut, maxPhiMassCut;

		// Optional: Useful utility functions.
		const DAnalysisUtilities* dAnalysisUtilities;

		//Store any histograms as member variables here
		TH1I *dEgamma;
		TH2I *dKplus_deltaInvBeta_P, *dKminus_deltaInvBeta_P;
		TH2I *dM2pi_M2k, *dM2pi_M2k_ProtonTag;
		TH2I *dMM2_M2k, *dProton_dEdx_P, *dProton_P_Theta, *dDeltaE_M2k;

		TH2I *dMM2_M2k_ProtonTag, *dDeltaE_M2k_ProtonTag, *dEgamma_M2k_ProtonTag;
		TH2I *dKplus_deltaInvBeta_P_ProtonTag, *dKminus_deltaInvBeta_P_ProtonTag;
		TH2I *dKplus_deltaInvBeta_P_PhiTag, *dKminus_deltaInvBeta_P_PhiTag;
		TH2I *dKplus_deltaInvBeta_P_RhoTag, *dKminus_deltaInvBeta_P_RhoTag;
		TH2I *dKplus_Beta_P_PhiTag, *dKminus_Beta_P_PhiTag;
		TH2I *dKplus_Beta_P_RhoTag, *dKminus_Beta_P_RhoTag;
		TH2I *dKplus_P_Theta_PhiTag, *dKminus_P_Theta_PhiTag;
		TH2I *dKplus_P_Theta_RhoTag, *dKminus_P_Theta_RhoTag;
};

#endif // _DCustomAction_p2k_hists_

