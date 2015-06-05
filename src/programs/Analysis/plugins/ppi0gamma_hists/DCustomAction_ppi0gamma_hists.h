// $Id$
//
//    File: DCustomAction_ppi0gamma_hists.h
// Created: Wed Jan 21 16:53:41 EST 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _DCustomAction_ppi0gamma_hists_
#define _DCustomAction_ppi0gamma_hists_

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

class DCustomAction_ppi0gamma_hists : public DAnalysisAction
{
	public:

		DCustomAction_ppi0gamma_hists(const DReaction* locReaction, bool locUseKinFitResultsFlag, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Custom_ppi0gamma_hists", locUseKinFitResultsFlag, locActionUniqueString) {}

		void Initialize(JEventLoop* locEventLoop);

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		// Optional: Useful utility functions.
		const DAnalysisUtilities* dAnalysisUtilities;

		//Store any histograms as member variables here
		TH1I *dEgamma;

		TH2I *dMM2_MPi0, *dMM2_MOmega;
		TH2I *dDeltaPhi_MOmega, *dPhiOmega_PhiP;
		TH2I *dMM2_MOmegaCoplanarTag, *dDeltaE_MOmegaCoplanarTag, *dMM2_DeltaE_CoplanarTag;
		TH2I *dMM2_MOmegaProtonTag, *dDeltaE_MOmegaProtonTag, *dMM2_DeltaE_ProtonTag, *dEgamma_MOmegaProtonTag, *dEgamma_MOmegaKinTag;
		TH2I *dProton_dEdx_P, *dProton_P_Theta;
		TH2I *dDeltaE_MOmega;
		
};

#endif // _DCustomAction_ppi0gamma_hists_

