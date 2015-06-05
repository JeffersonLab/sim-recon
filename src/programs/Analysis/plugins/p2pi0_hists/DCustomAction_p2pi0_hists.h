// $Id$
//
//    File: DCustomAction_p2pi0_hists.h
// Created: Wed Jan 21 16:53:41 EST 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _DCustomAction_p2pi0_hists_
#define _DCustomAction_p2pi0_hists_

#include <string>
#include <iostream>

#include "TH1.h"
#include "TMath.h"
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Plane3D.h"

#include "JANA/JEventLoop.h"
#include "JANA/JApplication.h"

#include "ANALYSIS/DAnalysisAction.h"
#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DAnalysisUtilities.h"

using namespace std;
using namespace jana;
using namespace ROOT::Math;

class DCustomAction_p2pi0_hists : public DAnalysisAction
{
	public:

		DCustomAction_p2pi0_hists(const DReaction* locReaction, bool locUseKinFitResultsFlag, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Custom_p2pi0_hists", locUseKinFitResultsFlag, locActionUniqueString) {}

		void Initialize(JEventLoop* locEventLoop);

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		// Optional: Useful utility functions.
		const DAnalysisUtilities* dAnalysisUtilities;

		//Store any histograms as member variables here
		TH1I *dEgamma;

		TH2I *dMM2_M2pi0, *dDeltaE_M2pi0, *dMpi0_corr;
		TH2I *dDeltaPhi_M2pi0, *dPhi2pi0_PhiP;
                TH2I *dMM2_M2pi0_CoplanarTag, *dDeltaE_M2pi0_CoplanarTag, *dMpi0_corr_CoplanarTag;
		TH2I *dMM2_M2pi0_ProtonTag, *dDeltaE_M2pi0_ProtonTag, *dMpi0_corr_ProtonTag;
		TH2I *dMpi0_corr_MMTag;
		TH2I *dMM2_M2pi0_Pi0Tag, *dDeltaE_M2pi0_Pi0Tag, *dMM2_DeltaE_Pi0Tag, *dEgamma_M2pi0_Pi0Tag;
		TH2I *dDalitz1_Pi0Tag, *dDalitz2_Pi0Tag, *dDalitz3_Pi0Tag, *dDalitz4_Pi0Tag,;

		TH2I *dProton_dEdx_P, *dProton_P_Theta;

};

#endif // _DCustomAction_p2pi0_hists_

