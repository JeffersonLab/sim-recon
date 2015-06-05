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

		//Store any histograms as member variables here
		TH1I *dEgamma;
		TH2I *dMM2_M3pi, *dProton_dEdx_P, *dProton_P_Theta;
		TH2I *dDeltaE_M3pi;
		TH2I *dMM2_M3pi_ProtonTag, *dDeltaE_M3pi_ProtonTag, *dMM2_DeltaE_ProtonTag;
                TH2I *dMM2_M3pi_CoplanarTag, *dDeltaE_M3pi_CoplanarTag, *dMM2_DeltaE_CoplanarTag;
                TH2I *dDeltaPhi_M3pi, *dPhi3pi_PhiP;
		TH2I *dEgamma_M3pi_ProtonTag;
		TH2I *dMppizero_M2pi, *dMppiplus_M2pi, *dMppiminus_M2pi;

};

#endif // _DCustomAction_p3pi_hists_

