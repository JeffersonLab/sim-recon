// $Id$
//
//    File: DCustomAction_CLcomp.h
// Created: Tue Jan 30 11:16:36 EST 2018
// Creator: aebarnes (on Linux egbert 2.6.32-696.13.2.el6.x86_64 x86_64)
//

#ifndef _DCustomAction_CLcomp_
#define _DCustomAction_CLcomp_

#include <string>
#include <iostream>

#include "JANA/JEventLoop.h"
#include "JANA/JApplication.h"

#include "ANALYSIS/DAnalysisAction.h"
#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DAnalysisUtilities.h"

#include "KINFITTER/DKinFitter.h"
#include "ANALYSIS/DKinFitUtils_GlueX.h"
#include "PID/DKinematicData.h"

#include <TH1.h>
#include <TH2.h>

using namespace std;
using namespace jana;

class DCustomAction_CLcomp : public DAnalysisAction
{
	public:

		DCustomAction_CLcomp(const DReaction* locReaction, bool locUseKinFitResultsFlag, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Custom_CLcomp", locUseKinFitResultsFlag, locActionUniqueString) {}

		void Initialize(JEventLoop* locEventLoop);
		void Reset_NewEvent(void){}; //RESET HISTOGRAM DUPLICATE-CHECK TRACKING HERE!!
	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		DKinFitter* dKinFitter;
                DKinFitUtils_GlueX* dKinFitUtils;
                DAnalysisUtilities* dAnalysisUtils;

		//Store any histograms as member variables here
		TH1I *dHist_IM;
		TH1I *dHist_IM_kept;
		TH1I *dHist_IM_cut;
		TH1I *dHist_CL_KK;
		TH1I *dHist_CL_KK_kept;
		TH1I *dHist_CL_KK_cut;
		TH1I *dHist_CL_pipi;
		TH1I *dHist_CL_pipi_kept;
		TH1I *dHist_CL_pipi_cut;
		TH2I *dHist_log10_ratio_vs_E;
		TH2I *dHist_log10_ratio_vs_mass;
		TH2I *dHist_mass_vs_E;
		TH2I *dHist_mass_vs_E_kept;
		TH2I *dHist_mass_vs_E_cut;
};

#endif // _DCustomAction_CLcomp_

