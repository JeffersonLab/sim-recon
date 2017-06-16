// $Id$
//
//    File: DCustomAction_ee_convert_cuts.h
// Created: Wed Jun 14 06:26:48 EDT 2017
// Creator: jrsteven (on Linux ifarm1402.jlab.org 3.10.0-327.el7.x86_64 x86_64)
//

#ifndef _DCustomAction_ee_convert_cuts_
#define _DCustomAction_ee_convert_cuts_

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

class DCustomAction_ee_convert_cuts : public DAnalysisAction
{
	public:

		DCustomAction_ee_convert_cuts(const DReaction* locReaction, bool locUseKinFitResultsFlag, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Custom_ee_convert_cuts", locUseKinFitResultsFlag, locActionUniqueString) {}

		void Initialize(JEventLoop* locEventLoop);
		void Reset_NewEvent(void){}; //RESET HISTOGRAM DUPLICATE-CHECK TRACKING HERE!!
	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		//Store any histograms as member variables here
		TH2I *dHist_Vertex, *dHist_Vertex_DOCA, *dHist_DOCAVertex_DOCA, *dHist_Vertex_LowMass, *dHist_Vertex_Final;
		TH2I *dHist_DOCA_DOCAVertex, *dHist_VertexZ_DOCAVertexZ, *dHist_VertexR_DOCAVertexR;
		TH2I *dHist_DOCA_DeltaPhi, *dHist_DOCA_DeltaTheta, *dHist_DeltaPhi_DeltaTheta, *dHist_DOCA_PhiV;
		TH1I *dHist_Mee;
		TH2I *dHist_EOverP, *dHist_EOverP_FCALBCAL, *dHist_EOverP_BothBCAL, *dHist_EOverP_BothFCAL;

		set<set<pair<const JObject*, unsigned int> > > dPreviousSourceObjects;
};

#endif // _DCustomAction_ee_convert_cuts_

