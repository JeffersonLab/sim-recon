// $Id$
// DCustomAction_Z2pi_trees.h. Modeled after DCustomAction_p2pi_hists.h
//
//    File: DCustomAction_p2pi_hists.h
// Created: Wed Jan 21 16:53:41 EST 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _DCustomAction_Z2pi_trees_
#define _DCustomAction_Z2pi_trees_

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

class DCustomAction_Z2pi_trees : public DAnalysisAction
{
	public:

		DCustomAction_Z2pi_trees(const DReaction* locReaction, bool locUseKinFitResultsFlag, string locActionUniqueString = "") : 
	        DAnalysisAction(locReaction, "Custom_Z2pi_trees", locUseKinFitResultsFlag, locActionUniqueString) {}

		void Initialize(JEventLoop* locEventLoop);

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		// Parameters for event selection to fill histograms
		int endpoint_energy_bins;
		double cohmin_energy, cohedge_energy, endpoint_energy;
		double dEdxCut, minMMCut, maxMMCut, minMM2Cut, maxMM2Cut, missingEnergyCut, min2piMassCut, max2piMassCut;
		double Mtarg;

		// Optional: Useful utility functions.
		const DAnalysisUtilities* dAnalysisUtilities;

		// need PID algos for SC matching
                const DParticleID* dParticleID;

		Particle_t dMissingPID;

		//Store any histograms as member variables here
		TH1I *dEgamma;
		TH2I *dMM_M2pi;
		TH2I *dMM2_M2pi, *dDeltaE_M2pi;
		TH2I *dEgamma_M2pi;
		TH2I *dPiPlusPsi_t;
		TH2I *dPiPlusPsi_Egamma, *dPiPlusPhi_Egamma;
};

#endif // _DCustomAction_Z2pi_trees_

