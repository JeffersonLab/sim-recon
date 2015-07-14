// $Id$
//
//    File: DCustomAction_CutExtraPi0.h
// Created: Sun Jun 28 15:10:48 EDT 2015
// Creator: pmatt (on Darwin Pauls-MacBook-Pro-2.local 13.4.0 i386)
//

#ifndef _DCustomAction_CutExtraPi0_
#define _DCustomAction_CutExtraPi0_

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

class DCustomAction_CutExtraPi0 : public DAnalysisAction
{
	public:

		DCustomAction_CutExtraPi0(const DReaction* locReaction, double locLowMassCut, double locHighMassCut, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Custom_CutExtraPi0", false, locActionUniqueString),
		dLowMassCut(locLowMassCut), dHighMassCut(locHighMassCut) {}

		void Initialize(JEventLoop* locEventLoop);

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		// Cut away combos with pi0 invariant mass BETWEEN these bounds
		double dLowMassCut;
		double dHighMassCut;

		const DAnalysisUtilities* dAnalysisUtilities;

		//Store any histograms as member variables here
		TH1I* dHist_Pi0InvariantMass;

		//To check for double counting
		set<set<const DNeutralParticleHypothesis*> > dPreviousSourceObjects;
};

#endif // _DCustomAction_CutExtraPi0_

