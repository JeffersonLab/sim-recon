// $Id$
//
//    File: DCustomAction_dEdxCut.h
// Created: Thu Oct  1 11:18:05 EDT 2015
// Creator: pmatt (on Darwin Pauls-MacBook-Pro-2.local 13.4.0 i386)
//

#ifndef _DCustomAction_dEdxCut_
#define _DCustomAction_dEdxCut_

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

class DCustomAction_dEdxCut : public DAnalysisAction
{
	public:

		DCustomAction_dEdxCut(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Custom_dEdxCut", locActionUniqueString) {}

		void Initialize(JEventLoop* locEventLoop);

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		// Optional: Useful utility functions.
		// const DAnalysisUtilities* dAnalysisUtilities;

		//Store any histograms as member variables here
		TF1* dFunc_dEdxCut_SelectHeavy; //e.g. proton
		TF1* dFunc_dEdxCut_SelectLight; //e.g. pion, kaon
};

#endif // _DCustomAction_dEdxCut_

