// $Id$
//
//    File: DCustomAction_Pi0Cuts.h
// Created: Thu Jan 22 11:19:46 EST 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _DCustomAction_Pi0Cuts_
#define _DCustomAction_Pi0Cuts_

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

class DCustomAction_Pi0Cuts : public DAnalysisAction
{
	public:

                DCustomAction_Pi0Cuts(const DReaction* locReaction, bool locUseKinFitResultsFlag, double locMinFCAL, string locActionUniqueString = "") : 
	        DAnalysisAction(locReaction, "Custom_Pi0Cuts", locUseKinFitResultsFlag, locActionUniqueString), dMinFCAL(locMinFCAL){}

		void Initialize(JEventLoop* locEventLoop);

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		//Store any histograms as member variables here
		double dMinFCAL;
		
};

#endif // _DCustomAction_Pi0Cuts_

