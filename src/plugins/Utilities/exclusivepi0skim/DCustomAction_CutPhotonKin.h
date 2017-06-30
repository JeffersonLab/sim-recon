// $Id$
//
//    File: DCustomAction_CutPhotonKin.h
// Created: Fri Jul 18 12:51:03 EDT 2014
// Creator: jrsteven (on Linux ifarm1102 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#ifndef _DCustomAction_CutPhotonKin_
#define _DCustomAction_CutPhotonKin_

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

class DCustomAction_CutPhotonKin : public DAnalysisAction
{
	public:

		DCustomAction_CutPhotonKin(const DReaction* locReaction, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Custom_CutPhotonKin", false, locActionUniqueString) {}

		void Initialize(JEventLoop* locEventLoop);

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		// Optional: Useful utility functions.
		// const DAnalysisUtilities* dAnalysisUtilities;

		//Store any histograms as member variables here
};

#endif // _DCustomAction_CutPhotonKin_

