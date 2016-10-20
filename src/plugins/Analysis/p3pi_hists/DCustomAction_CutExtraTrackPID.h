// $Id$
//
//    File: DCustomAction_CutExtraTrackPID.h
// Created: Thu Oct  1 21:41:32 EDT 2015
// Creator: pmatt (on Darwin Pauls-MacBook-Pro-2.local 13.4.0 i386)
//

#ifndef _DCustomAction_CutExtraTrackPID_
#define _DCustomAction_CutExtraTrackPID_

#include <string>
#include <iostream>

#include "JANA/JEventLoop.h"
#include "JANA/JApplication.h"

#include "ANALYSIS/DAnalysisAction.h"
#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DAnalysisUtilities.h"

#include "DCustomAction_dEdxCut_p3pi.h"

using namespace std;
using namespace jana;

class DCustomAction_CutExtraTrackPID : public DAnalysisAction
{
	public:

		DCustomAction_CutExtraTrackPID(const DReaction* locReaction, Particle_t locExtraTrackTargetPID, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Custom_CutExtraTrackPID", false, locActionUniqueString), 
		dExtraTrackTargetPID(locExtraTrackTargetPID) {}

		void Initialize(JEventLoop* locEventLoop);

	private:

		Particle_t dExtraTrackTargetPID;
		map<DetectorSystem_t, double> dPIDCuts;
		DCustomAction_dEdxCut_p3pi* ddEdxCutAction;

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		// Optional: Useful utility functions.
		const DAnalysisUtilities* dAnalysisUtilities;

		//Store any histograms as member variables here
};

#endif // _DCustomAction_CutExtraTrackPID_

