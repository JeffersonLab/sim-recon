// $Id$
//
//    File: DCustomAction_TrackingEfficiency.h
// Created: Wed Feb 25 09:38:06 EST 2015
// Creator: pmatt (on Linux pmattdesktop.jlab.org 2.6.32-504.8.1.el6.x86_64 x86_64)
//

#ifndef _DCustomAction_TrackingEfficiency_
#define _DCustomAction_TrackingEfficiency_

#include <string>
#include <iostream>

#include "JANA/JEventLoop.h"
#include "JANA/JApplication.h"

#include "BCAL/DBCALShower.h"
#include "TRACKING/DTrackTimeBased.h"
#include "PID/DParticleID.h"
#include "PID/DDetectorMatches.h"

#include "ANALYSIS/DAnalysisAction.h"
#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DAnalysisUtilities.h"
#include "ANALYSIS/DTreeInterface.h"

using namespace std;
using namespace jana;

class DCustomAction_TrackingEfficiency : public DAnalysisAction
{
	public:

		DCustomAction_TrackingEfficiency(const DReaction* locReaction, bool locUseKinFitResultsFlag, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Custom_TrackingEfficiency", locUseKinFitResultsFlag, locActionUniqueString)
		{
			dTreeInterface = NULL;
		}

		void Initialize(JEventLoop* locEventLoop);

		~DCustomAction_TrackingEfficiency(void)
		{
			if(dTreeInterface != NULL)
				delete dTreeInterface; //SAVES FILES, TREES
		}

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);
		double Calc_MatchFOM(const DVector3& locDeltaP3, TMatrixFSym locInverse3x3Matrix) const;

		Particle_t dMissingPID;

		// Optional: Useful utility functions.
		const DAnalysisUtilities* dAnalysisUtilities;
		const DParticleID* dParticleID;

		DTreeFillData dTreeFillData;
		DTreeInterface* dTreeInterface;
};

#endif // _DCustomAction_TrackingEfficiency_

