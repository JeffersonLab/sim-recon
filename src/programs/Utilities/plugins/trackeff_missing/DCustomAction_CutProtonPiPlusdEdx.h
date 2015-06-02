// $Id$
//
//    File: DCustomAction_CutProtonPiPlusdEdx.h
// Created: Tue Jun  2 18:25:06 EDT 2015
// Creator: pmatt (on Linux pmattdesktop.jlab.org 2.6.32-504.16.2.el6.x86_64 x86_64)
//

#ifndef _DCustomAction_CutProtonPiPlusdEdx_
#define _DCustomAction_CutProtonPiPlusdEdx_

#include <string>
#include <iostream>

#include "JANA/JEventLoop.h"
#include "JANA/JApplication.h"

#include "ANALYSIS/DAnalysisAction.h"
#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DParticleCombo.h"

using namespace std;
using namespace jana;

class DCustomAction_CutProtonPiPlusdEdx : public DAnalysisAction
{
	public:

		DCustomAction_CutProtonPiPlusdEdx(const DReaction* locReaction, double locTrackdEdxCut_InKeV, bool locCutProtonsInOverlapRegionFlag = false, string locActionUniqueString = "") :
		DAnalysisAction(locReaction, "Custom_CutdEdx", false, locActionUniqueString),
		dTrackdEdxCut_InKeV(locTrackdEdxCut_InKeV), dCutProtonsInOverlapRegionFlag(locCutProtonsInOverlapRegionFlag), dOverlapRegionMinP(1.0) {}

		void Initialize(JEventLoop* locEventLoop){};
		string Get_ActionName(void) const;

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		double dTrackdEdxCut_InKeV;
		bool dCutProtonsInOverlapRegionFlag;
		double dOverlapRegionMinP;
};

#endif // _DCustomAction_CutProtonPiPlusdEdx_
