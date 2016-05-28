// $Id$
//
//    File: DCustomAction_CutExtraShowers.h
// Created: Sun Jun 28 15:10:48 EDT 2015
// Creator: pmatt (on Darwin Pauls-MacBook-Pro-2.local 13.4.0 i386)
//

#ifndef _DCustomAction_CutExtraShowers_
#define _DCustomAction_CutExtraShowers_

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

class DCustomAction_CutExtraShowers : public DAnalysisAction
{
	public:

		DCustomAction_CutExtraShowers(const DReaction* locReaction, double locMaxEnergyCut, string locActionUniqueString = "") : 
		DAnalysisAction(locReaction, "Custom_CutExtraShowers", false, locActionUniqueString),
		dMaxEnergyCut(locMaxEnergyCut) {}

		void Initialize(JEventLoop* locEventLoop);

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		double dMaxEnergyCut;

		const DAnalysisUtilities* dAnalysisUtilities;
};

#endif // _DCustomAction_CutExtraShowers_

