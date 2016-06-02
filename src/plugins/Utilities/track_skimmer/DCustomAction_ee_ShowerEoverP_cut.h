// $Id$
//
//    File: DCustomAction_ee_ShowerEoverP_cut.h
// Created: Mon Mar  9 18:27:49 EDT 2015
// Creator: sdobbs (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _DCustomAction_ee_ShowerEoverP_cut_
#define _DCustomAction_ee_ShowerEoverP_cut_

#include <string>
#include <iostream>

#include "JANA/JEventLoop.h"
#include "JANA/JApplication.h"

#include "ANALYSIS/DAnalysisAction.h"
#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DAnalysisUtilities.h"

#include "TH1D.h"
#include "TH2D.h"

using namespace std;
using namespace jana;

class DCustomAction_ee_ShowerEoverP_cut : public DAnalysisAction
{
	public:

        DCustomAction_ee_ShowerEoverP_cut(const DReaction* locReaction, bool locUseKinFitResultsFlag, 
                                          double in_BCAL_EP_min=0., double in_BCAL_EP_max=100., 
                                          double in_FCAL_EP_min=0., double in_FCAL_EP_max=100., string locActionUniqueString = ""):
    DAnalysisAction(locReaction, "CustomAction_ShowerEoverP_cut", locUseKinFitResultsFlag, locActionUniqueString),
                dBCAL_EP_min(in_BCAL_EP_min),dBCAL_EP_max(in_BCAL_EP_max),dFCAL_EP_min(in_FCAL_EP_min),dFCAL_EP_max(in_FCAL_EP_max) {}


		void Initialize(JEventLoop* locEventLoop);

	private:

		bool Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo);

		double dBCAL_EP_min, dBCAL_EP_max;
		double dFCAL_EP_min, dFCAL_EP_max;
};

#endif // _DCustomAction_ee_ShowerEoverP_cut_

