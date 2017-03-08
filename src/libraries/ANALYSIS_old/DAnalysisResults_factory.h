// $Id$
//
//    File: DAnalysisResults_factory.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DAnalysisResults_factory_
#define _DAnalysisResults_factory_

#include <map>
#include <deque>
#include <vector>

#include "TH1D.h"
#include "TH2D.h"
#include "TDirectoryFile.h"

#include "JANA/JFactory.h"
#include "DANA/DApplication.h"

#include "TRACKING/DMCThrown.h"

#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DCutActions.h"
#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DAnalysisAction.h"
#include "ANALYSIS/DAnalysisUtilities.h"
#include "ANALYSIS/DAnalysisResults.h"
#include "ANALYSIS/DHistogramActions.h"

using namespace jana;
using namespace std;

class DAnalysisResults_factory : public jana::JFactory<DAnalysisResults>
{
	public:
		DAnalysisResults_factory(){};
		~DAnalysisResults_factory(){};

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(JEventLoop *locEventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *locEventLoop, uint64_t eventnumber);	///< Called every event.

		void Get_Reactions(JEventLoop* locEventLoop, vector<const DReaction*>& locReactions) const;
		void Make_ControlHistograms(vector<const DReaction*>& locReactions);

		unsigned int dDebugLevel;
		DApplication* dApplication;
		double dMinThrownMatchFOM;

		map<const DReaction*, bool> dMCReactionExactMatchFlags;
		map<const DReaction*, DCutAction_TrueCombo*> dTrueComboCuts;

		map<const DReaction*, TH1D*> dHistMap_NumParticleCombos;
		map<const DReaction*, TH1D*> dHistMap_NumEventsSurvivedAction_All;
		map<const DReaction*, TH1D*> dHistMap_NumEventsWhereTrueComboSurvivedAction;
		map<const DReaction*, TH2D*> dHistMap_NumCombosSurvivedAction;
		map<const DReaction*, TH1D*> dHistMap_NumCombosSurvivedAction1D;
};

#endif // _DAnalysisResults_factory_

