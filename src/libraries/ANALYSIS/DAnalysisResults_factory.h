// $Id$
//
//    File: DAnalysisResults_factory.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DAnalysisResults_factory_
#define _DAnalysisResults_factory_

#include <unordered_map>
#include <map>
#include <set>
#include <vector>

#include "TH1D.h"
#include "TH2D.h"
#include "TDirectoryFile.h"
#include "TROOT.h"

#include "JANA/JFactory.h"
#include "DANA/DApplication.h"

#include "TRACKING/DMCThrown.h"
#include "TRIGGER/DTrigger.h"

#include "KINFITTER/DKinFitter.h"
#include "ANALYSIS/DKinFitResults.h"
#include "ANALYSIS/DKinFitUtils_GlueX.h"

#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DReactionVertexInfo.h"
#include "ANALYSIS/DCutActions.h"
#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DAnalysisAction.h"
#include "ANALYSIS/DAnalysisUtilities.h"
#include "ANALYSIS/DAnalysisResults.h"
#include "ANALYSIS/DHistogramActions.h"
#include "ANALYSIS/DSourceComboer.h"
#include "ANALYSIS/DParticleComboCreator.h"

using namespace jana;
using namespace std;

class DAnalysisResults_factory : public jana::JFactory<DAnalysisResults>
{
	public:
		~DAnalysisResults_factory(void){delete dSourceComboer;}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(JEventLoop *locEventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *locEventLoop, uint64_t eventnumber);	///< Called every event.

		void Make_ControlHistograms(vector<const DReaction*>& locReactions);
		void Check_ReactionNames(vector<const DReaction*>& locReactions) const;
		const DParticleCombo* Find_TrueCombo(JEventLoop *locEventLoop, const DReaction* locReaction, const vector<const DParticleCombo*>& locCombos);

		bool Execute_Actions(JEventLoop* locEventLoop, bool locIsKinFit, const DParticleCombo* locCombo, const DParticleCombo* locTrueCombo, bool locPreKinFitFlag, const vector<DAnalysisAction*>& locActions, size_t& locActionIndex, vector<size_t>& locNumCombosSurvived, int& locLastActionTrueComboSurvives);

		const DParticleCombo* Handle_ComboFit(const DReactionVertexInfo* locReactionVertexInfo, const DParticleCombo* locParticleCombo, const DReaction* locReaction);
		pair<shared_ptr<const DKinFitChain>, const DKinFitResults*> Fit_Kinematics(const DReactionVertexInfo* locReactionVertexInfo, const DReaction* locReaction, const DParticleCombo* locParticleCombo, DKinFitType locKinFitType, bool locUpdateCovMatricesFlag);
		DKinFitResults* Build_KinFitResults(const DParticleCombo* locParticleCombo, DKinFitType locKinFitType, const shared_ptr<const DKinFitChain>& locKinFitChain);

		unsigned int dDebugLevel = 0;
		DApplication* dApplication;
		double dMinThrownMatchFOM;
		DSourceComboer* dSourceComboer;
		DParticleComboCreator* dParticleComboCreator;
		bool dIsMCFlag = false;

		bool dRequireKinFitConvergence = true;
		unsigned int dKinFitDebugLevel = 0;
		DKinFitter* dKinFitter;
		DKinFitUtils_GlueX* dKinFitUtils;
		map<pair<set<shared_ptr<DKinFitConstraint>>, bool>, DKinFitResults*> dConstraintResultsMap; //used for determining if kinfit results will be identical //bool: update cov matrix flag
		map<tuple<const DParticleCombo*, DKinFitType, bool>, const DParticleCombo*> dPreToPostKinFitComboMap;

		DResourcePool<DKinFitResults> dResourcePool_KinFitResults;
		vector<DKinFitResults*> dCreatedKinFitResults;
		DKinFitResults* Get_KinFitResultsResource(void)
		{
			auto locKinFitResults = dResourcePool_KinFitResults.Get_Resource();
			locKinFitResults->Reset();
			dCreatedKinFitResults.push_back(locKinFitResults);
			return locKinFitResults;
		}

		unordered_map<const DReaction*, bool> dMCReactionExactMatchFlags;
		unordered_map<const DReaction*, DCutAction_TrueCombo*> dTrueComboCuts;

		unordered_map<const DReaction*, TH1*> dHistMap_NumParticleCombos;
		unordered_map<const DReaction*, TH1*> dHistMap_NumEventsSurvivedAction_All;
		unordered_map<const DReaction*, TH1*> dHistMap_NumEventsWhereTrueComboSurvivedAction;
		unordered_map<const DReaction*, TH2*> dHistMap_NumCombosSurvivedAction;
		unordered_map<const DReaction*, TH1*> dHistMap_NumCombosSurvivedAction1D;
};

#endif // _DAnalysisResults_factory_

