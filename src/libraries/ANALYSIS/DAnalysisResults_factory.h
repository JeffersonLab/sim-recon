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
		size_t Get_KinFitParticlePoolSize(void) const{return dKinFitUtils->Get_KinFitParticlePoolSize();};
		size_t Get_KinFitParticlePoolSize_Shared(void) const{return dKinFitUtils->Get_KinFitParticlePoolSize_Shared();};
		size_t Get_KinFitConstraintVertexPoolSize(void) const{return dKinFitUtils->Get_KinFitConstraintVertexPoolSize();};
		size_t Get_KinFitConstraintSpacetimePoolSize(void) const{return dKinFitUtils->Get_KinFitConstraintSpacetimePoolSize();};
		size_t Get_KinFitConstraintP4PoolSize(void) const{return dKinFitUtils->Get_KinFitConstraintP4PoolSize();};
		size_t Get_KinFitConstraintMassPoolSize(void) const{return dKinFitUtils->Get_KinFitConstraintMassPoolSize();};
		size_t Get_KinFitChainPoolSize(void) const{return dKinFitUtils->Get_KinFitChainPoolSize();};
		size_t Get_KinFitChainStepPoolSize(void) const{return dKinFitUtils->Get_KinFitChainStepPoolSize();};
		size_t Get_SymMatrixPoolSize(void) const{return dKinFitUtils->Get_SymMatrixPoolSize();};

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(JEventLoop *locEventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *locEventLoop, uint64_t eventnumber);	///< Called every event.

		void Make_ControlHistograms(vector<const DReaction*>& locReactions, bool locIsMCFlag);
		void Check_ReactionNames(vector<const DReaction*>& locReactions) const;
		const DParticleCombo* Find_TrueCombo(JEventLoop *locEventLoop, const DReaction* locReaction, const vector<const DParticleCombo*>& locCombos);

		bool Execute_Actions(JEventLoop* locEventLoop, const DParticleCombo* locCombo, const DParticleCombo* locTrueCombo, bool locPreKinFitFlag, const vector<DAnalysisAction*>& locActions, size_t& locActionIndex, vector<size_t>& locNumCombosSurvived, int& locLastActionTrueComboSurvives);

		const DParticleCombo* Handle_ComboFit(const DReactionVertexInfo* locReactionVertexInfo, const DParticleCombo* locParticleCombo, const DReaction* locReaction);
		pair<const DKinFitChain*, const DKinFitResults*> Fit_Kinematics(const DReactionVertexInfo* locReactionVertexInfo, const DParticleCombo* locParticleCombo, DKinFitType locKinFitType, bool locUpdateCovMatricesFlag);
		DKinFitResults* Build_KinFitResults(const DParticleCombo* locParticleCombo, DKinFitType locKinFitType, const DKinFitChain* locKinFitChain);

		unsigned int dDebugLevel = 0;
		DApplication* dApplication;
		double dMinThrownMatchFOM;
		DSourceComboer* dSourceComboer;
		DParticleComboCreator* dParticleComboCreator;

		unsigned int dKinFitDebugLevel = 0;
		DKinFitter* dKinFitter;
		DKinFitUtils_GlueX* dKinFitUtils;
		map<pair<set<DKinFitConstraint*>, bool>, DKinFitResults*> dConstraintResultsMap; //used for determining if kinfit results will be identical //bool: update cov matrix flag
		map<tuple<const DParticleCombo*, DKinFitType, bool>, const DParticleCombo*> dPreToPostKinFitComboMap;

		unordered_map<const DReaction*, bool> dMCReactionExactMatchFlags;
		unordered_map<const DReaction*, DCutAction_TrueCombo*> dTrueComboCuts;

		unordered_map<const DReaction*, TH1*> dHistMap_NumParticleCombos;
		unordered_map<const DReaction*, TH1*> dHistMap_NumEventsSurvivedAction_All;
		unordered_map<const DReaction*, TH1*> dHistMap_NumEventsWhereTrueComboSurvivedAction;
		unordered_map<const DReaction*, TH2*> dHistMap_NumCombosSurvivedAction;
		unordered_map<const DReaction*, TH1*> dHistMap_NumCombosSurvivedAction1D;
};

#endif // _DAnalysisResults_factory_

