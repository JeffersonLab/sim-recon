#ifndef _DReaction_
#define _DReaction_

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

#include "JANA/JObject.h"
#include "JANA/JEventLoop.h"
#include "JANA/JFactory.h"
#include "particleType.h"
#include "ANALYSIS/DReactionStep.h"

using namespace std;
using namespace jana;
using namespace DAnalysis;

namespace DAnalysis
{

class DAnalysisAction;

enum DKinFitType
{
	d_NoFit = 0, 
	d_P4Fit, //also includes invariant mass constraints
	d_VertexFit,
	d_SpacetimeFit,
	d_P4AndVertexFit, //also includes invariant mass constraints
	d_P4AndSpacetimeFit //also includes invariant mass constraints
};

class DReaction : public JObject
{
	public:
		JOBJECT_PUBLIC(DReaction);

		// CONSTRUCTOR: //User must specify a unique reaction name upon construction
		DReaction(void) = delete;
		DReaction(string locReactionName, vector<const DReactionStep*> locSteps = {}, DKinFitType locKinFitType = d_NoFit, string locTreeFileName = "");

		// DESTRUCTOR:
		virtual ~DReaction(void);

		// SET OBJECT DATA:
		void Add_ReactionStep(const DReactionStep* locReactionStep){dReactionSteps.push_back(locReactionStep);}
		void Clear_ReactionSteps(void){dReactionSteps.clear();}
		void Add_AnalysisAction(DAnalysisAction* locAnalysisAction){dAnalysisActions.push_back(locAnalysisAction);}

		// SET KINFIT CONTROL
		void Set_KinFitType(DKinFitType locKinFitType){dKinFitType = locKinFitType;}
		void Set_KinFitUpdateCovarianceMatricesFlag(bool locUpdateFlag){dKinFitUpdateCovarianceMatricesFlag = locUpdateFlag;}

		// SET PRE-DPARTICLECOMBO CUT VALUES //Command-line values will override these values
		void Set_NumPlusMinusRFBunches(size_t locNumPlusMinusRFBunches){dNumPlusMinusRFBunches = locNumPlusMinusRFBunches;}
		void Set_MaxExtraGoodTracks(size_t locMaxExtraGoodTracks){dMaxExtraGoodTracks = pair<bool, size_t>(true, locMaxExtraGoodTracks);}
		void Set_MaxNumBeamPhotonsInBunch(size_t locMaxNumBeamPhotonsInBunch){dMaxNumBeamPhotonsInBunch = pair<bool, size_t>(true, locMaxNumBeamPhotonsInBunch);}

		// SET EventStore SKIMS //comma-separated list expected
		void Set_EventStoreSkims(string locEventStoreSkims){dEventStoreSkims = locEventStoreSkims;}

		// GET CONTROL INFO:
		string Get_ReactionName(void) const{return dReactionName;}
		bool Get_IsInclusiveFlag(void) const;

		// GET KINFIT CONTROL
		DKinFitType Get_KinFitType(void) const{return dKinFitType;}
		bool Get_KinFitUpdateCovarianceMatricesFlag(void) const{return dKinFitUpdateCovarianceMatricesFlag;}

		// GET REACTION STEPS:
		size_t Get_NumReactionSteps(void) const{return dReactionSteps.size();}
		const DReactionStep* Get_ReactionStep(size_t locStepIndex) const{return dReactionSteps.at(locStepIndex);}
		vector<const DReactionStep*> Get_ReactionSteps(void) const{return dReactionSteps;}

		// GET PIDS //step index of -1: All steps
		vector<Particle_t> Get_FinalPIDs(int locStepIndex = -1, bool locIncludeMissingFlag = true, bool locIncludeDecayingFlag = true,
				Charge_t locCharge = d_AllCharges, bool locIncludeDuplicatesFlag = true) const;
		vector<Particle_t> Get_MissingPIDs(int locStepIndex = -1, Charge_t locCharge = d_AllCharges, bool locIncludeDuplicatesFlag = true) const;

		// GET ANALYSIS ACTIONS:
		size_t Get_NumAnalysisActions(void) const{return dAnalysisActions.size();}
		vector<DAnalysisAction*> Get_AnalysisActions(void) const{return dAnalysisActions;}

		// GET PRE-DPARTICLECOMBO CUT VALUES //Command-line values will override these values
		size_t Get_NumPlusMinusRFBunches(void) const{return dNumPlusMinusRFBunches;}
		pair<bool, size_t> Get_MaxExtraGoodTracks(void) const{return dMaxExtraGoodTracks;}
		pair<bool, size_t> Get_MaxNumBeamPhotonsInBunch(void) const{return dMaxNumBeamPhotonsInBunch;}

		// GET EventStore SKIMS //comma-separated list expected
		string Get_EventStoreSkims(void) const{return dEventStoreSkims;}

		// ROOT OUTPUT:
		void Enable_TTreeOutput(string locTTreeOutputFileName, bool locSaveUnusedFlag = false);
		string Get_TTreeOutputFileName(void) const{return dTTreeOutputFileName;}
		bool Get_SaveUnusedFlag(void) const{return dSaveUnusedFlag;}
		bool Get_EnableTTreeOutputFlag(void) const{return dEnableTTreeOutputFlag;}

		// BUILD ANY FLAGS
		//Default false. If true: Once one is built, don't bother making others. 
		bool Get_AnyComboFlag(void) const{return dAnyComboFlag;}
		void Set_AnyComboFlag(bool locAnyComboFlag){dAnyComboFlag = locAnyComboFlag;}

	private:
		// PRIVATE METHODS:

		// CONTROL MEMBERS:
		string dReactionName; //must be unique
		DKinFitType dKinFitType = d_NoFit; //defined in ANALYSIS/DKinFitResults.h
		bool dKinFitUpdateCovarianceMatricesFlag = false; //true to create new error matrices post-kinfit, false to keep the old ones

		// ROOT TTREE OUTPUT:
		bool dEnableTTreeOutputFlag = false;
		bool dSaveUnusedFlag = false;
		string dTTreeOutputFileName = "";

		// REACTION AND ANALYSIS MEMBERS:
		vector<const DReactionStep*> dReactionSteps;
		vector<DAnalysisAction*> dAnalysisActions;

		// PRE-DPARTICLECOMBO CONTROL-CUT VALUES
			//bool = true/false for cut enabled/disabled, double = cut value
			//Command-line values (variable names are below in all-caps) will override these values
		size_t dNumPlusMinusRFBunches = 99999; //COMBO:NUM_PLUSMINUS_RF_BUNCHES //e.g. if 0 then only center bunch, if 2 then +/- 2 bunches
		pair<bool, size_t> dMaxExtraGoodTracks = make_pair(false, size_t(0)); //COMBO:MAX_EXTRA_GOOD_TRACKS - "good" defined by PreSelect factory
		pair<bool, int> dMaxNumBeamPhotonsInBunch = make_pair(false, 0); //COMBO:MAX_NUM_BEAM_PHOTONS cut out combos with more than this # of beam photons surviving the RF delta-t cut

		// EVENT STORE QUERY
		string dEventStoreSkims = ""; // First is skim name (default = "all"), second is additional query (default = "")

		// BUILD ANY FLAGS
		//Default false. If true: Once one is built, don't bother making others. 
		bool dAnyComboFlag = false;
};

/****************************************************** NAMESPACE-SCOPE NON-INLINE FUNCTION DECLARATIONS *******************************************************/

int Get_DecayStepIndex(const DReaction* locReaction, size_t locStepIndex, size_t locParticleIndex);
pair<int, int> Get_InitialParticleDecayFromIndices(const DReaction* locReaction, int locStepIndex);
vector<Particle_t> Get_ChainPIDs(const DReaction* locReaction, size_t locStepIndex, int locUpToStepIndex, vector<Particle_t> locUpThroughPIDs, bool locExpandDecayingFlag);
int Get_DefinedParticleStepIndex(const DReaction* locReaction);
vector<const DReaction*> Get_Reactions(JEventLoop* locEventLoop);
size_t Get_ParticleInstanceIndex(const DReactionStep* locStep, size_t locParticleIndex);

/****************************************************** CONSTRUCTORS AND DESTRUCTORS *******************************************************/

inline DReaction::DReaction(string locReactionName, vector<const DReactionStep*> locSteps, DKinFitType locKinFitType, string locTreeFileName) :
dReactionName(locReactionName), dKinFitType(locKinFitType), dKinFitUpdateCovarianceMatricesFlag(false), dEnableTTreeOutputFlag(locTreeFileName != ""),
dSaveUnusedFlag(false), dTTreeOutputFileName(locTreeFileName), dReactionSteps(locSteps) {}

/************************************************************** DREACTION INLINE FUNCTIONS ***************************************************************/

inline bool DReaction::Get_IsInclusiveFlag(void) const
{
	auto locInclusiveSearcher = [](const DReactionStep* locStep) -> bool{return locStep->Get_IsInclusiveFlag();};
	return (std::find_if(dReactionSteps.begin(), dReactionSteps.end(), locInclusiveSearcher) != dReactionSteps.end());
}

inline void DReaction::Enable_TTreeOutput(string locTTreeOutputFileName, bool locSaveUnusedFlag)
{
	dEnableTTreeOutputFlag = true;
	dSaveUnusedFlag = locSaveUnusedFlag;
	dTTreeOutputFileName = locTTreeOutputFileName;
}

/****************************************************** NAMESPACE-SCOPE INLINE FUNCTIONS: MISC *******************************************************/

inline bool Get_IsFirstStepBeam(const DReaction* locReaction)
{
	//impossible for first step to be rescattering: makes no sense: if has target, treat as beam. else treat as decaying & don't care about production mechanism
	auto locFirstStep = locReaction->Get_ReactionStep(0);
	return ((locFirstStep->Get_TargetPID() != Unknown) || (locFirstStep->Get_SecondBeamPID() != Unknown));
}

inline bool Check_ChannelEquality(const DReaction* lhs, const DReaction* rhs, bool locSameOrderFlag = true, bool locRightSubsetOfLeftFlag = false)
{
	//assume for now that the steps have to be in the same order
	auto locSteps_lhs = lhs->Get_ReactionSteps();
	auto locSteps_rhs = rhs->Get_ReactionSteps();
	if(locSteps_lhs.size() != locSteps_rhs.size())
		return false;

	auto Equality_Checker = [&locSameOrderFlag, &locRightSubsetOfLeftFlag](const DReactionStep* lhs, const DReactionStep* rhs) -> bool
		{return DAnalysis::Check_ChannelEquality(lhs, rhs, locSameOrderFlag, locRightSubsetOfLeftFlag);};
	return std::equal(locSteps_lhs.begin(), locSteps_lhs.end(), locSteps_rhs.begin(), Equality_Checker);
}

/****************************************************** NAMESPACE-SCOPE INLINE FUNCTIONS: PIDS *******************************************************/

inline vector<Particle_t> Get_ChainPIDs(const DReaction* locReaction, Particle_t locInitialPID, bool locExpandDecayingFlag)
{
	//if multiple decay steps have locInitialPID as the parent, only the first listed is used
	return Get_ChainPIDs(locReaction, locInitialPID, -1, vector<Particle_t>(), locExpandDecayingFlag);
}

inline vector<Particle_t> Get_ChainPIDs(const DReaction* locReaction, Particle_t locInitialPID, int locUpToStepIndex, vector<Particle_t> locUpThroughPIDs, bool locExpandDecayingFlag)
{
	//if multiple decay steps have locInitialPID as the parent, only the first listed is used
	auto locReactionSteps = locReaction->Get_ReactionSteps();
	auto locPIDSearcher = [&locInitialPID](const DReactionStep* locStep) -> bool{return (locStep->Get_InitialPID() == locInitialPID);};
	auto locStepIterator = std::find_if(locReactionSteps.begin(), locReactionSteps.end(), locPIDSearcher);
	if(locStepIterator == locReactionSteps.end())
		return {};

	size_t locStepIndex = std::distance(locReactionSteps.begin(), locStepIterator);
	return Get_ChainPIDs(locReaction, locStepIndex, locUpToStepIndex, locUpThroughPIDs, locExpandDecayingFlag);
}

inline string Convert_PIDsToROOTName(const vector<Particle_t>& locPIDs)
{
	auto locPIDTransformer = [](Particle_t locPID) -> string {return ParticleName_ROOT(locPID);};
	vector<string> locParticleNames;
	std::transform(locPIDs.begin(), locPIDs.end(), std::back_inserter(locParticleNames), locPIDTransformer);
	return std::accumulate(locParticleNames.begin(), locParticleNames.end(), string(""));
}

inline bool Check_IfMissingDecayProduct(const DReaction* locReaction, size_t locStepIndex)
{
	const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(locStepIndex);
	if(locReactionStep->Get_IsInclusiveFlag())
		return true;
	if(locReactionStep->Get_MissingPID() != Unknown)
		return true;

	for(size_t loc_i = 0; loc_i < locReactionStep->Get_NumFinalPIDs(); ++loc_i)
	{
		int locDecayStepIndex = Get_DecayStepIndex(locReaction, locStepIndex, loc_i);
		if((locDecayStepIndex > 0) && Check_IfMissingDecayProduct(locReaction, locDecayStepIndex))
			return true;
	}
	return false;
}

template <typename DType> inline map<DType, size_t> Convert_VectorToCountMap(const vector<DType>& locVector)
{
	map<DType, size_t> locMap;
	for(const auto& locObject : locVector)
	{
		auto locIterator = locMap.find(locObject);
		if(locIterator != locMap.end())
			++(locIterator->second);
		else
			locMap.emplace(locObject, 1);
	}
	return locMap;
}

} //end DAnalysis namespace

#endif // _DReaction_

