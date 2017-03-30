#ifndef _DReaction_
#define _DReaction_

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

#include "JANA/JObject.h"
#include "particleType.h"
#include "ANALYSIS/DReactionStep.h"

using namespace std;
using namespace jana;

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

class DReactionBase : public JObject
{
	public:
		JOBJECT_PUBLIC(DReactionBase);

		// CONSTRUCTOR:
		DReactionBase(string locReactionName); //User must specify a unique reaction name upon construction
		// DESTRUCTOR:
		virtual ~DReactionBase(void);

		// SET OBJECT DATA:
		void Add_ReactionStep(const DReactionStep* locReactionStep){dReactionSteps.push_back(locReactionStep);}
		void Clear_ReactionSteps(void){dReactionSteps.clear();}

		// SET KINFIT CONTROL
		void Set_KinFitType(DKinFitType locKinFitType){dKinFitType = locKinFitType;}

		// GET CONTROL INFO:
		string Get_ReactionName(void) const{return dReactionName;}
		bool Get_IsInclusiveFlag(void) const;

		// GET KINFIT CONTROL
		DKinFitType Get_KinFitType(void) const{return dKinFitType;}

		// GET REACTION STEPS:
		size_t Get_NumReactionSteps(void) const{return dReactionSteps.size();}
		const DReactionStep* Get_ReactionStep(size_t locStepIndex) const{return dReactionSteps.at(locStepIndex);}
		vector<const DReactionStep*> Get_ReactionSteps(void) const{return dReactionSteps;}

		// GET PIDS
		vector<Particle_t> Get_FinalPIDs(bool locIncludeMissingFlag = true, bool locIncludeDecayingFlag = true, Charge_t locCharge = d_AllCharges, bool locIncludeDuplicatesFlag = true) const;

	private:
		// PRIVATE METHODS:
		DReactionBase(void); //make default constructor private. MUST set a name upon construction (and must be unique!)

		// CONTROL MEMBERS:
		string dReactionName; //must be unique
		DKinFitType dKinFitType; //defined in ANALYSIS/DKinFitResults.h

		// REACTION AND ANALYSIS MEMBERS:
		vector<const DReactionStep*> dReactionSteps;
};

class DReaction : public DReactionBase
{
	public:
		JOBJECT_PUBLIC(DReaction);

		// CONSTRUCTOR: //User must specify a unique reaction name upon construction
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
		void Set_MaxPhotonRFDeltaT(double locMaxPhotonRFDeltaT){dMaxPhotonRFDeltaT = pair<bool, double>(true, locMaxPhotonRFDeltaT);}
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

		// GET PIDS
		vector<Particle_t> Get_FinalPIDs(bool locIncludeMissingFlag = true, bool locIncludeDecayingFlag = true, Charge_t locCharge = d_AllCharges, bool locIncludeDuplicatesFlag = true) const;

		// GET ANALYSIS ACTIONS:
		size_t Get_NumAnalysisActions(void) const{return dAnalysisActions.size();}
		DAnalysisAction* Get_AnalysisAction(size_t locIndex) const{return dAnalysisActions.at(locIndex);}

		// GET PRE-DPARTICLECOMBO CUT VALUES //Command-line values will override these values
		pair<bool, double> Get_MaxPhotonRFDeltaT(void) const{return dMaxPhotonRFDeltaT;}
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
		DReaction(void); //make default constructor private. MUST set a name upon construction (and must be unique!)

		// CONTROL MEMBERS:
		string dReactionName; //must be unique
		DKinFitType dKinFitType; //defined in ANALYSIS/DKinFitResults.h
		bool dKinFitUpdateCovarianceMatricesFlag; //true to create new error matrices post-kinfit, false to keep the old ones

		// ROOT TTREE OUTPUT:
		bool dEnableTTreeOutputFlag; //default is false
		bool dSaveUnusedFlag; //default is false
		string dTTreeOutputFileName;

		// REACTION AND ANALYSIS MEMBERS:
		vector<const DReactionStep*> dReactionSteps;
		vector<DAnalysisAction*> dAnalysisActions;

		// PRE-DPARTICLECOMBO CONTROL-CUT VALUES
			//bool = true/false for cut enabled/disabled, double = cut value
			//Command-line values (variable names are below in all-caps) will override these values
		pair<bool, double> dMaxPhotonRFDeltaT; //COMBO:MAX_PHOTON_RF_DELTAT - the maximum photon-rf time difference: used for photon selection
		pair<bool, size_t> dMaxExtraGoodTracks; //COMBO:MAX_EXTRA_GOOD_TRACKS - "good" defined by PreSelect factory
		pair<bool, int> dMaxNumBeamPhotonsInBunch; //COMBO:MAX_NUM_BEAM_PHOTONS cut out combos with more than this # of beam photons surviving the RF delta-t cut

		// EVENT STORE QUERY
		string dEventStoreSkims; // First is skim name (default = "all"), second is additional query (default = "")

		// BUILD ANY FLAGS
		//Default false. If true: Once one is built, don't bother making others. 
		bool dAnyComboFlag;
};

/****************************************************** CONSTRUCTORS AND DESTRUCTORS *******************************************************/

DReaction::DReaction(string locReactionName, vector<const DReactionStep*> locSteps, DKinFitType locKinFitType, string locTreeFileName) :
dReactionName(locReactionName), dKinFitType(locKinFitType), dKinFitUpdateCovarianceMatricesFlag(false), dEnableTTreeOutputFlag(locTreeFileName != ""),
dSaveUnusedFlag(false), dTTreeOutputFileName(locTreeFileName), dReactionSteps(locSteps), dAnalysisActions{}, dMaxPhotonRFDeltaT(false, 0.0),
dMaxExtraGoodTracks(false, 0), dMaxNumBeamPhotonsInBunch(false, 0), dEventStoreSkims(""), dAnyComboFlag(false) {}

DReaction::~DReaction(void)
{
	//DO NOT DELETE REACTION STEPS: MIGHT BE SHARED BETWEEN DIFFERENT DREACTIONS
	for(size_t loc_i = 0; loc_i < dAnalysisActions.size(); ++loc_i)
		delete dAnalysisActions[loc_i];
}


/****************************************************** NAMESPACE-SCOPE NON-INLINE FUNCTION DECLARATIONS *******************************************************/

int Get_DecayStepIndex(const DReaction* locReaction, size_t locStepIndex, size_t locParticleIndex);
pair<int, int> Get_InitialParticleDecayFromIndices(const DReaction* locReaction, int locStepIndex);

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


inline bool DReaction::Get_IsFirstStepBeam(void) const
{
	//impossible for first step to be rescattering: makes no sense: if has target, treat as beam. else treat as decaying & don't care about production mechanism
	return (dReactionSteps[0]->Get_TargetPID() != Unknown);
}

/****************************************************** NAMESPACE-SCOPE INLINE FUNCTIONS: NAMES *******************************************************/

inline string Get_DecayChainFinalParticlesROOTNames(const DReaction* locReaction, Particle_t locInitialPID, bool locKinFitResultsFlag)
{
	//if multiple decay steps have locInitialPID as the parent, only the first listed is used
	return Get_DecayChainFinalParticlesROOTNames(locReaction, locInitialPID, -1, vector<Particle_t>(), locKinFitResultsFlag);
}

inline string Get_DecayChainFinalParticlesROOTNames(const DReaction* locReaction, Particle_t locInitialPID, int locUpToStepIndex, vector<Particle_t> locUpThroughPIDs, bool locKinFitResultsFlag)
{
	//if multiple decay steps have locInitialPID as the parent, only the first listed is used
	auto locReactionSteps = locReaction->Get_ReactionSteps();
	auto locPIDSearcher = [](const DReactionStep* locStep) -> bool{return (locStep->Get_InitialPID() == locInitialPID);};
	auto locStepIterator = std::find_if(locReactionSteps.begin(), locReactionSteps.end(), locPIDSearcher);
	if(locStepIterator == locReactionSteps.end());
		return string("");

	size_t locStepIndex = std::distance(locReactionSteps.begin(), locStepIterator);
	return Get_DecayChainFinalParticlesROOTNames(locReaction, locStepIndex, locUpToStepIndex, locUpThroughPIDs, locKinFitResultsFlag, false);
}

string Get_DecayChainFinalParticlesROOTNames(const DReaction* locReaction, size_t locStepIndex, int locUpToStepIndex, vector<Particle_t> locUpThroughPIDs, bool locKinFitResultsFlag, bool locExpandDecayingParticlesFlag)
{
	//excludes missing!!!
	//if locKinFitResultsFlag = true: don't expand decaying particles (through decay chain) that were included in the kinfit (still expand resonances)
	string locName = "";
	auto locPIDs = locReaction->Get_ReactionStep(locStepIndex)->Get_FinalPIDs();
	int locMissingParticleIndex = dReactionSteps[locStepIndex]->Get_MissingParticleIndex();
	bool locSearchPIDsFlag = !locUpThroughPIDs.empty();
	for(size_t loc_j = 0; loc_j < locPIDs.size(); ++loc_j)
	{
		if(int(loc_j) == locMissingParticleIndex)
			continue; //exclude missing!

		if(locSearchPIDsFlag && (int(locStepIndex) == locUpToStepIndex))
		{
			bool locPIDFoundFlag = false;
			for(vector<Particle_t>::iterator locIterator = locUpThroughPIDs.begin(); locIterator != locUpThroughPIDs.end(); ++locIterator)
			{
				if((*locIterator) != locPIDs[loc_j])
					continue;
				locUpThroughPIDs.erase(locIterator);
				locPIDFoundFlag = true;
				break;
			}
			if(!locPIDFoundFlag)
				continue; //skip it: don't want to include it
		}

		//find all future steps in which this pid is a parent
		int locDecayingStepIndex = -1;
		//expand if: not-kinfitting, or not a fixed mass
		if((!locKinFitResultsFlag) || (!IsFixedMass(locPIDs[loc_j])) || locExpandDecayingParticlesFlag)
			locDecayingStepIndex = Get_DecayStepIndex(locReaction, locStepIndex, loc_j);

		if(locDecayingStepIndex == -1)
			locName += ParticleName_ROOT(locPIDs[loc_j]);
		else
			locName += Get_DecayChainFinalParticlesROOTNames(locDecayingStepIndex, locUpToStepIndex, locUpThroughPIDs, locKinFitResultsFlag, locExpandDecayingParticlesFlag);
	}
	return locName;
}

inline bool DReaction::Check_IfMissingDecayProduct(size_t locStepIndex) const
{
	const DReactionStep* locReactionStep = Get_ReactionStep(locStepIndex);
	Particle_t locMissingPID;
	if(locReactionStep->Get_MissingPID(locMissingPID))
		return true;
	for(size_t loc_j = 0; loc_j < locReactionStep->Get_NumFinalParticleIDs(); ++loc_j)
	{
		int locDecayStepIndex = Get_DecayStepIndex(locStepIndex, loc_j);
		if(locDecayStepIndex <= 0)
			continue; //doesn't decay
		if(Check_IfMissingDecayProduct(locDecayStepIndex))
			return true;
	}
	return false;
}

string DReaction::Get_DetectedParticlesROOTName(void) const
{
	string locDetectedParticlesROOTName;

	Particle_t locPID;
	for(size_t loc_i = 0; loc_i < dReactionSteps.size(); ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < dReactionSteps[loc_i]->Get_NumFinalParticleIDs(); ++loc_j)
		{
			locPID = dReactionSteps[loc_i]->Get_FinalParticleID(loc_j);
			if(dReactionSteps[loc_i]->Get_MissingParticleIndex() == int(loc_j))
				continue; //missing particle

			//see if this pid is a parent in a future step
			int locDecayStepIndex = Get_DecayStepIndex(loc_i, loc_j);
			if(locDecayStepIndex >= 0)
				continue;

			locDetectedParticlesROOTName += ParticleName_ROOT(locPID);
		}
	}
	return locDetectedParticlesROOTName;
}

} //end DAnalysis namespace

#endif // _DReaction_

