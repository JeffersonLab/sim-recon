#ifndef _DReaction_
#define _DReaction_

#include <deque>
#include <vector>
#include <string>
#include <iostream>

#include "JANA/JObject.h"
#include "particleType.h"
#include "ANALYSIS/DReactionStep.h"

using namespace std;
using namespace jana;

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

		// CONSTRUCTOR:
		DReaction(string locReactionName); //User must specify a unique reaction name upon construction
		// DESTRUCTOR:
		virtual ~DReaction(void);

		// SET OBJECT DATA:
		void Add_ReactionStep(const DReactionStep* locReactionStep){dReactionSteps.push_back(locReactionStep);}
		void Clear_ReactionSteps(void){dReactionSteps.clear();}
		void Add_AnalysisAction(DAnalysisAction* locAnalysisAction){dAnalysisActions.push_back(locAnalysisAction);}

		// SET KINFIT CONTROL
		void Set_KinFitType(DKinFitType locKinFitType){dKinFitType = locKinFitType;}
		void Set_KinFitUpdateCovarianceMatricesFlag(bool locKinFitUpdateCovarianceMatricesFlag){dKinFitUpdateCovarianceMatricesFlag = locKinFitUpdateCovarianceMatricesFlag;}

		// SET PRE-DPARTICLECOMBO CUT VALUES //Command-line values will override these values
		void Set_MinChargedPIDFOM(double locMinChargedPIDFOM){dMinChargedPIDFOM = pair<bool, double>(true, locMinChargedPIDFOM);}
		void Set_MinPhotonPIDFOM(double locMinPhotonPIDFOM){dMinPhotonPIDFOM = pair<bool, double>(true, locMinPhotonPIDFOM);}
		void Set_MaxPhotonRFDeltaT(double locMaxPhotonRFDeltaT){dMaxPhotonRFDeltaT = pair<bool, double>(true, locMaxPhotonRFDeltaT);}
		void Set_MinProtonMomentum(double locMinProtonMomentum){dMinProtonMomentum = pair<bool, double>(true, locMinProtonMomentum);}
		void Set_MaxExtraGoodTracks(size_t locMaxExtraGoodTracks){dMaxExtraGoodTracks = pair<bool, size_t>(true, locMaxExtraGoodTracks);}
		void Set_MaxNumBeamPhotonsInBunch(size_t locMaxNumBeamPhotonsInBunch){dMaxNumBeamPhotonsInBunch = pair<bool, size_t>(true, locMaxNumBeamPhotonsInBunch);}

		// SET PRE-COMBO-BLUEPRINT INVARIANT MASS CUTS
		void Set_InvariantMassCut(Particle_t locStepInitialPID, double locMinInvariantMass, double locMaxInvariantMass);

		// SET EventStore SKIMS //comma-separated list expected
		void Set_EventStoreSkims(string locEventStoreSkims){dEventStoreSkims = locEventStoreSkims;}

		// ADD COMBO PRE-SELECTION ACTION
		void Add_ComboPreSelectionAction(DAnalysisAction* locAction){dComboPreSelectionActions.push_back(locAction);}

		// GET CONTROL INFO:
		string Get_ReactionName(void) const{return dReactionName;}
		int Get_DecayStepIndex(int locStepIndex, int locParticleIndex) const;
		bool Check_IfMissingDecayProduct(size_t locStepIndex) const;
		pair<int, int> Get_InitialParticleDecayFromIndices(int locStepIndex) const; //1st is step index, 2nd is particle index
		int Get_DefinedParticleStepIndex(void) const; //-1 if none //defined: missing or open-ended-decaying
		bool Get_IsInclusiveChannelFlag(void) const;

		// GET KINFIT CONTROL
		DKinFitType Get_KinFitType(void) const{return dKinFitType;}
		bool Get_KinFitUpdateCovarianceMatricesFlag(void) const{return dKinFitUpdateCovarianceMatricesFlag;}

		// GET REACTION STEPS:
		size_t Get_NumReactionSteps(void) const{return dReactionSteps.size();}
		const DReactionStep* Get_ReactionStep(size_t locStepIndex) const;
		void Get_ReactionSteps(Particle_t locInitialPID, deque<const DReactionStep*>& locReactionSteps) const;

		// GET ANALYSIS ACTIONS:
		size_t Get_NumAnalysisActions(void) const{return dAnalysisActions.size();}
		DAnalysisAction* Get_AnalysisAction(size_t locActionIndex) const;
		DAnalysisAction* Get_LastAnalysisAction(void) const;

		// GET COMBO PRE-SELECTION CUTS:
		size_t Get_NumComboPreSelectionActions(void) const{return dComboPreSelectionActions.size();}
		DAnalysisAction* Get_ComboPreSelectionAction(size_t locActionIndex) const;

		// GET PIDs:
		//locChargeFlag: 0/1/2/3/4 for all, charged, neutral, q+, q- particles
		void Get_DetectedFinalPIDs(deque<Particle_t>& locDetectedPIDs, int locChargeFlag = 0, bool locIncludeDuplicatesFlag = false) const;
		void Get_DetectedFinalPIDs(deque<deque<Particle_t> >& locDetectedPIDs, int locChargeFlag = 0, bool locIncludeDuplicatesFlag = false) const;
		void Get_FinalStatePIDs(deque<Particle_t>& locFinalStatePIDs, bool locIncludeDuplicatesFlag = false) const;
		bool Get_MissingPID(Particle_t& locPID) const; //false if none missing

		// GET PARTICLE NAME STRINGS:
		string Get_DetectedParticlesROOTName(void) const;
		string Get_InitialParticlesROOTName(void) const;
		string Get_DecayChainFinalParticlesROOTNames(Particle_t locInitialPID, bool locKinFitResultsFlag) const;
		string Get_DecayChainFinalParticlesROOTNames(Particle_t locInitialPID, int locUpToStepIndex, deque<Particle_t> locUpThroughPIDs, bool locKinFitResultsFlag) const;
		string Get_DecayChainFinalParticlesROOTNames(size_t locStepIndex, int locUpToStepIndex, deque<Particle_t> locUpThroughPIDs, bool locKinFitResultsFlag, bool locExpandDecayingParticlesFlag) const;

		// GET PRE-DPARTICLECOMBO CUT VALUES //Command-line values will override these values
		pair<bool, double> Get_MinChargedPIDFOM(void) const{return dMinChargedPIDFOM;}
		pair<bool, double> Get_MinPhotonPIDFOM(void) const{return dMinPhotonPIDFOM;}
		pair<bool, double> Get_MaxPhotonRFDeltaT(void) const{return dMaxPhotonRFDeltaT;}
		pair<bool, double> Get_MinProtonMomentum(void) const{return dMinProtonMomentum;}
		pair<bool, size_t> Get_MaxExtraGoodTracks(void) const{return dMaxExtraGoodTracks;}
		pair<bool, size_t> Get_MaxNumBeamPhotonsInBunch(void) const{return dMaxNumBeamPhotonsInBunch;}

		// GET PRE-COMBO-BLUEPRINT MASS CUTS
		bool Get_InvariantMassCut(Particle_t locStepInitialPID, double& locMinInvariantMass, double& locMaxInvariantMass) const;

		// GET EventStore SKIMS //comma-separated list expected
		string Get_EventStoreSkims(void) const{return dEventStoreSkims;}

		// ROOT OUTPUT:
		void Enable_TTreeOutput(string locTTreeOutputFileName, bool locSaveUnusedFlag = false);
		string Get_TTreeOutputFileName(void) const{return dTTreeOutputFileName;}
		bool Get_SaveUnusedFlag(void) const{return dSaveUnusedFlag;}
		bool Get_EnableTTreeOutputFlag(void) const{return dEnableTTreeOutputFlag;}

		// BUILD ANY FLAGS
		//Default false. If true: Once one is built, don't bother making others. 
		bool Get_AnyBlueprintFlag(void) const{return dAnyBlueprintFlag;} //If true, no need to change dAnyComboFlag
		bool Get_AnyComboFlag(void) const{return dAnyComboFlag;}
		void Set_AnyBlueprintFlag(bool locAnyBlueprintFlag){dAnyBlueprintFlag = locAnyBlueprintFlag;} //If true, no need to change dAnyComboFlag
		void Set_AnyComboFlag(bool locAnyComboFlag){dAnyComboFlag = locAnyComboFlag;}

		// OTHER:
		bool Check_AreStepsIdentical(const DReaction* locReaction) const;

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
			//all cuts are disabled by default except dMinProtonMomentum: 300 MeV/c (value used during track reconstruction)
			//note: tracks with no PID information are not cut-by/included-in the PID cuts
		pair<bool, double> dMinChargedPIDFOM; //COMBO:MIN_CHARGED_PID_FOM - the minimum PID FOM for a particle used for this DReaction
		pair<bool, double> dMinPhotonPIDFOM; //COMBO:MIN_PHOTON_PID_FOM - the minimum PID FOM for a neutral particle used for this DReaction
		pair<bool, double> dMaxPhotonRFDeltaT; //COMBO:MAX_PHOTON_RF_DELTAT - the maximum photon-rf time difference: used for photon selection
		pair<bool, double> dMinProtonMomentum; //COMBO:MIN_PROTON_MOMENTUM - when testing whether a non-proton DChargedTrackHypothesis could be a proton, this is the minimum momentum it can have
		pair<bool, size_t> dMaxExtraGoodTracks; //COMBO:MAX_EXTRA_GOOD_TRACKS - "good" defined by PreSelect factory
		pair<bool, int> dMaxNumBeamPhotonsInBunch; //COMBO:MAX_NUM_BEAM_PHOTONS cut out combos with more than this # of beam photons surviving the RF delta-t cut

		// PRE-COMBO-BLUEPRINT MASS CUTS
		map<Particle_t, pair<double, double> > dInvariantMassCuts;

		// COMBO PRE-SELECTION CUTS
		vector<DAnalysisAction*> dComboPreSelectionActions;

		// EVENT STORE QUERY
		string dEventStoreSkims; // First is skim name (default = "all"), second is additional query (default = "")

		// BUILD ANY FLAGS
		//Default false. If true: Once one is built, don't bother making others. 
		bool dAnyBlueprintFlag; //If true, don't need to bother changing dAnyComboFlag
		bool dAnyComboFlag;
};

inline void DReaction::Set_InvariantMassCut(Particle_t locStepInitialPID, double locMinInvariantMass, double locMaxInvariantMass)
{
	dInvariantMassCuts[locStepInitialPID] = pair<double, double>(locMinInvariantMass, locMaxInvariantMass);
}

inline const DReactionStep* DReaction::Get_ReactionStep(size_t locStepIndex) const
{
	if(locStepIndex >= dReactionSteps.size())
		return NULL;
	return dReactionSteps[locStepIndex];
}

inline string DReaction::Get_InitialParticlesROOTName(void) const
{
	if(dReactionSteps.empty())
		return (string());
	return dReactionSteps[0]->Get_InitialParticlesROOTName();
}

inline string DReaction::Get_DecayChainFinalParticlesROOTNames(Particle_t locInitialPID, bool locKinFitResultsFlag) const
{
	//if multiple decay steps have locInitialPID as the parent, only the first listed is used
	return Get_DecayChainFinalParticlesROOTNames(locInitialPID, -1, deque<Particle_t>(), locKinFitResultsFlag);
}

inline string DReaction::Get_DecayChainFinalParticlesROOTNames(Particle_t locInitialPID, int locUpToStepIndex, deque<Particle_t> locUpThroughPIDs, bool locKinFitResultsFlag) const
{
	//if multiple decay steps have locInitialPID as the parent, only the first listed is used
	deque<Particle_t> locPIDs;
	string locName = "";
	for(size_t loc_i = 0; loc_i < dReactionSteps.size(); ++loc_i)
	{
		if(dReactionSteps[loc_i]->Get_InitialParticleID() != locInitialPID)
			continue;
		return Get_DecayChainFinalParticlesROOTNames(loc_i, locUpToStepIndex, locUpThroughPIDs, locKinFitResultsFlag, false);
	}
	return string("");
}

inline void DReaction::Get_ReactionSteps(Particle_t locInitialPID, deque<const DReactionStep*>& locReactionSteps) const
{
	locReactionSteps.clear();
	for(size_t loc_i = 0; loc_i < dReactionSteps.size(); ++loc_i)
	{
		if(dReactionSteps[loc_i]->Get_InitialParticleID() == locInitialPID)
			locReactionSteps.push_back(dReactionSteps[loc_i]);
	}
}

inline DAnalysisAction* DReaction::Get_AnalysisAction(size_t locActionIndex) const
{
	if(locActionIndex >= dAnalysisActions.size())
		return NULL;
	return dAnalysisActions[locActionIndex];
}

inline DAnalysisAction* DReaction::Get_LastAnalysisAction(void) const
{
	if(dAnalysisActions.empty())
		return NULL;
	return dAnalysisActions.back();
}

inline DAnalysisAction* DReaction::Get_ComboPreSelectionAction(size_t locActionIndex) const
{
	if(locActionIndex >= dComboPreSelectionActions.size())
		return NULL;
	return dComboPreSelectionActions[locActionIndex];
}

inline bool DReaction::Get_MissingPID(Particle_t& locPID) const
{
	for(size_t loc_i = 0; loc_i < Get_NumReactionSteps(); ++loc_i)
	{
		if(Get_ReactionStep(loc_i)->Get_MissingPID(locPID))
			return true;
	}
	return false;
}

inline bool DReaction::Check_AreStepsIdentical(const DReaction* locReaction) const
{
	if(locReaction->Get_NumReactionSteps() != dReactionSteps.size())
		return false;
	for(size_t loc_i = 0; loc_i < Get_NumReactionSteps(); ++loc_i)
	{
		if(locReaction->Get_ReactionStep(loc_i) != dReactionSteps[loc_i])
			return false;
	}
	return true;
}

inline void DReaction::Enable_TTreeOutput(string locTTreeOutputFileName, bool locSaveUnusedFlag)
{
	dEnableTTreeOutputFlag = true;
	dSaveUnusedFlag = locSaveUnusedFlag;
	dTTreeOutputFileName = locTTreeOutputFileName;
}

inline bool DReaction::Get_InvariantMassCut(Particle_t locStepInitialPID, double& locMinInvariantMass, double& locMaxInvariantMass) const
{
	map<Particle_t, pair<double, double> >::const_iterator locIterator = dInvariantMassCuts.find(locStepInitialPID);
	if(locIterator == dInvariantMassCuts.end())
		return false;
	locMinInvariantMass = locIterator->second.first;
	locMaxInvariantMass = locIterator->second.second;
	return true;
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

#endif // _DReaction_

