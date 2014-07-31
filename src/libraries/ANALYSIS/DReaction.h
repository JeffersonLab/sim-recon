#ifndef _DReaction_
#define _DReaction_

#include <deque>
#include <string>
#include <iostream>

#include "JANA/JObject.h"
#include "particleType.h"
#include "ANALYSIS/DReactionStep.h"
#include "ANALYSIS/DKinFitResults.h"

using namespace std;
using namespace jana;

class DAnalysisAction;

class DReaction : public JObject
{
	public:
		JOBJECT_PUBLIC(DReaction);

		// CONSTRUCTOR:
		DReaction(string locReactionName); //User must specify a unique reaction name upon construction

		// SET OBJECT DATA:
		inline void Set_KinFitType(DKinFitType locKinFitType){dKinFitType = locKinFitType;}
		inline void Add_ReactionStep(const DReactionStep* locReactionStep){dReactionSteps.push_back(locReactionStep);}
		inline void Add_AnalysisAction(DAnalysisAction* locAnalysisAction){dAnalysisActions.push_back(locAnalysisAction);}

		// SET TRACK SELECTION FACTORIES //Command-line values will override these values
		inline void Set_ChargedTrackFactoryTag(string locChargedTrackFactoryTag){dChargedTrackFactoryTag = locChargedTrackFactoryTag;}
		inline void Set_NeutralShowerFactoryTag(string locNeutralShowerFactoryTag){dNeutralShowerFactoryTag = locNeutralShowerFactoryTag;}

		// SET PRE-DPARTICLECOMBO CUT VALUES //Command-line values will override these values
		inline void Set_MinChargedPIDFOM(double locMinChargedPIDFOM){dMinChargedPIDFOM = pair<bool, double>(true, locMinChargedPIDFOM);}
		inline void Set_MinPhotonPIDFOM(double locMinPhotonPIDFOM){dMinPhotonPIDFOM = pair<bool, double>(true, locMinPhotonPIDFOM);}
		inline void Set_MinCombinedPIDFOM(double locMinCombinedPIDFOM){dMinCombinedPIDFOM = pair<bool, double>(true, locMinCombinedPIDFOM);}
		inline void Set_MinTrackingFOM(double locMinTrackingFOM){dMinTrackingFOM = pair<bool, double>(true, locMinTrackingFOM);}
		inline void Set_MinCombinedTrackingFOM(double locMinCombinedTrackingFOM){dMinCombinedTrackingFOM = pair<bool, double>(true, locMinCombinedTrackingFOM);}
		inline void Set_MaxPhotonRFDeltaT(double locMaxPhotonRFDeltaT){dMaxPhotonRFDeltaT = pair<bool, double>(true, locMaxPhotonRFDeltaT);}
		inline void Set_MinProtonMomentum(double locMinProtonMomentum){dMinProtonMomentum = pair<bool, double>(true, locMinProtonMomentum);}
		inline void Set_HasDetectorMatchFlag(bool locHasDetectorMatchFlag){dHasDetectorMatchFlag = pair<bool, bool>(true, locHasDetectorMatchFlag);}

		// GET CONTROL MEMBERS:
		inline string Get_ReactionName(void) const{return dReactionName;}
		inline DKinFitType Get_KinFitType(void) const{return dKinFitType;}

		// GET REACTION STEPS:
		inline size_t Get_NumReactionSteps(void) const{return dReactionSteps.size();}
		const DReactionStep* Get_ReactionStep(size_t locStepIndex) const;
		void Get_ReactionSteps(Particle_t locInitialPID, deque<const DReactionStep*>& locReactionSteps) const;

		// GET ANALYSIS ACTIONS:
		inline size_t Get_NumAnalysisActions(void) const{return dAnalysisActions.size();}
		DAnalysisAction* Get_AnalysisAction(size_t locActionIndex) const;
		DAnalysisAction* Get_LastAnalysisAction(void) const;

		// GET PIDs:
		void Get_DetectedFinalPIDs(deque<Particle_t>& locDetectedPIDs, bool locIncludeDuplicatesFlag = false) const;
		void Get_DetectedFinalPIDs(deque<deque<Particle_t> >& locDetectedPIDs, bool locIncludeDuplicatesFlag = false) const;
		void Get_DetectedFinalChargedPIDs(deque<Particle_t>& locDetectedChargedPIDs, bool locIncludeDuplicatesFlag = false) const;
		void Get_DetectedFinalChargedPIDs(deque<deque<Particle_t> >& locDetectedChargedPIDs, bool locIncludeDuplicatesFlag = false) const;
		void Get_FinalStatePIDs(deque<Particle_t>& locFinalStatePIDs, bool locIncludeDuplicatesFlag = false) const;
		bool Get_MissingPID(Particle_t& locPID) const; //false if none missing

		// GET PARTICLE NAME STRINGS:
		string Get_DetectedParticlesROOTName(void) const;
		string Get_InitialParticlesROOTName(void) const;
		string Get_DecayChainFinalParticlesROOTNames(Particle_t locInitialPID, bool locKinFitResultsFlag) const;
		string Get_DecayChainFinalParticlesROOTNames(Particle_t locInitialPID, int locUpToStepIndex, deque<Particle_t> locUpThroughPIDs, bool locKinFitResultsFlag) const;
		string Get_DecayChainFinalParticlesROOTNames(size_t locStepIndex, int locUpToStepIndex, deque<Particle_t> locUpThroughPIDs, bool locKinFitResultsFlag, bool locExpandDecayingParticlesFlag) const;

		// GET TRACK SELECTION FACTORIES //Command-line values will override these values
		inline string Get_ChargedTrackFactoryTag(void) const{return dChargedTrackFactoryTag;}
		inline string Get_NeutralShowerFactoryTag(void) const{return dNeutralShowerFactoryTag;}

		// GET PRE-DPARTICLECOMBO CUT VALUES //Command-line values will override these values
		inline pair<bool, double> Get_MinChargedPIDFOM(void) const{return dMinChargedPIDFOM;}
		inline pair<bool, double> Get_MinPhotonPIDFOM(void) const{return dMinPhotonPIDFOM;}
		inline pair<bool, double> Get_MinCombinedPIDFOM(void) const{return dMinCombinedPIDFOM;}
		inline pair<bool, double> Get_MinTrackingFOM(void) const{return dMinTrackingFOM;}
		inline pair<bool, double> Get_MinCombinedTrackingFOM(void) const{return dMinCombinedTrackingFOM;}
		inline pair<bool, double> Get_MaxPhotonRFDeltaT(void) const{return dMaxPhotonRFDeltaT;}
		inline pair<bool, double> Get_MinProtonMomentum(void) const{return dMinProtonMomentum;}
		inline pair<bool, bool> Get_HasDetectorMatchFlag(void) const{return dHasDetectorMatchFlag;}

		// ROOT OUTPUT:
		void Enable_TTreeOutput(string locTTreeOutputFileName);
		inline string Get_TTreeOutputFileName(void) const{return dTTreeOutputFileName;}
		inline bool Get_EnableTTreeOutputFlag(void) const{return dEnableTTreeOutputFlag;}
		inline void Set_MinThrownMatchFOMForROOT(double locMinThrownMatchFOMForROOT){dMinThrownMatchFOMForROOT = locMinThrownMatchFOMForROOT;}
		inline double Get_MinThrownMatchFOMForROOT(void) const{return dMinThrownMatchFOMForROOT;}

		// OTHER:
		bool Check_IsDecayingParticle(Particle_t locPID, size_t locSearchStartIndex = 1) const;
		bool Check_AreStepsIdentical(const DReaction* locReaction) const;

	private:
		// PRIVATE METHODS:
		DReaction(void); //make default constructor private. MUST set a name upon construction (and must be unique!)

		// CONTROL MEMBERS:
		string dReactionName; //must be unique
		DKinFitType dKinFitType; //defined in ANALYSIS/DKinFitResults.h

		// ROOT TTREE OUTPUT:
		bool dEnableTTreeOutputFlag; //default is false
		string dTTreeOutputFileName;
		double dMinThrownMatchFOMForROOT; //cut applied when setting matching when filling ROOT tree //default is -1 (always passes)

		// REACTION AND ANALYSIS MEMBERS:
		deque<const DReactionStep*> dReactionSteps;
		deque<DAnalysisAction*> dAnalysisActions;

		// TRACK SELECTION FACTORIES:
			//Command-line values will override these values
		string dChargedTrackFactoryTag; //default is ""
		string dNeutralShowerFactoryTag; //default is ""

		// PRE-DPARTICLECOMBO CUT VALUES
			//bool = true/false for cut enabled/disabled, double = cut value
			//Command-line values (variable names are below in all-caps) will override these values
			//all cuts are disabled by default except dMinProtonMomentum: 300 MeV/c (value used during track reconstruction)
			//note: tracks with no PID information are not cut-by/included-in the PID cuts
		pair<bool, double> dMinChargedPIDFOM; //COMBO:MIN_CHARGED_PID_FOM - the minimum PID FOM for a particle used for this DReaction
		pair<bool, double> dMinPhotonPIDFOM; //COMBO:MIN_PHOTON_PID_FOM - the minimum PID FOM for a neutral particle used for this DReaction
		pair<bool, double> dMinCombinedPIDFOM; //COMBO:MIN_COMBINED_PID_FOM - the minimum combined PID FOM for all particles used for this DReaction
		pair<bool, double> dMinTrackingFOM; //COMBO:MIN_TRACKING_FOM - the minimum Tracking FOM for a charged track used for this DReaction
		pair<bool, double> dMinCombinedTrackingFOM; //COMBO:MIN_COMBINED_TRACKING_FOM - the minimum combined Tracking FOM for all charged tracks used for this DReaction
		pair<bool, double> dMaxPhotonRFDeltaT; //COMBO:MAX_PHOTON_RF_DELTAT - the maximum photon-rf time difference: used for photon selection
		pair<bool, double> dMinProtonMomentum; //COMBO:MIN_PROTON_MOMENTUM - when testing whether a non-proton DChargedTrackHypothesis could be a proton, this is the minimum momentum it can have
		pair<bool, bool> dHasDetectorMatchFlag; //COMBO:HAS_DETECTOR_MATCH_FLAG - if both are true, require tracks to have a detector match
};

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

inline bool DReaction::Check_IsDecayingParticle(Particle_t locPID, size_t locSearchStartIndex) const
{
	//see if this pid is a parent in a future step
	for(size_t loc_k = locSearchStartIndex; loc_k < dReactionSteps.size(); ++loc_k)
	{
		if(dReactionSteps[loc_k]->Get_InitialParticleID() == locPID)
			return true;
	}
	return false;
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

inline void DReaction::Enable_TTreeOutput(string locTTreeOutputFileName)
{
	dEnableTTreeOutputFlag = true;
	dTTreeOutputFileName = locTTreeOutputFileName;
}

#endif // _DReaction_

