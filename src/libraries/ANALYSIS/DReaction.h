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
		inline void Exclude_DecayingParticleFromP4KinFit(size_t locStepIndex){dDecayingParticlesExcludedFromP4Kinfit.push_back(locStepIndex);}
		inline void Add_AnalysisAction(DAnalysisAction* locAnalysisAction){dAnalysisActions.push_back(locAnalysisAction);}

		// SET TRACK SELECTION FACTORIES //Command-line values will override these values
		inline void Set_ChargedTrackFactoryTag(string locChargedTrackFactoryTag){dChargedTrackFactoryTag = locChargedTrackFactoryTag;}
		inline void Set_NeutralShowerFactoryTag(string locNeutralShowerFactoryTag){dNeutralShowerFactoryTag = locNeutralShowerFactoryTag;}

		// SET PRE-DPARTICLECOMBO CUT VALUES //Command-line values will override these values
		inline void Set_MinIndividualChargedPIDFOM(double locMinIndividualChargedPIDFOM){dMinIndividualChargedPIDFOM = pair<bool, double>(true, locMinIndividualChargedPIDFOM);}
		inline void Set_MinCombinedChargedPIDFOM(double locMinCombinedChargedPIDFOM){dMinCombinedChargedPIDFOM = pair<bool, double>(true, locMinCombinedChargedPIDFOM);}
		inline void Set_MinIndividualTrackingFOM(double locMinIndividualTrackingFOM){dMinIndividualTrackingFOM = pair<bool, double>(true, locMinIndividualTrackingFOM);}
		inline void Set_MinCombinedTrackingFOM(double locMinCombinedTrackingFOM){dMinCombinedTrackingFOM = pair<bool, double>(true, locMinCombinedTrackingFOM);}
		inline void Set_MaxPhotonRFDeltaT(double locMaxPhotonRFDeltaT){dMaxPhotonRFDeltaT = pair<bool, double>(true, locMaxPhotonRFDeltaT);}
		inline void Set_MinProtonMomentum(double locMinProtonMomentum){dMinProtonMomentum = pair<bool, double>(true, locMinProtonMomentum);}


		// GET CONTROL MEMBERS:
		inline string Get_ReactionName(void) const{return dReactionName;}
		inline DKinFitType Get_KinFitType(void) const{return dKinFitType;}
		inline void Get_DecayingParticlesExcludedFromP4Kinfit(deque<size_t>& locExcludedParticleStepIndices) const{locExcludedParticleStepIndices = dDecayingParticlesExcludedFromP4Kinfit;}

		// GET REACTION STEPS:
		inline size_t Get_NumReactionSteps(void) const{return dReactionSteps.size();}
		const DReactionStep* Get_ReactionStep(size_t locStepIndex) const;
		void Get_ReactionSteps(Particle_t locInitialPID, deque<const DReactionStep*>& locReactionSteps) const;

		// GET ANALYSIS ACTIONS:
		inline size_t Get_NumAnalysisActions(void) const{return dAnalysisActions.size();}
		DAnalysisAction* Get_AnalysisAction(size_t locActionIndex) const;

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
		void Get_DecayChainFinalParticlesROOTNames(Particle_t locInitialPID, deque<string>& locNames, bool locKinFitResultsFlag = false) const;
		void Get_DecayChainFinalParticlesROOTNames(Particle_t locInitialPID, deque<deque<string> >& locParticleNames, deque<string>& locNames, bool locKinFitResultsFlag = false) const;

		// GET TRACK SELECTION FACTORIES //Command-line values will override these values
		inline string Get_ChargedTrackFactoryTag(void) const{return dChargedTrackFactoryTag;}
		inline string Get_NeutralShowerFactoryTag(void) const{return dNeutralShowerFactoryTag;}

		// GET PRE-DPARTICLECOMBO CUT VALUES //Command-line values will override these values
		inline pair<bool, double> Get_MinIndividualChargedPIDFOM(void) const{return dMinIndividualChargedPIDFOM;}
		inline pair<bool, double> Get_MinCombinedChargedPIDFOM(void) const{return dMinCombinedChargedPIDFOM;}
		inline pair<bool, double> Get_MinIndividualTrackingFOM(void) const{return dMinIndividualTrackingFOM;}
		inline pair<bool, double> Get_MinCombinedTrackingFOM(void) const{return dMinCombinedTrackingFOM;}
		inline pair<bool, double> Get_MaxPhotonRFDeltaT(void) const{return dMaxPhotonRFDeltaT;}
		inline pair<bool, double> Get_MinProtonMomentum(void) const{return dMinProtonMomentum;}

		// ROOT OUTPUT:
		inline void Enable_TTreeOutput(string locTTreeOutputFileName)
		{
			dEnableTTreeOutputFlag = true;
			dTTreeOutputFileName = locTTreeOutputFileName;
		}
		inline string Get_TTreeOutputFileName(void) const{return dTTreeOutputFileName;}
		inline bool Get_EnableTTreeOutputFlag(void) const{return dEnableTTreeOutputFlag;}

		// OTHER:
		bool Check_IsDecayingParticle(Particle_t locPID, size_t locSearchStartIndex = 1) const;
		bool Check_IfDecayingParticleExcludedFromP4KinFit(size_t locStepIndex) const;

	private:
		// PRIVATE METHODS:
		DReaction(void); //make default constructor private. MUST set a name upon construction (and must be unique!)
		void Get_DecayChainFinalParticlesROOTNames(size_t locStepIndex, deque<deque<string> >& locNames, bool locKinFitResultsFlag = false) const;

		// CONTROL MEMBERS:
		string dReactionName; //must be unique
		DKinFitType dKinFitType; //defined in ANALYSIS/DKinFitResults.h
		deque<size_t> dDecayingParticlesExcludedFromP4Kinfit; //to exclude decaying particles from the kinematic fit (resonances are automatically excluded) //size_t is step index where it is a parent

		// ROOT TTREE OUTPUT:
		bool dEnableTTreeOutputFlag; //default is false
		string dTTreeOutputFileName;

		// REACTION AND ANALYSIS MEMBERS:
		deque<const DReactionStep*> dReactionSteps;
		deque<DAnalysisAction*> dAnalysisActions;

		// TRACK SELECTION FACTORIES
			//Command-line values will override these values
		string dChargedTrackFactoryTag; //default is ""
		string dNeutralShowerFactoryTag; //default is ""

		// PRE-DPARTICLECOMBO CUT VALUES
			//bool = true/false for cut enabled/disabled, double = cut value
			//Command-line values (variable names are below in all-caps) will override these values
			//all cuts are disabled by default except dMinProtonMomentum: 300 MeV/c (value used during track reconstruction)
			//note: tracks with no PID information are not cut-by/included-in the PID cuts
		pair<bool, double> dMinIndividualChargedPIDFOM; //COMBO:MININDIVIDUALCHARGEDPIDFOM - the minimum PID FOM for a charged track used for this DReaction
		pair<bool, double> dMinCombinedChargedPIDFOM; //COMBO:MINCOMBINEDCHARGEDPIDFOM - the minimum combined PID FOM for all charged tracks used for this DReaction
		pair<bool, double> dMinIndividualTrackingFOM; //COMBO:MININDIVIDUALTRACKINGFOM - the minimum Tracking FOM for a charged track used for this DReaction
		pair<bool, double> dMinCombinedTrackingFOM; //COMBO:MINCOMBINEDTRACKINGFOM - the minimum combined Tracking FOM for all charged tracks used for this DReaction
		pair<bool, double> dMaxPhotonRFDeltaT; //COMBO:PHOTONRFDELTAT - the maximum photon-rf time difference: used for photon selection
		pair<bool, double> dMinProtonMomentum; //COMBO:MINPROTONMOMENTUM - when testing whether a non-proton DChargedTrackHypothesis could be a proton, this is the minimum momentum it can have
};

#endif // _DReaction_

