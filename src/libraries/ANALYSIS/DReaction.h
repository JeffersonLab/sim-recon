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
		inline void Set_MinIndividualPIDFOM(double locMinIndividualPIDFOM){dMinIndividualPIDFOM = pair<bool, double>(true, locMinIndividualPIDFOM);}
		inline void Set_MinCombinedPIDFOM(double locMinCombinedPIDFOM){dMinCombinedPIDFOM = pair<bool, double>(true, locMinCombinedPIDFOM);}
		inline void Set_MinIndividualTrackingFOM(double locMinIndividualTrackingFOM){dMinIndividualTrackingFOM = pair<bool, double>(true, locMinIndividualTrackingFOM);}
		inline void Set_MinCombinedTrackingFOM(double locMinCombinedTrackingFOM){dMinCombinedTrackingFOM = pair<bool, double>(true, locMinCombinedTrackingFOM);}
		inline void Set_MaxPhotonRFDeltaT(double locMaxPhotonRFDeltaT){dMaxPhotonRFDeltaT = pair<bool, double>(true, locMaxPhotonRFDeltaT);}
		inline void Set_MinProtonMomentum(double locMinProtonMomentum){dMinProtonMomentum = pair<bool, double>(true, locMinProtonMomentum);}

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
		inline pair<bool, double> Get_MinIndividualPIDFOM(void) const{return dMinIndividualPIDFOM;}
		inline pair<bool, double> Get_MinCombinedPIDFOM(void) const{return dMinCombinedPIDFOM;}
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
		pair<bool, double> dMinIndividualPIDFOM; //COMBO:MIN_INDIVIDUAL_PID_FOM - the minimum PID FOM for a charged track used for this DReaction
		pair<bool, double> dMinCombinedPIDFOM; //COMBO:MIN_COMBINED_PID_FOM - the minimum combined PID FOM for all charged tracks used for this DReaction
		pair<bool, double> dMinIndividualTrackingFOM; //COMBO:MIN_INDIVIDUAL_TRACKING_FOM - the minimum Tracking FOM for a charged track used for this DReaction
		pair<bool, double> dMinCombinedTrackingFOM; //COMBO:MIN_COMBINED_TRACKING_FOM - the minimum combined Tracking FOM for all charged tracks used for this DReaction
		pair<bool, double> dMaxPhotonRFDeltaT; //COMBO:MAX_PHOTON_RF_DELTAT - the maximum photon-rf time difference: used for photon selection
		pair<bool, double> dMinProtonMomentum; //COMBO:MIN_PROTON_MOMENTUM - when testing whether a non-proton DChargedTrackHypothesis could be a proton, this is the minimum momentum it can have
};

#endif // _DReaction_

