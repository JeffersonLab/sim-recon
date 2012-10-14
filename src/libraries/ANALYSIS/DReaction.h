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
		DReaction(string locReactionName) : dReactionName(locReactionName){} //User must specify a unique reaction name upon construction

		// SET OBJECT DATA:
		inline void Set_KinFitType(DKinFitType locKinFitType){dKinFitType = locKinFitType;}
		inline void Add_ReactionStep(const DReactionStep* locReactionStep){dReactionSteps.push_back(locReactionStep);}
		inline void Exclude_DecayingParticleFromP4KinFit(size_t locStepIndex){dDecayingParticlesExcludedFromP4Kinfit.push_back(locStepIndex);}
		inline void Add_AnalysisAction(DAnalysisAction* locAnalysisAction){dAnalysisActions.push_back(locAnalysisAction);}

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
		void Get_DetectedFinalPIDs(deque<Particle_t>& locDetectedPIDs) const;
		void Get_DetectedFinalPIDs(deque<deque<Particle_t> >& locDetectedPIDs) const;
		void Get_DetectedFinalChargedPIDs(deque<Particle_t>& locDetectedChargedPIDs) const;
		void Get_DetectedFinalChargedPIDs(deque<deque<Particle_t> >& locDetectedChargedPIDs) const;
		void Get_FinalStatePIDs(deque<Particle_t>& locFinalStatePIDs) const;

		// GET PARTICLE NAME STRINGS:
		string Get_DetectedParticlesROOTName(void) const;
		string Get_InitialParticlesROOTName(void) const;
		void Get_DecayChainFinalParticlesROOTNames(Particle_t locInitialPID, deque<string>& locNames, bool locKinFitResultsFlag = false) const;
		void Get_DecayChainFinalParticlesROOTNames(Particle_t locInitialPID, deque<deque<string> >& locParticleNames, deque<string>& locNames, bool locKinFitResultsFlag = false) const;

		// OTHER:
		bool Check_IsDecayingParticle(Particle_t locPID, size_t locSearchStartIndex = 1) const;
		bool Check_IfDecayingParticleExcludedFromP4KinFit(size_t locStepIndex) const;

	private:
		// PRIVATE METHODS:
		DReaction(void); //make default constructor private. MUST set a name upon construction (and must be unique!)
		void Get_DecayChainFinalParticlesROOTNames(size_t locStepIndex, deque<deque<string> >& locNames, bool locKinFitResultsFlag = false) const;

		// REACTION AND ANALYSIS MEMBERS:
		deque<const DReactionStep*> dReactionSteps;
		deque<DAnalysisAction*> dAnalysisActions;

		// CONTROL MEMBERS:
		string dReactionName; //must be unique
		DKinFitType dKinFitType; //defined in ANALYSIS/DKinFitResults.h
		deque<size_t> dDecayingParticlesExcludedFromP4Kinfit; //to exclude decaying particles from the kinematic fit (resonances are automatically excluded) //size_t is step index where it is a parent
};

#endif // _DReaction_

