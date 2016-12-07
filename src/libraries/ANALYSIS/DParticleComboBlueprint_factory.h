#ifndef _DParticleComboBlueprint_factory_
#define _DParticleComboBlueprint_factory_

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>
#include <ANALYSIS/DParticleComboBlueprint.h>
#include <ANALYSIS/DReaction.h>
#include "ANALYSIS/DTrackTimeBased_factory_Combo.h"
#include <EVENTSTORE/DESSkimData.h>

#include <JANA/JObject.h>
#include <particleType.h>
#include <SplitString.h>
#include <PID/DChargedTrack.h>
#include <PID/DNeutralShower.h>
#include <PID/DVertex.h>
#include <PID/DDetectorMatches.h>
#include <TRIGGER/DTrigger.h>

#include <deque>
#include <map>
#include <set>
#include <vector>

using namespace std;
using namespace jana;

class DParticleComboBlueprint_factory : public jana::JFactory<DParticleComboBlueprint>
{
	public:
		DParticleComboBlueprint_factory(){};
		~DParticleComboBlueprint_factory(){};

		size_t Get_ParticleComboBlueprintPoolSize(void) const{return dParticleComboBlueprintPool_Acquired.size();};
		size_t Get_ParticleComboBlueprintStepPoolSize(void const{return dParticleComboBlueprintStepPool_Acquired.size();};
		size_t Get_ParticleComboBlueprintPoolSize_Shared(void) const;
		size_t Get_ParticleComboBlueprintStepPoolSize_Shared(void) const;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop* locEventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop* locEventLoop, uint64_t eventnumber);	///< Called every event.

		/********************************************************** CONTROL SETTINGS ************************************************************/

		unsigned int dDebugLevel;
		size_t dMaxNumNeutralShowers;

		string dTrackSelectionTag;
		string dShowerSelectionTag;

		// Cuts
			//(first) bool = true/false for cut enabled/disabled, double = cut value
			//Command-line values will override those set in the DReaction
		pair<bool, double> dMinProtonMomentum; //when testing whether another charged hypo could be a proton, this is the minimum momentum it can have
		pair<bool, size_t> dMaxExtraGoodTracks; //"good" defined as: at least one PID hypothesis passes other implemented "good" cuts

		//PIDs that should have been created during track reconstruction
		set<Particle_t> dAvailablePIDs;

		DTrackTimeBased_factory_Combo* dTrackTimeBasedFactory_Combo;

		/****************************************************** INPUT DATA FOR THE EVENT ********************************************************/

		const DVertex* dVertex;
		const DDetectorMatches* dDetectorMatches;

		vector<const DChargedTrack*> dChargedTracks;
		vector<const DChargedTrack*> dPositiveChargedTracks;
		vector<const DChargedTrack*> dNegativeChargedTracks;
		vector<const DNeutralShower*> dNeutralShowers;

		/********************************************************* RESOURCE MANAGEMENT **********************************************************/

		//The resource pools need to be shared between threads, to balance out the memory load
			//They should be accessible in the header file (inline functions)
		//However, you cannot make them global/extern/static/static-member variables in the header file:
			//They would be in the header file, and the header file is included in the ANALYSIS library AND in each plugin that uses it
				//When a header file is included in a src file, it's contents are essentially copied directly into it
			//Thus there are two instances of each static variable: one in each translation unit (library)
			//Supposedly(?) they are linked together during runtime when loading, so there is (supposedly) no undefined behavior.
			//However, this causes a double free (double-deletion) when these libraries are closed at the end of the program, crashing it.
		//Thus the variables must be in a single source file that is compiled into a single library
		//However, you (somehow?) cannot make them global/extern variables in the .cc file
			//This also (somehow?) causes the double-free problem above for (at least) stl containers
			//It works for pointers-to-stl-containers and fundamental types, but I dunno why.
			//It's not good encapsulation anyway though.
		//THE SOLUTION:
			//Define the variables as static, in the source file, WITHIN A PRIVATE MEMBER FUNCTION.
			//Thus the static (shared between threads) variables themselves only have function scope.
			//Access is only available via the private member function, thus access is fully controlled.
			//They are shared amongst threads, so locks are necessary, but since they are private this class can handle it internally

		//Get resources
		DParticleComboBlueprint* Get_ParticleComboBlueprintResource(void);
		DParticleComboBlueprintStep* Get_ParticleComboBlueprintStepResource(void);

		//Recycle / reset
		void Reset_Memory(void); //call for each new event
		void Recycle_ParticleComboBlueprintStep(DParticleComboBlueprintStep* locStep){dParticleComboBlueprintStepPool_Available.push_back(locStep);}

		//Return static variables by reference
		deque<DParticleComboBlueprint*>& Get_AvailableComboDeque(void) const;
		deque<DParticleComboBlueprintStep*>& Get_AvailableComboSteps(void) const;

		//Acquire resources from the shared pool
		void Acquire_Combos(size_t locNumRequestedCombos);
		void Acquire_Steps(size_t locNumRequestedSteps);

		//the target max-size for the available pools
		size_t dTargetMaxNumAvailableCombos;
		size_t dTargetMaxNumAvailableComboSteps;

		//when no resources are at hand, the number of objects to request from the shared pool
		size_t dNumFillBufferCombos;
		size_t dNumFillBufferSteps;

		//acquired from the shared pool for this event
		deque<DParticleComboBlueprintStep*> dParticleComboBlueprintStepPool_Acquired;
		deque<DParticleComboBlueprint*> dParticleComboBlueprintPool_Acquired;

		//used as a local buffers for this thread //all recycling goes to the buffer for this thread
		deque<DParticleComboBlueprintStep*> dParticleComboBlueprintStepPool_Available;
		deque<DParticleComboBlueprint*> dParticleComboBlueprintPool_Available; //available buffer from those acquired from the pool

		/******************************************************* BUILD COMBOS FUNCTIONS *********************************************************/

		void Get_Reactions(JEventLoop *locEventLoop, vector<const DReaction*>& locReactions) const;
		void Check_ReactionNames(vector<const DReaction*>& locReactions) const;

		void Build_ParticleComboBlueprints(const DReaction* locReaction);
		bool Setup_ComboLoop(const DReaction* locReaction, deque<deque<int> >& locResumeAtIndexDeque, deque<deque<int> >& locNumPossibilitiesDeque, map<int, int>& locInitialParticleStepFromIndexMap, map<pair<int, int>, int>& locFinalStateDecayStepIndexMap);
		void Find_Combos(const DReaction* locReaction, deque<deque<int> >& locResumeAtIndexDeque, const deque<deque<int> >& locNumPossibilitiesDeque, map<int, int>& locInitialParticleStepFromIndexMap, map<pair<int, int>, int>& locFinalStateDecayStepIndexMap);

		bool Handle_EndOfReactionStep(const DReaction* locReaction, DParticleComboBlueprint*& locParticleComboBlueprint, DParticleComboBlueprintStep*& locParticleComboBlueprintStep, int& locStepIndex, int& locParticleIndex, deque<deque<int> >& locResumeAtIndexDeque, const deque<deque<int> >& locNumPossibilitiesDeque, map<int, int>& locInitialParticleStepFromIndexMap);
		bool Handle_Decursion(DParticleComboBlueprint* locParticleComboBlueprint, deque<deque<int> >& locResumeAtIndexDeque, const deque<deque<int> >& locNumPossibilitiesDeque, int& locParticleIndex, int& locStepIndex, DParticleComboBlueprintStep*& locParticleComboBlueprintStep);

		bool Check_IfDuplicateStepCombo(const DParticleComboBlueprint* locParticleComboBlueprint, const DParticleComboBlueprintStep* locCurrentStep, int locStepIndex, deque<deque<int> >& locResumeAtIndexDeque, const deque<deque<int> >& locNumPossibilitiesDeque) const;
		bool Check_IfStepsAreIdentical(const DParticleComboBlueprint* locParticleComboBlueprint, const DParticleComboBlueprintStep* locCurrentStep, const DParticleComboBlueprintStep* locPreviousStep) const;

		const JObject* Grab_DetectedParticle(const DReaction* locReaction, Particle_t locAnalysisPID, int& locResumeAtIndex);
		const JObject* Grab_NeutralShower(vector<const DNeutralShower*>& locNeutralShowers, int& locResumeAtIndex);
		const JObject* Grab_DetectedTrack(const DReaction* locReaction, Particle_t locAnalysisPID, vector<const DChargedTrack*>& locTracks, int& locResumeAtIndex) const;

		const DChargedTrackHypothesis* Get_ChargedHypothesisToUse(const DChargedTrack* locChargedTrack, Particle_t locAnalysisPID) const;

		bool Calc_FinalStateP4(size_t locTotalNumSteps, const DParticleComboBlueprint* locParticleComboBlueprint, const DParticleComboBlueprintStep* locNewParticleComboBlueprintStep, int locStepIndexToGrab, DLorentzVector& locFinalStateP4) const;

		/******************************************************** BUILD COMBOS MEMORY ***********************************************************/

		set<const JObject*> dCurrentComboSourceObjects;
		set<const DParticleComboBlueprintStep*> dSavedBlueprintSteps;

		//used to see if can resuse memory with an identical, previously-created step
			//uses a custom comparator to sort/find by step == operator, rather than pointer
		set<DParticleComboBlueprintStep*, DParticleComboBlueprintStep_Comparator> dBlueprintStepSet;
};

#endif // _DParticleComboBlueprint_factory_

