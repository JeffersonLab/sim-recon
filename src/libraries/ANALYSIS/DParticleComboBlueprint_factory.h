#ifndef _DParticleComboBlueprint_factory_
#define _DParticleComboBlueprint_factory_

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>
#include <ANALYSIS/DParticleComboBlueprint.h>
#include <ANALYSIS/DReaction.h>
#include "ANALYSIS/DTrackTimeBased_factory_Combo.h"

#include <JANA/JObject.h>
#include <particleType.h>
#include <PID/DChargedTrack.h>
#include <PID/DNeutralShower.h>
#include <PID/DVertex.h>
#include <PID/DDetectorMatches.h>

#include <deque>
#include <map>
#include <vector>

using namespace std;
using namespace jana;

class DParticleComboBlueprint_factory : public jana::JFactory<DParticleComboBlueprint>
{
	public:
		DParticleComboBlueprint_factory(){};
		~DParticleComboBlueprint_factory(){};

		void Reset_Pools(void);

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop* locEventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop* locEventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		void Get_Reactions(JEventLoop *locEventLoop, vector<const DReaction*>& locReactions) const;

		void Build_ParticleComboBlueprints(const DReaction* locReaction);
		bool Setup_ComboLoop(const DReaction* locReaction, deque<deque<int> >& locResumeAtIndexDeque, deque<deque<int> >& locNumPossibilitiesDeque, map<int, int>& locInitialParticleStepFromIndexMap, map<pair<int, int>, int>& locFinalStateDecayStepIndexMap);
		void Find_Combos(const DReaction* locReaction, deque<deque<int> >& locResumeAtIndexDeque, const deque<deque<int> >& locNumPossibilitiesDeque, map<int, int>& locInitialParticleStepFromIndexMap, map<pair<int, int>, int>& locFinalStateDecayStepIndexMap);

		bool Handle_EndOfReactionStep(const DReaction* locReaction, DParticleComboBlueprint*& locParticleComboBlueprint, DParticleComboBlueprintStep*& locParticleComboBlueprintStep, int& locStepIndex, int& locParticleIndex, deque<deque<int> >& locResumeAtIndexDeque, const deque<deque<int> >& locNumPossibilitiesDeque, map<int, int>& locInitialParticleStepFromIndexMap);
		bool Handle_Decursion(DParticleComboBlueprint* locParticleComboBlueprint, deque<deque<int> >& locResumeAtIndexDeque, const deque<deque<int> >& locNumPossibilitiesDeque, int& locParticleIndex, int& locStepIndex, DParticleComboBlueprintStep*& locParticleComboBlueprintStep);

		bool Check_IfDuplicateStepCombo(const DParticleComboBlueprint* locParticleComboBlueprint, const DParticleComboBlueprintStep* locCurrentStep, int locStepIndex, deque<deque<int> >& locResumeAtIndexDeque, const deque<deque<int> >& locNumPossibilitiesDeque) const;
		bool Check_IfStepsAreIdentical(const DParticleComboBlueprint* locParticleComboBlueprint, const DParticleComboBlueprintStep* locCurrentStep, const DParticleComboBlueprintStep* locPreviousStep) const;
		int Grab_DecayingParticle(Particle_t locAnalysisPID, int& locResumeAtIndex, const DReaction* locReaction, int locStepIndex, int locParticleIndex);

		const JObject* Grab_DetectedParticle(const DReaction* locReaction, Particle_t locAnalysisPID, int& locResumeAtIndex);
		const JObject* Grab_NeutralShower(vector<const DNeutralShower*>& locNeutralShowers, int& locResumeAtIndex);
		const JObject* Grab_DetectedTrack(const DReaction* locReaction, Particle_t locAnalysisPID, vector<const DChargedTrack*>& locTracks, int& locResumeAtIndex) const;

		const DChargedTrackHypothesis* Get_ChargedHypothesisToUse(const DChargedTrack* locChargedTrack, Particle_t locAnalysisPID) const;

		bool Calc_FinalStateP4(size_t locTotalNumSteps, const DParticleComboBlueprint* locParticleComboBlueprint, const DParticleComboBlueprintStep* locNewParticleComboBlueprintStep, int locStepIndexToGrab, DLorentzVector& locFinalStateP4) const;

		DParticleComboBlueprintStep* Get_ParticleComboBlueprintStepResource(void);
		inline void Recycle_ParticleComboBlueprintStep(DParticleComboBlueprintStep* locParticleComboBlueprintStep){dParticleComboBlueprintStepPool_Available.push_back(locParticleComboBlueprintStep);}

		deque<DParticleComboBlueprintStep*> dParticleComboBlueprintStepPool_All;
		deque<DParticleComboBlueprintStep*> dParticleComboBlueprintStepPool_Available;

		set<const JObject*> dCurrentComboSourceObjects;
		set<const DParticleComboBlueprintStep*> dSavedBlueprintSteps;

		//Objects for this event
		vector<const DChargedTrack*> dChargedTracks;
		vector<const DChargedTrack*> dPositiveChargedTracks;
		vector<const DChargedTrack*> dNegativeChargedTracks;
		vector<const DNeutralShower*> dNeutralShowers;

		unsigned int dDebugLevel;
		size_t MAX_DParticleComboBlueprintStepPoolSize;
		const DVertex* dVertex;
		const DDetectorMatches* dDetectorMatches;

		// PRE-DPARTICLECOMBO TRACK SELECTION FACTORY TAGS
		string dTrackSelectionTag;
		string dShowerSelectionTag;

		// PRE-DPARTICLECOMBO CUT VALUES
			//(first) bool = true/false for cut enabled/disabled, double = cut value
			//Command-line values will override those set in the DReaction
		pair<bool, double> dMinProtonMomentum; //when testing whether a non-proton DChargedTrackHypothesis could be a proton, this is the minimum momentum it can have
		pair<bool, size_t> dMaxExtraGoodTracks; //"good" defined as: at least one PID hypothesis passes other implemented "good" cuts

		//PIDs that should have been created during track reconstruction
		set<Particle_t> dAvailablePIDs;

		//used to see if can resuse memory with an identical, previously-created step
			//a map is used instead of a loop over previous combos because map access is significantly faster if #combos is very large
		map<DParticleComboBlueprintStep, DParticleComboBlueprintStep*> dBlueprintStepMap;

		DTrackTimeBased_factory_Combo* dTrackTimeBasedFactory_Combo;
};

#endif // _DParticleComboBlueprint_factory_

