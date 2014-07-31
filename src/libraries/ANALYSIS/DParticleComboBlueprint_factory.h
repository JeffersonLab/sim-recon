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
#include <HDGEOMETRY/DGeometry.h>
#include <DANA/DApplication.h>

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

		jerror_t Build_ParticleComboBlueprints(JEventLoop* locEventLoop, const DReaction* locReaction);
		bool Setup_ComboLoop(const DReaction* locReaction, int locNumDetectedNeutralParticles, int locNumDetectedChargedParticles, int locNumDetectedPositiveParticles, int locNumDetectedNegativeParticles, deque<deque<int> >& locResumeAtIndexDeque, deque<deque<int> >& locNumPossibilitiesDeque);
		void Find_Combos(const DReaction* locReaction, deque<const JObject*>& locNeutralShowerDeque, deque<const JObject*>& locChargedTrackDeque_Positive, deque<const JObject*>& locChargedTrackDeque_Negative, deque<deque<int> >& locResumeAtIndexDeque, const deque<deque<int> >& locNumPossibilitiesDeque, vector<DParticleComboBlueprint*>& locParticleComboBlueprints);

		bool Handle_EndOfReactionStep(const DReaction* locReaction, DParticleComboBlueprint*& locParticleComboBlueprint, DParticleComboBlueprintStep*& locParticleComboBlueprintStep, int& locStepIndex, int& locParticleIndex, deque<deque<int> >& locResumeAtIndexDeque, const deque<deque<int> >& locNumPossibilitiesDeque, vector<DParticleComboBlueprint*>& locParticleComboBlueprints);
		bool Handle_Decursion(DParticleComboBlueprint* locParticleComboBlueprint, deque<deque<int> >& locResumeAtIndexDeque, const deque<deque<int> >& locNumPossibilitiesDeque, int& locParticleIndex, int& locStepIndex, DParticleComboBlueprintStep*& locParticleComboBlueprintStep);

		bool Check_IfDuplicateStepCombo(const DParticleComboBlueprint* locParticleComboBlueprint, int locStepIndex, deque<deque<int> >& locResumeAtIndexDeque, const deque<deque<int> >& locNumPossibilitiesDeque) const;
		bool Check_IfStepsAreIdentical(const DParticleComboBlueprint* locParticleComboBlueprint, const DParticleComboBlueprintStep* locCurrentStep, const DParticleComboBlueprintStep* locPreviousStep) const;
		int Grab_DecayingParticle(Particle_t locAnalysisPID, int& locResumeAtIndex, const DReaction* locReaction, int locStepIndex, int locParticleIndex);
		const JObject* Grab_DetectedTrack(DParticleComboBlueprint* locParticleComboBlueprint, Particle_t locAnalysisPID, int& locResumeAtIndex, deque<const JObject*>& locNeutralShowerDeque, deque<const JObject*>& locChargedTrackDeque_Positive, deque<const JObject*>& locChargedTrackDeque_Negative);
		const JObject* Choose_SourceObject(const DReaction* locReaction, Particle_t locAnalysisPID, DParticleComboBlueprint* locParticleComboBlueprint, deque<const JObject*>& locSourceObjects, int& locResumeAtIndex) const;

		bool Cut_TrackingFOM(const DReaction* locReaction, const DChargedTrackHypothesis* locChargedTrackHypothesis) const;
		bool Cut_HasDetectorMatch(const DReaction* locReaction, const DChargedTrackHypothesis* locChargedTrackHypothesis) const;

		DParticleComboBlueprintStep* Get_ParticleComboBlueprintStepResource(void);
		inline void Recycle_ParticleComboBlueprintStep(DParticleComboBlueprintStep* locParticleComboBlueprintStep){dParticleComboBlueprintStepPool_Available.push_back(locParticleComboBlueprintStep);}

		deque<DParticleComboBlueprintStep*> dParticleComboBlueprintStepPool_All;
		deque<DParticleComboBlueprintStep*> dParticleComboBlueprintStepPool_Available;

		set<const JObject*> dCurrentComboSourceObjects;

		unsigned int dDebugLevel;
		size_t MAX_DParticleComboBlueprintStepPoolSize;

		// PRE-DPARTICLECOMBO TRACK SELECTION FACTORY TAGS
			//bool = true to get tracks from specified factory, false otherwise
			//Command-line values will override those set in the DReaction
		pair<bool, string> dReactionTrackSelectionTag;
		pair<bool, string> dReactionShowerSelectionTag;

		// PRE-DPARTICLECOMBO CUT VALUES
			//(first) bool = true/false for cut enabled/disabled, double = cut value
			//Command-line values will override those set in the DReaction
		pair<bool, double> dMinTrackingFOM; //the minimum Tracking FOM for a charged track used for this DReaction
		pair<bool, double> dMinProtonMomentum; //when testing whether a non-proton DChargedTrackHypothesis could be a proton, this is the minimum momentum it can have
		pair<bool, bool> dHasDetectorMatchFlag; //if both are true, require tracks to have a detector match

		//used to see if can resuse memory with an identical, previously-created step
			//a map is used instead of a loop over previous combos because map access is significantly faster if #combos is very large
		map<DParticleComboBlueprintStep, DParticleComboBlueprintStep*> dBlueprintStepMap;

		DTrackTimeBased_factory_Combo* dTrackTimeBasedFactory_Combo;
};

#endif // _DParticleComboBlueprint_factory_

