#ifndef _DParticleComboBlueprint_factory_
#define _DParticleComboBlueprint_factory_

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>
#include <ANALYSIS/DParticleComboBlueprint.h>
#include <ANALYSIS/DReaction.h>

#include <JANA/JObject.h>
#include <particleType.h>
#include <PID/DChargedTrack.h>
#include <PID/DNeutralShower.h>

#include <deque>
#include <vector>

using namespace std;
using namespace jana;

class DParticleComboBlueprint_factory : public jana::JFactory<DParticleComboBlueprint>
{
	public:
		DParticleComboBlueprint_factory(){};
		~DParticleComboBlueprint_factory(){};

		void Reset_Pools(void);
		inline double Get_MinimumProtonMomentum(void) const{return dMinimumProtonMomentum;}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop* locEventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop* locEventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		jerror_t Build_ParticleComboBlueprints(JEventLoop* locEventLoop, const DReaction* locReaction);
		void Find_Combos(const DReaction* locReaction, deque<const JObject*>& locNeutralShowerDeque, deque<const JObject*>& locChargedTrackDeque_Positive, deque<const JObject*>& locChargedTrackDeque_Negative, deque<deque<int> >& locResumeAtIndexDeque, deque<deque<int> >& locNumPossibilitiesDeque, vector<DParticleComboBlueprint*>& locParticleComboBlueprints);

		DParticleComboBlueprint* Clone_ParticleComboBlueprint(const DParticleComboBlueprint* locParticleComboBlueprint);
		inline void Recycle_ParticleComboBlueprintStep(DParticleComboBlueprintStep* locParticleComboBlueprintStep){dParticleComboBlueprintStepPool_Available.push_back(locParticleComboBlueprintStep);}

		bool Handle_EndOfReactionStep(const DReaction* locReaction, DParticleComboBlueprint*& locParticleComboBlueprint, DParticleComboBlueprintStep*& locParticleComboBlueprintStep, int& locStepIndex, int& locParticleIndex, deque<deque<int> >& locResumeAtIndexDeque, deque<deque<int> >& locNumPossibilitiesDeque, vector<DParticleComboBlueprint*>& locParticleComboBlueprints);
		bool Handle_Decursion(DParticleComboBlueprint* locParticleComboBlueprint, deque<deque<int> >& locResumeAtIndexDeque, const deque<deque<int> >& locNumPossibilitiesDeque, int& locParticleIndex, int& locStepIndex, DParticleComboBlueprintStep*& locParticleComboBlueprintStep);
		void Setup_ComboLoop(const DReaction* locReaction, int locNumDetectedNeutralParticles, int locNumDetectedPositiveParticles, int locNumDetectedNegativeParticles, deque<deque<int> >& locResumeAtIndexDeque, deque<deque<int> >& locNumPossibilitiesDeque);

		int Grab_DecayingParticle(DParticleComboBlueprint* locParticleComboBlueprint, Particle_t locAnalysisPID, int& locResumeAtIndex, const DReaction* locReaction, int locStepIndex, int locParticleIndex, DParticleComboBlueprintStep* locParticleComboBlueprintStep);
		const JObject* Grab_DetectedTrack(DParticleComboBlueprint* locParticleComboBlueprint, Particle_t locAnalysisPID, int& locResumeAtIndex, deque<const JObject*>& locNeutralShowerDeque, deque<const JObject*>& locChargedTrackDeque_Positive, deque<const JObject*>& locChargedTrackDeque_Negative);
		const JObject* Choose_SourceObject(Particle_t locAnalysisPID, DParticleComboBlueprint* locParticleComboBlueprint, deque<const JObject*>& locSourceObjects, int& locResumeAtIndex) const;

		bool Cut_TrackingChiSqPerDF(const DChargedTrackHypothesis* locChargedTrackHypothesis) const;
		bool Cut_PIDFOM(const DChargedTrackHypothesis* locChargedTrackHypothesis) const;
		bool Cut_VertexZ(const DChargedTrackHypothesis* locChargedTrackHypothesis) const;

		DParticleComboBlueprintStep* Get_ParticleComboBlueprintStepResource(void);

		deque<DParticleComboBlueprintStep*> dParticleComboBlueprintStepPool_All;
		deque<DParticleComboBlueprintStep*> dParticleComboBlueprintStepPool_Available;

		unsigned int dDebugLevel;
		size_t MAX_DParticleComboBlueprintStepPoolSize;
		double dMinimumProtonMomentum; //if detected track has momentum less than this, don't attempt to assign a Proton (or any heaveier) PID to it

		bool dVertexZCutFlag;
		double dMinVertexZ;
		double dMaxVertexZ;
		int dMaximumNumTracks;

		double dMinChargedPIDFOM;
		double dMaxTrackingChiSqPerDF;
};

#endif // _DParticleComboBlueprint_factory_

