#ifndef _DParticleComboBlueprint_
#define _DParticleComboBlueprint_

#include <deque>

#include "JANA/JObject.h"
#include "particleType.h"
#include "PID/DKinematicData.h"
#include "PID/DChargedTrack.h"
#include "PID/DNeutralShower.h"
#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DParticleComboBlueprintStep.h"

using namespace std;
using namespace jana;

class DParticleComboBlueprint : public JObject
{
	public:

		JOBJECT_PUBLIC(DParticleComboBlueprint);

		const DParticleComboBlueprintStep* Get_ParticleComboBlueprintStep(size_t locStepIndex) const;
		const DParticleComboBlueprintStep* Pop_ParticleComboBlueprintStep(void);
		inline void Add_ParticleComboBlueprintStep(const DParticleComboBlueprintStep* locParticleComboBlueprintStep){dParticleComboBlueprintSteps.push_back(locParticleComboBlueprintStep);}
		inline size_t Get_NumParticleComboBlueprintSteps(void) const{return dParticleComboBlueprintSteps.size();}
		void Set_ParticleComboBlueprintStep(const DParticleComboBlueprintStep* locParticleComboBlueprintStep, size_t locStepIndex);

		inline const DReaction* Get_Reaction(void) const{return dReaction;}
		inline void Set_Reaction(const DReaction* locReaction){dReaction = locReaction;}

	private:
		const DReaction* dReaction;
		deque<const DParticleComboBlueprintStep*> dParticleComboBlueprintSteps; //must be in order you want to evaluate them
};

#endif // _DParticleComboBlueprint_

