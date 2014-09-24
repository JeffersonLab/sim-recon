#ifndef _DAnalysisResults_
#define _DAnalysisResults_

#include <map>
#include <deque>

#include "JANA/JObject.h"
#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DParticleCombo.h"

using namespace std;
using namespace jana;

class DAnalysisResults : public JObject
{
	public:
		JOBJECT_PUBLIC(DAnalysisResults);

		inline void Set_Reaction(const DReaction* locReaction){dReaction = locReaction;}
		inline const DReaction* Get_Reaction(void) const{return dReaction;}

		inline size_t Get_NumPassedParticleCombos(void) const{return dPassedParticleCombos.size();}
		inline void Get_PassedParticleCombos(deque<const DParticleCombo*>& locPassedParticleCombos) const{locPassedParticleCombos = dPassedParticleCombos;}
		inline void Add_PassedParticleCombo(const DParticleCombo* locPassedParticleCombo){dPassedParticleCombos.push_back(locPassedParticleCombo);}

	private:
		const DReaction* dReaction;
		deque<const DParticleCombo*> dPassedParticleCombos; //DParticleCombo objects that passed all DAnalysisAction cuts. 
};

#endif // _DAnalysisResults_

