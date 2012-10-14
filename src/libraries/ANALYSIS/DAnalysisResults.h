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

		inline void Add_FailedParticleCombo(const DParticleCombo* locParticleCombo, size_t locAnalysisActionIndex){dFailedParticleComboMap[locParticleCombo] = locAnalysisActionIndex;}
		inline void Get_FailedParticleComboMap(map<const DParticleCombo*, size_t>& locFailedParticleComboMap) const{locFailedParticleComboMap = dFailedParticleComboMap;}

		inline void Get_PassedParticleCombos(deque<const DParticleCombo*>& locPassedParticleCombos) const{locPassedParticleCombos = dPassedParticleCombos;}
		inline void Add_PassedParticleCombo(const DParticleCombo* locPassedParticleCombo){dPassedParticleCombos.push_back(locPassedParticleCombo);}

	private:
		const DReaction* dReaction;
		map<const DParticleCombo*, size_t> dFailedParticleComboMap; //indicates at which action these DParticleCombo failed to pass the cut
		deque<const DParticleCombo*> dPassedParticleCombos; //DParticleCombo objects that passed all DAnalysisAction cuts. 
};

#endif // _DAnalysisResults_

