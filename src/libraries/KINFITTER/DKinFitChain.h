#ifndef _DKinFitChain_
#define _DKinFitChain_

#include <vector>
#include <set>
#include <algorithm>

#include "DKinFitParticle.h"
#include "DKinFitChainStep.h"

//This class is not necessary to use the kinematic fitter, but it is necessary to use some of the setup help functions in DKinFitUtils
	//Is mostly useful when coding for the generic situation of ANY possible decay chain (rather than handling a specific one)

using namespace std;

class DKinFitChain
{
	public:

		DKinFitChain(void) : dDefinedParticleStepIndex(-1), dIsInclusiveChannelFlag(false) {}
		void Reset(void);

		//GET, ADD STEPS
		const DKinFitChainStep* Get_KinFitChainStep(size_t locStepIndex) const;
		void Add_KinFitChainStep(DKinFitChainStep* locKinFitChainStep){dKinFitChainSteps.push_back(locKinFitChainStep);}
		size_t Get_NumKinFitChainSteps(void) const{return dKinFitChainSteps.size();}

		//GET ALL PARTICLES
		set<DKinFitParticle*> Get_AllParticles(void) const;

		//GET CONTROL INFO
		char Get_DefinedParticleStepIndex(void) const{return dDefinedParticleStepIndex;}
		bool Get_IsInclusiveChannelFlag(void) const{return dIsInclusiveChannelFlag;}
		char Get_DecayStepIndex(DKinFitParticle* locKinFitParticle) const;

		//SET CONTROL INFO
		void Set_DefinedParticleStepIndex(int locDefinedParticleStepIndex){dDefinedParticleStepIndex = locDefinedParticleStepIndex;}
		void Set_IsInclusiveChannelFlag(bool locIsInclusiveChannelFlag){dIsInclusiveChannelFlag = locIsInclusiveChannelFlag;}
		void Set_DecayStepIndex(DKinFitParticle* locKinFitParticle, int locDecayStepIndex){dDecayStepIndices[locKinFitParticle] = locDecayStepIndex;}

		//PRINT INFO
		void Print_InfoToScreen(void) const;

	private:

		vector<DKinFitChainStep*> dKinFitChainSteps;
		map<DKinFitParticle*, char> dDecayStepIndices; //key is decaying particle, value is the step representing the particle decay
		char dDefinedParticleStepIndex; //step containing the missing or open-ended-decaying particle, -1 if none
		bool dIsInclusiveChannelFlag; //i.e. does the missing particle have PID 0 (unknown)
};

inline void DKinFitChain::Reset(void)
{
	dDefinedParticleStepIndex = -1;
	dIsInclusiveChannelFlag = false;
	dKinFitChainSteps.clear();
	dDecayStepIndices.clear();
}

inline char DKinFitChain::Get_DecayStepIndex(DKinFitParticle* locKinFitParticle) const
{
	map<DKinFitParticle*, char>::const_iterator locIterator = dDecayStepIndices.find(locKinFitParticle);
	return ((locIterator != dDecayStepIndices.end()) ? locIterator->second : -1);
}

inline const DKinFitChainStep* DKinFitChain::Get_KinFitChainStep(size_t locStepIndex) const
{
	return ((locStepIndex < dKinFitChainSteps.size()) ? dKinFitChainSteps[locStepIndex] : NULL);
}

inline set<DKinFitParticle*> DKinFitChain::Get_AllParticles(void) const
{
	set<DKinFitParticle*> locAllParticles;
	for(size_t loc_i = 0; loc_i < dKinFitChainSteps.size(); ++loc_i)
	{
		set<DKinFitParticle*> locStepParticles = dKinFitChainSteps[loc_i]->Get_AllParticles();
		locAllParticles.insert(locStepParticles.begin(), locStepParticles.end());
	}
	return locAllParticles;
}

inline void DKinFitChain::Print_InfoToScreen(void) const
{
	for(size_t loc_i = 0; loc_i < dKinFitChainSteps.size(); ++loc_i)
	{
		cout << "DKinFitChain: Printing step " << loc_i << endl;
		dKinFitChainSteps[loc_i]->Print_InfoToScreen();
	}

	cout << "DKinFitChain: PID, Pointer, decay-step indices:" << endl;
	map<DKinFitParticle*, char>::const_iterator locIterator = dDecayStepIndices.begin();
	for(; locIterator != dDecayStepIndices.end(); ++locIterator)
		cout << locIterator->first->Get_PID() << ", " << locIterator->first << ", " << locIterator->second << endl;

	cout << "DKinFitChain: defined particle step index, inclusive channel flag = " << dDefinedParticleStepIndex << ", " << dIsInclusiveChannelFlag << endl;
}

#endif // _DKinFitChain_
