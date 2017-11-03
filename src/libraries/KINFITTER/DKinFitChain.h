#ifndef _DKinFitChain_
#define _DKinFitChain_

#include <vector>
#include <algorithm>
#include <memory>

#include "DKinFitParticle.h"
#include "DKinFitChainStep.h"
#include "DResettable.h"

//This class is not necessary to use the kinematic fitter, but it is necessary to use some of the setup help functions in DKinFitUtils
	//Is mostly useful when coding for the generic situation of ANY possible decay chain (rather than handling a specific one)

using namespace std;

class DKinFitChain : public DResettable
{
	public:

		void Reset(void);
		void Release(void);

		//GET, ADD STEPS
		shared_ptr<const DKinFitChainStep> Get_KinFitChainStep(size_t locStepIndex) const;
		void Add_KinFitChainStep(const shared_ptr<DKinFitChainStep>& locKinFitChainStep){dKinFitChainSteps.push_back(locKinFitChainStep);}
		size_t Get_NumKinFitChainSteps(void) const{return dKinFitChainSteps.size();}

		//GET ALL PARTICLES
		vector<shared_ptr<DKinFitParticle>> Get_AllParticles(void) const;

		//GET CONTROL INFO
		signed char Get_DefinedParticleStepIndex(void) const{return dDefinedParticleStepIndex;}
		bool Get_IsInclusiveChannelFlag(void) const{return dIsInclusiveChannelFlag;}
		signed char Get_DecayStepIndex(const shared_ptr<DKinFitParticle>& locKinFitParticle) const;

		//SET CONTROL INFO
		void Set_DefinedParticleStepIndex(signed char locDefinedParticleStepIndex){dDefinedParticleStepIndex = locDefinedParticleStepIndex;}
		void Set_IsInclusiveChannelFlag(bool locIsInclusiveChannelFlag){dIsInclusiveChannelFlag = locIsInclusiveChannelFlag;}
		void Set_DecayStepIndex(const shared_ptr<DKinFitParticle>& locKinFitParticle, int locDecayStepIndex){dDecayStepIndices[locKinFitParticle] = locDecayStepIndex;}

		//PRINT INFO
		void Print_InfoToScreen(void) const;

	private:

		vector<shared_ptr<DKinFitChainStep>> dKinFitChainSteps;
		map<shared_ptr<DKinFitParticle>, char> dDecayStepIndices; //key is decaying particle, value is the step representing the particle decay
		signed char dDefinedParticleStepIndex = -1; //step containing the missing or open-ended-decaying particle, -1 if none
		bool dIsInclusiveChannelFlag = false; //i.e. does the missing particle have PID 0 (unknown)
};

inline void DKinFitChain::Reset(void)
{
	dDefinedParticleStepIndex = -1;
	dIsInclusiveChannelFlag = false;
	dKinFitChainSteps.clear();
	dDecayStepIndices.clear();
}

inline void DKinFitChain::Release(void)
{
	dKinFitChainSteps.clear();
	dDecayStepIndices.clear();
}

inline signed char DKinFitChain::Get_DecayStepIndex(const shared_ptr<DKinFitParticle>& locKinFitParticle) const
{
	auto locIterator = dDecayStepIndices.find(locKinFitParticle);
	return ((locIterator != dDecayStepIndices.end()) ? locIterator->second : -1);
}

inline shared_ptr<const DKinFitChainStep> DKinFitChain::Get_KinFitChainStep(size_t locStepIndex) const
{
	return ((locStepIndex < dKinFitChainSteps.size()) ? std::const_pointer_cast<const DKinFitChainStep>(dKinFitChainSteps[locStepIndex]) : nullptr);
}

inline vector<shared_ptr<DKinFitParticle>> DKinFitChain::Get_AllParticles(void) const
{
	vector<shared_ptr<DKinFitParticle>> locAllParticles;
	for(size_t loc_i = 0; loc_i < dKinFitChainSteps.size(); ++loc_i)
	{
		auto locStepParticles = dKinFitChainSteps[loc_i]->Get_AllParticles();
		locAllParticles.insert(locAllParticles.end(), locStepParticles.begin(), locStepParticles.end());
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
	auto locIterator = dDecayStepIndices.begin();
	for(; locIterator != dDecayStepIndices.end(); ++locIterator)
		cout << locIterator->first->Get_PID() << ", " << locIterator->first << ", " << int(locIterator->second) << endl;

	cout << "DKinFitChain: defined particle step index, inclusive channel flag = " << int(dDefinedParticleStepIndex) << ", " << dIsInclusiveChannelFlag << endl;
}

#endif // _DKinFitChain_
