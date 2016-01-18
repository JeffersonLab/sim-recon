#ifndef _DKinFitChain_
#define _DKinFitChain_

#include <deque>
#include <set>
#include <algorithm>

#include "DKinFitParticle.h"

//This class is not necessary to use the kinematic fitter, but it is necessary to use some of the setup help functions in DKinFitUtils
	//Is mostly useful when coding for the generic situation of ANY possible decay chain (rather than handling a specific one)

using namespace std;

class DKinFitChainStep
{
	public:

		DKinFitChainStep(void) : dInitialParticleDecayFromStepIndex(-1), dConstrainDecayingMassFlag(false) {}
		void Reset(void);

		//GET PARTICLES
		set<DKinFitParticle*> Get_InitialParticles(void) const{return dInitialParticles;}
		set<DKinFitParticle*> Get_FinalParticles(void) const{return dFinalParticles;}
		set<DKinFitParticle*> Get_AllParticles(void) const;

		//GET CONTROL INFO
		int Get_InitialParticleDecayFromStepIndex(void) const{return dInitialParticleDecayFromStepIndex;}
		bool Get_ConstrainDecayingMassFlag(void) const{return dConstrainDecayingMassFlag;}

		//ADD PARTICLES
		void Add_InitialParticle(DKinFitParticle* locInitialParticle){dInitialParticles.insert(locInitialParticle);}
		void Add_FinalParticle(DKinFitParticle* locFinalParticle){dFinalParticles.insert(locFinalParticle);}

		//SET CONTROL INFO
		void Set_InitialParticleDecayFromStepIndex(int locDecayFromStepIndex){dInitialParticleDecayFromStepIndex = locDecayFromStepIndex;}
		void Set_ConstrainDecayingMassFlag(bool locConstrainDecayingMassFlag){dConstrainDecayingMassFlag = locConstrainDecayingMassFlag;}

	private:

		//refers to the decaying particle in dInitialParticles //-1 if none, else index points to step index it is produced at
		int dInitialParticleDecayFromStepIndex;
		bool dConstrainDecayingMassFlag; //true to constrain mass of the initial state particle

		set<DKinFitParticle*> dInitialParticles;
		set<DKinFitParticle*> dFinalParticles;
};

inline void DKinFitChainStep::Reset(void)
{
	dInitialParticleDecayFromStepIndex = -1;
	dConstrainDecayingMassFlag = false;
	dInitialParticles.clear();
	dFinalParticles.clear();
}

inline set<DKinFitParticle*> DKinFitChainStep::Get_AllParticles(void) const
{
	set<DKinFitParticle*> locAllParticles;
	set_union(dInitialParticles.begin(), dInitialParticles.end(), dFinalParticles.begin(), dFinalParticles.end(), inserter(locAllParticles, locAllParticles.begin()));
	return locAllParticles;
}

class DKinFitChain
{
	public:

		DKinFitChain(void) : dDefinedParticleStepIndex(-1), dIsInclusiveChannelFlag(false) {}
		void Reset(void);

		//GET, ADD STEPS
		const DKinFitChainStep* Get_KinFitChainStep(size_t locStepIndex) const;
		void Add_KinFitChainStep(DKinFitChainStep* locKinFitChainStep){dKinFitChainSteps.push_back(locKinFitChainStep);}
		size_t Get_NumKinFitChainSteps(void) const{return dKinFitChainSteps.size();}

		//GET CONTROL INFO
		int Get_DefinedParticleStepIndex(void) const{return dDefinedParticleStepIndex;}
		bool Get_IsInclusiveChannelFlag(void) const{return dIsInclusiveChannelFlag;}
		int Get_DecayStepIndex(DKinFitParticle* locKinFitParticle) const;

		//SET CONTROL INFO
		void Set_DefinedParticleStepIndex(int locDefinedParticleStepIndex){dDefinedParticleStepIndex = locDefinedParticleStepIndex;}
		void Set_IsInclusiveChannelFlag(bool locIsInclusiveChannelFlag){dIsInclusiveChannelFlag = locIsInclusiveChannelFlag;}
		void Set_DecayStepIndex(DKinFitParticle* locKinFitParticle, int locDecayStepIndex){dDecayStepIndices[locKinFitParticle] = locDecayStepIndex;}

	private:

		deque<DKinFitChainStep*> dKinFitChainSteps;
		map<DKinFitParticle*, int> dDecayStepIndices; //key is decaying particle, value is the step representing the particle decay
		int dDefinedParticleStepIndex; //step containing the missing or open-ended-decaying particle, -1 if none
		bool dIsInclusiveChannelFlag; //i.e. does the missing particle have PID 0 (unknown)
};

inline void DKinFitChain::Reset(void)
{
	dDefinedParticleStepIndex = -1;
	dIsInclusiveChannelFlag = false;
	dKinFitChainSteps.clear();
	dDecayStepIndices.clear();
}

inline int DKinFitChain::Get_DecayStepIndex(DKinFitParticle* locKinFitParticle) const
{
	map<DKinFitParticle*, int>::const_iterator locIterator = dDecayStepIndices.find(locKinFitParticle);
	return ((locIterator != dDecayStepIndices.end()) ? locIterator->second : -1);
}

inline const DKinFitChainStep* DKinFitChain::Get_KinFitChainStep(size_t locStepIndex) const
{
	return ((locStepIndex < dKinFitChainSteps.size()) ? dKinFitChainSteps[locStepIndex] : NULL);
}

#endif // _DKinFitChain_

