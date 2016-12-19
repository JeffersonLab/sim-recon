#ifndef _DKinFitChainStep_
#define _DKinFitChainStep_

#include <vector>
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
		char Get_InitialParticleDecayFromStepIndex(void) const{return dInitialParticleDecayFromStepIndex;}
		bool Get_ConstrainDecayingMassFlag(void) const{return dConstrainDecayingMassFlag;}

		//ADD PARTICLES
		void Add_InitialParticle(DKinFitParticle* locInitialParticle){dInitialParticles.insert(locInitialParticle);}
		void Add_FinalParticle(DKinFitParticle* locFinalParticle){dFinalParticles.insert(locFinalParticle);}

		//SET CONTROL INFO
		void Set_InitialParticleDecayFromStepIndex(int locDecayFromStepIndex){dInitialParticleDecayFromStepIndex = locDecayFromStepIndex;}
		void Set_ConstrainDecayingMassFlag(bool locConstrainDecayingMassFlag){dConstrainDecayingMassFlag = locConstrainDecayingMassFlag;}

		//PRINT INFO
		void Print_InfoToScreen(void) const;

	private:

		//refers to the decaying particle in dInitialParticles //-1 if none, else index points to step index it is produced at
		char dInitialParticleDecayFromStepIndex;
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

inline void DKinFitChainStep::Print_InfoToScreen(void) const
{
	cout << "DKinFitChainStep decay from, constrain mass flags = " << dInitialParticleDecayFromStepIndex << ", " << dConstrainDecayingMassFlag << endl;

	cout << "DKinFitChainStep init particles: PIDs, pointers:" << endl;
	set<DKinFitParticle*>::const_iterator locIterator = dInitialParticles.begin();
	for(; locIterator != dInitialParticles.end(); ++locIterator)
		cout << (*locIterator)->Get_PID() << ", " << *locIterator << endl;

	cout << "DKinFitChainStep final particles: PIDs, pointers:" << endl;
	for(locIterator = dFinalParticles.begin(); locIterator != dFinalParticles.end(); ++locIterator)
		cout << (*locIterator)->Get_PID() << ", " << *locIterator << endl;
}

#endif // _DKinFitChainStep_
