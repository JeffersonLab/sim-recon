#ifndef _DKinFitChainStep_
#define _DKinFitChainStep_

#include <vector>

#include "DKinFitParticle.h"
#include "DResettable.h"

//This class is not necessary to use the kinematic fitter, but it is necessary to use some of the setup help functions in DKinFitUtils
	//Is mostly useful when coding for the generic situation of ANY possible decay chain (rather than handling a specific one)

using namespace std;

class DKinFitChainStep : public DResettable
{
	public:

		void Reset(void);
		void Release(void);

		//GET PARTICLES
		vector<shared_ptr<DKinFitParticle>> Get_InitialParticles(void) const{return dInitialParticles;}
		vector<shared_ptr<DKinFitParticle>> Get_FinalParticles(void) const{return dFinalParticles;}
		vector<shared_ptr<DKinFitParticle>> Get_AllParticles(void) const;
		shared_ptr<DKinFitParticle> Get_InitialParticle(size_t locIndex) const{return dInitialParticles[locIndex];}
		shared_ptr<DKinFitParticle> Get_FinalParticle(size_t locIndex) const{return dFinalParticles[locIndex];}

		//GET CONTROL INFO
		signed char Get_InitialParticleDecayFromStepIndex(void) const{return dInitialParticleDecayFromStepIndex;}
		bool Get_ConstrainDecayingMassFlag(void) const{return dConstrainDecayingMassFlag;}

		//ADD PARTICLES
		void Add_InitialParticle(const shared_ptr<DKinFitParticle>& locInitialParticle){dInitialParticles.push_back(locInitialParticle);}
		void Add_FinalParticle(const shared_ptr<DKinFitParticle>& locFinalParticle){dFinalParticles.push_back(locFinalParticle);}
		void Set_InitialParticle(const shared_ptr<DKinFitParticle>& locInitialParticle, size_t locIndex){dInitialParticles[locIndex] = locInitialParticle;}
		void Set_FinalParticle(const shared_ptr<DKinFitParticle>& locFinalParticle, size_t locIndex){dFinalParticles[locIndex] = locFinalParticle;}

		//SET CONTROL INFO
		void Set_InitialParticleDecayFromStepIndex(int locDecayFromStepIndex){dInitialParticleDecayFromStepIndex = locDecayFromStepIndex;}
		void Set_ConstrainDecayingMassFlag(bool locConstrainDecayingMassFlag){dConstrainDecayingMassFlag = locConstrainDecayingMassFlag;}

		//PRINT INFO
		void Print_InfoToScreen(void) const;

	private:

		//refers to the decaying particle in dInitialParticles //-1 if none, else index points to step index it is produced at
		signed char dInitialParticleDecayFromStepIndex = -1;
		bool dConstrainDecayingMassFlag = false; //true to constrain mass of the initial state particle

		vector<shared_ptr<DKinFitParticle>> dInitialParticles;
		vector<shared_ptr<DKinFitParticle>> dFinalParticles;
};

inline void DKinFitChainStep::Reset(void)
{
	dInitialParticleDecayFromStepIndex = -1;
	dConstrainDecayingMassFlag = false;
	dInitialParticles.clear();
	dFinalParticles.clear();
}

inline void DKinFitChainStep::Release(void)
{
	dInitialParticles.clear();
	dFinalParticles.clear();
}

inline vector<shared_ptr<DKinFitParticle>> DKinFitChainStep::Get_AllParticles(void) const
{
	auto locAllParticles = dInitialParticles;
	locAllParticles.insert(locAllParticles.end(), dFinalParticles.begin(), dFinalParticles.end());
	locAllParticles.erase(std::remove(locAllParticles.begin(), locAllParticles.end(), nullptr), locAllParticles.end());
	return locAllParticles;
}

inline void DKinFitChainStep::Print_InfoToScreen(void) const
{
	cout << "DKinFitChainStep decay from, constrain mass flags = " << int(dInitialParticleDecayFromStepIndex) << ", " << dConstrainDecayingMassFlag << endl;

	cout << "DKinFitChainStep init particles: PIDs, pointers:" << endl;
	for(auto& locParticle : dInitialParticles)
	{
		if(locParticle == nullptr)
			cout << "X, nullptr" << endl;
		else
			cout << locParticle->Get_PID() << ", " << locParticle << endl;
	}

	cout << "DKinFitChainStep final particles: PIDs, pointers:" << endl;
	for(auto& locParticle : dFinalParticles)
	{
		if(locParticle == nullptr)
			cout << "X, nullptr" << endl;
		else
			cout << locParticle->Get_PID() << ", " << locParticle << endl;
	}
}

#endif // _DKinFitChainStep_
