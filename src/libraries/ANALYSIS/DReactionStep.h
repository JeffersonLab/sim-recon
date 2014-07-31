#ifndef _DReactionStep_
#define _DReactionStep_

#include <deque>
#include <string>
#include <iostream>
#include <stdlib.h>

#include "particleType.h"

using namespace std;

class DReactionStep
{
	public:
		// CONSTRUCTOR
		DReactionStep(void) : dInitialParticleID(Unknown), dTargetParticleID(Unknown), dMissingParticleIndex(-1), dApplyKinFitMassConstraintOnInitialParticleFlag(true) { };

		bool Are_ParticlesIdentical(const DReactionStep* locReactionStep) const; //note order can be re-arranged!

		// SET OBJECT DATA:
		void Set_InitialParticleID(Particle_t locPID, bool locIsMissingFlag = false); //TRUE IS NOT SUPPORTED YET!
		inline void Set_TargetParticleID(Particle_t locPID){dTargetParticleID = locPID;}
		void Add_FinalParticleID(Particle_t locPID, bool locIsMissingFlag = false);
		inline void Set_ApplyKinFitMassConstraintOnInitialParticleFlag(bool locFlag){dApplyKinFitMassConstraintOnInitialParticleFlag = locFlag;}
		void Reset(void);

		// GET INITIAL, TARGET, AND MISSING DATA:
		inline Particle_t Get_InitialParticleID(void) const{return dInitialParticleID;}
		inline Particle_t Get_TargetParticleID(void) const{return dTargetParticleID;}
		inline int Get_MissingParticleIndex(void) const{return dMissingParticleIndex;}
		inline bool Get_ApplyKinFitMassConstraintOnInitialParticleFlag(void) const{return dApplyKinFitMassConstraintOnInitialParticleFlag;}

		// GET FINAL PARTICLE PIDs:
		Particle_t Get_FinalParticleID(size_t locFinalParticleIndex) const;
		inline size_t Get_NumFinalParticleIDs(void) const{return dFinalParticleIDs.size();}
		inline void Get_FinalParticleIDs(deque<Particle_t>& locPIDs) const{locPIDs = dFinalParticleIDs;}

		// GET PARTICLE NAME STRINGS:
		string Get_StepName(void) const; //e.g. Lambda_->_Pi-_(Proton)
		string Get_StepROOTName(void) const; //e.g. #it{#Lambda}#rightarrow#it{#pi}^{-}(#it{p})
		string Get_FinalParticlesROOTName(void) const; //e.g. #it{#pi}^{-}(#it{p})
		void Get_FinalParticlesROOTName(deque<string>& locParticleNames) const; //e.g. #it{#pi}^{-}, #it{p}
		string Get_InitialParticlesROOTName(void) const; //e.g. #it{#Lambda}
		string Get_FinalNonMissingParticlesROOTName(void) const; //e.g. #it{#pi}^{-}

		// GET FINAL PARTICLE PIDs - BY TRAIT:
		void Get_NonMissingFinalChargedPIDs(deque<Particle_t>& locPIDs) const;
		void Get_NonMissingFinalNeutralPIDs(deque<Particle_t>& locPIDs) const;
		void Get_NonMissingFinalPIDs(deque<Particle_t>& locPIDs) const;
		bool Get_MissingPID(Particle_t& locPID) const; //false if none missing

	private:
		// PID MEMBERS:
		Particle_t dInitialParticleID; //e.g. lambda, gamma
		Particle_t dTargetParticleID; //Unknown for no target (default)
		deque<Particle_t> dFinalParticleIDs;

		// CONTROL MEMBERS:
		int dMissingParticleIndex; //-1 for no missing particles, -2 for missing init (beam) particle (not yet supported!), else final state particle at this index is missing (0 -> x)
		bool dApplyKinFitMassConstraintOnInitialParticleFlag; //default true, is ignored when not applicable (e.g. init is non-decaying (beam) particle)
};

inline void DReactionStep::Reset(void)
{
	dInitialParticleID = Unknown;
	dTargetParticleID = Unknown;
	dMissingParticleIndex = -1;
	dApplyKinFitMassConstraintOnInitialParticleFlag = true;
	dFinalParticleIDs.clear();
}

inline void DReactionStep::Set_InitialParticleID(Particle_t locPID, bool locIsMissingFlag)
{
	if(IsResonance(locPID))
	{
		cout << "ERROR: CANNOT SET RESONANCE PID. ABORTING." << endl;
		abort();
	}

	if(locIsMissingFlag)
	{
//		dMissingParticleIndex = -2;
		cout << "ERROR: MISSING BEAM PARTICLE IS NOT YET SUPPORTED! ABORTING." << endl;
		abort();
	}

	dInitialParticleID = locPID;
}

inline void DReactionStep::Add_FinalParticleID(Particle_t locPID, bool locIsMissingFlag)
{
	if(IsResonance(locPID))
	{
		cout << "ERROR: CANNOT SET RESONANCE PID. ABORTING." << endl;
		abort();
	}

	if(locIsMissingFlag)
	{
		if(dMissingParticleIndex != -1)
		{
			cout << "ERROR: MORE THAN ONE MISSING PARTICLE. ABORTING." << endl;
			abort();
		}
		dFinalParticleIDs.push_back(locPID);
		dMissingParticleIndex = dFinalParticleIDs.size() - 1;
	}
	else if(locPID == Unknown)
	{
		cout << "ERROR: CANNOT SET UNKNOWN PID AS NON-MISSING FINAL PARTICLE. ABORTING." << endl;
		abort();
	}
	else
		dFinalParticleIDs.push_back(locPID);
}

inline Particle_t DReactionStep::Get_FinalParticleID(size_t locFinalParticleIndex) const
{
	if(locFinalParticleIndex >= dFinalParticleIDs.size())
		return Unknown;
	return dFinalParticleIDs[locFinalParticleIndex];
}

inline string DReactionStep::Get_InitialParticlesROOTName(void) const
{
	string locStepROOTName = ParticleName_ROOT(dInitialParticleID);
	if(dTargetParticleID != Unknown)
		locStepROOTName += ParticleName_ROOT(dTargetParticleID);
	return locStepROOTName;
}

inline string DReactionStep::Get_FinalParticlesROOTName(void) const
{
	string locStepROOTName;
	for(size_t loc_i = 0; loc_i < dFinalParticleIDs.size(); ++loc_i)
	{
		if(int(loc_i) == dMissingParticleIndex)
			continue;
		locStepROOTName += ParticleName_ROOT(dFinalParticleIDs[loc_i]);
	}
	if(dMissingParticleIndex >= 0)
		locStepROOTName += string("(") + ParticleName_ROOT(dFinalParticleIDs[dMissingParticleIndex]) + string(")");
	return locStepROOTName;
}

inline void DReactionStep::Get_FinalParticlesROOTName(deque<string>& locParticleNames) const
{
	locParticleNames.clear();
	for(size_t loc_i = 0; loc_i < dFinalParticleIDs.size(); ++loc_i)
	{
		if(int(loc_i) == dMissingParticleIndex)
			continue;
		locParticleNames.push_back(ParticleName_ROOT(dFinalParticleIDs[loc_i]));
	}
	if(dMissingParticleIndex >= 0)
		locParticleNames.push_back(string("(") + ParticleName_ROOT(dFinalParticleIDs[dMissingParticleIndex]) + string(")"));
}

inline string DReactionStep::Get_FinalNonMissingParticlesROOTName(void) const
{
	string locStepROOTName;
	for(size_t loc_i = 0; loc_i < dFinalParticleIDs.size(); ++loc_i)
	{
		if(int(loc_i) == dMissingParticleIndex)
			continue;
		locStepROOTName += ParticleName_ROOT(dFinalParticleIDs[loc_i]);
	}
	return locStepROOTName;
}

inline string DReactionStep::Get_StepROOTName(void) const
{
	string locStepROOTName = Get_InitialParticlesROOTName();
	locStepROOTName += "#rightarrow";
	locStepROOTName += Get_FinalParticlesROOTName();
	return locStepROOTName;
}

inline string DReactionStep::Get_StepName(void) const
{
	string locStepName = ParticleType(dInitialParticleID);
	if(dTargetParticleID != Unknown)
		locStepName += string("_") + ParticleType(dTargetParticleID);
	locStepName += "_->";
	for(size_t loc_i = 0; loc_i < dFinalParticleIDs.size(); ++loc_i)
	{
		if(int(loc_i) == dMissingParticleIndex)
			continue;
		locStepName += string("_") + ParticleType(dFinalParticleIDs[loc_i]);
	}
	if(dMissingParticleIndex >= 0)
		locStepName += string("_(") + ParticleType(dFinalParticleIDs[dMissingParticleIndex]) + string(")");
	return locStepName;
}

inline void DReactionStep::Get_NonMissingFinalChargedPIDs(deque<Particle_t>& locPIDs) const
{
	locPIDs.clear();
	for(size_t loc_i = 0; loc_i < dFinalParticleIDs.size(); ++loc_i)
	{
		if(int(loc_i) == dMissingParticleIndex)
			continue;
		if(ParticleCharge(dFinalParticleIDs[loc_i]) != 0)
			locPIDs.push_back(dFinalParticleIDs[loc_i]);
	}
}

inline void DReactionStep::Get_NonMissingFinalNeutralPIDs(deque<Particle_t>& locPIDs) const
{
	locPIDs.clear();
	for(size_t loc_i = 0; loc_i < dFinalParticleIDs.size(); ++loc_i)
	{
		if(int(loc_i) == dMissingParticleIndex)
			continue;
		if(ParticleCharge(dFinalParticleIDs[loc_i]) == 0)
			locPIDs.push_back(dFinalParticleIDs[loc_i]);
	}
}

inline void DReactionStep::Get_NonMissingFinalPIDs(deque<Particle_t>& locPIDs) const
{
	locPIDs.clear();
	for(size_t loc_i = 0; loc_i < dFinalParticleIDs.size(); ++loc_i)
	{
		if(int(loc_i) == dMissingParticleIndex)
			continue;
		locPIDs.push_back(dFinalParticleIDs[loc_i]);
	}
}

inline bool DReactionStep::Get_MissingPID(Particle_t& locPID) const
{
	if((dMissingParticleIndex == -1) || (dMissingParticleIndex >= int(dFinalParticleIDs.size())))
		return false;
	locPID = dFinalParticleIDs[dMissingParticleIndex];
	return true;
}

#endif // _DReactionStep_

