#ifndef _DReactionStep_
#define _DReactionStep_

#include <deque>
#include <string>
#include <iostream>

#include "particleType.h"

using namespace std;

class DReactionStep
{
	public:
		// CONSTRUCTOR
		DReactionStep(void) : dInitialParticleID(Unknown), dTargetParticleID(Unknown), dMissingParticleIndex(-1) { };

		// SET OBJECT DATA:
		inline void Set_InitialParticleID(Particle_t locPID){dInitialParticleID = locPID;}
		inline void Set_TargetParticleID(Particle_t locPID){dTargetParticleID = locPID;}
		inline void Set_MissingParticleIndex(int locMissingParticleIndex){dMissingParticleIndex = locMissingParticleIndex;}
		inline void Add_FinalParticleID(Particle_t locPID){dFinalParticleIDs.push_back(locPID);}

		// GET INITIAL, TARGET, AND MISSING DATA:
		inline Particle_t Get_InitialParticleID(void) const{return dInitialParticleID;}
		inline Particle_t Get_TargetParticleID(void) const{return dTargetParticleID;}
		inline int Get_MissingParticleIndex(void) const{return dMissingParticleIndex;}

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
		Particle_t dTargetParticleID; //Unknown for no target
		deque<Particle_t> dFinalParticleIDs;

		// CONTROL MEMBERS:
		int dMissingParticleIndex; //-1 for no missing particles, else final state particle at this index is missing (0 -> x)
};

#endif // _DReactionStep_

