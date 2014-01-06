#include "DReactionStep.h"

void DReactionStep::Reset(void)
{
	dInitialParticleID = Unknown;
	dTargetParticleID = Unknown;
	dMissingParticleIndex = -1;
	dFinalParticleIDs.clear();
}

void DReactionStep::Set_InitialParticleID(Particle_t locPID, bool locIsMissingFlag)
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

void DReactionStep::Add_FinalParticleID(Particle_t locPID, bool locIsMissingFlag)
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

Particle_t DReactionStep::Get_FinalParticleID(size_t locFinalParticleIndex) const
{
	if(locFinalParticleIndex >= dFinalParticleIDs.size())
		return Unknown;
	return dFinalParticleIDs[locFinalParticleIndex];
}

string DReactionStep::Get_StepName(void) const
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

string DReactionStep::Get_InitialParticlesROOTName(void) const
{
	string locStepROOTName = ParticleName_ROOT(dInitialParticleID);
	if(dTargetParticleID != Unknown)
		locStepROOTName += ParticleName_ROOT(dTargetParticleID);
	return locStepROOTName;
}

string DReactionStep::Get_FinalParticlesROOTName(void) const
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

void DReactionStep::Get_FinalParticlesROOTName(deque<string>& locParticleNames) const
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

string DReactionStep::Get_FinalNonMissingParticlesROOTName(void) const
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

string DReactionStep::Get_StepROOTName(void) const
{
	string locStepROOTName = Get_InitialParticlesROOTName();
	locStepROOTName += "#rightarrow";
	locStepROOTName += Get_FinalParticlesROOTName();
	return locStepROOTName;
}

void DReactionStep::Get_NonMissingFinalChargedPIDs(deque<Particle_t>& locPIDs) const
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

void DReactionStep::Get_NonMissingFinalNeutralPIDs(deque<Particle_t>& locPIDs) const
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

void DReactionStep::Get_NonMissingFinalPIDs(deque<Particle_t>& locPIDs) const
{
	locPIDs.clear();
	for(size_t loc_i = 0; loc_i < dFinalParticleIDs.size(); ++loc_i)
	{
		if(int(loc_i) == dMissingParticleIndex)
			continue;
		locPIDs.push_back(dFinalParticleIDs[loc_i]);
	}
}

bool DReactionStep::Get_MissingPID(Particle_t& locPID) const
{
	if((dMissingParticleIndex == -1) || (dMissingParticleIndex >= int(dFinalParticleIDs.size())))
		return false;
	locPID = dFinalParticleIDs[dMissingParticleIndex];
	return true;
}

