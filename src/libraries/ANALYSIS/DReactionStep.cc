#include "DReactionStep.h"

bool DReactionStep::Are_ParticlesIdentical(const DReactionStep* locReactionStep, bool locExceptMissingUnknownInInputFlag) const
{
	if(dInitialParticleID != locReactionStep->dInitialParticleID)
		return false;
	if(dTargetParticleID != locReactionStep->dTargetParticleID)
		return false;

	int locSizeChange = 0;
	if(locExceptMissingUnknownInInputFlag)
	{
		Particle_t locMissingPID = Unknown;
		bool locIsMissingPID = locReactionStep->Get_MissingPID(locMissingPID);
		if(locIsMissingPID && (locMissingPID == Unknown))
			locSizeChange = 1;
	}

	if(dFinalParticleIDs.size() != (locReactionStep->dFinalParticleIDs.size() - locSizeChange))
		return false;

	//note order can be re-arranged!
	deque<Particle_t> locFinalParticleIDsCopy = locReactionStep->dFinalParticleIDs;
	for(size_t loc_i = 0; loc_i < dFinalParticleIDs.size(); ++loc_i)
	{
		bool locMatchFoundFlag = false;
		deque<Particle_t>::iterator locIterator = locFinalParticleIDsCopy.begin();
		for(; locIterator != locFinalParticleIDsCopy.end(); ++locIterator)
		{
			if(dFinalParticleIDs[loc_i] != (*locIterator))
				continue;
			locMatchFoundFlag = true;
			locFinalParticleIDsCopy.erase(locIterator);
			break;
		}
		if(!locMatchFoundFlag)
			return false;
	}
	return true;
}

