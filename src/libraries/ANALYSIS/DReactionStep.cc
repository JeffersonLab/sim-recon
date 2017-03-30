#include "DReactionStep.h"

namespace DAnalysis
{

/************************************************************** DREACTIONSTEP FUNCTIONS ***************************************************************/

int DReactionStep::Prepare_InfoArguments(vector<Particle_t>& locFinalPIDs, Particle_t locMissingFinalPID, bool locInclusiveFlag, bool locBeamMissingFlag) const
{
	if((locInclusiveFlag || locBeamMissingFlag) && (locMissingFinalPID != Unknown))
	{
		cout << "ERROR: CANNOT HAVE MISSING PID + MISSING BEAM OR INCLUSIVE. ABORTING." << endl;
		abort();
	}

	if(locMissingFinalPID != Unknown)
	{
		locFinalPIDs.push_back(locMissingFinalPID);
		return locFinalPIDs.size() - 1;
	}
	else if(locInclusiveFlag)
		return DReactionStepInfo::Get_MissingIndex_Inclusive();
	else if(locBeamMissingFlag)
		return DReactionStepInfo::Get_MissingIndex_Beam();
	else
		return DReactionStepInfo::Get_MissingIndex_None();
}

void DReactionStep::Add_FinalParticleID(Particle_t locPID, bool locIsMissingFlag)
{
	Check_IsResonance(locPID);

	if(!locIsMissingFlag)
	{
		if(locPID == Unknown)
		{
			cout << "ERROR: CANNOT SET UNKNOWN PID AS NON-MISSING FINAL PARTICLE. ABORTING." << endl;
			abort();
		}
		dReactionStepInfo->dFinalPIDs.push_back(locPID);
		return;
	}

	//now is missing, check to make sure none others are missing
	if(dReactionStepInfo->dMissingParticleIndex != DReactionStepInfo::Get_MissingIndex_None())
	{
		cout << "ERROR: CANNOT SET MORE THAN ONE MISSING PARTICLE IN A STEP. ABORTING." << endl;
		abort();
	}

	//if unknown, instead set as inclusive
	if(locPID == Unknown)
		dReactionStepInfo->dMissingParticleIndex = DReactionStepInfo::Get_MissingIndex_Inclusive();
	else
	{
		dReactionStepInfo->dMissingParticleIndex = dReactionStepInfo->dFinalPIDs.size();
		dReactionStepInfo->dFinalPIDs.push_back(locPID);
	}
}

/************************************************************** NAMESPACE-SCOPE FUNCTIONS ***************************************************************/

string Get_InitialParticlesName(const DReactionStep* locStep, bool locTLatexFlag)
{
    std::function<char*(Particle_t)> locGetNameFunc = locTLatexFlag ? ParticleName_ROOT : ParticleType;
    string locInitPIDString = locGetNameFunc(locStep->Get_InitialPID());
	string locStepName = locStep->Get_IsMissingBeamFlag() ? string("(") + locInitPIDString + string(")") : locInitPIDString;

	Particle_t locSecondBeamPID = locStep->Get_SecondBeamPID();
	Particle_t locTargetPID = locStep->Get_TargetPID();
	Particle_t locSecondPID = (locTargetPID != Unknown) ? locTargetPID : locSecondBeamPID;
	if(locSecondPID != Unknown)
	{
		if(!locTLatexFlag)
			locStepName += "_";
		locStepName += locGetNameFunc(locSecondPID);
	}

	return locStepName;
}

vector<string> Get_FinalParticleNames(const DReactionStep* locStep, bool locIncludeMissingFlag, bool locTLatexFlag)
{
	vector<string> locParticleNames;
    std::function<char*(Particle_t)> locGetNameFunc = locTLatexFlag ? ParticleName_ROOT : ParticleType;

	vector<Particle_t> locFinalPIDs = locStep->Get_FinalPIDs(false);
	auto locTransformer = [](Particle_t& locPID) -> string {return locGetNameFunc(locPID);};
	std::transform(locFinalPIDs.begin(), locFinalPIDs.end(), std::back_inserter(locParticleNames), locTransformer);

	if(!locIncludeMissingFlag)
		return locParticleNames;

	if(locStep->Get_IsInclusiveFlag())
	{
		locParticleNames.push_back("(X)");
		return locParticleNames;
	}

	Particle_t locMissingPID = locStep->Get_MissingPID();
	if(locMissingPID != Unknown)
		locParticleNames.push_back(string("(") + locGetNameFunc(locMissingPID) + string(")"));

	return locParticleNames;
}

bool Are_ParticlesIdentical(const DReactionStep* locStep1, const DReactionStep* locStep2, bool locExceptMissingUnknownInInputFlag)
{
	if(locStep1->Get_InitialPID() != locStep2->Get_InitialPID())
		return false;
	if(locStep1->Get_SecondBeamPID() != locStep2->Get_SecondBeamPID())
		return false;
	if(locStep1->Get_TargetPID() != locStep2->Get_TargetPID())
		return false;

	//FiX THIS
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
	vector<Particle_t> locFinalParticleIDsCopy = locReactionStep->dFinalParticleIDs;
	for(size_t loc_i = 0; loc_i < dFinalParticleIDs.size(); ++loc_i)
	{
		bool locMatchFoundFlag = false;
		vector<Particle_t>::iterator locIterator = locFinalParticleIDsCopy.begin();
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

} //end DAnalysis namespace
