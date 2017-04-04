#include "DReactionStep.h"

namespace DAnalysis
{

/************************************************************** MEMBER FUNCTIONS ***************************************************************/

int DReactionStep::Prepare_InfoArguments(vector<Particle_t>& locFinalPIDs, Particle_t locMissingFinalPID, bool locInclusiveFlag, bool locBeamMissingFlag, bool locSecondBeamMissingFlag) const
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
		return DReactionStep::Get_ParticleIndex_Inclusive();
	else if(locBeamMissingFlag)
		return DReactionStep::Get_ParticleIndex_Initial();
	else if(locSecondBeamMissingFlag)
		return DReactionStep::Get_ParticleIndex_SecondBeam();
	else
		return DReactionStep::Get_ParticleIndex_None();
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
	if(dReactionStepInfo->dMissingParticleIndex != DReactionStep::Get_ParticleIndex_None())
	{
		cout << "ERROR: CANNOT SET MORE THAN ONE MISSING PARTICLE IN A STEP. ABORTING." << endl;
		abort();
	}

	//if unknown, instead set as inclusive
	if(locPID == Unknown)
		dReactionStepInfo->dMissingParticleIndex = DReactionStep::Get_ParticleIndex_Inclusive();
	else
	{
		dReactionStepInfo->dMissingParticleIndex = dReactionStepInfo->dFinalPIDs.size();
		dReactionStepInfo->dFinalPIDs.push_back(locPID);
	}
}

vector<Particle_t> DReactionStep::Get_FinalPIDs(bool locIncludeMissingFlag, Charge_t locCharge, bool locIncludeDuplicatesFlag) const
{
	vector<Particle_t> locFinalPIDs = dReactionStepInfo->dFinalPIDs;
	if(!locIncludeMissingFlag && (dReactionStepInfo->dMissingParticleIndex >= 0))
		locFinalPIDs.erase(locFinalPIDs.begin() + dReactionStepInfo->dMissingParticleIndex);

	auto locComparator = [&](Particle_t& locPID) -> bool {return !Is_CorrectCharge(locPID, locCharge);};
	locFinalPIDs.erase(std::remove_if(locFinalPIDs.begin(), locFinalPIDs.end(), locComparator), locFinalPIDs.end());

	if(!locIncludeDuplicatesFlag)
	{
		std::sort(locFinalPIDs.begin(), locFinalPIDs.end());
		locFinalPIDs.erase(std::unique(locFinalPIDs.begin(), locFinalPIDs.end()), locFinalPIDs.end());
	}
	return locFinalPIDs;
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
	auto locPIDTransformer = [](Particle_t& locPID) -> string {return locGetNameFunc(locPID);};
	std::transform(locFinalPIDs.begin(), locFinalPIDs.end(), std::back_inserter(locParticleNames), locPIDTransformer);

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
	if(locStep1->Get_MissingPID() != locStep2->Get_MissingPID())
		return false;

	auto locFinalPIDs1 = locStep1->Get_FinalPIDs();
	auto locFinalPIDs2 = locStep2->Get_FinalPIDs();
	if(locFinalPIDs1.size() != locFinalPIDs2.size())
		return false;

	//note order can be re-arranged!
	return std::is_permutation(locFinalPIDs1.begin(), locFinalPIDs1.end(), locFinalPIDs2.begin());
}

} //end DAnalysis namespace
