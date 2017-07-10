#include "ANALYSIS/DReactionStepVertexInfo.h"

namespace DAnalysis
{

/************************************************************** MEMBER FUNCTIONS ***************************************************************/

void DReactionStepVertexInfo::Register_DecayingParticleConstraints(bool locFitFlag, const vector<pair<int, int>>& locNoConstrainDecayingParticles,
		const map<pair<int, int>, const DReactionStepVertexInfo*>& locFullConstrainDecayingParticles = {})
{
	dDecayingParticles_FullConstrain[locFitFlag] = locFullConstrainDecayingParticles;
	for(auto locMapPair : locFullConstrainDecayingParticles)
		dFullConstrainParticles[locFitFlag].emplace_back(locMapPair.first);
	for(auto locParticlePair : locNoConstrainDecayingParticles)
	{
		dDecayingParticles_NoConstrain[locFitFlag].emplace(locParticlePair, nullptr);
		dNoConstrainParticles[locFitFlag].emplace_back(locParticlePair);
	}
	std::sort(dFullConstrainParticles[locFitFlag].begin(), dFullConstrainParticles[locFitFlag].end());
	std::sort(dNoConstrainParticles[locFitFlag].begin(), dNoConstrainParticles[locFitFlag].end());
}

void DReactionStepVertexInfo::Set_ParticleIndices(bool locFitFlag, const vector<pair<int, int>>& locFullConstrainParticles, const vector<pair<int, int>>& locDecayingParticles,
		const vector<pair<int, int>>& locOnlyConstrainTimeParticles, const vector<pair<int, int>>& locNoConstrainParticles)
{
	//store
	dFullConstrainParticles[locFitFlag] = locFullConstrainParticles;
	dDecayingParticles = locDecayingParticles;
	dOnlyConstrainTimeParticles = locOnlyConstrainTimeParticles;
	dNoConstrainParticles[locFitFlag] = locNoConstrainParticles;

	//sort
	std::sort(dFullConstrainParticles[locFitFlag].begin(), dFullConstrainParticles[locFitFlag].end());
	std::sort(dDecayingParticles.begin(), dDecayingParticles.end());
	std::sort(dOnlyConstrainTimeParticles.begin(), dOnlyConstrainTimeParticles.end());
	std::sort(dNoConstrainParticles[locFitFlag].begin(), dNoConstrainParticles[locFitFlag].end());
}

//if you want beam, say initial state not target or decaying
vector<pair<int, int>> DReactionStepVertexInfo::Filter_Particles(vector<pair<int, int>> locParticles, DReactionState_t locState, Charge_t locCharge,
		bool locIncludeDecayingFlag, bool locIncludeMissingFlag, bool locIncludeTargetFlag) const
{
	if(locState != d_EitherState)
	{
		auto Check_State = [&locState](const pair<int, int>& locIndices) -> bool
		{return ((locIndices.second >= 0) != (locState == d_FinalState));};
		locParticles.erase(std::remove_if(locParticles.begin(), locParticles.end(), Check_State), locParticles.end());
	}
	if(locCharge != d_AllCharges)
	{
		auto Check_Charge = [this, &locCharge](const pair<int, int>& locIndices) -> bool
		{
			Particle_t locPID = dReaction->Get_ReactionStep(locIndices.first)->Get_PID(locIndices.second);
			return !Is_CorrectCharge(locPID, locCharge);
		};
		locParticles.erase(std::remove_if(locParticles.begin(), locParticles.end(), Check_Charge), locParticles.end());
	}
	if(!locIncludeDecayingFlag)
	{
		auto Check_Decaying = [this](const pair<int, int>& locIndices) -> bool
		{return std::binary_search(dDecayingParticles.begin(), dDecayingParticles.end(), locIndices);};
		locParticles.erase(std::remove_if(locParticles.begin(), locParticles.end(), Check_Decaying), locParticles.end());
	}
	if(!locIncludeMissingFlag)
	{
		auto Check_Missing = [this](const pair<int, int>& locIndices) -> bool
		{return (locIndices.second == dReaction->Get_ReactionStep(locIndices.first)->Get_MissingParticleIndex());};
		locParticles.erase(std::remove_if(locParticles.begin(), locParticles.end(), Check_Missing), locParticles.end());
	}
	if(!locIncludeTargetFlag)
	{
		auto Check_Target = [this](const pair<int, int>& locIndices) -> bool
		{return (locIndices.second == DReactionStep::Get_ParticleIndex_Target());};
		locParticles.erase(std::remove_if(locParticles.begin(), locParticles.end(), Check_Target), locParticles.end());
	}
	return locParticles;
}

} //end namespace DAnalysis
