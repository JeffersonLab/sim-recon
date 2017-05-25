#include "ANALYSIS/DReactionVertexInfo_factory.h"

using namespace std;
using namespace jana;
using namespace DAnalysis;

jerror_t DReactionVertexInfo_factory::init(void)
{
	// This is so that the created objects persist throughout the life of the program instead of being cleared each event.
	SetFactoryFlag(PERSISTANT);
	return NOERROR;
}

jerror_t DReactionVertexInfo_factory::evnt(JEventLoop *locEventLoop, uint64_t locEventNumber)
{
	auto locReactions = DAnalysis::Get_Reactions(locEventLoop);
	for(auto locReactionIterator = locReactions.begin(); locReactionIterator != locReactions.end(); ++locReactionIterator)
	{
		auto locReaction = *locReactionIterator;

		//see if a reaction that we've already done has an identical channel to this
		auto Equality_Checker = [&locReaction](const pair<const DReaction*, DReactionVertexInfo*>& locReactionPair) -> bool
			{return DAnalysis::Check_ChannelEquality(locReaction, locReactionPair.first, true);};

		//do the check
		auto locSearchIterator = std::find_if(dVertexInfoMap.begin(), dVertexInfoMap.end(), Equality_Checker);
		if(locSearchIterator != dVertexInfoMap.end())
		{
			//not unique! register reaction and move on
			auto locVertexInfo = locSearchIterator->second;
			locVertexInfo->Add_Reaction(locReaction);
			continue;
		}

		//unique channel: create vertex info
		auto locVertexInfo = Build_VertexInfo(locReaction);
		_data.push_back(locVertexInfo);
		dVertexInfoMap.emplace(locReaction, locVertexInfo);
	}

	return NOERROR;
}

DReactionVertexInfo* DReactionVertexInfo_factory::Build_VertexInfo(const DReaction* locReaction) const
{
	set<size_t> locStepIndexSet;
	vector<shared_ptr<DReactionStepVertexInfo>> locVertexInfos;
	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		if(locStepIndexSet.find(loc_i) != locStepIndexSet.end())
			continue; //already did this step

		//Start a new vertex, and save when done
		auto locVertexInfo = Setup_VertexInfo(locReaction, loc_i, nullptr);
		locVertexInfos.push_back(locVertexInfo);

		auto locStepIndexVector = locVertexInfo->Get_StepIndices();
		std::copy(locStepIndexVector.begin(), locStepIndexVector.end(), std::inserter(locStepIndexSet, locStepIndexSet.end()));

		Group_VertexParticles(locVertexInfo.get());
	}

	locVertexInfos = Link_Vertices(locReaction, locVertexInfos); //sorted by dependency order
	return new DReactionVertexInfo(locReaction, locVertexInfos);
}

shared_ptr<DReactionStepVertexInfo> DReactionVertexInfo_factory::Setup_VertexInfo(const DReaction* locReaction, size_t locStepIndex, DReactionStepVertexInfo* locVertexInfo) const
{
	//create/update vertex info
	if(locVertexInfo == nullptr)
		locVertexInfo = make_shared<DReactionStepVertexInfo>(locReaction, locStepIndex);
	else
		locVertexInfo->Add_ReactionStep(locStepIndex);

	//loop over final particles: add to the vertex constraint, dive through decaying particles that decay in-place
		//if decaying in-place: don't add (would add to constraint, but not here (purely internal))
	const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(locStepIndex);
	for(size_t loc_i = 0; loc_i < locReactionStep->Get_NumFinalPIDs(); ++loc_i)
	{
		//check if particle decays, and if so, if in-place
		int locDecayStepIndex = Get_DecayStepIndex(locReaction, locStepIndex, loc_i);
		Particle_t locPID = locReactionStep->Get_FinalPID(loc_i);
		if((locDecayStepIndex >= 0) && !IsDetachedVertex(locPID)) //yes: combine with the decay products
			Setup_VertexInfo(locReaction, locDecayStepIndex, locVertexInfo);
	}

	return locVertexInfo;
}

void DReactionVertexInfo_factory::Group_VertexParticles(DReactionStepVertexInfo* locVertexInfo)
{
	auto locReaction = locVertexInfo->Get_Reaction();
	auto locStepIndices = locVertexInfo->Get_StepIndices();
	vector<pair<int, int>> locFullConstrainParticles, locDecayingParticles, locOnlyConstrainTimeParticles, locNoConstrainParticles;

	//Decaying: Only those that can conceivably be used to constrain: All unless dLinkVerticesFlag disabled (then none)
	bool locFirstStepFlag = true;
	for(auto locStepIndex : locStepIndices)
	{
		auto locStep = locReaction->Get_ReactionStep(locStepIndex);

		//beam
		if((locStepIndex == 0) && Get_IsFirstStepBeam(locReaction)) //production
		{
			//beam 1
			if(locStep->Get_IsBeamMissingFlag())
				locNoConstrainParticles.emplace_back(locStepIndex, DReactionStep::Get_ParticleIndex_Initial());
			else
				locFullConstrainParticles.emplace_back(locStepIndex, DReactionStep::Get_ParticleIndex_Initial());

			//beam 2
			if(locStep->Get_SecondBeamPID() != Unknown)
			{
				if(locStep->Get_IsSecondBeamMissingFlag())
					locNoConstrainParticles.emplace_back(locStepIndex, DReactionStep::Get_ParticleIndex_SecondBeam());
				else
					locFullConstrainParticles.emplace_back(locStepIndex, DReactionStep::Get_ParticleIndex_SecondBeam());
			}
		}
		else if(locFirstStepFlag) //decaying //if false: not detached: only save once at this vertex (in final state), not twice!!
			locDecayingParticles.emplace_back(locStepIndex, DReactionStep::Get_ParticleIndex_Initial());

		//target
		if(locStep->Get_TargetPID() != Unknown) //target
			locNoConstrainParticles.emplace_back(locStepIndex, DReactionStep::Get_ParticleIndex_Target());

		//final state
		auto locFinalPIDs = locStep->Get_FinalPIDs();
		for(size_t loc_i = 0; loc_i < locFinalPIDs.size(); ++loc_i)
		{
			int locDecayStepIndex = Get_DecayStepIndex(locReaction, locStepIndex, loc_i);

			if(locDecayStepIndex >= 0) //decaying
				locDecayingParticles.emplace_back(locStepIndex, loc_i);
			else if(loc_i == locStep->Get_MissingParticleIndex()) //missing
				locNoConstrainParticles.emplace_back(locStepIndex, loc_i);
			else if(ParticleCharge(locFinalPIDs[loc_i]) != 0) //detected charged
				locFullConstrainParticles.emplace_back(locStepIndex, loc_i);
			else if(locFinalPIDs[loc_i] == Gamma) //photon
				locOnlyConstrainTimeParticles.emplace_back(locStepIndex, loc_i);
			else //massive neutrals can't constrain position or time!
				locNoConstrainParticles.emplace_back(locStepIndex, loc_i);
		}

		locFirstStepFlag = false;
	}

	locVertexInfo->Set_ParticleIndices(locFullConstrainParticles, locDecayingParticles, locOnlyConstrainTimeParticles, locNoConstrainParticles);
}

vector<shared_ptr<DReactionStepVertexInfo>> DReactionVertexInfo_factory::Link_Vertices(const DReaction* locReaction, vector<shared_ptr<DReactionStepVertexInfo>> locVertexInfos) const
{
	//loop over vertex-constraints-to-sort:
		//find which constraints decaying particles should be defined-by/constrained-to
		//find order in which constraints need to be constrained

	bool locProgessMadeFlag = false;
	auto locVertexIterator = locVertexInfos.begin();
	map<pair<int, int>, shared_ptr<DReactionStepVertexInfo>> locDefinedDecayingParticles;
	vector<shared_ptr<DReactionStepVertexInfo>> locSortedVertexInfos; //sorted by dependency
	while(!locVertexInfos.empty())
	{
		if(locVertexIterator == locVertexInfos.end())
		{
			//made a full loop through
			if(!locProgessMadeFlag)
				break; //no progress made: cannot constrain remaining vertices

			//reset for next pass through
			locVertexIterator = locVertexInfos.begin();
			locProgessMadeFlag = false;
			continue;
		}
		auto locVertexInfo = *locVertexIterator;

		locProgessMadeFlag = Associate_DecayingParticles(true, locVertexInfo, locDefinedDecayingParticles);
		if(!locProgessMadeFlag)
			++locVertexIterator; //try again later
		else //Erase this vertex from future consideration
		{
			locVertexIterator = locVertexInfos.erase(locVertexIterator);
			locSortedVertexInfos.push_back(locVertexInfo);
		}
	}

	//the remaining vertices are dangling: not enough info to constrain
	//save dangling info
	for(auto locVertexInfo : locVertexInfos)
	{
		Associate_DecayingParticles(false, locVertexInfo, locDefinedDecayingParticles);
		locSortedVertexInfos.push_back(locVertexInfo);
	}

	//set parent vertex infos
	unordered_map<size_t, shared_ptr<DReactionStepVertexInfo>> locVertexInfoMap;
	for(auto locVertexInfo : locSortedVertexInfos)
	{
		for(auto locStepIndex : locVertexInfo->Get_StepIndices())
			locVertexInfoMap[locStepIndex] = locVertexInfo;
	}
	for(auto locVertexInfo : locSortedVertexInfos)
	{
		auto locPrimaryStepIndex = locVertexInfo->Get_StepIndices().front();
		if(locPrimaryStepIndex == 0)
			continue; //for 0 is nullptr, but is nullptr by default
		auto locParentStepIndex = DAnalysis::Get_InitialParticleDecayFromIndices(locReaction, locPrimaryStepIndex).first;
		locVertexInfo->Set_ParentVertexInfo(locVertexInfoMap[locParentStepIndex]);
	}
	return locSortedVertexInfos;
}

bool DReactionVertexInfo_factory::Associate_DecayingParticles(bool locLinkingFlag, shared_ptr<DReactionStepVertexInfo>& locVertexInfo,
		map<pair<int, int>, shared_ptr<DReactionStepVertexInfo>>& locDefinedDecayingParticles) const
{
	//find which decaying particles at this vertex have/haven't been previously defined
	vector<pair<int, int>> locNoConstrainDecayingParticles;
	map<pair<int, int>, shared_ptr<DReactionStepVertexInfo>> locConstrainingDecayingParticles;
	auto locDecayingParticles = locVertexInfo->Get_DecayingParticles();
	for(auto locParticlePair : locDecayingParticles)
	{
		auto locIterator = locDefinedDecayingParticles.find(locParticlePair);
		if(locIterator != locDefinedDecayingParticles.end())
		{
			locConstrainingDecayingParticles.emplace(*locIterator);
			locDefinedDecayingParticles.erase(locIterator); //can't use the same one twice
		}
		else //not found
			locNoConstrainDecayingParticles.emplace_back(locParticlePair);
	}

	//see if enough tracks //if not linking, then don't need "enough": will register as dangling vertex instead
	bool locEnoughTracksFlag = (locConstrainingDecayingParticles.size() + locVertexInfo->Get_FullConstrainParticles().size() >= 2);
	if(locLinkingFlag && !locEnoughTracksFlag)
		return false; //trying to link vertices, not enough tracks yet. try again later

	//Save results
	locVertexInfo->Register_DecayingParticleConstraints(locNoConstrainDecayingParticles, locConstrainingDecayingParticles);

	//each of the constraining decaying particles was a no-constrain at another vertex. set their vertex-info pointer to the new one
	for(const auto& locMapPair : locConstrainingDecayingParticles)
	{
		auto& locDefiningVertexInfo = locMapPair.second;
		locDefiningVertexInfo->Register_DecayingNoConstrainUseVertex(locMapPair.first, locVertexInfo);
	}

	if(!locLinkingFlag)
	{
		locVertexInfo->Set_DanglingVertexFlag(true);
		return false; //this is a dangling vertex, nothing further to do
	}

	//The positions of these decaying particles are now defined: Can use to constrain vertices in later constraints
	//since we need to match with particles in other constraints, save the OTHER index for the particle
		//if was in initial state, save final-state pair. and vice versa
	auto locReaction = locVertexInfo->Get_Reaction();
	for(auto locParticlePair : locNoConstrainDecayingParticles)
	{
		if(locParticlePair.second < 0) //was in initial state: save final state
			locDefinedDecayingParticles.emplace(DAnalysis::Get_InitialParticleDecayFromIndices(locReaction, locParticlePair.first), locVertexInfo);
		else //was in final state: save initial state
		{
			int locDecayStepIndex = DAnalysis::Get_DecayStepIndex(locReaction, locParticlePair.first, locParticlePair.second);
			locDefinedDecayingParticles.emplace(std::make_pair(locDecayStepIndex, DReactionStep::Get_ParticleIndex_Initial()), locVertexInfo);
		}
	}

	return true;
}

