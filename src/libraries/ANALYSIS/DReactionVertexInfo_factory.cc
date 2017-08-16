#include "ANALYSIS/DReactionVertexInfo_factory.h"

using namespace std;
using namespace jana;
using namespace DAnalysis;

jerror_t DReactionVertexInfo_factory::init(void)
{
	// This is so that the created objects persist throughout the life of the program instead of being cleared each event.
	SetFactoryFlag(PERSISTANT);
	dResourcePool_ReactionStepVertexInfo = new DResourcePool<DReactionStepVertexInfo>();

	gPARMS->SetDefaultParameter("VERTEXINFO:DEBUG_LEVEL", dDebugLevel);

	return NOERROR;
}

jerror_t DReactionVertexInfo_factory::evnt(JEventLoop *locEventLoop, uint64_t locEventNumber)
{
	dKinFitUtils = new DKinFitUtils_GlueX(locEventLoop);

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
		if(dDebugLevel > 0)
		{
			cout << "CREATED REACTION VERTEX INFO:" << endl;
			DAnalysis::Print_ReactionVertexInfo(locVertexInfo);
		}

		_data.push_back(locVertexInfo);
		dVertexInfoMap.emplace(locReaction, locVertexInfo);
	}

	return NOERROR;
}

DReactionVertexInfo* DReactionVertexInfo_factory::Build_VertexInfo(const DReaction* locReaction)
{
	set<size_t> locStepIndexSet;
	vector<DReactionStepVertexInfo*> locVertexInfos;
	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		if(locStepIndexSet.find(loc_i) != locStepIndexSet.end())
			continue; //already did this step

		//Start a new vertex, and save when done
		auto locVertexInfo = Setup_VertexInfo(locReaction, loc_i, nullptr);
		locVertexInfos.push_back(locVertexInfo);

		auto locStepIndexVector = locVertexInfo->Get_StepIndices();
		std::copy(locStepIndexVector.begin(), locStepIndexVector.end(), std::inserter(locStepIndexSet, locStepIndexSet.end()));

		Group_VertexParticles(locVertexInfo);
		if(dDebugLevel > 0)
		{
			cout << "Grouped vertex info: " << endl;
			DAnalysis::Print_ReactionStepVertexInfo(locVertexInfo);
		}
	}

	//work in reverse: try to build decaying particles out of decay products first, rather than missing mass
	std::reverse(std::begin(locVertexInfos), std::end(locVertexInfos));

	//results are sorted by dependency order
	locVertexInfos = Link_Vertices(locReaction, locVertexInfos, false);
	if(dDebugLevel > 0)
	{
		cout << "Linked vertex infos, fitflag = false: " << endl;
		for(auto& locVertexInfo : locVertexInfos)
			DAnalysis::Print_ReactionStepVertexInfo(locVertexInfo);
	}
	locVertexInfos = Link_Vertices(locReaction, locVertexInfos, true);
	return new DReactionVertexInfo(locReaction, locVertexInfos);
}

DReactionStepVertexInfo* DReactionVertexInfo_factory::Setup_VertexInfo(const DReaction* locReaction, size_t locStepIndex, DReactionStepVertexInfo* locVertexInfo)
{
	//create/update vertex info
	if(locVertexInfo == nullptr)
	{
		locVertexInfo = dResourcePool_ReactionStepVertexInfo->Get_Resource();
		locVertexInfo->Reset();
		locVertexInfo->Set_Members(locReaction, locStepIndex);
	}
	else
		locVertexInfo->Add_ReactionStep(locStepIndex);

	//loop over final particles: add to the vertex constraint, dive through decaying particles that decay in-place
		//if decaying in-place: don't add (would add to constraint, but not here (purely internal))
	auto locReactionStep = locReaction->Get_ReactionStep(locStepIndex);
	if(locReactionStep->Get_IsInclusiveFlag())
		locVertexInfo->Set_IsInclusiveVertexFlag(true);
	for(size_t loc_i = 0; loc_i < locReactionStep->Get_NumFinalPIDs(); ++loc_i)
	{
		//check if particle decays, and if so, if in-place
		auto locDecayStepIndex = Get_DecayStepIndex(locReaction, locStepIndex, loc_i);
		auto locPID = locReactionStep->Get_FinalPID(loc_i);
		if((locDecayStepIndex > 0) && !IsDetachedVertex(locPID)) //yes: combine with the decay products
			Setup_VertexInfo(locReaction, locDecayStepIndex, locVertexInfo);
	}

	return locVertexInfo;
}

void DReactionVertexInfo_factory::Group_VertexParticles(DReactionStepVertexInfo* locVertexInfo)
{
	auto locReaction = locVertexInfo->Get_Reaction();
	auto locStepIndices = locVertexInfo->Get_StepIndices();
	vector<pair<int, int>> locDecayingParticles, locOnlyConstrainTimeParticles;
	vector<pair<int, int>> locFullConstrainParticles_Fit, locNoConstrainParticles_Fit, locFullConstrainParticles_Recon, locNoConstrainParticles_Recon;

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
			{
				locNoConstrainParticles_Fit.emplace_back(locStepIndex, DReactionStep::Get_ParticleIndex_Initial());
				locNoConstrainParticles_Recon.emplace_back(locStepIndex, DReactionStep::Get_ParticleIndex_Initial());
			}
			else
			{
				locFullConstrainParticles_Recon.emplace_back(locStepIndex, DReactionStep::Get_ParticleIndex_Initial());
				if(dKinFitUtils->Get_IncludeBeamlineInVertexFitFlag())
					locFullConstrainParticles_Fit.emplace_back(locStepIndex, DReactionStep::Get_ParticleIndex_Initial());
				else //not in fit!
					locNoConstrainParticles_Fit.emplace_back(locStepIndex, DReactionStep::Get_ParticleIndex_Initial());
			}

			//beam 2
			if(locStep->Get_SecondBeamPID() != Unknown)
			{
				if(locStep->Get_IsBeamMissingFlag())
				{
					locNoConstrainParticles_Fit.emplace_back(locStepIndex, DReactionStep::Get_ParticleIndex_SecondBeam());
					locNoConstrainParticles_Recon.emplace_back(locStepIndex, DReactionStep::Get_ParticleIndex_SecondBeam());
				}
				else
				{
					locFullConstrainParticles_Recon.emplace_back(locStepIndex, DReactionStep::Get_ParticleIndex_SecondBeam());
					if(dKinFitUtils->Get_IncludeBeamlineInVertexFitFlag())
						locFullConstrainParticles_Fit.emplace_back(locStepIndex, DReactionStep::Get_ParticleIndex_SecondBeam());
					else //not in fit!
						locNoConstrainParticles_Fit.emplace_back(locStepIndex, DReactionStep::Get_ParticleIndex_SecondBeam());
				}
			}
		}
		else if(locFirstStepFlag) //decaying //if false: not detached: only save once at this vertex (in final state), not twice!!
			locDecayingParticles.emplace_back(locStepIndex, DReactionStep::Get_ParticleIndex_Initial());

		//target
		if(locStep->Get_TargetPID() != Unknown) //target
		{
			locNoConstrainParticles_Fit.emplace_back(locStepIndex, DReactionStep::Get_ParticleIndex_Target());
			locNoConstrainParticles_Recon.emplace_back(locStepIndex, DReactionStep::Get_ParticleIndex_Target());
		}

		//final state
		auto locFinalPIDs = locStep->Get_FinalPIDs();
		for(size_t loc_i = 0; loc_i < locFinalPIDs.size(); ++loc_i)
		{
			int locDecayStepIndex = Get_DecayStepIndex(locReaction, locStepIndex, loc_i);

			if(locDecayStepIndex >= 0) //decaying
				locDecayingParticles.emplace_back(locStepIndex, loc_i);
			else if(int(loc_i) == locStep->Get_MissingParticleIndex()) //missing
			{
				locNoConstrainParticles_Fit.emplace_back(locStepIndex, loc_i);
				locNoConstrainParticles_Recon.emplace_back(locStepIndex, loc_i);
			}
			else if(ParticleCharge(locFinalPIDs[loc_i]) != 0) //detected charged
			{
				locFullConstrainParticles_Fit.emplace_back(locStepIndex, loc_i);
				locFullConstrainParticles_Recon.emplace_back(locStepIndex, loc_i);
			}
			else if(locFinalPIDs[loc_i] == Gamma) //photon
				locOnlyConstrainTimeParticles.emplace_back(locStepIndex, loc_i);
			else //massive neutrals can't constrain position or time!
			{
				locNoConstrainParticles_Fit.emplace_back(locStepIndex, loc_i);
				locNoConstrainParticles_Recon.emplace_back(locStepIndex, loc_i);
			}
		}

		locFirstStepFlag = false;
	}

	locVertexInfo->Set_ParticleIndices(true, locFullConstrainParticles_Fit, locDecayingParticles, locOnlyConstrainTimeParticles, locNoConstrainParticles_Fit);
	locVertexInfo->Set_ParticleIndices(false, locFullConstrainParticles_Recon, locDecayingParticles, locOnlyConstrainTimeParticles, locNoConstrainParticles_Recon);
}

vector<DReactionStepVertexInfo*> DReactionVertexInfo_factory::Link_Vertices(const DReaction* locReaction, vector<DReactionStepVertexInfo*> locVertexInfos, bool locFitFlag) const
{
	//loop over vertex-constraints-to-sort:
		//find which constraints decaying particles should be defined-by/constrained-to
		//find order in which constraints need to be constrained

	bool locProgessMadeFlag = false;
	bool locTryMissingParticleVertexFlag = false;
	auto locVertexIterator = locVertexInfos.begin();
	map<pair<int, int>, DReactionStepVertexInfo*> locDefinedDecayingParticles;
	vector<DReactionStepVertexInfo*> locSortedVertexInfos; //sorted by dependency
	while(!locVertexInfos.empty())
	{
		if(locVertexIterator == locVertexInfos.end())
		{
			//made a full loop through
			if(!locProgessMadeFlag)
			{
				if(locTryMissingParticleVertexFlag)
					break; //no progress made: cannot constrain remaining vertices
				locTryMissingParticleVertexFlag = true; //try this now
				if(dDebugLevel > 0)
					cout << "try another loop, using missing mass" << endl;
			}

			//reset for next pass through
			locVertexIterator = locVertexInfos.begin();
			locProgessMadeFlag = false;
			continue;
		}

		//Try a vertex info
		auto locVertexInfo = *locVertexIterator;
		if(!locTryMissingParticleVertexFlag && !locVertexInfo->Get_MissingParticles().empty())
		{
			++locVertexIterator; //try again later
			continue; //first try to reconstruct decaying particles without using the missing mass
		}

		locProgessMadeFlag = Associate_DecayingParticles(locFitFlag, true, locVertexInfo, locDefinedDecayingParticles);
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
		Associate_DecayingParticles(locFitFlag, false, locVertexInfo, locDefinedDecayingParticles);
		locSortedVertexInfos.push_back(locVertexInfo);
	}

	//set parent vertex infos
	unordered_map<size_t, DReactionStepVertexInfo*> locVertexInfoMap;
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

bool DReactionVertexInfo_factory::Associate_DecayingParticles(bool locFitFlag, bool locLinkingFlag, DReactionStepVertexInfo* locVertexInfo, map<pair<int, int>, DReactionStepVertexInfo*>& locDefinedDecayingParticles) const
{
	if(dDebugLevel > 0)
	{
		cout << "Associate_DecayingParticles: fit flag, link flag, info, #defined decaying: " << locFitFlag << ", " << locLinkingFlag << ", " << locVertexInfo << ", " << locDefinedDecayingParticles.size() << endl;
		DAnalysis::Print_ReactionStepVertexInfo(locVertexInfo);
	}

	//find which decaying particles at this vertex have/haven't been previously defined
	vector<pair<int, int>> locNoConstrainDecayingParticles;
	map<pair<int, int>, const DReactionStepVertexInfo*> locConstrainingDecayingParticles;
	auto locDecayingParticles = locVertexInfo->Get_DecayingParticles();
	for(auto locParticlePair : locDecayingParticles)
	{
		auto locIterator = locDefinedDecayingParticles.find(locParticlePair);
		if(locIterator != locDefinedDecayingParticles.end())
			locConstrainingDecayingParticles.emplace(*locIterator);
		else //not found
			locNoConstrainDecayingParticles.emplace_back(locParticlePair);
	}
	if(dDebugLevel > 0)
		cout << "# decaying particles, constrain/no-constrain: " << locConstrainingDecayingParticles.size() << ", " << locNoConstrainDecayingParticles.size() << endl;

	//tricky: beamline can be used to find vertex, but not in a kinfit (because the errors are zero)
	//so, if a production vertex has only one charged/known-decaying track, then:
		//it's vertex position can be found: not dangling
		//it's vertex position cannot be fit (not enough constraints)
		//have a separate flag:  
	//enough tracks for kinfit vs. enough tracks to determine separate vertex

	//see if enough tracks //if not linking, then don't need "enough": will register as dangling vertex instead
	auto locNumConstrainingParticles = locConstrainingDecayingParticles.size() + locVertexInfo->Get_FullConstrainParticles(locFitFlag).size();
	if(locFitFlag)
		locVertexInfo->Set_FittableVertexFlag((locNumConstrainingParticles >= 2));
	if(dDebugLevel > 0)
		cout << "#constraining: " << locNumConstrainingParticles << endl;

	bool locEnoughTracksFlag = (locNumConstrainingParticles >= 2);
	if(locLinkingFlag && !locEnoughTracksFlag)
		return false; //trying to link vertices, not enough tracks yet. try again later

	//Save results
	locVertexInfo->Register_DecayingParticleConstraints(locFitFlag, locNoConstrainDecayingParticles, locConstrainingDecayingParticles);
	//each of the constraining decaying particles was a no-constrain at another vertex. set their vertex-info pointer to the new one
	for(const auto& locMapPair : locConstrainingDecayingParticles)
	{
		auto locDefiningVertexInfo = const_cast<DReactionStepVertexInfo*>(locMapPair.second); //easier this way
		locDefiningVertexInfo->Register_DecayingNoConstrainUseVertex(locFitFlag, locMapPair.first, locVertexInfo);
		locDefinedDecayingParticles.erase(locMapPair.first); //can't use the same one twice
	}

	if(!locLinkingFlag)
	{
		if(!locFitFlag)
			locVertexInfo->Set_DanglingVertexFlag(true);
		return false; //this is a dangling vertex, nothing further to do
	}

	//The positions of these decaying particles are now defined: Can use to constrain vertices in later constraints
		//However, we can only do so if their p4 is defined as well!  It won't be if there's also a missing particle at this vertex
	//since we need to match with particles in other constraints, save the OTHER index for the particle
		//if was in initial state, save final-state pair. and vice versa
	if(locVertexInfo->Get_IsInclusiveVertexFlag() || !locVertexInfo->Get_MissingParticles().empty())
		return true; //missing particle: decaying p4 not defined! (have already tried defining it the other way (missing vs invariant) and couldn't)
	auto locReaction = locVertexInfo->Get_Reaction();
	for(auto locParticlePair : locNoConstrainDecayingParticles)
	{
		if(dDebugLevel > 0)
			cout << "defined decaying particle indices: " << locParticlePair.first << ", " << locParticlePair.second << endl;
		if(locParticlePair.second < 0) //was in initial state: save final state
		{
			auto locParticlePairToRegister = DAnalysis::Get_InitialParticleDecayFromIndices(locReaction, locParticlePair.first);
			if(dDebugLevel > 0)
				cout << "registering decaying particle indices: " << locParticlePairToRegister.first << ", " << locParticlePairToRegister.second << endl;
			locDefinedDecayingParticles.emplace(locParticlePairToRegister, locVertexInfo);
		}
		else //was in final state: save initial state
		{
			int locDecayStepIndex = DAnalysis::Get_DecayStepIndex(locReaction, locParticlePair.first, locParticlePair.second);
			auto locParticlePairToRegister = std::make_pair(locDecayStepIndex, DReactionStep::Get_ParticleIndex_Initial());
			if(dDebugLevel > 0)
				cout << "registering decaying particle indices: " << locParticlePairToRegister.first << ", " << locParticlePairToRegister.second << endl;
			locDefinedDecayingParticles.emplace(locParticlePairToRegister, locVertexInfo);
		}
	}

	return true;
}

