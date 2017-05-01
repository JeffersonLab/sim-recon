#include "ANALYSIS/DSourceComboVertexer.h"
#include "ANALYSIS/DSourceComboer.h"

namespace DAnalysis
{

DSourceComboVertexer::DSourceComboVertexer(JEventLoop* locEventLoop, DSourceComboer* locSourceComboer, DSourceComboP4Handler* locSourceComboP4Handler) :
dSourceComboer(locSourceComboer), dSourceComboP4Handler(locSourceComboP4Handler)
{
	locEventLoop->GetSingle(dAnalysisUtilities);

	//GET THE GEOMETRY
	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry* locGeometry = locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());

	//TARGET INFORMATION
	double locTargetCenterZ = 65.0;
	locGeometry->GetTargetZ(locTargetCenterZ);
	dTargetCenter.SetXYZ(0.0, 0.0, locTargetCenterZ);
}

signed char DSourceComboVertexer::Get_VertexZBin(bool locIsProductionVertex, const DSourceCombo* locSourceCombo) const
{
	auto locConstrainingParticles = Get_ConstrainingParticles(locIsProductionVertex, locSourceCombo);
	if(locConstrainingParticles.empty())
		return DSourceComboInfo::Get_VertexZIndex_Unknown();
	auto locVertexZ = Get_Vertex(locIsProductionVertex, locConstrainingParticles).Z();
	return dSourceComboer->Get_PhotonVertexZBin(locVertexZ);
}

void DSourceComboVertexer::Calc_VertexTimeOffsets_WithCharged(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionChargedCombo)
{
	auto locIsPrimaryProductionVertex = locReactionVertexInfo->Get_StepVertexInfos().front()->Get_ProductionVertexFlag();
	if(dTimeOffsets.find(std::make_pair(locIsPrimaryProductionVertex, locReactionChargedCombo)) != dTimeOffsets.end())
		return; //already done!! //e.g. the same channel used for 2 different DReactions

	//loop through vertices
	map<pair<int, int>, const DKinematicData*> locReconDecayParticleMap; //decaying particle indices -> kinematic data //indices: when the decaying particle is in the INITIAL state
	unordered_map<const DReactionStepVertexInfo*, const DSourceCombo*> locVertexPrimaryComboMap;
	for(const auto& locStepVertexInfo : locReactionVertexInfo->Get_StepVertexInfos())
	{
		/***************************************** CHECK IF VERTEX POSITION IS INDETERMINATE A THIS STAGE ********************************************/

		auto locVertexPrimaryCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionChargedCombo, locStepVertexInfo);
		auto locIsProductionVertexFlag = locStepVertexInfo->Get_ProductionVertexFlag();
		auto locComboProductionPair = std::make_pair(locIsProductionVertexFlag, locVertexPrimaryCombo);

		if(locStepVertexInfo->Get_DanglingVertexFlag())
		{
			//is forever indeterminate, even with neutrals and beam energy
			//If this is the production vertex, choose the center of the target.  If not, choose the vertex where the decay parent was produced.
			auto locParentVertexInfo = locStepVertexInfo->Get_ParentVertexInfo();
			if(locParentVertexInfo == nullptr) //initial vertex
				dConstrainingParticlesByCombo.emplace(locComboProductionPair, {}); //empty: target center
			else //decay products
			{
				auto locParentCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionChargedCombo, locParentVertexInfo);
				auto locIsParentProductionVertex = locParentVertexInfo->Get_ProductionVertexFlag();
				auto& locVertexParticles = dConstrainingParticlesByCombo[std::make_pair(locIsParentProductionVertex, locParentCombo)];
				dConstrainingParticlesByCombo.emplace(locComboProductionPair, locVertexParticles);
			}
			continue;
		}

		//get combo & info
		locVertexPrimaryComboMap.emplace(locStepVertexInfo, locVertexPrimaryCombo);

		//get particles
		auto locChargedSourceParticles = DAnalysis::Get_SourceParticles_ThisVertex(locVertexPrimaryCombo);
		auto locDecayingParticles = Get_FullConstrainDecayingParticles(locStepVertexInfo, locReconDecayParticleMap);

		//now, just because it isn't dangling, it doesn't mean we have enough information to actually determine the vertex
		//e.g. we may need neutrals or beam energy to define constraining decay particle p4/trajectory
		auto locDeterminableIterator = dVertexDeterminableWithChargedMap.find(locComboProductionPair);
		if(locDeterminableIterator == dVertexDeterminableWithChargedMap.end())
		{
			//we don't know yet if it is determinable or not.  figure it out
			auto locNumConstrainingParticles = locChargedSourceParticles.size() + locDecayingParticles.size();

			//determine it, save result, and continue if can't
			bool locVertexDeterminableFlag = locIsProductionVertexFlag ? (locNumConstrainingParticles > 0) : (locNumConstrainingParticles > 1);
			dVertexDeterminableWithChargedMap.emplace(locComboProductionPair, locVertexDeterminableFlag);
			if(!locVertexDeterminableFlag)
				continue; //can't determine yet, but will be able to in the future
		}
		if(locDeterminableIterator->second == false)
			continue;

		vector<const DKinematicData*> locVertexParticles;
		auto locVertex = Calc_Vertex(locIsProductionVertexFlag, locChargedSourceParticles, locDecayingParticles, locVertexParticles);
		dConstrainingParticlesByCombo.emplace(locComboProductionPair, locVertexParticles);

		//set the vertices (really, the vertex particles) for the other combos for faster lookup during neutral comboing
		auto locVertexCombos = DAnalysis::Get_SourceCombos_ThisVertex(locVertexPrimaryCombo);
		for(auto& locVertexCombo : locVertexCombos)
		{
			auto locVertexParticles = dConstrainingParticlesByCombo[std::make_pair(locIsProductionVertexFlag, locVertexPrimaryCombo)];
			dConstrainingParticlesByCombo.emplace(std::make_pair(locIsProductionVertexFlag, locVertexCombo), locVertexParticles);
		}

		Construct_DecayingParticle_InvariantMass(locStepVertexInfo, locVertexPrimaryCombo, nullptr, locVertex, locReconDecayParticleMap);
	}

	//do time offsets once all the vertices have been found
	Calc_TimeOffsets(locReactionVertexInfo, locReactionChargedCombo, nullptr);
}

void DSourceComboVertexer::Calc_VertexTimeOffsets_WithPhotons(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionChargedCombo, const DSourceCombo* locReactionFullCombo)
{
	//Calculate vertex positions & time offsets using photons
	//not likely to have any effect, but it's necessary sometimes (but rarely)
	//E.g. g, p ->  K0, Sigma+    K0 -> 3pi: The selected pi0 photons could help define the production vertex

	auto locIsPrimaryProductionVertex = locReactionVertexInfo->Get_StepVertexInfos().front()->Get_ProductionVertexFlag();
	if(dTimeOffsets.find(std::make_pair(locIsPrimaryProductionVertex, locReactionFullCombo)) != dTimeOffsets.end())
		return; //already done!! //e.g. the same channel used for 2 different DReactions

	//loop over vertices in dependency order
	map<pair<int, int>, const DKinematicData*> locReconDecayParticleMap; //decaying particle indices -> kinematic data //indices: when the decaying particle is in the INITIAL state
	for(const auto& locStepVertexInfo : locReactionVertexInfo->Get_StepVertexInfos())
	{
		if(locStepVertexInfo->Get_DanglingVertexFlag())
			continue; //is forever indeterminate, even with neutrals and beam energy

		auto locIsProductionVertexFlag = locStepVertexInfo->Get_ProductionVertexFlag();
		auto locVertexPrimaryChargedCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionChargedCombo, locStepVertexInfo);
		auto locVertexPrimaryFullCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionFullCombo, locStepVertexInfo);
		auto locComboProductionPair = std::make_pair(locIsProductionVertexFlag, locVertexPrimaryFullCombo);

		//see if vertex has already been found
		auto locChargedComboProductionPair = std::make_pair(locIsProductionVertexFlag, locVertexPrimaryChargedCombo);
		auto locFullComboProductionPair = std::make_pair(locIsProductionVertexFlag, locVertexPrimaryFullCombo);
		if(dConstrainingParticlesByCombo.find(locChargedComboProductionPair) != dConstrainingParticlesByCombo.end())
		{
			//already done! get/construct any recon decaying particles
			auto locVertex = Get_Vertex(locIsProductionVertexFlag, locVertexPrimaryChargedCombo);
			auto locVertexParticles = dConstrainingParticlesByCombo[locChargedComboProductionPair];
			dConstrainingParticlesByCombo.emplace(locFullComboProductionPair, locVertexParticles); //so that the vertex can be retrieved by either the charged or full combo
			Construct_DecayingParticle_InvariantMass(locStepVertexInfo, locVertexPrimaryChargedCombo, locVertexPrimaryFullCombo, locVertex, locReconDecayParticleMap);
			continue;
		}

		//get particles
		auto locChargedSourceParticles = DAnalysis::Get_SourceParticles_ThisVertex(locVertexPrimaryChargedCombo);
		auto locDecayingParticles = Get_FullConstrainDecayingParticles(locStepVertexInfo, locReconDecayParticleMap);

		//now, just because it isn't dangling, it doesn't mean we have enough information to actually determine the vertex
		//e.g. we may need beam energy to define constraining decay particle p4/trajectory
		auto locDeterminableIterator = dVertexDeterminableWithPhotonsMap.find(locComboProductionPair);
		if(locDeterminableIterator == dVertexDeterminableWithPhotonsMap.end())
		{
			//we don't know yet if it is determinable or not.  figure it out
			auto locNumConstrainingParticles = locChargedSourceParticles.size() + locDecayingParticles.size();

			//determine it, save result, and continue if can't
			bool locVertexDeterminableFlag = locIsProductionVertexFlag ? (locNumConstrainingParticles > 0) : (locNumConstrainingParticles > 1);
			dVertexDeterminableWithPhotonsMap.emplace(locComboProductionPair, locVertexDeterminableFlag);
			if(!locVertexDeterminableFlag)
				continue; //can't determine yet, but will be able to in the future
		}
		if(locDeterminableIterator->second == false)
			continue;

		//find the vertex, save the results
		vector<const DKinematicData*> locVertexParticles;
		auto locVertex = Calc_Vertex(locIsProductionVertexFlag, locChargedSourceParticles, locDecayingParticles, locVertexParticles);
		dConstrainingParticlesByCombo.emplace(locChargedComboProductionPair, locVertexParticles); //so that the vertex can be retrieved via the charged combo
		dConstrainingParticlesByCombo.emplace(locFullComboProductionPair, locVertexParticles);

		Construct_DecayingParticle_InvariantMass(locStepVertexInfo, locVertexPrimaryChargedCombo, locVertexPrimaryFullCombo, locVertex, locReconDecayParticleMap);
	}

	//CALC TIME OFFSETS
	Calc_TimeOffsets(locReactionVertexInfo, locReactionChargedCombo, locReactionFullCombo);
}

void DSourceComboVertexer::Calc_VertexTimeOffsets_WithBeam(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionChargedCombo, const DSourceCombo* locReactionFullCombo, const DKinematicData* locBeamParticle)
{
	//IF PRIMARY VERTEX IS UNKNOWN THEN DO NOTHING!
	if(!Get_IsVertexFoundFlag(true, locReactionChargedCombo, locReactionFullCombo))
		return;
	auto locPrimaryVertexZ = Get_PrimaryVertex(locReactionVertexInfo, locReactionChargedCombo).Z();

	int locRFBunch = dSourceComboTimeHandler->Calc_RFBunchShift(locBeamParticle->time());

	//uses the beam to define remaining vertices
	//loop over vertices in dependency order
	map<pair<int, int>, const DKinematicData*> locReconDecayParticleMap; //decaying particle indices -> kinematic data //indices: when the decaying particle is in the INITIAL state
	for(const auto& locStepVertexInfo : locReactionVertexInfo->Get_StepVertexInfos())
	{
		if(locStepVertexInfo->Get_DanglingVertexFlag())
			continue; //is forever indeterminate, even with beam energy

		auto locIsProductionVertexFlag = locStepVertexInfo->Get_ProductionVertexFlag();
		auto locVertexPrimaryChargedCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionChargedCombo, locStepVertexInfo);
		auto locVertexPrimaryFullCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionFullCombo, locStepVertexInfo);

		//see if vertex has already been found //can search with either charged or full
		auto locChargedComboProductionPair = std::make_pair(locIsProductionVertexFlag, locVertexPrimaryChargedCombo);
		if(dConstrainingParticlesByCombo.find(locChargedComboProductionPair) != dConstrainingParticlesByCombo.end())
		{
			//already done! get/construct any recon decaying particles
			auto locVertex = Get_Vertex(locIsProductionVertexFlag, locVertexPrimaryChargedCombo);
			Construct_DecayingParticle_InvariantMass(locStepVertexInfo, locVertexPrimaryChargedCombo, locVertexPrimaryFullCombo, locVertex, locReconDecayParticleMap);
			if(locIsProductionVertexFlag) //do missing if needed
			{
				auto locRFVertexTime = dSourceComboTimeHandler->Calc_PropagatedRFTime(locPrimaryVertexZ, locRFBunch, 0.0);
				Construct_DecayingParticle_MissingMass(locStepVertexInfo, locReactionChargedCombo, locReactionFullCombo, locVertexPrimaryFullCombo, locBeamParticle, locVertex, locRFBunch, locRFVertexTime, locReconDecayParticleMap);
			}
			continue;
		}

		auto locFullComboProductionPair = std::make_pair(locIsProductionVertexFlag, locVertexPrimaryFullCombo);

		//get particles
		auto locChargedSourceParticles = DAnalysis::Get_SourceParticles_ThisVertex(locVertexPrimaryChargedCombo);
		auto locDecayingParticles = Get_FullConstrainDecayingParticles(locStepVertexInfo, locReconDecayParticleMap);

		//find the vertex, save the results
		vector<const DKinematicData*> locVertexParticles;
		auto locVertex = Calc_Vertex(locIsProductionVertexFlag, locChargedSourceParticles, locDecayingParticles, locVertexParticles);
		dConstrainingParticlesByCombo.emplace(locChargedComboProductionPair, locVertexParticles); //so that the vertex can be retrieved via the charged combo
		dConstrainingParticlesByCombo.emplace(locFullComboProductionPair, locVertexParticles);

		//CALC AND STORE TIME OFFSET
		auto locChargedReactionPair_TimeOffset = std::make_pair(false, locReactionChargedCombo);
		auto locFullReactionPair_TimeOffset = std::make_pair(false, locReactionFullCombo);

		//get parent
		auto locParentVertexInfo = locStepVertexInfo->Get_ParentVertexInfo();
		auto locParentFullCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionFullCombo, locParentVertexInfo);
		auto locParentTimeOffset = dTimeOffsets[locFullReactionPair_TimeOffset][locVertexPrimaryFullCombo];

		//get vertices & path length
		auto locParentProductionVertex = Get_Vertex(locParentVertexInfo->Get_ProductionVertexFlag(), locParentFullCombo);
		auto locPathLength = (locVertex - locParentProductionVertex).Mag();

		//compute and save result
		auto locP4 = locDecayingParticles.back()->lorentzMomentum(); //when locDecayingParticles retrieved, this particle was saved to the back!
		auto locTimeOffset = locPathLength/(locP4.Beta()*SPEED_OF_LIGHT) + locParentTimeOffset;

		dTimeOffsets[locChargedReactionPair_TimeOffset].emplace(locVertexPrimaryChargedCombo, locTimeOffset);
		dTimeOffsets[locFullReactionPair_TimeOffset].emplace(locVertexPrimaryFullCombo, locTimeOffset);

		//construct decaying particles by missing mass
		auto locRFVertexTime = dSourceComboTimeHandler->Calc_PropagatedRFTime(locPrimaryVertexZ, locRFBunch, locTimeOffset);
		Construct_DecayingParticle_MissingMass(locStepVertexInfo, locReactionChargedCombo, locReactionFullCombo, locVertexPrimaryFullCombo, locBeamParticle, locVertex, locRFBunch, locRFVertexTime, locReconDecayParticleMap);
	}
}

DVector3 DSourceComboVertexer::Calc_Vertex(bool locIsProductionVertexFlag, const vector<pair<Particle_t, const JObject*>>& locChargedSourceParticles, const vector<const DKinematicData*>& locDecayingParticles, vector<const DKinematicData*>& locVertexParticles)
{

	/************************************************** CALCULATING VERTEX POSITIONS ***********************************************
	*
	* Production vertex:
	* 1) If there is at least one charged track at the production vertex with a theta > 30 degrees:
	*    The production vertex is the POCA to the beamline of the track with the largest theta.
	* 2) If not, then for each detached vertex:
	*    a) If there are any neutral or missing particles, or there are < 2 detected charged tracks at the detached vertex: ignore it
	*    b) Otherwise:
	*       i) The detached vertex is at the center of the lines between the POCAs of the two closest tracks.
	*       ii) Calculate the p3 of the decaying particles at this point and then find the POCA of the decaying particle to the beamline.
	* 3) Now, the production vertex is the POCA to the beamline of the track/decaying-particle with the largest theta.
	* 4) Otherwise, the production vertex is the center of the target.
	*
	* Detached vertices:
	* 1) If at least 2 decay products have defined trajectories (detected & charged or decaying & reconstructed):
	*    The detached vertex is at the center of the lines between the POCAs of the two closest particles.
	* 2) If one decay product is detected & charged, and the decaying particle production vertex was well defined (i.e. not forced to be center of target):
	*    a) Determine the decaying particle trajectory from the missing mass of the system (must be done after beam photon selection!!)
	*    b) The detached vertex is at the POCA of the charged decay product and the decaying particle
	* 3) Otherwise use the decaying particle production vertex for its position.
	*
	*******************************************************************************************************************************/

	//Get detected charged track hypotheses at this vertex
	locVertexParticles.clear();
	locVertexParticles.reserve(locChargedSourceParticles.size());

	auto Get_Hypothesis = [](const pair<Particle_t, const JObject*>& locPair) -> const DChargedTrackHypothesis*
		{return static_cast<const DChargedTrack*>(locPair.second)->Get_Hypothesis(locPair.first);};
	std::transform(locChargedSourceParticles.begin(), locChargedSourceParticles.end(), std::back_inserter(locVertexParticles), Get_Hypothesis);
	locVertexParticles.insert(locVertexParticles.end(), locDecayingParticles.begin(), locDecayingParticles.end());

	if(locIsProductionVertexFlag)
	{
		//use track with theta nearest 90 degrees
		auto locThetaNearest90Iterator = Get_ThetaNearest90Iterator(locVertexParticles);
		double locThetaNearest90 = (*locThetaNearest90Iterator)->momentum().Theta()*180.0/TMath::Pi();
		if(locThetaNearest90 < dMinThetaForVertex)
		{
			//try decaying particles instead
			auto locThetaNearest90Iterator_Decaying = Get_ThetaNearest90Iterator(locDecayingParticles);
			double locLargestTheta_Decaying = (*locThetaNearest90Iterator_Decaying)->momentum().Theta()*180.0/TMath::Pi();
			if(locLargestTheta_Decaying > locThetaNearest90)
				locThetaNearest90Iterator = locThetaNearest90Iterator_Decaying;
		}
		locVertexParticles = {*locThetaNearest90Iterator};
		auto locVertex = locVertexParticles[0]->position();
		dVertexMap.emplace(std::make_pair(locIsProductionVertexFlag, locVertexParticles), locVertex);
		return locVertex;
	}

	//detached vertex
	if(locVertexParticles.size() < 2)
		locVertexParticles.insert(locVertexParticles.end(), locDecayingParticles.begin(), locDecayingParticles.end());
	std::sort(locVertexParticles.begin(), locVertexParticles.end());

	//if vertex already computed, don't bother
	auto locVertexIterator = dVertexMap.find(std::make_pair(locIsProductionVertexFlag, locVertexParticles));
	if(locVertexIterator == dVertexMap.end())
	{
		auto locVertex = dAnalysisUtilities->Calc_CrudeVertex(locVertexParticles);
		dVertexMap.emplace(std::make_pair(locIsProductionVertexFlag, locVertexParticles), locVertex);
		return locVertex;
	}

	return locVertexIterator->second;
}

vector<const DKinematicData*> DSourceComboVertexer::Get_FullConstrainDecayingParticles(const DReactionStepVertexInfo* locStepVertexInfo, const map<pair<int, int>, const DKinematicData*>& locReconDecayParticleMap)
{
	auto locConstrainingDecayingParticles = locStepVertexInfo->Get_DecayingParticles_FullConstrain();
	vector<const DKinematicData*> locDecayingParticles;
	pair<int, int> locDecayMissingMassReconPair(-99, -99); //if there is one (and there will be at most one), save for last (easy to retrieve)
	for(auto locDecayPair : locConstrainingDecayingParticles)
	{
		//the particle pair in the map is where the decaying particle was when it was created (in the initial/final state for invariant/missing mass)
		//however, the pairs in locConstrainingDecayingParticles are when then the decaying particles are USED (the other state)
		//so, we have to convert from final state to initial state (or vice versa) to do the map lookup
		pair<int, int> locParticlePair(-99, -99);
		if(locDecayPair.first.second >= 0) //need it in the final state, was created in the initial state: convert from final to initial for map lookup
		{
			int locDecayStepIndex = DAnalysis::Get_DecayStepIndex(locStepVertexInfo->Get_Reaction(), locDecayPair.first.first, locDecayPair.first.second);
			locParticlePair = std::make_pair(locDecayStepIndex, DReactionStep::Get_ParticleIndex_Initial());
		}
		else //need it in the initial state, was created in the final state: convert: convert from initial to final for map lookup
		{
			locDecayMissingMassReconPair = DAnalysis::Get_InitialParticleDecayFromIndices(locStepVertexInfo->Get_Reaction(), locDecayPair.first.first);
			continue;
		}

		auto locReconIterator = locReconDecayParticleMap.find(locParticlePair);
		if(locReconIterator != locReconDecayParticleMap.end())
			locDecayingParticles.push_back(locReconIterator->second);
	}

	//retrieve the decay missing-mass recon particle, if there was one
	auto locReconIterator = locReconDecayParticleMap.find(locDecayMissingMassReconPair);
	if(locReconIterator != locReconDecayParticleMap.end())
		locDecayingParticles.push_back(locReconIterator->second);

	return locDecayingParticles;
}

void DSourceComboVertexer::Construct_DecayingParticle_InvariantMass(const DReactionStepVertexInfo* locReactionStepVertexInfo, const DSourceCombo* locChargedVertexCombo, const DSourceCombo* locFullVertexCombo, DVector3 locVertex, map<pair<int, int>, const DKinematicData*>& locReconDecayParticleMap)
{
	//we also can't compute the p4 yet if the decay products contain neutrals: this is done before comboing them!
	auto locSourceComboUse = dSourceComboer->Get_SourceComboUse(locReactionStepVertexInfo);
	if(locFullVertexCombo == nullptr) //charged only
	{
		if(dSourceComboer->Get_ChargeContent(std::get<2>(locSourceComboUse)) != d_Charged)
			return;
		locFullVertexCombo = locChargedVertexCombo;
	}
	else if(dSourceComboer->Get_HasMassiveNeutrals(std::get<2>(locSourceComboUse)))
	{
		//can't compute the p4 yet if the decay products contain massive neutrals: time offset is needed for massive neutral p4, but isn't known yet because decay p4 unknown!
		//will have to be computed via missing mass instead (although it technically isn't needed for it)
		return;
	}

	auto locReaction = locReactionStepVertexInfo->Get_Reaction();
	auto locIsProductionVertexFlag = locReactionStepVertexInfo->Get_ProductionVertexFlag();
	auto locChargedComboProductionPair = std::make_pair(locIsProductionVertexFlag, locChargedVertexCombo);

	//loop over decaying no-constrain decaying particles
	map<pair<int, int>, weak_ptr<DReactionStepVertexInfo>> locNoConstrainDecayingParticles = locReactionStepVertexInfo->Get_DecayingParticles_NoConstrain();
	for(const auto& locNoConstrainPair : locNoConstrainDecayingParticles)
	{
		//we cannot define decaying p4 via missing mass, because the beam is not chosen yet
		//therefore, if the no-constrain decaying particle is a final-state particle at this vertex, we cannot define its p4
		const auto& locParticlePair = locNoConstrainPair.first;
		if(locParticlePair.second >= 0)
			continue; //momentum must be determined by missing mass, cannot do yet! (or isn't detached, for which we don't care)

		//only create kinematic data for detached PIDs
		auto locDecayPID = std::get<0>(locSourceComboUse);
		if(!IsDetachedVertex(locDecayPID))
			continue; //either not detached (we don't care!), or is Unknown (missing decay product: can't compute)

		//if the particle has already been created, reuse it!
		auto locIterator = dReconDecayParticles_FromProducts.find(locSourceComboUse);
		if(locIterator != dReconDecayParticles_FromProducts.end())
		{
			locReconDecayParticleMap.emplace(locSourceComboUse, locIterator->second);
			continue;
		}

		//create a new one
		auto locVertexZBin = dSourceComboer->Get_PhotonVertexZBin(locVertex.Z());
		auto locP4 = dSourceComboP4Handler->Calc_P4_NoMassiveNeutrals(locFullVertexCombo, locVertexZBin);
		auto locKinematicData = new DKinematicData(locDecayPID, locP4.Vect(), locVertex, 0.0);

		//register it
		locReconDecayParticleMap.emplace(locParticlePair, locKinematicData);
		dReconDecayParticles_FromProducts.emplace(locSourceComboUse, locKinematicData);
	}
}

void DSourceComboVertexer::Construct_DecayingParticle_MissingMass(const DReactionStepVertexInfo* locReactionStepVertexInfo, const DSourceCombo* locReactionChargedCombo, const DSourceCombo* locReactionFullCombo, const DSourceCombo* locFullVertexCombo, const DKinematicData* locBeamParticle, DVector3 locVertex, int locRFBunch, double locRFVertexTime, map<pair<int, int>, const DKinematicData*>& locReconDecayParticleMap)
{
	//we can only calculate up to 1 at a time
	//the input full combo contains the decaying particle for which the missing mass is to be computed
	//if there is more than one, then this is impossible

	auto locSourceComboUse = dSourceComboer->Get_SourceComboUse(locReactionStepVertexInfo);
	auto locReaction = locReactionStepVertexInfo->Get_Reaction();
	auto locIsProductionVertexFlag = locReactionStepVertexInfo->Get_ProductionVertexFlag();
	auto locFullComboProductionPair = std::make_pair(locIsProductionVertexFlag, locFullVertexCombo);

	//ASSUMES FIXED TARGET EXPERIMENT!
	DLorentzVector locInitialStateP4; //lookup or is beam + target
	if(locIsProductionVertexFlag)
	{
		Particle_t locPID = locReaction->Get_ReactionStep(0)->Get_TargetPID();
		locInitialStateP4 = locBeamParticle->lorentzMomentum() + DLorentzVector(0.0, 0.0, 0.0, ParticleMass(locPID));
	}
	else //already computed, look it up!
		locInitialStateP4 = dReconDecayParticles_FromMissing[std::make_pair(locReactionFullCombo, locSourceComboUse)]->lorentzMomentum();

	//loop over decaying no-constrain decaying particles, figure out how many we have to do this for (can't do more than 1!)
	map<pair<int, int>, weak_ptr<DReactionStepVertexInfo>> locNoConstrainDecayingParticles = locReactionStepVertexInfo->Get_DecayingParticles_NoConstrain();
	pair<int, int> locToReconParticleIndices(-99, -99);
	for(const auto& locNoConstrainPair : locNoConstrainDecayingParticles)
	{
		//we are defining p4 through missing mass. therefore, the decaying particle must be in the final state at this vertex
		//therefore, if the no-constrain decaying particle is a final-state particle at this vertex, we cannot define its p4
		const auto& locParticlePair = locNoConstrainPair.first;
		if(locParticlePair.second < 0)
			continue; //momentum already determined by invariant mass!

		//the particle must be a constraining particle at another vertex, or else we aren't interested (i.e. detached)
		auto locReactionStep = locReaction->Get_ReactionStep(locParticlePair.first);
		auto locDecayPID = locReactionStep->Get_PID(locParticlePair.second);
		if(!IsDetachedVertex(locDecayPID))
			continue; //not detached: we don't care!

		if(locToReconParticleIndices.first >= 0)
			return; //can't recon more than 1 decaying particle via missing mass!
		locToReconParticleIndices = locParticlePair;
	}

	auto locDecayStepIndex = Get_DecayStepIndex(locReaction, locToReconParticleIndices.first, locToReconParticleIndices.second);
	auto locDecayUse = dSourceComboer->Get_SourceComboUse(locReaction, locDecayStepIndex);

	//if the particle has already been created, reuse it!
	auto locMissingPair = std::make_pair(locReactionFullCombo, locDecayUse);
	auto locIterator = dReconDecayParticles_FromMissing.find(locMissingPair);
	if(locIterator != dReconDecayParticles_FromMissing.end())
	{
		locReconDecayParticleMap.emplace(locToReconParticleIndices, locIterator->second);
		return;
	}

	//make sure there isn't another missing particle
	for(auto locStepIndex : locReactionStepVertexInfo->Get_StepIndices())
	{
		if(locStepIndex == locDecayStepIndex)
			continue;
		if(DAnalysis::Get_HasMissingParticle_FinalState(locReaction->Get_ReactionStep(locStepIndex)))
			return; //missing decay product! can't compute 2 missing p4's at once!
	}

	//create a new one
	//calc final state p4
	auto locVertexZBin = dSourceComboer->Get_PhotonVertexZBin(locVertex.Z());
	DLorentzVector locFinalStateP4;
	if(!dSourceComboP4Handler->Calc_P4_HasMassiveNeutrals(locIsProductionVertexFlag, locReactionChargedCombo, locFullVertexCombo, locVertexZBin, locVertex, locRFBunch, locRFVertexTime, locDecayUse, locFinalStateP4))
		return; //invalid somehow

	auto locReactionStep = locReaction->Get_ReactionStep(locToReconParticleIndices.first);
	auto locDecayPID = locReactionStep->Get_PID(locToReconParticleIndices.second);

	auto locP4 = locInitialStateP4 - locFinalStateP4;
	auto locKinematicData = new DKinematicData(locDecayPID, locP4.Vect(), locVertex, 0.0);

	//register it
	locReconDecayParticleMap.emplace(locToReconParticleIndices, locKinematicData);
	dReconDecayParticles_FromMissing.emplace(locMissingPair, locKinematicData);
}

void DSourceComboVertexer::Calc_TimeOffsets(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locChargedReactionCombo, const DSourceCombo* locFullReactionCombo)
{
	auto locIsPrimaryProductionVertex = locReactionVertexInfo->Get_StepVertexInfos().front()->Get_ProductionVertexFlag();
	auto locActiveReactionCombo = (locFullReactionCombo != nullptr) ? locFullReactionCombo : locChargedReactionCombo;

	auto locChargedReactionPair = std::make_pair(locIsPrimaryProductionVertex, locChargedReactionCombo);
	auto locNeutralReactionPair = std::make_pair(locIsPrimaryProductionVertex, locFullReactionCombo);

	auto& locChargedTimeOffsetMap = dTimeOffsets[locChargedReactionPair];
	auto& locFullTimeOffsetMap = dTimeOffsets[locNeutralReactionPair];

	//loop over vertices
	//MUST GO IN STEP ORDER!!
	for(auto locStepVertexInfo : locReactionVertexInfo->Get_StepVertexInfos())
	{
		auto locChargedVertexCombo = dSourceComboer->Get_VertexPrimaryCombo(locChargedReactionCombo, locStepVertexInfo);
		auto locFullVertexCombo = (locFullReactionCombo != nullptr) ? dSourceComboer->Get_VertexPrimaryCombo(locFullReactionCombo, locStepVertexInfo) : nullptr;
		auto locActiveVertexCombo = (locFullReactionCombo != nullptr) ? locFullVertexCombo : locChargedVertexCombo;

		//see if already determined
		auto locChargedIterator = locChargedTimeOffsetMap.find(locChargedVertexCombo);
		if(locChargedIterator != locChargedTimeOffsetMap.end())
		{
			if(locFullReactionCombo != nullptr) //save in full map
				locFullTimeOffsetMap.emplace(locFullVertexCombo, locChargedIterator->second);
			continue; //previously determined!
		}

		if(locStepVertexInfo->Get_ProductionVertexFlag())
		{
			locChargedTimeOffsetMap.emplace(locChargedReactionCombo, 0.0);
			continue;
		}

		//get parent information
		auto locParentVertexInfo = locStepVertexInfo->Get_ParentVertexInfo();
		auto locParentCombo = dSourceComboer->Get_VertexPrimaryCombo(locActiveReactionCombo, locParentVertexInfo);
		auto locParentTimeOffset = locChargedTimeOffsetMap[locChargedVertexCombo];

		//get vertices & path length
		auto locVertex = Get_Vertex(false, locActiveVertexCombo);
		auto locParentProductionVertex = Get_Vertex(locParentVertexInfo->Get_ProductionVertexFlag(), locParentCombo);
		auto locPathLength = (locVertex - locParentProductionVertex).Mag();

		//compute and save result
		auto locVertexZBin = Get_VertexZBin(false, locChargedVertexCombo);
		auto locP4 = dSourceComboP4Handler->Calc_P4_NoMassiveNeutrals(locActiveVertexCombo, locVertexZBin);
		auto locTimeOffset = locPathLength/(locP4.Beta()*SPEED_OF_LIGHT) + locParentTimeOffset;
		locChargedTimeOffsetMap.emplace(locChargedVertexCombo, locTimeOffset);
		if(locFullVertexCombo != nullptr)
			locFullTimeOffsetMap.emplace(locFullVertexCombo, locTimeOffset);
	}
}

} //end DAnalysis namespace

