#include "ANALYSIS/DSourceComboVertexer.h"
#include "ANALYSIS/DSourceComboer.h"

namespace DAnalysis
{

DSourceComboVertexer::DSourceComboVertexer(JEventLoop* locEventLoop, const DSourceComboer* locSourceComboer, DSourceComboP4Handler* locSourceComboP4Handler) :
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

signed char DSourceComboVertexer::Get_VertexZBin(bool locIsProductionVertex, const DSourceCombo* locChargedCombo) const
{
	if(!Get_VertexDeterminableWithCharged(locIsProductionVertex, locChargedCombo))
		return DSourceComboInfo::Get_VertexZIndex_Unknown();
	auto locVertexZ = Get_Vertex(locIsProductionVertex, locChargedCombo).Z();
	return dSourceComboer->Get_PhotonVertexZBin(locVertexZ);
}

DVector3 DSourceComboVertexer::Get_PrimaryVertex(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionChargedCombo) const
{
	auto locIsProductionVertex = locReactionVertexInfo->Get_StepVertexInfo(0)->Get_ProductionVertexFlag();
	return Get_Vertex(locIsProductionVertex, locReactionChargedCombo);
}

void DSourceComboVertexer::Calc_VertexTimeOffsets(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionChargedCombo)
{
	auto locIsPrimaryProductionVertex = locReactionVertexInfo->Get_StepVertexInfos().front()->Get_ProductionVertexFlag();
	if(dTimeOffsets.find(std::make_pair(locIsPrimaryProductionVertex, locReactionChargedCombo)) != dTimeOffsets.end())
		return; //already done!! //e.g. the same channel used for 2 different DReactions

	//loop through vertices, determining initial guesses
	map<pair<int, int>, const DKinematicData*> locReconDecayParticleMap; //decaying particle indices -> kinematic data //indices: when the decaying particle is in the INITIAL state
	unordered_map<const DReactionStepVertexInfo*, const DSourceCombo*> locVertexPrimaryComboMap;
	for(const auto& locStepVertexInfo : locReactionVertexInfo->Get_StepVertexInfos())
	{
		/***************************************** CHECK IF VERTEX POSITION IS INDETERMINATE A THIS STAGE ********************************************/

		auto locSourceCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionChargedCombo, locStepVertexInfo);

		if(locStepVertexInfo->Get_DanglingVertexFlag())
		{
			//is forever indeterminate, even with neutrals and beam energy
			//If this is the production vertex, choose the center of the target.  If not, choose the vertex where the decay parent was produced.
			auto locParentVertexInfo = locStepVertexInfo->Get_ParentVertexInfo();
			if(locParentVertexInfo == nullptr) //production vertex
				dVertexParticlesByCombo.emplace(std::make_pair(true, locSourceCombo), {}); //empty: target center
			else //decay products
			{
				auto locParentCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionChargedCombo, locParentVertexInfo);
				auto locIsParentProductionVertex = locParentVertexInfo->Get_ProductionVertexFlag();
				dVertexParticlesByCombo.emplace(std::make_pair(false, locSourceCombo), dVertexParticlesByCombo[std::make_pair(locIsParentProductionVertex, locParentCombo)]);
			}
			continue;
		}

		//get combo & info
		bool locIsProductionVertexFlag = locStepVertexInfo->Get_ProductionVertexFlag();
		locVertexPrimaryComboMap.emplace(locStepVertexInfo, locSourceCombo);

		//get particles
		auto locSourceParticles = DAnalysis::Get_SourceParticles_ThisVertex(locSourceCombo);
		auto locDecayingParticles = Get_FullConstrainDecayingParticles(locStepVertexInfo, locReconDecayParticleMap);

		//now, just because it isn't dangling, it doesn't mean we have enough information to actually determine the vertex
		//e.g. we may need neutrals or beam energy to define constraining decay particle p4/trajectory
		auto locComboProductionPair = std::make_pair(locIsProductionVertexFlag, locSourceCombo);
		auto locDeterminableIterator = dVertexDeterminableWithChargedMap.find(locComboProductionPair);
		if(locDeterminableIterator == dVertexDeterminableWithChargedMap.end())
		{
			//we don't know yet if it is determinable or not.  figure it out
			auto locNumConstrainingParticles = locSourceParticles.size() + locDecayingParticles.size();

			//determine it, save result, and continue if can't
			bool locVertexDeterminableFlag = locIsProductionVertexFlag ? (locNumConstrainingParticles > 0) : (locNumConstrainingParticles > 1);
			dVertexDeterminableWithChargedMap.emplace(locComboProductionPair, locVertexDeterminableFlag);
			if(!locVertexDeterminableFlag)
				continue; //can't determine yet, but will be able to in the future
		}
		if(locDeterminableIterator->second == false)
			continue;

		//Get detected charged track hypotheses at this vertex
		vector<const DKinematicData*> locVertexParticles;
		locVertexParticles.reserve(locSourceParticles.size());

		auto Get_Hypothesis = [](const pair<Particle_t, const JObject*>& locPair) -> bool
			{return static_cast<const DChargedTrack*>(locPair.second)->Get_Hypothesis(locPair.first);};
		std::transform(locVertexParticles.begin(), locVertexParticles.end(), std::back_inserter(locVertexParticles), Get_Hypothesis);
		locVertexParticles.insert(locVertexParticles.end(), locDecayingParticles.begin(), locDecayingParticles.end());

		/*********************************************************** INITIAL VERTEX GUESS ************************************************************/

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

		DVector3 locVertex;
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
			locVertex = locVertexParticles[0]->position();
		}
		else //detached vertex
		{
			if(locVertexParticles.size() < 2)
				locVertexParticles.insert(locVertexParticles.end(), locDecayingParticles.begin(), locDecayingParticles.end());
			std::sort(locVertexParticles.begin(), locVertexParticles.end());

			//if vertex already computed, don't bother
			auto locVertexIterator = dVertexMap.find(std::make_pair(locIsProductionVertexFlag, locVertexParticles));
			if(locVertexIterator == dVertexMap.end())
				locVertex = dAnalysisUtilities->Calc_CrudeVertex(locVertexParticles);
			else
				locVertex = locVertexIterator->second;
		}

		/************************************************************ UPDATE WITH RESULTS ************************************************************/

		dVertexParticlesByCombo.emplace(locComboProductionPair, locVertexParticles);
		dVertexMap.emplace(std::make_pair(locIsProductionVertexFlag, locVertexParticles), locVertex);
		Construct_DecayingParticle(locStepVertexInfo, locSourceCombo, locVertex, locReconDecayParticleMap);
	}

	//do time offsets once all the vertices have been found
	Calc_TimeOffsets(locReactionVertexInfo, locReactionChargedCombo);

	//set the vertices (really, the vertex particles) for the other combos for faster lookup
	for(const auto& locVertexPrimaryComboPair : locVertexPrimaryComboMap)
	{
		bool locIsProductionVertexFlag = locVertexPrimaryComboPair.first->Get_ProductionVertexFlag();
		auto locVertexCombos = DAnalysis::Get_SourceCombos_ThisVertex(locVertexPrimaryComboPair.second);
		for(auto& locVertexCombo : locVertexCombos)
			dVertexParticlesByCombo.emplace(std::make_pair(locIsProductionVertexFlag, locVertexCombo), dVertexParticlesByCombo[std::make_pair(locIsProductionVertexFlag, locVertexPrimaryComboPair.second)]);
	}
}

vector<const DKinematicData*> DSourceComboVertexer::Get_FullConstrainDecayingParticles(const DReactionStepVertexInfo* locStepVertexInfo, const map<pair<int, int>, const DKinematicData*>& locReconDecayParticleMap)
{
	auto locConstrainingDecayingParticles = locStepVertexInfo->Get_DecayingParticles_FullConstrain();
	vector<const DKinematicData*> locDecayingParticles;
	for(auto locDecayPair : locConstrainingDecayingParticles)
	{
		//the particle pair in the map is where the decaying particle was in the initial state
		//however, the pairs in locConstrainingDecayingParticles are when then the decaying particles were in the final state
		//so, we have to convert from final state to initial state to do the map lookup
		int locDecayStepIndex = DAnalysis::Get_DecayStepIndex(locStepVertexInfo->Get_Reaction(), locDecayPair.first.first, locDecayPair.first.second);
		auto locParticlePair = std::make_pair(locDecayStepIndex, DReactionStep::Get_ParticleIndex_Initial());

		auto locReconIterator = locReconDecayParticleMap.find(locParticlePair);
		if(locReconIterator != locReconDecayParticleMap.end())
			locDecayingParticles.push_back(locReconIterator->second);
	}
	return locDecayingParticles;
}

void DSourceComboVertexer::Construct_DecayingParticle(const DReactionStepVertexInfo* locReactionStepVertexInfo, const DSourceCombo* locSourceCombo, DVector3 locVertex, map<pair<int, int>, const DKinematicData*>& locReconDecayParticleMap)
{
	//we also can't compute the p4 yet if the decay products contain neutrals: this is done before comboing them!
	auto locSourceComboUse = dSourceComboer->Get_SourceComboUse(locReactionStepVertexInfo);
	if(dSourceComboer->Get_ChargeContent(std::get<2>(locSourceComboUse)) != d_Charged)
		return;

	//get "decaying" no-constrain decaying particles
	map<pair<int, int>, weak_ptr<DReactionStepVertexInfo>> locNoConstrainDecayingParticles = locReactionStepVertexInfo->Get_DecayingParticles_NoConstrain();
	auto locReaction = locReactionStepVertexInfo->Get_Reaction();
	for(const auto& locNoConstrainPair : locNoConstrainDecayingParticles)
	{
		//we cannot define decaying p4 via missing mass, because the beam is not chosen yet
		//therefore, if the no-constrain decaying particle is a final-state particle at this vertex, we cannot define its p4
		const auto& locParticlePair = locNoConstrainPair.first;
		if(locParticlePair.second >= 0)
			continue; //momentum must be determined by missing mass, cannot do yet! (or isn't detached, for which we don't care)

		//only create kinematic data for detached PIDs
		auto locDecayPID = locReaction->Get_ReactionStep(locParticlePair.first)->Get_PID(locParticlePair.second);
		if(!IsDetachedVertex(locDecayPID))
			continue; //not detached: we don't care!

		//if the particle has already been created, reuse it!
		auto locDecayParticleTuple = std::make_tuple(locDecayPID, dVertexParticlesByCombo[locSourceCombo]);
		auto locIterator = dReconDecayParticles.find(locDecayParticleTuple);
		if(locIterator != dReconDecayParticles.end())
		{
			locReconDecayParticleMap.emplace(locParticlePair, locIterator->second);
			continue;
		}

		//create a new one
		auto locP4 = dSourceComboP4Handler->Get_P4(locSourceCombo);
//switch to resource pooL?
		auto locKinematicData = new DKinematicData(locDecayPID, locP4.Vect(), locVertex, 0.0);

		//register it
		locReconDecayParticleMap.emplace(locParticlePair, locKinematicData);
		dReconDecayParticles.emplace(locDecayParticleTuple, locKinematicData);
	}
}

void DSourceComboVertexer::Calc_TimeOffsets(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionCombo)
{
	auto locIsPrimaryProductionVertex = locReactionVertexInfo->Get_StepVertexInfos().front()->Get_ProductionVertexFlag();
	auto locReactionPair = std::make_pair(locIsPrimaryProductionVertex, locReactionCombo);

	//loop over vertices
	for(auto locStepVertexInfo : locReactionVertexInfo->Get_StepVertexInfos())
	{
		if(locStepVertexInfo->Get_ProductionVertexFlag())
		{
			dTimeOffsets[locReactionPair].emplace(locReactionCombo, 0.0);
			continue;
		}
		auto locParentVertexInfo = locStepVertexInfo->Get_ParentVertexInfo();

		//get combos
		auto locVertexCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionCombo, locStepVertexInfo);
		auto locParentCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionCombo, locParentVertexInfo);

		//get vertices & path length
		auto locVertex = Get_Vertex(false, locVertexCombo);
		auto locParentProductionVertex = Get_Vertex(locParentVertexInfo->Get_ProductionVertexFlag(), locParentCombo);
		auto locPathLength = (locVertex - locParentProductionVertex).Mag();

		//compute and save result
		auto locP4 = dSourceComboP4Handler->Get_P4(locVertexCombo);
		auto locTimeOffset = locPathLength/(locP4.Beta()*SPEED_OF_LIGHT);
		dTimeOffsets[locReactionPair].emplace(locVertexCombo, locTimeOffset);
	}
}

} //end DAnalysis namespace

