#include "ANALYSIS/DSourceComboVertexer.h"
#include "ANALYSIS/DSourceComboer.h"
#include "ANALYSIS/DSourceComboP4Handler.h"
#include "ANALYSIS/DSourceComboTimeHandler.h"

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

	gPARMS->SetDefaultParameter("COMBO:DEBUG_LEVEL", dDebugLevel);
}

vector<signed char> DSourceComboVertexer::Get_VertexZBins(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionCombo, const DKinematicData* locBeamParticle) const
{
	if(locReactionCombo == nullptr)
		return {};

	vector<signed char> locVertexZBins;
	for(auto locStepInfo : locReactionVertexInfo->Get_StepVertexInfos())
	{
		auto locVertexPrimaryCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionCombo, locStepInfo);
		locVertexZBins.emplace_back(Get_VertexZBin(locStepInfo->Get_ProductionVertexFlag(), locVertexPrimaryCombo, locBeamParticle));
	}
	return locVertexZBins;
}

signed char DSourceComboVertexer::Get_VertexZBin(const DReactionStepVertexInfo* locStepVertexInfo, const DSourceCombo* locReactionCombo, const DKinematicData* locBeamParticle) const
{
	auto locVertexPrimaryCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionCombo, locStepVertexInfo);
	return Get_VertexZBin(locStepVertexInfo->Get_ProductionVertexFlag(), locVertexPrimaryCombo, locBeamParticle);
}

signed char DSourceComboVertexer::Get_VertexZBin(bool locIsProductionVertex, const DSourceCombo* locSourceCombo, const DKinematicData* locBeamParticle) const
{
	if(locSourceCombo == nullptr)
		return DSourceComboInfo::Get_VertexZIndex_Unknown();

	auto locConstrainingParticles = Get_ConstrainingParticles(locIsProductionVertex, locSourceCombo, locBeamParticle);
	if(locConstrainingParticles.empty())
		return DSourceComboInfo::Get_VertexZIndex_Unknown();
	auto locVertexZ = Get_Vertex(locIsProductionVertex, locConstrainingParticles).Z();
	return dSourceComboTimeHandler->Get_PhotonVertexZBin(locVertexZ);
}

void DSourceComboVertexer::Calc_VertexTimeOffsets_WithCharged(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionChargedCombo)
{
	if(dDebugLevel >= 10)
		cout << "DSourceComboVertexer::Calc_VertexTimeOffsets_WithCharged()" << endl;

	auto locIsPrimaryProductionVertex = locReactionVertexInfo->Get_StepVertexInfo(0)->Get_ProductionVertexFlag();
	//even if below is true, we still need to register step vertex infos
	auto locEverythingFoundFlag = (dTimeOffsets.find(std::make_tuple(locIsPrimaryProductionVertex, locReactionChargedCombo, (const DKinematicData*)nullptr)) != dTimeOffsets.end());

	//loop through vertices
	map<pair<int, int>, const DKinematicData*> locReconDecayParticleMap; //decaying particle indices -> kinematic data //indices: when the decaying particle is in the INITIAL state
	unordered_map<const DReactionStepVertexInfo*, const DSourceCombo*> locVertexPrimaryComboMap;
	for(const auto& locStepVertexInfo : locReactionVertexInfo->Get_StepVertexInfos())
	{
		if(dDebugLevel >= 10)
			cout << "Step: " << locStepVertexInfo->Get_StepIndices().front() << endl;

		/***************************************** CHECK IF VERTEX POSITION IS INDETERMINATE A THIS STAGE ********************************************/

		auto locVertexPrimaryCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionChargedCombo, locStepVertexInfo);
		auto locIsProductionVertexFlag = locStepVertexInfo->Get_ProductionVertexFlag();
		auto locComboProductionTuple = std::make_tuple(locIsProductionVertexFlag, locVertexPrimaryCombo, (const DKinematicData*)nullptr);

		if(locStepVertexInfo->Get_DanglingVertexFlag())
		{
			//is forever indeterminate, even with neutrals and beam energy
			//If this is the production vertex, choose the center of the target.  If not, choose the vertex where the decay parent was produced.
			auto locParentVertexInfo = locStepVertexInfo->Get_ParentVertexInfo();
			if(dDebugLevel >= 10)
				cout << "Dangling, parent info = " << locParentVertexInfo << endl;
			if(locParentVertexInfo == nullptr) //initial vertex
				dConstrainingParticlesByCombo.emplace(locComboProductionTuple, vector<const DKinematicData*>{}); //empty: target center
			else //decay products
			{
				auto locParentCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionChargedCombo, locParentVertexInfo);
				auto locIsParentProductionVertex = locParentVertexInfo->Get_ProductionVertexFlag();
				auto& locVertexParticles = dConstrainingParticlesByCombo[std::make_tuple(locIsParentProductionVertex, locParentCombo, (const DKinematicData*)nullptr)];
				dConstrainingParticlesByCombo.emplace(locComboProductionTuple, locVertexParticles);
			}
			continue;
		}

		//get combo & info
		locVertexPrimaryComboMap.emplace(locStepVertexInfo, locVertexPrimaryCombo);

		//get particles
		auto locChargedSourceParticles = DAnalysis::Get_SourceParticles_ThisVertex(locVertexPrimaryCombo);
		auto locDecayingParticles = Get_FullConstrainDecayingParticles(locStepVertexInfo, locReconDecayParticleMap);
		if(dDebugLevel >= 10)
			cout << "# charged/decaying particles at vertex = " << locChargedSourceParticles.size() << ", " << locDecayingParticles.size() << endl;

		//now, just because it isn't dangling, it doesn't mean we have enough information to actually determine the vertex
		//e.g. we may need neutrals or beam energy to define constraining decay particle p4/trajectory
		auto locDeterminableIterator = dVertexDeterminableWithChargedMap.find(locStepVertexInfo);
		if(locDeterminableIterator == dVertexDeterminableWithChargedMap.end())
		{
			//we don't know yet if it is determinable or not.  figure it out
			auto locNumConstrainingParticles = locChargedSourceParticles.size() + locDecayingParticles.size();

			//determine it, save result, and continue if can't
			bool locVertexDeterminableFlag = locIsProductionVertexFlag ? (locNumConstrainingParticles > 0) : (locNumConstrainingParticles > 1);
			if(dDebugLevel >= 10)
				cout << "vertex determinable flag = " << locVertexDeterminableFlag << endl;
			dVertexDeterminableWithChargedMap.emplace(locStepVertexInfo, locVertexDeterminableFlag);
			if(!locVertexDeterminableFlag)
				continue; //can't determine yet, but will be able to in the future
		}
		else if(!locDeterminableIterator->second)
			continue;
		if(locEverythingFoundFlag)
			continue; //already done

		vector<const DKinematicData*> locVertexParticles;
		auto locVertex = Calc_Vertex(locIsProductionVertexFlag, locChargedSourceParticles, locDecayingParticles, locVertexParticles);
		dConstrainingParticlesByCombo.emplace(locComboProductionTuple, locVertexParticles);

		//set the vertices (really, the vertex particles) for the other combos for faster lookup during neutral comboing
		auto locVertexCombos = DAnalysis::Get_SourceCombos_ThisVertex(locVertexPrimaryCombo);
		for(const auto& locVertexCombo : locVertexCombos)
		{
			auto locVertexParticles = dConstrainingParticlesByCombo[std::make_tuple(locIsProductionVertexFlag, locVertexPrimaryCombo, (const DKinematicData*)nullptr)];
			dConstrainingParticlesByCombo.emplace(std::make_tuple(locIsProductionVertexFlag, locVertexCombo, (const DKinematicData*)nullptr), locVertexParticles);
		}
		Construct_DecayingParticle_InvariantMass(locStepVertexInfo, locVertexPrimaryCombo, locVertex, locReconDecayParticleMap);
	}
	if(locEverythingFoundFlag)
		return; //already done

	//do time offsets once all the vertices have been found
	Calc_TimeOffsets(locReactionVertexInfo, locReactionChargedCombo, nullptr);
}

void DSourceComboVertexer::Calc_VertexTimeOffsets_WithPhotons(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionChargedCombo, const DSourceCombo* locReactionFullCombo)
{
	//Calculate vertex positions & time offsets using photons
	//not likely to have any effect, but it's necessary sometimes (but rarely)
	//E.g. g, p ->  K0, Sigma+    K0 -> 3pi: The selected pi0 photons could help define the production vertex

	auto locIsPrimaryProductionVertex = locReactionVertexInfo->Get_StepVertexInfo(0)->Get_ProductionVertexFlag();
	//even if below is true, we still need to register step vertex infos
	auto locEverythingFoundFlag = (dTimeOffsets.find(std::make_tuple(locIsPrimaryProductionVertex, locReactionFullCombo, (const DKinematicData*)nullptr)) != dTimeOffsets.end());

	//loop over vertices in dependency order
	map<pair<int, int>, const DKinematicData*> locReconDecayParticleMap; //decaying particle indices -> kinematic data //indices: when the decaying particle is in the INITIAL state
	for(const auto& locStepVertexInfo : locReactionVertexInfo->Get_StepVertexInfos())
	{
		if(locStepVertexInfo->Get_DanglingVertexFlag())
			continue; //is forever indeterminate, even with neutrals and beam energy

		auto locIsProductionVertexFlag = locStepVertexInfo->Get_ProductionVertexFlag();
		auto locVertexPrimaryChargedCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionChargedCombo, locStepVertexInfo);
		auto locVertexPrimaryFullCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionFullCombo, locStepVertexInfo);

		//see if vertex has already been found
		auto locChargedComboProductionTuple = std::make_tuple(locIsProductionVertexFlag, locVertexPrimaryChargedCombo, (const DKinematicData*)nullptr);
		auto locFullComboProductionTuple = std::make_tuple(locIsProductionVertexFlag, locVertexPrimaryFullCombo, (const DKinematicData*)nullptr);
		if(dConstrainingParticlesByCombo.find(locChargedComboProductionTuple) != dConstrainingParticlesByCombo.end())
		{
			//already done! get/construct any recon decaying particles
			auto locVertex = Get_Vertex(locIsProductionVertexFlag, locVertexPrimaryChargedCombo, nullptr);
			auto locVertexParticles = dConstrainingParticlesByCombo[locChargedComboProductionTuple];
			dConstrainingParticlesByCombo.emplace(locFullComboProductionTuple, locVertexParticles); //so that the vertex can be retrieved by either the charged or full combo
			Construct_DecayingParticle_InvariantMass(locStepVertexInfo, locVertexPrimaryFullCombo, locVertex, locReconDecayParticleMap);
			continue;
		}

		//get particles
		auto locChargedSourceParticles = DAnalysis::Get_SourceParticles_ThisVertex(locVertexPrimaryChargedCombo);
		auto locDecayingParticles = Get_FullConstrainDecayingParticles(locStepVertexInfo, locReconDecayParticleMap);

		//now, just because it isn't dangling, it doesn't mean we have enough information to actually determine the vertex
		//e.g. we may need beam energy to define constraining decay particle p4/trajectory
		auto locDeterminableIterator = dVertexDeterminableWithPhotonsMap.find(locStepVertexInfo);
		if(locDeterminableIterator == dVertexDeterminableWithPhotonsMap.end())
		{
			//we don't know yet if it is determinable or not.  figure it out
			auto locNumConstrainingParticles = locChargedSourceParticles.size() + locDecayingParticles.size();

			//determine it, save result, and continue if can't
			bool locVertexDeterminableFlag = locIsProductionVertexFlag ? (locNumConstrainingParticles > 0) : (locNumConstrainingParticles > 1);
			dVertexDeterminableWithPhotonsMap.emplace(locStepVertexInfo, locVertexDeterminableFlag);
			if(!locVertexDeterminableFlag)
				continue; //can't determine yet, but will be able to in the future
		}
		else if(!locDeterminableIterator->second)
			continue;
		if(locEverythingFoundFlag)
			continue; //already done

		//find the vertex, save the results
		vector<const DKinematicData*> locVertexParticles;
		auto locVertex = Calc_Vertex(locIsProductionVertexFlag, locChargedSourceParticles, locDecayingParticles, locVertexParticles);
		dConstrainingParticlesByCombo.emplace(locFullComboProductionTuple, locVertexParticles);

		Construct_DecayingParticle_InvariantMass(locStepVertexInfo, locVertexPrimaryFullCombo, locVertex, locReconDecayParticleMap);
	}
	if(locEverythingFoundFlag)
		return; //already done

	//CALC TIME OFFSETS
	Calc_TimeOffsets(locReactionVertexInfo, locReactionChargedCombo, locReactionFullCombo);
}

void DSourceComboVertexer::Calc_VertexTimeOffsets_WithBeam(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, const DKinematicData* locBeamParticle)
{
	if(dDebugLevel >= 10)
		cout << "DSourceComboVertexer::Calc_VertexTimeOffsets_WithBeam()" << endl;

	//IF PRIMARY VERTEX IS UNKNOWN THEN DO NOTHING!
	if(!Get_IsVertexFoundFlag(true, locReactionFullCombo, locBeamParticle))
	{
		if(dDebugLevel >= 10)
			cout << "primary vertex unknown, nothing we can do." << endl;
		return;
	}

	auto locIsPrimaryProductionVertex = locReactionVertexInfo->Get_StepVertexInfo(0)->Get_ProductionVertexFlag();
	auto locPrimaryVertexZ = Get_PrimaryVertex(locReactionVertexInfo, locReactionFullCombo, locBeamParticle).Z();
	int locRFBunch = dSourceComboTimeHandler->Calc_RFBunchShift(locBeamParticle->time());
	if(dDebugLevel >= 10)
		cout << "primary vertex-z rf bunch: " << locPrimaryVertexZ << ", " << locRFBunch << endl;

	//uses the beam to define remaining vertices
	//loop over vertices in dependency order
	map<pair<int, int>, const DKinematicData*> locReconDecayParticleMap; //decaying particle indices -> kinematic data //indices: when the decaying particle is in the INITIAL state
	for(const auto& locStepVertexInfo : locReactionVertexInfo->Get_StepVertexInfos())
	{
		if(dDebugLevel >= 10)
			cout << "Step: " << locStepVertexInfo->Get_StepIndices().front() << ", dangling-flag = " << locStepVertexInfo->Get_DanglingVertexFlag() << endl;
		if(locStepVertexInfo->Get_DanglingVertexFlag())
			continue; //is forever indeterminate, even with beam energy

		auto locIsProductionVertexFlag = locStepVertexInfo->Get_ProductionVertexFlag();
		auto locVertexPrimaryFullCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionFullCombo, locStepVertexInfo);

		//see if vertex has already been found //can search with either charged or full
		auto locFullComboProductionTuple = std::make_tuple(locIsProductionVertexFlag, locVertexPrimaryFullCombo, (const DKinematicData*)nullptr);
		if(dConstrainingParticlesByCombo.find(locFullComboProductionTuple) != dConstrainingParticlesByCombo.end())
		{
			if(dDebugLevel >= 10)
				cout << "vertex already found" << endl;

			//already done! get/construct any recon decaying particles
			auto locVertex = Get_Vertex(locIsProductionVertexFlag, locVertexPrimaryFullCombo, locBeamParticle);
			Construct_DecayingParticle_InvariantMass(locStepVertexInfo, locVertexPrimaryFullCombo, locVertex, locReconDecayParticleMap);
			if(locIsProductionVertexFlag) //do missing if needed
			{
				auto locRFVertexTime = dSourceComboTimeHandler->Calc_PropagatedRFTime(locPrimaryVertexZ, locRFBunch, 0.0);
				Construct_DecayingParticle_MissingMass(locStepVertexInfo, locReactionFullCombo, locVertexPrimaryFullCombo, locBeamParticle, locVertex, locRFBunch, locRFVertexTime, locReconDecayParticleMap);
			}
			continue;
		}

		if(dDebugLevel >= 10)
			cout << "vertex not found yet" << endl;

		//get particles
		auto locChargedSourceParticles = DAnalysis::Get_SourceParticles_ThisVertex(locVertexPrimaryFullCombo, d_Charged);
		auto locDecayingParticles = Get_FullConstrainDecayingParticles(locStepVertexInfo, locReconDecayParticleMap);

		//find the vertex, save the results
		vector<const DKinematicData*> locVertexParticles;
		auto locVertex = Calc_Vertex(locIsProductionVertexFlag, locChargedSourceParticles, locDecayingParticles, locVertexParticles);
		dConstrainingParticlesByCombo.emplace(locFullComboProductionTuple, locVertexParticles);

		//CALC AND STORE TIME OFFSET
		auto locFullReactionTuple_TimeOffset = std::make_tuple(locIsPrimaryProductionVertex, locReactionFullCombo, locBeamParticle);

		//get parent
		auto locParentVertexInfo = locStepVertexInfo->Get_ParentVertexInfo();
		auto locParentVertexFullCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionFullCombo, locParentVertexInfo);
		auto locParentTimeOffset = dTimeOffsets[locFullReactionTuple_TimeOffset][locParentVertexFullCombo];

		//get vertices & path length
		auto locParentProductionVertex = Get_Vertex(locParentVertexInfo->Get_ProductionVertexFlag(), locParentVertexFullCombo, locBeamParticle);
		auto locPathLength = (locVertex - locParentProductionVertex).Mag();

		//compute and save result
		auto locP4 = locDecayingParticles.back()->lorentzMomentum(); //when locDecayingParticles retrieved, this particle was saved to the back!
		auto locTimeOffset = locPathLength/(locP4.Beta()*SPEED_OF_LIGHT) + locParentTimeOffset;
		dTimeOffsets[locFullReactionTuple_TimeOffset].emplace(locVertexPrimaryFullCombo, locTimeOffset);

		//construct decaying particles by missing mass
		auto locRFVertexTime = dSourceComboTimeHandler->Calc_PropagatedRFTime(locPrimaryVertexZ, locRFBunch, locTimeOffset);
		if(dDebugLevel >= 10)
			cout << "parent time offset, beta, path length, time offset, prop rf time: " << locParentTimeOffset << ", " << locP4.Beta() << ", " << locPathLength << ", " << locTimeOffset << ", " << locRFVertexTime << endl;
		Construct_DecayingParticle_MissingMass(locStepVertexInfo, locReactionFullCombo, locVertexPrimaryFullCombo, locBeamParticle, locVertex, locRFBunch, locRFVertexTime, locReconDecayParticleMap);
	}
}

DVector3 DSourceComboVertexer::Calc_Vertex(bool locIsProductionVertexFlag, const vector<pair<Particle_t, const JObject*>>& locChargedSourceParticles, const vector<const DKinematicData*>& locDecayingParticles, vector<const DKinematicData*>& locVertexParticles)
{
	if(dDebugLevel >= 10)
		cout << "DSourceComboVertexer::Calc_Vertex()" << endl;

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

	if(dDebugLevel >= 10)
		cout << "locIsProductionVertexFlag = " << locIsProductionVertexFlag << endl;
	if(locIsProductionVertexFlag)
	{
		//use track with theta nearest 90 degrees
		auto locThetaNearest90Iterator = Get_ThetaNearest90Iterator(locVertexParticles);
		double locThetaNearest90 = (locThetaNearest90Iterator != locVertexParticles.end()) ? (*locThetaNearest90Iterator)->momentum().Theta()*180.0/TMath::Pi() : 0.0;
		if(locThetaNearest90 < dMinThetaForVertex)
		{
			//try decaying particles instead
			auto locThetaNearest90Iterator_Decaying = Get_ThetaNearest90Iterator(locDecayingParticles);
			double locLargestTheta_Decaying = (locThetaNearest90Iterator_Decaying != locDecayingParticles.end()) ? (*locThetaNearest90Iterator_Decaying)->momentum().Theta()*180.0/TMath::Pi() : 0.0;
			if(locLargestTheta_Decaying > locThetaNearest90)
				locThetaNearest90Iterator = locThetaNearest90Iterator_Decaying;
		}
		locVertexParticles = {*locThetaNearest90Iterator};
		//vertex is 1/2-way between track POCA to beamline and the beamline itself: if POCA not on beamline, likely due to resolution issues, 
		auto locTrackPosition = locVertexParticles[0]->position();
		auto locVertex = DVector3(0.5*locTrackPosition.X(), 0.5*locTrackPosition.Y(), locTrackPosition.Z());
//		if(false) //COMPARE: Comparison-to-old mode
			locVertex = dVertex->dSpacetimeVertex.Vect();
		dVertexMap.emplace(std::make_pair(locIsProductionVertexFlag, locVertexParticles), locVertex);
		if(dDebugLevel >= 10)
			cout << "particle PID, theta, production vertex = " << locVertexParticles[0]->PID() << ", " << locVertexParticles[0]->momentum().Theta()*180.0/TMath::Pi() << ", " << locVertex.X() << ", " << locVertex.Y() << ", " << locVertex.Z() << endl;
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
		if(dDebugLevel >= 10)
			cout << "crude vertex = " << locVertex.X() << ", " << locVertex.Y() << ", " << locVertex.Z() << endl;
		dVertexMap.emplace(std::make_pair(locIsProductionVertexFlag, locVertexParticles), locVertex);
		return locVertex;
	}

	return locVertexIterator->second;
}

vector<const DKinematicData*> DSourceComboVertexer::Get_FullConstrainDecayingParticles(const DReactionStepVertexInfo* locStepVertexInfo, const map<pair<int, int>, const DKinematicData*>& locReconDecayParticleMap)
{
	auto locConstrainingDecayingParticles = locStepVertexInfo->Get_DecayingParticles_FullConstrain(false);
	vector<const DKinematicData*> locDecayingParticles;
	pair<int, int> locDecayMissingMassReconPair(-99, -99); //if there is one (and there will be at most one), save for last (easy to retrieve)
	for(const auto& locDecayPair : locConstrainingDecayingParticles)
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

void DSourceComboVertexer::Construct_DecayingParticle_InvariantMass(const DReactionStepVertexInfo* locStepVertexInfo, const DSourceCombo* locVertexCombo, DVector3 locVertex, map<pair<int, int>, const DKinematicData*>& locReconDecayParticleMap)
{
	if(dDebugLevel >= 10)
		cout << "DSourceComboVertexer::Construct_DecayingParticle_InvariantMass()" << endl;
	//we also can't compute the p4 yet if the decay products contain neutrals: this is done before comboing them!
	auto locSourceComboUse = dSourceComboer->Get_SourceComboUse(locStepVertexInfo);
	if(dSourceComboer->Get_HasMassiveNeutrals(std::get<2>(locSourceComboUse)))
	{
		//can't compute the p4 yet if the decay products contain massive neutrals: time offset is needed for massive neutral p4, but isn't known yet because decay p4 unknown!
		//will have to be computed via missing mass instead (although it technically isn't needed for it)
		return;
	}

	auto locIsProductionVertexFlag = locStepVertexInfo->Get_ProductionVertexFlag();
	if(dDebugLevel >= 10)
		cout << "locIsProductionVertexFlag = " << locIsProductionVertexFlag << endl;

	//loop over decaying no-constrain decaying particles
	map<pair<int, int>, const DReactionStepVertexInfo*> locNoConstrainDecayingParticles = locStepVertexInfo->Get_DecayingParticles_NoConstrain(false);
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

		if(dDebugLevel >= 10)
			cout << "detached decaying pid = " << locDecayPID << endl;

		//if the particle has already been created, reuse it!
		auto locDecayParticleTuple = std::make_tuple(locDecayPID, locIsProductionVertexFlag, locVertexCombo, (const DKinematicData*)nullptr);
		auto locIterator = dReconDecayParticles_FromProducts.find(locDecayParticleTuple);
		if(locIterator != dReconDecayParticles_FromProducts.end())
		{
			locReconDecayParticleMap.emplace(locParticlePair, locIterator->second);
			continue;
		}

		//create a new one
		auto locP4 = dSourceComboP4Handler->Calc_P4_NoMassiveNeutrals(locVertexCombo, locVertex, std::get<1>(locSourceComboUse), nullptr, false);
		auto locKinematicData = dResourcePool_KinematicData.Get_Resource();
		locKinematicData->Reset();
		locKinematicData->Set_Members(locDecayPID, locP4.Vect(), locVertex, 0.0);

		if(dDebugLevel >= 10)
			cout << "decaying particle created, pxyzE = " << locP4.X() << ", " << locP4.Y() << ", " << locP4.Z() << ", " << locP4.E() << endl;

		//register it
		locReconDecayParticleMap.emplace(locParticlePair, locKinematicData);
		dReconDecayParticles_FromProducts.emplace(locDecayParticleTuple, locKinematicData);
	}
}

void DSourceComboVertexer::Construct_DecayingParticle_MissingMass(const DReactionStepVertexInfo* locStepVertexInfo, const DSourceCombo* locReactionFullCombo, const DSourceCombo* locFullVertexCombo, const DKinematicData* locBeamParticle, DVector3 locVertex, int locRFBunch, double locRFVertexTime, map<pair<int, int>, const DKinematicData*>& locReconDecayParticleMap)
{
	if(dDebugLevel >= 10)
		cout << "DSourceComboVertexer::Construct_DecayingParticle_MissingMass()" << endl;

	if(locStepVertexInfo->Get_IsInclusiveVertexFlag() || !locStepVertexInfo->Get_MissingParticles().empty())
		return; //decaying particles are not defined!

	//we can only calculate up to 1 at a time
	//the input full combo contains the decaying particle for which the missing mass is to be computed
	//if there is more than one, then this is impossible

	auto locSourceComboUse = dSourceComboer->Get_SourceComboUse(locStepVertexInfo);
	auto locReaction = locStepVertexInfo->Get_Reaction();
	auto locIsProductionVertexFlag = locStepVertexInfo->Get_ProductionVertexFlag();

	//loop over decaying no-constrain decaying particles, figure out how many we have to do this for (can't do more than 1!)
	map<pair<int, int>, const DReactionStepVertexInfo*> locNoConstrainDecayingParticles = locStepVertexInfo->Get_DecayingParticles_NoConstrain(false);
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

	if(dDebugLevel >= 10)
		cout << "to-recon decay particle indices: " << locToReconParticleIndices.first << ", " << locToReconParticleIndices.second << endl;
	if(locToReconParticleIndices.first < 0)
		return;

	auto locDecayStepIndex = Get_DecayStepIndex(locReaction, locToReconParticleIndices.first, locToReconParticleIndices.second);
	auto locDecayUse = dSourceComboer->Get_SourceComboUse(locReaction, locDecayStepIndex);

	auto locReactionStep = locReaction->Get_ReactionStep(locToReconParticleIndices.first);
	auto locDecayPID = locReactionStep->Get_PID(locToReconParticleIndices.second);

	//if the particle has already been created, reuse it!
	auto locDecayParticleTuple = std::make_tuple(locDecayPID, locReactionFullCombo, locIsProductionVertexFlag, locFullVertexCombo, locBeamParticle);
	auto locIterator = dReconDecayParticles_FromMissing.find(locDecayParticleTuple);
	if(locIterator != dReconDecayParticles_FromMissing.end())
	{
		locReconDecayParticleMap.emplace(locToReconParticleIndices, locIterator->second);
		return;
	}

	//make sure there isn't another missing particle
	for(const auto& locStepIndex : locStepVertexInfo->Get_StepIndices())
	{
		if(int(locStepIndex) == locDecayStepIndex)
			continue;
		if(DAnalysis::Get_HasMissingParticle_FinalState(locReaction->Get_ReactionStep(locStepIndex)))
			return; //missing decay product! can't compute 2 missing p4's at once!
	}

	//create a new one
	//calc final state p4
	DLorentzVector locFinalStateP4;
	if(!dSourceComboP4Handler->Calc_P4_HasMassiveNeutrals(locIsProductionVertexFlag, true, locReactionFullCombo, locFullVertexCombo, locVertex, locRFBunch, locRFVertexTime, locDecayUse, locFinalStateP4, locBeamParticle, true))
		return; //invalid somehow

	//ASSUMES FIXED TARGET EXPERIMENT!
	DLorentzVector locInitialStateP4; //lookup or is beam + target
	if(locIsProductionVertexFlag)
	{
		Particle_t locPID = locReaction->Get_ReactionStep(0)->Get_TargetPID();
		locInitialStateP4 = locBeamParticle->lorentzMomentum() + DLorentzVector(0.0, 0.0, 0.0, ParticleMass(locPID));
	}
	else //already computed, look it up!
	{
		auto locPreviousVertexInfo = (locStepVertexInfo->Get_ParentVertexInfo() == nullptr) ? locStepVertexInfo : locStepVertexInfo->Get_ParentVertexInfo();
		auto locPreviousIsProdVertexFlag = locPreviousVertexInfo->Get_ProductionVertexFlag();
		auto locPreviousFullVertexCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionFullCombo, locPreviousVertexInfo);
		auto locPreviousDecayPID = std::get<0>(locSourceComboUse);
		auto locPreviousDecayParticleTuple = std::make_tuple(locPreviousDecayPID, locReactionFullCombo, locPreviousIsProdVertexFlag, locPreviousFullVertexCombo, locBeamParticle);
		locInitialStateP4 = dReconDecayParticles_FromMissing[locPreviousDecayParticleTuple]->lorentzMomentum();
		if(dDebugLevel >= 10)
			cout << "retrieved decaying pid, p4: " << locPreviousDecayPID << ", " << locInitialStateP4.Px() << ", " << locInitialStateP4.Py() << ", " << locInitialStateP4.Pz() << ", " << locInitialStateP4.E() << endl;
	}

	auto locP4 = locInitialStateP4 - locFinalStateP4;
	if(dDebugLevel >= 10)
		cout << "pid, missing p4: " << locDecayPID << ", " << locP4.Px() << ", " << locP4.Py() << ", " << locP4.Pz() << ", " << locP4.E() << endl;
	auto locKinematicData = dResourcePool_KinematicData.Get_Resource();
	locKinematicData->Reset();
	locKinematicData->Set_Members(locDecayPID, locP4.Vect(), locVertex, 0.0);

	//register it
	locReconDecayParticleMap.emplace(locToReconParticleIndices, locKinematicData);
	dReconDecayParticles_FromMissing.emplace(locDecayParticleTuple, locKinematicData);
}

void DSourceComboVertexer::Calc_TimeOffsets(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locChargedReactionCombo, const DSourceCombo* locFullReactionCombo)
{
	if(dDebugLevel >= 10)
		cout << "DSourceComboVertexer::Calc_TimeOffsets()" << endl;

	//only for calculating from invariant mass, not missing mass vertices
	auto locIsPrimaryProductionVertex = locReactionVertexInfo->Get_StepVertexInfo(0)->Get_ProductionVertexFlag();
	auto locActiveReactionCombo = (locFullReactionCombo != nullptr) ? locFullReactionCombo : locChargedReactionCombo;

	auto locChargedReactionTuple = std::make_tuple(locIsPrimaryProductionVertex, locChargedReactionCombo, (const DKinematicData*)nullptr);
	auto locNeutralReactionTuple = std::make_tuple(locIsPrimaryProductionVertex, locFullReactionCombo, (const DKinematicData*)nullptr);

	auto& locChargedTimeOffsetMap = dTimeOffsets[locChargedReactionTuple];
	auto& locFullTimeOffsetMap = dTimeOffsets[locNeutralReactionTuple];
	auto& locActiveTimeOffsetMap = (locFullReactionCombo != nullptr) ? locFullTimeOffsetMap : locChargedTimeOffsetMap;

	//loop over vertices
	//MUST GO IN STEP ORDER!!
	for(const auto& locStepVertexInfo : DAnalysis::Get_StepVertexInfos_OrderByStep(locReactionVertexInfo))
	{
		if(dDebugLevel >= 10)
			cout << "Step: " << locStepVertexInfo->Get_StepIndices().front() << endl;

		auto locChargedVertexCombo = dSourceComboer->Get_VertexPrimaryCombo(locChargedReactionCombo, locStepVertexInfo);
		auto locFullVertexCombo = (locFullReactionCombo != nullptr) ? dSourceComboer->Get_VertexPrimaryCombo(locFullReactionCombo, locStepVertexInfo) : nullptr;
		auto locActiveVertexCombo = (locFullReactionCombo != nullptr) ? locFullVertexCombo : locChargedVertexCombo;

		//see if already determined
		auto locActiveIterator = locActiveTimeOffsetMap.find(locFullVertexCombo);
		if(locActiveIterator != locActiveTimeOffsetMap.end())
			continue; //previously determined!

		//see if at full-combo stage (not charged only).  If so, and offset already found at charged stage, just copy the result
		if(locFullReactionCombo != nullptr) //full stage
		{
			auto locChargedIterator = locChargedTimeOffsetMap.find(locChargedVertexCombo);
			if(locChargedIterator != locChargedTimeOffsetMap.end())
			{
				//save in full map
				locFullTimeOffsetMap.emplace(locFullVertexCombo, locChargedIterator->second);
				continue; //previously determined!
			}
		}

		auto locParentVertexInfo = locStepVertexInfo->Get_ParentVertexInfo();
		if(locStepVertexInfo->Get_ProductionVertexFlag() || (locParentVertexInfo == nullptr))
		{
			if(dDebugLevel >= 10)
				cout << "Primary vertex: time offset = 0" << endl;
			locChargedTimeOffsetMap.emplace(locChargedReactionCombo, 0.0);
			continue;
		}

		//if this vertex has not been determined yet: save calcing of time offset for later
		if((locFullReactionCombo == nullptr) && !Get_VertexDeterminableWithCharged(locStepVertexInfo))
			continue; //don't know the vertex yet, try the next vertex
		if((locFullReactionCombo != nullptr) && !Get_VertexDeterminableWithPhotons(locStepVertexInfo))
			continue; //don't know the vertex yet, try the next vertex

		//get parent information
		auto locParentCombo = dSourceComboer->Get_VertexPrimaryCombo(locActiveReactionCombo, locParentVertexInfo);
		auto locParentTimeOffset = locActiveTimeOffsetMap[locParentCombo];

		if(dDebugLevel >= 10)
			cout << "Parent time offset = " << locParentTimeOffset << endl;

		//get vertices & path length
		auto locVertex = Get_Vertex(false, locActiveVertexCombo, nullptr);
		auto locParentProductionVertex = Get_Vertex(locParentVertexInfo->Get_ProductionVertexFlag(), locParentCombo, nullptr);
		auto locPathLength = (locVertex - locParentProductionVertex).Mag();
		if(dDebugLevel >= 10)
			cout << "Vertex, parent vertex, path length: " << locVertex.X() << ", " << locVertex.Y() << ", " << locVertex.Z() << ", " << locParentProductionVertex.X() << ", " << locParentProductionVertex.Y() << ", " << locParentProductionVertex.Z() << ", " << locPathLength << endl;

		//compute and save result
		auto locVertexZBin = Get_VertexZBin(false, locActiveVertexCombo, nullptr);
		auto locP4 = dSourceComboP4Handler->Calc_P4_NoMassiveNeutrals(locActiveVertexCombo, locVertex, locVertexZBin, nullptr, false);
		auto locTimeOffset = locPathLength/(locP4.Beta()*SPEED_OF_LIGHT) + locParentTimeOffset;

		if(dDebugLevel >= 10)
			cout << "Time offset = " << locTimeOffset << endl;

		locActiveTimeOffsetMap.emplace(locActiveVertexCombo, locTimeOffset);
	}
}

} //end DAnalysis namespace

