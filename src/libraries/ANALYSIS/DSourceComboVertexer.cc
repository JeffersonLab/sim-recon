#include "ANALYSIS/DSourceComboVertexer.h"
#include "ANALYSIS/DSourceComboer.h"

namespace DAnalysis
{

DSourceComboVertexer::DSourceComboVertexer(JEventLoop* locEventLoop, DSourceComboer* locSourceComboer, DSourceComboP4Handler* locSourceComboP4Handler) :
dSourceComboer(locSourceComboer), dSourceComboP4Handler(locSourceComboP4Handler)
{
	locEventLoop->GetSingle(dAnalysisUtilities);

	//INITIALIZE KINFITUTILS
	dKinFitUtils = new DKinFitUtils_GlueX(locEventLoop);
	dKinFitUtils->Set_IncludeBeamlineInVertexFitFlag(true);

	//GET THE GEOMETRY
	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry* locGeometry = locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());

	//TARGET INFORMATION
	double locTargetCenterZ = 65.0;
	locGeometry->GetTargetZ(locTargetCenterZ);
	dTargetCenter.SetXYZ(0.0, 0.0, locTargetCenterZ);
}


void DSourceComboVertexer::Calc_VertexTimeOffsets(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionCombo)
{
	//we can't calculate the time yet: we need the RF time.  however, we can compute the time offset FROM the RF time
	auto locStepVertexInfos = locReactionVertexInfo->Get_StepVertexInfos();
	auto locReaction = locReactionVertexInfo->Get_Reaction();

	//loop through vertices, determining initial guesses
	map<pair<int, int>, const DKinematicData*> locReconDecayParticleMap; //decaying particle indices -> kinematic data //indices: when the decaying particle is in the INITIAL state
	unordered_map<const DReactionStepVertexInfo*, const DSourceCombo*> locVertexPrimaryComboMap;
	for(const auto& locStepVertexInfo : locStepVertexInfos)
	{
		/***************************************** CHECK IF VERTEX POSITION IS INDETERMINATE A THIS STAGE ********************************************/

		if(locStepVertexInfo->Get_DanglingVertexFlag())
		{
			//is forever indeterminate, even with neutrals and beam energy
			//If this is the production vertex, choose the center of the target.  If not, choose the vertex where the decay parent was produced.
			auto locParentVertexInfo = locStepVertexInfo->Get_ParentVertexInfo();
			if(locParentVertexInfo == nullptr) //production vertex
				dVertexTimeOffsets[locReactionCombo].emplace(locStepVertexInfo, std::make_pair(dTargetCenter, 0.0));
			else //decay products
				dVertexTimeOffsets[locReactionCombo].emplace(locStepVertexInfo, dVertexTimeOffsets[locReactionCombo][locParentVertexInfo]);
			continue;
		}

		//get combo & info
		bool locIsProductionVertexFlag = locStepVertexInfo->Get_ProductionVertexFlag();
		auto locSourceCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionCombo, locStepVertexInfo);
		locVertexPrimaryComboMap.emplace(locStepVertexInfo, locSourceCombo);

		//now, just because it isn't dangling, it doesn't mean we have enough information to actually determine the vertex
		//e.g. we may need neutrals or beam energy to define constraining decay particle p4/trajectory
		auto locDeterminableIterator = dVertexDeterminableWithChargedMap.find(locStepVertexInfo);
		if(locDeterminableIterator == dVertexDeterminableWithChargedMap.end())
		{
			//we don't know yet if it is determinable or not.  figure it out

			//get particles
			auto locSourceParticles = DAnalysis::Get_SourceParticles_ThisVertex(locSourceCombo);
			auto locDefinedDecayingParticles = Get_FullConstrainDecayingParticles(locStepVertexInfo, locReconDecayParticleMap);
			auto locNumConstrainingParticles = locSourceParticles.size() + locDefinedDecayingParticles.size();

			//determine it, save result, and continue if can't
			bool locVertexDeterminableFlag = locIsProductionVertexFlag ? (locNumConstrainingParticles > 0) : (locNumConstrainingParticles > 1);
			dVertexDeterminableWithChargedMap.emplace(locStepVertexInfo, locVertexDeterminableFlag);
			if(!locVertexDeterminableFlag)
				continue; //can't determine yet, but will be able to in the future
		}
		if(locDeterminableIterator->second == false)
			continue;

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

		//Get detected charged track hypotheses at this vertex
		vector<const DKinematicData*> locVertexParticles;
		locVertexParticles.reserve(locSourceParticles.size());

		auto Get_Hypothesis = [](const pair<Particle_t, const JObject*>& locPair) -> bool
			{return static_cast<const DChargedTrack*>(locPair.second)->Get_Hypothesis(locPair.first);};
		std::transform(locVertexParticles.begin(), locVertexParticles.end(), std::back_inserter(locVertexParticles), Get_Hypothesis);

		DVector3 locVertex;
		if(locIsProductionVertexFlag)
		{
			//use track with theta nearest 90 degrees
			auto locThetaNearest90Iterator = Get_ThetaNearest90Iterator(locVertexParticles);
			double locThetaNearest90 = (*locThetaNearest90Iterator)->momentum().Theta()*180.0/TMath::Pi();
			if(locThetaNearest90 < dMinThetaForVertex)
			{
				//try decaying particles instead
				auto locDecayingParticles = Get_FullConstrainDecayingParticles(locStepVertexInfo, locReconDecayParticleMap);
				auto locThetaNearest90Iterator_Decaying = Get_ThetaNearest90Iterator(locDecayingParticles);
				double locLargestTheta_Decaying = (*locThetaNearest90Iterator_Decaying)->momentum().Theta()*180.0/TMath::Pi();
				if(locLargestTheta_Decaying > locThetaNearest90)
					locThetaNearest90Iterator = locThetaNearest90Iterator_Decaying;
			}
			locVertex = (*locThetaNearest90Iterator)->position();
		}
		else //detached vertex
		{
			locVertexParticles.insert(locVertexParticles.end(), locDefinedDecayingParticles.begin(), locDefinedDecayingParticles.end());
			locVertex = dAnalysisUtilities->Calc_CrudeVertex(locVertexParticles);
		}

		/************************************************************ UPDATE WITH RESULTS ************************************************************/

		dVertexTimeOffsets[locReactionCombo].emplace(locStepVertexInfo, std::make_pair(locVertex, 0.0));
		Construct_DecayingParticle(locStepVertexInfo, locSourceCombo, locVertex, locReconDecayParticleMap);
	}

	//do time offsets once all the vertices have been found
	Calc_TimeOffsets(locReactionCombo);

	//set time offsets by combo
	for(const auto& locVertexPrimaryComboPair : locVertexPrimaryComboMap)
	{
		auto locVertexTimeOffsetPair = dVertexTimeOffsets[locReactionCombo][locVertexPrimaryComboPair.first];
		Set_VertexTimeOffsets_ByCombo(locVertexPrimaryComboPair.second, locVertexTimeOffsetPair);
	}

	//free resources
	for(auto locDecayParticlePair : locReconDecayParticleMap)
		delete locDecayParticlePair.second;
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

		//create a new one
		//see if it's already been calculated
		auto locP4 = dSourceComboP4Handler->Get_P4(locSourceCombo);
		auto locKinematicData = new DKinematicData(locDecayPID, locP4.Vect(), locVertex, 0.0);
		//register it
		locReconDecayParticleMap[locParticlePair] = locKinematicData;
	}
}

void DSourceComboVertexer::Calc_TimeOffsets(const DSourceCombo* locReactionCombo)
{
	for(auto& locVertexPair : dVertexTimeOffsets[locReactionCombo])
	{
		auto locStepVertexInfo = locVertexPair.first;
		if(locStepVertexInfo->Get_ProductionVertexFlag())
			continue; //offset is already 0 by default
		auto locVertex = locVertexPair.second.first;

		auto locParentVertexInfo = locStepVertexInfo->Get_ParentVertexInfo();
		auto locParentProductionVertex = dVertexTimeOffsets[locReactionCombo][locParentVertexInfo].first;
		auto locPathLength = (locParentProductionVertex - locVertex).Mag();

		auto locSourceCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionCombo, locStepVertexInfo);
		auto locP4 = dSourceComboP4Handler->Get_P4(locSourceCombo);

		locVertexPair.second.second = locPathLength/(locP4.Beta()*SPEED_OF_LIGHT); //time offset
	}
}

} //end DAnalysis namespace

