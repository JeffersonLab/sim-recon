#ifndef DSourceComboVertexer_h
#define DSourceComboVertexer_h

#include <set>
#include <unordered_map>
#include <utility>
#include <memory>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <map>

#include "JANA/JEventLoop.h"

#include "particleType.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DNeutralShower.h"
#include "ANALYSIS/DSourceCombo.h"
#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DReactionStepVertexInfo.h"
#include "ANALYSIS/DKinFitUtils_GlueX.h"
#include "ANALYSIS/DAnalysisUtilities.h"

#include "PID/DVertex.h"
using namespace std;

namespace DAnalysis
{

class DSourceComboer;
class DSourceComboP4Handler;
class DSourceComboTimeHandler;

class DSourceComboVertexer
{
	public:

		//CONSTRUCTORS
		DSourceComboVertexer(void) = delete;
		DSourceComboVertexer(JEventLoop* locEventLoop, DSourceComboer* locSourceComboer, DSourceComboP4Handler* locSourceComboP4Handler);
		void Reset(void);
void Set_Vertex(const DVertex* locVertex){dVertex = locVertex;}
		//SETUP
		void Set_SourceComboTimeHandler(const DSourceComboTimeHandler* locSourceComboTimeHandler){dSourceComboTimeHandler = locSourceComboTimeHandler;}
		void Set_DebugLevel(int locDebugLevel){dDebugLevel = locDebugLevel;}

		//COMPUTE
		void Calc_VertexTimeOffsets_WithCharged(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionChargedCombo);
		void Calc_VertexTimeOffsets_WithPhotons(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionChargedCombo, const DSourceCombo* locReactionFullCombo);
		void Calc_VertexTimeOffsets_WithBeam(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, const DKinematicData* locBeamParticle);

		bool Get_VertexDeterminableWithCharged(const DReactionStepVertexInfo* locStepVertexInfo) const;
		bool Get_VertexDeterminableWithPhotons(const DReactionStepVertexInfo* locStepVertexInfo) const;

		//GET RESULTS
		bool Get_IsVertexKnown(bool locIsProductionVertex, const DSourceCombo* locSourceCombo, const DKinematicData* locBeamParticle) const;
		DVector3 Get_Vertex(bool locIsProductionVertex, const DSourceCombo* locSourceCombo, const DKinematicData* locBeamParticle) const;
		DVector3 Get_Vertex(bool locIsProductionVertex, const vector<const DKinematicData*>& locVertexParticles) const{return dVertexMap.find(std::make_pair(locIsProductionVertex, locVertexParticles))->second;}
		double Get_TimeOffset(bool locIsPrimaryProductionVertex, const DSourceCombo* locReactionCombo, const DSourceCombo* locVertexCombo, const DKinematicData* locBeamParticle) const;
		DVector3 Get_PrimaryVertex(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionCombo, const DKinematicData* locBeamParticle) const;
		vector<const DKinematicData*> Get_ConstrainingParticles(bool locIsProductionVertex, const DSourceCombo* locVertexCombo, const DKinematicData* locBeamParticle = nullptr) const; //beam = nullptr if not determined by missing mass! //this is empty if vertex not found yet!

		//GET VERTEX-Z BINS
		signed char Get_VertexZBin(bool locIsProductionVertex, const DSourceCombo* locSourceCombo, const DKinematicData* locBeamParticle) const;
		signed char Get_VertexZBin(const DReactionStepVertexInfo* locStepVertexInfo, const DSourceCombo* locReactionCombo, const DKinematicData* locBeamParticle) const;
		vector<signed char> Get_VertexZBins(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionCombo, const DKinematicData* locBeamParticle) const;

		bool Get_IsVertexFoundFlag(bool locIsProductionVertex, const DSourceCombo* locVertexPrimaryCombo, const DKinematicData* locBeamParticle) const;

	private:
		vector<const DKinematicData*>::const_iterator Get_ThetaNearest90Iterator(const vector<const DKinematicData*>& locParticles);
		vector<const DKinematicData*> Get_FullConstrainDecayingParticles(const DReactionStepVertexInfo* locStepVertexInfo, const map<pair<int, int>, const DKinematicData*>& locReconDecayParticleMap);

		DVector3 Calc_Vertex(bool locIsProductionVertexFlag, const vector<pair<Particle_t, const JObject*>>& locChargedSourceParticles, const vector<const DKinematicData*>& locDecayingParticles, vector<const DKinematicData*>& locVertexParticles);
		void Calc_TimeOffsets(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locChargedReactionCombo, const DSourceCombo* locFullReactionCombo = nullptr);

		void Construct_DecayingParticle_InvariantMass(const DReactionStepVertexInfo* locReactionStepVertexInfo, const DSourceCombo* locVertexCombo, DVector3 locVertex, map<pair<int, int>, const DKinematicData*>& locReconDecayParticleMap);
		void Construct_DecayingParticle_MissingMass(const DReactionStepVertexInfo* locReactionStepVertexInfo, const DSourceCombo* locReactionFullCombo, const DSourceCombo* locFullVertexCombo, const DKinematicData* locBeamParticle, DVector3 locVertex, int locRFBunch, double locRFVertexTime, map<pair<int, int>, const DKinematicData*>& locReconDecayParticleMap);

		//HANDLERS/ETC.
		DSourceComboer* dSourceComboer;
		DSourceComboP4Handler* dSourceComboP4Handler;
		const DSourceComboTimeHandler* dSourceComboTimeHandler = nullptr;
		const DAnalysisUtilities* dAnalysisUtilities;
		int dDebugLevel = 0;

		//EXPERIMENT INFORMATION
		DVector3 dTargetCenter;
const DVertex* dVertex;
		double dMinThetaForVertex = 30.0;

		//DETERMINABILITY
		unordered_map<const DReactionStepVertexInfo*, bool> dVertexDeterminableWithChargedMap; //excludes dangling vertex infos!! //only includes primary combos at each vertex
		unordered_map<const DReactionStepVertexInfo*, bool> dVertexDeterminableWithPhotonsMap; //excludes determinable-by-charged & dangling vertex infos!! //only includes primary combos at each vertex

		//time offsets depend on the ENTIRE reaction combo, not just the downstream ones! //time offset is from the RF time
		//bool: is the PRIMARY vertex a production vertex //why is bool used throughout here?
		//because the vertexing is different whether it's a production vertex or not, and a given combo of particles can be used either way
		map<tuple<bool, const DSourceCombo*, const DKinematicData*>, unordered_map<const DSourceCombo*, double>> dTimeOffsets; //first combo: primary reaction combo //kinematics: beam particle

		//VERTICES
		//Note that this only includes the particles used to find the vertex (+ rarely an extra or 2), not necessarily ALL of those at the vertex
		//bool: is production vertex, kinematicdata: beam particle
		map<tuple<bool, const DSourceCombo*, const DKinematicData*>, vector<const DKinematicData*>> dConstrainingParticlesByCombo;
		map<pair<bool, vector<const DKinematicData*>>, DVector3> dVertexMap; //vector: from dConstrainingParticlesByCombo

		//Reconstructed Decaying Particles
		map<tuple<Particle_t, bool, const DSourceCombo*, const DKinematicData*>, const DKinematicData*> dReconDecayParticles_FromProducts; //pid, is prod vertex, full vertex combo, beam particle
		map<tuple<Particle_t, const DSourceCombo*, bool, const DSourceCombo*, const DKinematicData*>, const DKinematicData*> dReconDecayParticles_FromMissing; //decay pid, full reaction combo, is prod vertex, full vertex combo, beam particle (this order)

		//RESOURCE POOL
		DResourcePool<DKinematicData> dResourcePool_KinematicData;
};

inline void DSourceComboVertexer::Reset(void)
{
	dConstrainingParticlesByCombo.clear();
	dVertexMap.clear();
	dTimeOffsets.clear();

	for(const auto& locParticlePair : dReconDecayParticles_FromProducts)
		dResourcePool_KinematicData.Recycle(locParticlePair.second);
	for(const auto& locParticlePair : dReconDecayParticles_FromMissing)
		dResourcePool_KinematicData.Recycle(locParticlePair.second);

	dReconDecayParticles_FromProducts.clear();
	dReconDecayParticles_FromMissing.clear();

	//undeterminable vertices
	dVertexMap.emplace(std::make_pair(false, vector<const DKinematicData*>()), dTargetCenter);
	dVertexMap.emplace(std::make_pair(true, vector<const DKinematicData*>()), dTargetCenter);
}

inline double DSourceComboVertexer::Get_TimeOffset(bool locIsPrimaryProductionVertex, const DSourceCombo* locReactionCombo, const DSourceCombo* locVertexCombo, const DKinematicData* locBeamParticle) const
{
	//the data member MAY be dependent on the beam particle, but it may not
	//so, first search with the beam particle; if not found then search without it
	auto locIterator = dTimeOffsets.find(std::make_tuple(locIsPrimaryProductionVertex, locReactionCombo, locBeamParticle));
	if((locBeamParticle != nullptr) && (locIterator == dTimeOffsets.end()))
		locIterator = dTimeOffsets.find(std::make_tuple(locIsPrimaryProductionVertex, locReactionCombo, nullptr));
	return ((locIterator != dTimeOffsets.end()) ? locIterator->second.find(locVertexCombo)->second : 0.0);
}

inline vector<const DKinematicData*> DSourceComboVertexer::Get_ConstrainingParticles(bool locIsProductionVertex, const DSourceCombo* locVertexCombo, const DKinematicData* locBeamParticle) const
{
	//the data member MAY be dependent on the beam particle, but it may not
	//so, first search with the beam particle; if not found then search without it
	auto locIterator = dConstrainingParticlesByCombo.find(std::make_tuple(locIsProductionVertex, locVertexCombo, locBeamParticle));
	if((locBeamParticle != nullptr) && (locIterator == dConstrainingParticlesByCombo.end()))
		locIterator = dConstrainingParticlesByCombo.find(std::make_tuple(locIsProductionVertex, locVertexCombo, (const DKinematicData*)nullptr));
	if(locIterator != dConstrainingParticlesByCombo.end())
		return locIterator->second;
	return {};
}

inline bool DSourceComboVertexer::Get_IsVertexKnown(bool locIsProductionVertex, const DSourceCombo* locSourceCombo, const DKinematicData* locBeamParticle) const
{
	return !Get_ConstrainingParticles(locIsProductionVertex, locSourceCombo, locBeamParticle).empty();
}

inline DVector3 DSourceComboVertexer::Get_Vertex(bool locIsProductionVertex, const DSourceCombo* locSourceCombo, const DKinematicData* locBeamParticle) const
{
	return Get_Vertex(locIsProductionVertex, Get_ConstrainingParticles(locIsProductionVertex, locSourceCombo, locBeamParticle));
}

inline DVector3 DSourceComboVertexer::Get_PrimaryVertex(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionCombo, const DKinematicData* locBeamParticle) const
{
	auto locIsProductionVertex = locReactionVertexInfo->Get_StepVertexInfo(0)->Get_ProductionVertexFlag();
	return Get_Vertex(locIsProductionVertex, locReactionCombo, locBeamParticle);
}

inline vector<const DKinematicData*>::const_iterator DSourceComboVertexer::Get_ThetaNearest90Iterator(const vector<const DKinematicData*>& locParticles)
{
	//true if first less than second
	auto Get_Nearer90Theta = [](const DKinematicData* lhs, const DKinematicData* rhs) -> bool
		{return fabs(rhs->momentum().Theta() - 0.5*TMath::Pi()) < fabs(lhs->momentum().Theta() - 0.5*TMath::Pi());};
	return std::max_element(locParticles.begin(), locParticles.end(), Get_Nearer90Theta);
}

inline bool DSourceComboVertexer::Get_IsVertexFoundFlag(bool locIsProductionVertex, const DSourceCombo* locVertexPrimaryCombo, const DKinematicData* locBeamParticle) const
{
	//the data member MAY be dependent on the beam particle, but it may not
	//so, first search with the beam particle; if not found then search without it
	auto locIterator = dConstrainingParticlesByCombo.find(std::make_tuple(locIsProductionVertex, locVertexPrimaryCombo, locBeamParticle));
	if((locBeamParticle != nullptr) && (locIterator == dConstrainingParticlesByCombo.end()))
		locIterator = dConstrainingParticlesByCombo.find(std::make_tuple(locIsProductionVertex, locVertexPrimaryCombo, nullptr));
	return (locIterator != dConstrainingParticlesByCombo.end());
}

inline bool DSourceComboVertexer::Get_VertexDeterminableWithCharged(const DReactionStepVertexInfo* locStepVertexInfo) const
{
	auto locIterator = dVertexDeterminableWithChargedMap.find(locStepVertexInfo);
	if(locIterator == dVertexDeterminableWithChargedMap.end())
		return false;
	return locIterator->second;
}

inline bool DSourceComboVertexer::Get_VertexDeterminableWithPhotons(const DReactionStepVertexInfo* locStepVertexInfo) const
{
	auto locIterator = dVertexDeterminableWithPhotonsMap.find(locStepVertexInfo);
	if(locIterator == dVertexDeterminableWithPhotonsMap.end())
		return false;
	return locIterator->second;
}

} //end DAnalysis namespace

#endif // DSourceComboVertexer_h
