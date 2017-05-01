#ifndef DSourceComboVertexer_h
#define DSourceComboVertexer_h

#include <deque>
#include <set>
#include <unordered_map>
#include <utility>
#include <memory>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>

#include "JANA/JEventLoop.h"

#include "particleType.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DChargedTrack.h"
#include "PID/DNeutralShower.h"
#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DReactionStepVertexInfo.h"
#include "ANALYSIS/DKinFitUtils_GlueX.h"
#include "ANALYSIS/DAnalysisUtilities.h"
#include "ANALYSIS/DSourceComboP4Handler.h"

using namespace std;

namespace DAnalysis
{

class DSourceComboer;

class DSourceComboVertexer
{
	public:

		//CONSTRUCTORS
		DSourceComboVertexer(void) = delete;
		DSourceComboVertexer(JEventLoop* locEventLoop, DSourceComboer* locSourceComboer, DSourceComboP4Handler* locSourceComboP4Handler);
		void Reset(void);

		//SETUP
		void Set_SourceComboTimeHandler(const DSourceComboTimeHandler* locSourceComboTimeHandler){dSourceComboTimeHandler = locSourceComboTimeHandler;}

		//COMPUTE
		void Calc_VertexTimeOffsets_WithCharged(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionChargedCombo);
		void Calc_VertexTimeOffsets_WithPhotons(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionChargedCombo, const DSourceCombo* locReactionFullCombo);
		void Calc_VertexTimeOffsets_WithBeam(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionChargedCombo, const DSourceCombo* locReactionFullCombo, const DKinematicData* locBeamParticle);

		bool Get_VertexDeterminableWithCharged(bool locIsProductionVertex, const DSourceCombo* locSourceCombo) const{return dVertexDeterminableWithChargedMap.find(std::make_pair(locIsProductionVertex, locSourceCombo))->second;}
		bool Get_VertexDeterminableWithPhotons(bool locIsProductionVertex, const DSourceCombo* locSourceCombo) const{return dVertexDeterminableWithPhotonsMap.find(std::make_pair(locIsProductionVertex, locSourceCombo))->second;}

		//GET RESULTS
		DVector3 Get_Vertex(bool locIsProductionVertex, const DSourceCombo* locSourceCombo) const;
		DVector3 Get_Vertex(bool locIsProductionVertex, const vector<const DKinematicData*>& locVertexParticles){return dVertexMap.find(std::make_pair(locIsProductionVertex, locVertexParticles))->second;}
		double Get_TimeOffset(bool locIsProductionVertex, const DSourceCombo* locReactionCombo, const DSourceCombo* locVertexCombo) const;
		DVector3 Get_PrimaryVertex(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionCombo) const;
		vector<const DKinematicData*> Get_ConstrainingParticles(bool locIsProductionVertex, const DSourceCombo* locVertexCombo) const; //this is empty if vertex not found yet!

		//GET VERTEX-Z BINS
		signed char Get_VertexZBin(bool locIsProductionVertex, const DSourceCombo* locSourceCombo) const;
		vector<signed char> Get_VertexZBins(bool locIsProductionVertex, const DSourceCombo* locReactionChargedCombo) const;

		bool Get_IsVertexFoundFlag(bool locIsProductionVertex, const DSourceCombo* locVertexPrimaryChargedCombo, const DSourceCombo* locVertexPrimaryFullCombo = nullptr) const;

	private:
		vector<const DKinematicData*>::const_iterator Get_ThetaNearest90Iterator(const vector<const DKinematicData*>& locParticles);
		vector<const DKinematicData*> Get_FullConstrainDecayingParticles(const DReactionStepVertexInfo* locStepVertexInfo, const map<pair<int, int>, const DKinematicData*>& locReconDecayParticleMap);

		DVector3 Calc_Vertex(bool locIsProductionVertexFlag, const vector<pair<Particle_t, const JObject*>>& locChargedSourceParticles, const vector<const DKinematicData*>& locDecayingParticles, vector<const DKinematicData*>& locVertexParticles);
		void Calc_TimeOffsets(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locChargedReactionCombo, const DSourceCombo* locFullReactionCombo = nullptr);

		void Construct_DecayingParticle_InvariantMass(const DReactionStepVertexInfo* locReactionStepVertexInfo, const DSourceCombo* locChargedVertexCombo, const DSourceCombo* locFullVertexCombo, DVector3 locVertex, map<pair<int, int>, const DKinematicData*>& locReconDecayParticleMap);
		void Construct_DecayingParticle_MissingMass(const DReactionStepVertexInfo* locReactionStepVertexInfo, const DSourceCombo* locReactionChargedCombo, const DSourceCombo* locReactionFullCombo, const DSourceCombo* locFullVertexCombo, const DKinematicData* locBeamParticle, DVector3 locVertex, int locRFBunch, double locRFVertexTime, map<pair<int, int>, const DKinematicData*>& locReconDecayParticleMap);

		DSourceComboer* dSourceComboer;
		DSourceComboP4Handler* dSourceComboP4Handler;
		const DSourceComboTimeHandler* dSourceComboTimeHandler = nullptr;
		const DAnalysisUtilities* dAnalysisUtilities;

		//EXPERIMENT INFORMATION
		DVector3 dTargetCenter;
		double dMinThetaForVertex = 30.0;

		//VERTICES AND TIME OFFSETS
		//bool: is production vertex
		unordered_map<pair<bool, const DSourceCombo*>, bool> dVertexDeterminableWithChargedMap; //excludes dangling vertex infos!! //only includes primary combos at each vertex
		unordered_map<pair<bool, const DSourceCombo*>, bool> dVertexDeterminableWithPhotonsMap; //excludes determinable-by-charged & dangling vertex infos!! //only includes primary combos at each vertex

		//time offsets depend on the ENTIRE reaction combo, not just the downstream ones! //time offset is from the RF time
		//bool: is production vertex
		unordered_map<pair<bool, const DSourceCombo*>, unordered_map<const DSourceCombo*, double>> dTimeOffsets; //first combo: primary reaction combo

		//Only works if beam not needed! (e.g. not by missing mass)
		//Note that this only includes the particles used to find the vertex (+ rarely an extra or 2), not necessarily ALL of those at the vertex
		//bool: is production vertex
		unordered_map<pair<bool, const DSourceCombo*>, vector<const DKinematicData*>> dConstrainingParticlesByCombo;
		unordered_map<pair<bool, vector<const DKinematicData*>>, DVector3> dVertexMap; //vector: from dConstrainingParticlesByCombo

		//Reconstructed Decaying Particles
		unordered_map<DSourceComboUse, const DKinematicData*> dReconDecayParticles_FromProducts; //use: where it decays
		unordered_map<pair<const DSourceCombo*, DSourceComboUse>, const DKinematicData*> dReconDecayParticles_FromMissing; //combo: full reaction combo //use: where it decays
};

inline void DSourceComboVertexer::Reset(void)
{
	dConstrainingParticlesByCombo.clear();
	dVertexMap.clear();
	dTimeOffsets.clear();

	//delete or recycle these!
	dReconDecayParticles_FromProducts.clear();
	dReconDecayParticles_FromMissing.clear();

	//undeterminable vertices
	dVertexMap.emplace({false, {}}, dTargetCenter);
	dVertexMap.emplace({true, {}}, dTargetCenter);
}

inline double DSourceComboVertexer::Get_TimeOffset(bool locIsProductionVertex, const DSourceCombo* locReactionCombo, const DSourceCombo* locVertexCombo) const
{
	return dTimeOffsets.find(std::make_pair(locIsProductionVertex, locReactionCombo))->second.find(locVertexCombo)->second;
}

inline vector<const DKinematicData*> DSourceComboVertexer::Get_ConstrainingParticles(bool locIsProductionVertex, const DSourceCombo* locVertexCombo) const
{
	auto locIterator = dConstrainingParticlesByCombo.find(std::make_pair(locIsProductionVertex, locVertexCombo));
	if(locIterator != dConstrainingParticlesByCombo.end())
		return locIterator->second;
	return {};
}

inline DVector3 DSourceComboVertexer::Get_Vertex(bool locIsProductionVertex, const DSourceCombo* locSourceCombo) const
{
	return Get_Vertex(locIsProductionVertex, Get_ConstrainingParticles(locIsProductionVertex, locSourceCombo));
}

inline DVector3 DSourceComboVertexer::Get_PrimaryVertex(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionCombo) const
{
	auto locIsProductionVertex = locReactionVertexInfo->Get_StepVertexInfo(0)->Get_ProductionVertexFlag();
	return Get_Vertex(locIsProductionVertex, locReactionCombo);
}

inline vector<signed char> DSourceComboVertexer::Get_VertexZBins(bool locIsProductionVertex, const DSourceCombo* locReactionChargedCombo) const
{
	vector<signed char> locVertexZBins;
	const auto& locReactionTimeOffsets = dTimeOffsets.find(std::make_pair(locIsProductionVertex, locReactionChargedCombo))->second;
	for(const auto& locComboTimeOffsetPair : locReactionTimeOffsets) //turns out this is fastest: we have the vertex combos!
	{
		locVertexZBins.emplace_back(Get_VertexZBin(locIsProductionVertex, locComboTimeOffsetPair.first));
		locIsProductionVertex = false;
	}
	return locVertexZBins;
}

inline vector<const DKinematicData*>::const_iterator DSourceComboVertexer::Get_ThetaNearest90Iterator(const vector<const DKinematicData*>& locParticles)
{
	auto Get_Nearer90Theta = [](const DKinematicData* lhs, const DKinematicData* rhs) -> bool
		{return fabs(lhs->momentum().Theta() - 0.5*TMath::Pi()) < fabs(rhs->momentum().Theta() - 0.5*TMath::Pi());};
	return std::max_element(locParticles.begin(), locParticles.end(), Get_Nearer90Theta);
}

inline bool DSourceComboVertexer::Get_IsVertexFoundFlag(bool locIsProductionVertex, const DSourceCombo* locVertexPrimaryChargedCombo, const DSourceCombo* locVertexPrimaryFullCombo) const
{
	auto locProductionPair = std::make_pair(locIsProductionVertex, locVertexPrimaryChargedCombo);
	if(dConstrainingParticlesByCombo.find(locProductionPair) != dConstrainingParticlesByCombo.end())
		return true;
	if(locVertexPrimaryFullCombo == nullptr)
		return false;

	locProductionPair = std::make_pair(locIsProductionVertex, locVertexPrimaryFullCombo);
	return (dConstrainingParticlesByCombo.find(locProductionPair) != dConstrainingParticlesByCombo.end());

}

} //end DAnalysis namespace

#endif // DSourceComboVertexer_h
