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
void Set_Vertex(const DVertex* locVertex){dVertex = locVertex;} //COMPARE
		//SETUP
		void Set_SourceComboTimeHandler(const DSourceComboTimeHandler* locSourceComboTimeHandler){dSourceComboTimeHandler = locSourceComboTimeHandler;}
		void Set_DebugLevel(int locDebugLevel){dDebugLevel = locDebugLevel;}

		//COMPUTE
		void Calc_VertexTimeOffsets_WithCharged(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionChargedCombo);
		void Calc_VertexTimeOffsets_WithPhotons(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionChargedCombo, const DSourceCombo* locReactionFullCombo);
		void Calc_VertexTimeOffsets_WithBeam(const DReactionVertexInfo* locReactionVertexInfo, const DSourceComboUse& locReactionFullComboUse, const DSourceCombo* locReactionFullCombo, const DKinematicData* locBeamParticle);

		bool Get_VertexDeterminableWithCharged(const DReactionStepVertexInfo* locStepVertexInfo) const;
		bool Get_VertexDeterminableWithPhotons(const DReactionStepVertexInfo* locStepVertexInfo) const;

		//GET IS-KNOWN
		bool Get_IsVertexKnown(bool locIsProductionVertex, const DSourceCombo* locReactionCombo, const DSourceCombo* locVertexCombo, const DKinematicData* locBeamParticle, bool locIsCombo2ndVertex) const;
		bool Get_IsVertexKnown_NoBeam(bool locIsProductionVertex, const DSourceCombo* locVertexCombo, bool locIsCombo2ndVertex) const;
		bool Get_IsTimeOffsetKnown(bool locIsPrimaryProductionVertex, const DSourceCombo* locReactionCombo, const DSourceCombo* locVertexCombo, const DKinematicData* locBeamParticle) const;

		//GET VERTEX
		DVector3 Get_Vertex_NoBeam(bool locIsProductionVertex, const DSourceCombo* locVertexCombo, bool locIsCombo2ndVertex) const;
		DVector3 Get_Vertex(bool locIsProductionVertex, const DSourceCombo* locReactionCombo, const DSourceCombo* locVertexCombo, const DKinematicData* locBeamParticle, bool locIsCombo2ndVertex) const;
		DVector3 Get_Vertex(bool locIsProductionVertex, const vector<const DKinematicData*>& locVertexParticles) const;
		DVector3 Get_Vertex(const DReactionStepVertexInfo* locStepVertexInfo, const DSourceCombo* locReactionCombo, const DKinematicData* locBeamParticle, bool locComboIsFullyCharged) const; //Only one that will handle dangling correctly!
		DVector3 Get_PrimaryVertex(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionCombo, const DKinematicData* locBeamParticle) const;

		//GET TIME OFFSET
		double Get_TimeOffset(bool locIsPrimaryProductionVertex, const DSourceCombo* locReactionCombo, const DSourceCombo* locVertexCombo, const DKinematicData* locBeamParticle) const;
		double Get_TimeOffset(const DReactionVertexInfo* locReactionVertexInfo, const DReactionStepVertexInfo* locStepVertexInfo, const DSourceCombo* locReactionCombo, const DKinematicData* locBeamParticle) const; //Only one that will handle dangling correctly!

		//GET CONSTRAINING PARTICLES (for vertex)
		vector<const DKinematicData*> Get_ConstrainingParticles(bool locIsProductionVertex, const DSourceCombo* locReactionCombo, const DSourceCombo* locVertexCombo, const DKinematicData* locBeamParticle, bool locIsCombo2ndVertex) const;
		vector<const DKinematicData*> Get_ConstrainingParticles_NoBeam(bool locIsProductionVertex, const DSourceCombo* locVertexCombo, bool locIsCombo2ndVertex) const;

		//GET VERTEX-Z BINS
		signed char Get_VertexZBin(bool locIsProductionVertex, const DSourceCombo* locReactionCombo, const DSourceCombo* locPrimaryVertexCombo, const DKinematicData* locBeamParticle, bool locIsCombo2ndVertex) const;
		signed char Get_VertexZBin_NoBeam(bool locIsProductionVertex, const DSourceCombo* locPrimaryVertexCombo, bool locIsCombo2ndVertex) const;
		signed char Get_VertexZBin(const DReactionStepVertexInfo* locStepVertexInfo, const DSourceCombo* locReactionCombo, const DKinematicData* locBeamParticle, bool locComboIsFullyCharged) const; //Only one that will handle dangling correctly!
		vector<signed char> Get_VertexZBins(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionCombo, const DKinematicData* locBeamParticle, bool locComboIsFullyCharged) const; //This will call the above

	private:
		vector<const DKinematicData*>::const_iterator Get_ThetaNearest90Iterator(const vector<const DKinematicData*>& locParticles);
		vector<const DKinematicData*> Get_FullConstrainDecayingParticles(const DReactionStepVertexInfo* locStepVertexInfo, const map<pair<int, int>, const DKinematicData*>& locReconDecayParticleMap);

		DVector3 Calc_Vertex(bool locIsProductionVertexFlag, const vector<pair<Particle_t, const JObject*>>& locChargedSourceParticles, const vector<const DKinematicData*>& locDecayingParticles, vector<const DKinematicData*>& locVertexParticles);
		void Calc_TimeOffsets(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locChargedReactionCombo, const DSourceCombo* locFullReactionCombo = nullptr);

		void Construct_DecayingParticle_InvariantMass(const DReactionStepVertexInfo* locReactionStepVertexInfo, const DSourceCombo* locVertexCombo, DVector3 locVertex, map<pair<int, int>, const DKinematicData*>& locReconDecayParticleMap);
		void Construct_DecayingParticle_MissingMass(const DReactionStepVertexInfo* locReactionStepVertexInfo, const DSourceComboUse& locReactionFullComboUse, const DSourceCombo* locReactionFullCombo, const DSourceCombo* locFullVertexCombo, const DKinematicData* locBeamParticle, DVector3 locVertex, int locRFBunch, double locRFVertexTime, map<pair<int, int>, const DKinematicData*>& locReconDecayParticleMap);

		//HANDLERS/ETC.
		DSourceComboer* dSourceComboer;
		DSourceComboP4Handler* dSourceComboP4Handler;
		const DSourceComboTimeHandler* dSourceComboTimeHandler = nullptr;
		const DAnalysisUtilities* dAnalysisUtilities;
		int dDebugLevel = 0;

		//EXPERIMENT INFORMATION
		DVector3 dTargetCenter;
const DVertex* dVertex; //COMPARE
		double dMinThetaForVertex = 30.0;

		//DETERMINABILITY
		unordered_map<const DReactionStepVertexInfo*, bool> dVertexDeterminableWithChargedMap; //excludes dangling vertex infos!! //only includes primary combos at each vertex
		unordered_map<const DReactionStepVertexInfo*, bool> dVertexDeterminableWithPhotonsMap; //excludes determinable-by-charged & dangling vertex infos!! //only includes primary combos at each vertex

		//TIME OFFSETS
		//time offsets & (sometimes) vertices depend on the ENTIRE reaction combo, not just the downstream ones! //time offset is from the RF time
		//bool: is the PRIMARY vertex a production vertex //why is bool used throughout here?
			//because the vertexing is different whether it's a production vertex or not, and a given combo of particles can be used either way
		map<tuple<bool, const DSourceCombo*, const DKinematicData*>, unordered_map<const DSourceCombo*, double>> dTimeOffsets; //first combo: primary reaction combo //kinematics: beam particle

		//VERTEX-CONSTRAINING PARTICLES
		//kinematic data: beam particle
		//first bool: is vertex combo production-vertex flag.  However, this flag is true for EVERY vertex that the beam is needed to find
		//second bool: is vertex combo a charged combo re-used in place of a has-neutrals combo (e.g. Xi0 -> pi0, Lambda: The lambda combo is used for Xi0 during charged-only stage)
			//false if actual slot (e.g. Lambda decay), true if re-used (e.g. Xi0 decay) //only true in charged-only stage!
			//when looping in dependency-order over vertex-infos, it will always be false first, then true //only true in charged-only stage!
			//worst-case?: g, p -> K0, Sigma+,  Sigma+ -> pi0, (p)
				//the pi+, pi- combo is used for the K0 decay, the K0 is used to find the production vertex
				//AND the production vertex is used for the Sigma+ decay vertex (is dangling, indeterminate): 3 vertices!
				//however, the 3rd vertex result is the same as the 2nd since copied: call with bool = true
		//If beam is NOT needed to find the vertex then the primary reaction combo (first combo) is nullptr!!! (also not needed)
		map<tuple<bool, const DSourceCombo*, const DSourceCombo*, const DKinematicData*, bool>, vector<const DKinematicData*>> dConstrainingParticlesByCombo; //first combo: primary reaction combo

		//VERTICES
		//Note that this only includes the particles used to find the vertex (+ rarely an extra or 2), not necessarily ALL of those at the vertex
		//bool: is production vertex
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

inline bool DSourceComboVertexer::Get_IsTimeOffsetKnown(bool locIsPrimaryProductionVertex, const DSourceCombo* locReactionCombo, const DSourceCombo* locVertexCombo, const DKinematicData* locBeamParticle) const
{
	//the data member MAY be dependent on the beam particle, but it may not
	//so, first search with the beam particle; if not found then search without it
	if(locBeamParticle == nullptr)
	{
		auto locIterator = dTimeOffsets.find(std::make_tuple(locIsPrimaryProductionVertex, locReactionCombo, nullptr));
		if(locIterator == dTimeOffsets.end())
			return false;
		auto& locComboMap = locIterator->second;
		auto locComboIterator = locComboMap.find(locVertexCombo);
		return (locComboIterator != locComboMap.end());
	}

	auto locIterator = dTimeOffsets.find(std::make_tuple(locIsPrimaryProductionVertex, locReactionCombo, locBeamParticle));
	if(locIterator == dTimeOffsets.end())
		return Get_TimeOffset(locIsPrimaryProductionVertex, locReactionCombo, locVertexCombo, nullptr); //try without beam

	auto& locComboMap = locIterator->second;
	auto locComboIterator = locComboMap.find(locVertexCombo);
	if(locComboIterator == locComboMap.end())
		return Get_TimeOffset(locIsPrimaryProductionVertex, locReactionCombo, locVertexCombo, nullptr); //try without beam

	return true;
}

inline double DSourceComboVertexer::Get_TimeOffset(bool locIsPrimaryProductionVertex, const DSourceCombo* locReactionCombo, const DSourceCombo* locVertexCombo, const DKinematicData* locBeamParticle) const
{
	//the data member MAY be dependent on the beam particle, but it may not
	//so, first search with the beam particle; if not found then search without it
	if(locBeamParticle == nullptr)
	{
		auto locIterator = dTimeOffsets.find(std::make_tuple(locIsPrimaryProductionVertex, locReactionCombo, nullptr));
		if(locIterator == dTimeOffsets.end())
			return 0.0;
		auto& locComboMap = locIterator->second;
		auto locComboIterator = locComboMap.find(locVertexCombo);
		return ((locComboIterator != locComboMap.end()) ? locComboIterator->second : 0.0);
	}

	auto locIterator = dTimeOffsets.find(std::make_tuple(locIsPrimaryProductionVertex, locReactionCombo, locBeamParticle));
	if(locIterator == dTimeOffsets.end())
		return Get_TimeOffset(locIsPrimaryProductionVertex, locReactionCombo, locVertexCombo, nullptr); //try without beam

	auto& locComboMap = locIterator->second;
	auto locComboIterator = locComboMap.find(locVertexCombo);
	if(locComboIterator == locComboMap.end())
		return Get_TimeOffset(locIsPrimaryProductionVertex, locReactionCombo, locVertexCombo, nullptr); //try without beam

	return locComboIterator->second;
}

inline vector<const DKinematicData*> DSourceComboVertexer::Get_ConstrainingParticles(bool locIsProductionVertex, const DSourceCombo* locReactionCombo, const DSourceCombo* locVertexCombo, const DKinematicData* locBeamParticle, bool locIsCombo2ndVertex) const
{
	//the data member MAY be dependent on the beam particle, but it may not
	//so, first search with the beam particle; if not found then search without it

//cout << "Query tuple: " << locIsProductionVertex << ", " << locReactionCombo << ", " << locVertexCombo << ", " << locBeamParticle << ", " << locIsCombo2ndVertex << endl;

	if(locBeamParticle == nullptr)
	{
		auto locIterator = dConstrainingParticlesByCombo.find(std::make_tuple(locIsProductionVertex, (const DSourceCombo*)nullptr, locVertexCombo, (const DKinematicData*)nullptr, locIsCombo2ndVertex));
		if(locIterator != dConstrainingParticlesByCombo.end())
			return locIterator->second;
		if(!locIsCombo2ndVertex)
			return {};

		//We MAY have the wrong flag for is-prod-vertex, due to dangling vertex on 2nd combo (e.g. worst-case scenario detailed above). try the other one
		locIterator = dConstrainingParticlesByCombo.find(std::make_tuple(!locIsProductionVertex, (const DSourceCombo*)nullptr, locVertexCombo, (const DKinematicData*)nullptr, locIsCombo2ndVertex));
		if(locIterator != dConstrainingParticlesByCombo.end())
			return locIterator->second;
		return {};
	}

	auto locIterator = dConstrainingParticlesByCombo.find(std::make_tuple(true, locReactionCombo, locVertexCombo, locBeamParticle, locIsCombo2ndVertex));
	if(locIterator == dConstrainingParticlesByCombo.end()) //see if beam is not needed
		locIterator = dConstrainingParticlesByCombo.find(std::make_tuple(locIsProductionVertex, (const DSourceCombo*)nullptr, locVertexCombo, (const DKinematicData*)nullptr, locIsCombo2ndVertex));
	if(locIterator != dConstrainingParticlesByCombo.end())
		return locIterator->second;
	return {};
}

inline DVector3 DSourceComboVertexer::Get_Vertex(bool locIsProductionVertex, const vector<const DKinematicData*>& locVertexParticles) const
{
	auto locIterator = dVertexMap.find(std::make_pair(locIsProductionVertex, locVertexParticles));
	if(locIterator != dVertexMap.end())
		return locIterator->second;
	return dTargetCenter;
}

inline signed char DSourceComboVertexer::Get_VertexZBin_NoBeam(bool locIsProductionVertex, const DSourceCombo* locPrimaryVertexCombo, bool locIsCombo2ndVertex) const
{
	return Get_VertexZBin(locIsProductionVertex, nullptr, locPrimaryVertexCombo, nullptr, locIsCombo2ndVertex);
}

inline vector<const DKinematicData*> DSourceComboVertexer::Get_ConstrainingParticles_NoBeam(bool locIsProductionVertex, const DSourceCombo* locVertexCombo, bool locIsCombo2ndVertex) const
{
	return Get_ConstrainingParticles(locIsProductionVertex, nullptr, locVertexCombo, nullptr, locIsCombo2ndVertex);
}

inline bool DSourceComboVertexer::Get_IsVertexKnown(bool locIsProductionVertex, const DSourceCombo* locReactionCombo, const DSourceCombo* locVertexCombo, const DKinematicData* locBeamParticle, bool locIsCombo2ndVertex) const
{
	return !Get_ConstrainingParticles(locIsProductionVertex, locReactionCombo, locVertexCombo, locBeamParticle, locIsCombo2ndVertex).empty();
}

inline bool DSourceComboVertexer::Get_IsVertexKnown_NoBeam(bool locIsProductionVertex, const DSourceCombo* locVertexCombo, bool locIsCombo2ndVertex) const
{
	return !Get_ConstrainingParticles_NoBeam(locIsProductionVertex, locVertexCombo, locIsCombo2ndVertex).empty();
}

//2 cases: you KNOW beam photon MUST be nullptr, and one where it MAY be nullptr
inline DVector3 DSourceComboVertexer::Get_Vertex_NoBeam(bool locIsProductionVertex, const DSourceCombo* locVertexCombo, bool locIsCombo2ndVertex) const
{
	return Get_Vertex(locIsProductionVertex, Get_ConstrainingParticles(locIsProductionVertex, nullptr, locVertexCombo, nullptr, locIsCombo2ndVertex));
}

inline DVector3 DSourceComboVertexer::Get_Vertex(bool locIsProductionVertex, const DSourceCombo* locReactionCombo, const DSourceCombo* locVertexCombo, const DKinematicData* locBeamParticle, bool locIsCombo2ndVertex) const
{
	//bool: for the vertex we want, not the primary
	auto locConstrainingParticles = Get_ConstrainingParticles(true, locReactionCombo, locVertexCombo, locBeamParticle, locIsCombo2ndVertex);
	if(locConstrainingParticles.empty())
		locConstrainingParticles = Get_ConstrainingParticles(locIsProductionVertex, nullptr, locVertexCombo, nullptr, locIsCombo2ndVertex);
	return Get_Vertex(locIsProductionVertex, locConstrainingParticles);
}

inline DVector3 DSourceComboVertexer::Get_PrimaryVertex(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionCombo, const DKinematicData* locBeamParticle) const
{
	auto locStepVertexInfo = locReactionVertexInfo->Get_StepVertexInfo(0);
	auto locIsProductionVertex = locStepVertexInfo->Get_ProductionVertexFlag();
	auto locComboIsFullyCharged = locReactionCombo->Get_SourceParticles(true, d_Neutral).empty();
	auto locIsCombo2ndVertex = (locComboIsFullyCharged && locStepVertexInfo->Get_FullConstrainParticles(false, d_FinalState, d_Charged, false).empty());
	return Get_Vertex(locIsProductionVertex, locReactionCombo, locReactionCombo, locBeamParticle, locIsCombo2ndVertex);
}

inline vector<const DKinematicData*>::const_iterator DSourceComboVertexer::Get_ThetaNearest90Iterator(const vector<const DKinematicData*>& locParticles)
{
	//true if first less than second
	auto Get_Nearer90Theta = [](const DKinematicData* lhs, const DKinematicData* rhs) -> bool
		{return fabs(rhs->momentum().Theta() - 0.5*TMath::Pi()) < fabs(lhs->momentum().Theta() - 0.5*TMath::Pi());};
	return std::max_element(locParticles.begin(), locParticles.end(), Get_Nearer90Theta);
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
