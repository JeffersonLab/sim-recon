#ifndef DSourceCombo_h
#define DSourceCombo_h

#include <map>
#include <vector>
#include <tuple>
#include <memory>
#include <algorithm>

#include "JANA/JObject.h"
#include "JANA/JEventLoop.h"

#include "particleType.h"

#include "PID/DNeutralShower.h"

using namespace std;
using namespace jana;

namespace DAnalysis
{

/****************************************************** DEFINE LAMBDAS, USING STATEMENTS *******************************************************/

//forward declarations
class DSourceComboInfo;
class DSourceCombo;

//DSourceComboUse is what the combo is USED for (the decay of Particle_t (if Unknown then is just a grouping)
using DSourceComboUse = tuple<Particle_t, signed char, const DSourceComboInfo*>; //e.g. Pi0, -> 2g //signed char: vertex-z bin of the final state (combo contents)
using DSourceCombosByUse_Small = map<DSourceComboUse, vector<const DSourceCombo*>>;

//Compare_SourceComboUses
inline auto Compare_SourceComboUses = [](const DSourceComboUse& lhs, const DSourceComboUse& rhs) -> bool
{
	//this puts mixed-charge first, then fully-neutral, then fully-charged
	auto locChargeContent_LHS = Get_ChargeContent(std::get<2>(lhs));
	auto locChargeContent_RHS = Get_ChargeContent(std::get<2>(rhs));
	if(locChargeContent_LHS != locChargeContent_RHS)
		return locChargeContent_LHS > locChargeContent_RHS;

	//within each of those, it puts the most-massive particles first
	if(std::get<0>(lhs) == std::get<0>(rhs))
	{
		if(std::get<1>(lhs) == std::get<1>(rhs))
			return *std::get<2>(lhs) > *std::get<2>(rhs);
		else
			return std::get<1>(lhs) > std::get<1>(rhs);
	}
	if(ParticleMass(std::get<0>(lhs)) == ParticleMass(std::get<0>(rhs)))
		return std::get<0>(lhs) > std::get<0>(rhs);
	return (ParticleMass(std::get<0>(lhs)) > ParticleMass(std::get<0>(rhs)));
};

/************************************************************** DEFINE CLASSES ***************************************************************/

//In theory, for safety, dynamically-allocated objects should be stored in a shared_ptr
//However, the combos take a TON of memory, and a shared_ptr<T*> takes 3x the memory of a regular T* pointer
//The info objects will exist for the life of the program, and so don't need to be recycled to a resource pool.
//The combo objects will be recycled after every event into a resource pool.

//If we REALLY need the memory, we can store these things in std::array instead of vector
//These classes would need to have template parameters for the array sizes, and have a common base class
//However, the base class would have to return either vectors (need to convert) or raw pointers (need bounds checking) instead of std::arrays (since type unknown to base class)
//So it would require more CPU.

//THE MOST NUMBER OF PARTICLES OF A GIVEN TYPE IS 255 (# stored in unsigned char)
class DSourceComboInfo
{
	public:

		//CONSTRUCTOR
		DSourceComboInfo(void) = delete;
		DSourceComboInfo(const vector<pair<Particle_t, unsigned char>>& locNumParticles, const vector<pair<DSourceComboUse, unsigned char>>& locFurtherDecays = {});

		//OPERATORS
		bool operator< (const DSourceComboInfo& rhs) const;

		//GET MEMBERS
		vector<pair<Particle_t, unsigned char>> Get_NumParticles(bool locEntireChainFlag = false) const{return dNumParticles;}
		vector<pair<DSourceComboUse, unsigned char>> Get_FurtherDecays(void) const{return dFurtherDecays;}

		//definitions of negative values for any particle index //in order such that operator< returns order expected for string (e.g. gp->...)
		static signed char Get_VertexZIndex_FCAL(void){return -2;}
		static signed char Get_VertexZIndex_Unknown(void){return -1;}

	private:

		//FORWARD DECLARE COMPARISON STRUCT
		struct DCompare_ParticlePairPIDs;

		//don't have decaying PID a direct member of this combo info
		//e.g. for a 2g pair, it has no idea whether or not it came from a Pi0, an eta, through direct production, etc.
		//this way, e.g., the 2g combos can be used for ANY of those possibilities, without creating new objects
		//it is the responsibility of the containing object to know what the combos are used for: DSourceComboUse

		vector<pair<Particle_t, unsigned char>> dNumParticles;

		//this will be sorted with Compare_SourceComboUses
		vector<pair<DSourceComboUse, unsigned char>> dFurtherDecays; //unsigned char: # of (e.g.) pi0s, etc.
};

class DSourceCombo
{
	public:

		//CONSTRUCTOR
		DSourceCombo(const vector<pair<Particle_t, const JObject*>>& locSourceParticles, const DSourceCombosByUse_Small& locFurtherDecayCombos, bool locIsZIndependent = false);

		//SET MEMBERS
		void Set_Members(const vector<pair<Particle_t, const JObject*>>& locSourceParticles, const DSourceCombosByUse_Small& locFurtherDecayCombos, bool locIsZIndependent = false){};

		//GET MEMBERS
		vector<pair<Particle_t, const JObject*>> Get_SourceParticles(bool locEntireChainFlag = false) const;
		DSourceCombosByUse_Small Get_FurtherDecayCombos(void) const{return dFurtherDecayCombos;}
		bool Get_IsComboingZIndependent(void) const{return dIsComboingZIndependent;}

	private:

		//particles & decays
		vector<pair<Particle_t, const JObject*>> dSourceParticles; //original DNeutralShower or DChargedTrack
		//FYI, if use PID is unknown: Vector is guaranteed to be size 1, and none of ITS further decays can have an unknown use PID
		//Note: Either:
		//1) The entire content is one PID
		//2) There is at most one particle of each PID
		//Otherwise, groupings of 2+ (e.g. 2pi-) will be in the further-decay section (X -> 2pi+)
		DSourceCombosByUse_Small dFurtherDecayCombos; //vector is e.g. size 3 if 3 pi0s needed

		//actually, everything is z-dependent for massive neutrals.
		//however the momentum is SO z-dependent, that we can't cut on it until the end when we have the final vertex, AFTER comboing
		//a mere z-bin is not enough.
		//So, as far as COMBOING is concerned, massive neutrals are Z-INDEPENDENT
		bool dIsComboingZIndependent; //is false for BCAL photons
};

struct DSourceComboInfo::DCompare_ParticlePairPIDs
{
	bool operator()(const pair<Particle_t, unsigned char>& lhs, const pair<Particle_t, unsigned char>& rhs) const{return lhs.first < rhs.first;} //sort
	bool operator()(const pair<Particle_t, unsigned char>& lhs, Particle_t rhs) const{return lhs.first < rhs;} //lookup
	bool operator()(Particle_t lhs, const pair<Particle_t, unsigned char>& rhs) const{return lhs < rhs.first;} //lookup
};

/*********************************************************** INLINE MEMBER FUNCTION DEFINITIONS ************************************************************/

inline DSourceComboInfo::DSourceComboInfo(const vector<pair<Particle_t, unsigned char>>& locNumParticles, const vector<pair<DSourceComboUse, unsigned char>>& locFurtherDecays) :
		dNumParticles(locNumParticles), dFurtherDecays(locFurtherDecays)
{
	std::sort(dNumParticles.begin(), dNumParticles.end());
	std::sort(dFurtherDecays.begin(), dFurtherDecays.end(), Compare_SourceComboUses);
}

inline bool DSourceComboInfo::operator< (const DSourceComboInfo& rhs) const
{
	if(dNumParticles != rhs.dNumParticles)
		return dNumParticles < rhs.dNumParticles;

	//check if maps have different sizes
	if(dFurtherDecays.size() != rhs.dFurtherDecays.size())
		return dFurtherDecays.size() < rhs.dFurtherDecays.size();

	//check if there's a mismatch between the maps
	auto locMismachIterators = std::mismatch(dFurtherDecays.begin(), dFurtherDecays.end(), rhs.dFurtherDecays.begin());
	if(locMismachIterators.first == dFurtherDecays.end())
		return false; //maps are identical

	//check if keys are equal
	if(locMismachIterators.first->first == locMismachIterators.second->first)
		return locMismachIterators.first->second < locMismachIterators.second->second; //compare values
	else
		return locMismachIterators.first->first < locMismachIterators.second->first; //compare keys
}

inline vector<pair<Particle_t, unsigned char>> DSourceComboInfo::Get_NumParticles(bool locEntireChainFlag) const
{
	if(!locEntireChainFlag || dFurtherDecays.empty())
		return dNumParticles;

	vector<pair<Particle_t, unsigned char>> locToReturnNumParticles = dNumParticles;
	for(const auto& locDecayPair : dFurtherDecays)
	{
		const auto& locDecayComboInfo = std::get<2>(locDecayPair.first);
		auto locNumDecayParticles = locDecayComboInfo->Get_NumParticles(true);

		for(const auto& locParticlePair : locNumDecayParticles)
		{
			//search through locToReturnNumParticles to retrieve the iterator corresponding to this PID if it's already present
			auto locIteratorPair = std::equal_range(locToReturnNumParticles.begin(), locToReturnNumParticles.end(), locParticlePair.first, DCompare_ParticlePairPIDs());
			if(locIteratorPair.first != locIteratorPair.second)
				(*locIteratorPair.first)->second += locParticlePair.second; //it exists: increase it
			else //doesn't exist. insert it
				locToReturnNumParticles.insert(locIteratorPair.first, locParticlePair);
		}
	}

	return locToReturnNumParticles;
}

inline DSourceCombo::DSourceCombo(const vector<pair<Particle_t, const JObject*>>& locSourceParticles, const DSourceCombosByUse_Small& locFurtherDecayCombos, bool locIsZIndependent) :
		dSourceParticles(locSourceParticles), dFurtherDecayCombos(locFurtherDecayCombos), dIsComboingZIndependent(locIsZIndependent)
{
	std::sort(dSourceParticles.begin(), dSourceParticles.end());
}

inline void DSourceCombo::Set_Members(const vector<pair<Particle_t, const JObject*>>& locSourceParticles, const DSourceCombosByUse_Small& locFurtherDecayCombos, bool locIsZIndependent)
{
	dSourceParticles = locSourceParticles;
	dFurtherDecayCombos = locFurtherDecayCombos;
	dIsComboingZIndependent = locIsZIndependent;
	std::sort(dSourceParticles.begin(), dSourceParticles.end());
}

inline vector<pair<Particle_t, const JObject*>> DSourceCombo::Get_SourceParticles(bool locEntireChainFlag) const
{
	if(!locEntireChainFlag || dFurtherDecayCombos.empty())
		return dSourceParticles;

	vector<pair<Particle_t, const JObject*>> locToReturnParticles = dSourceParticles;
	for(const auto& locDecayPair : dFurtherDecayCombos)
	{
		const auto& locDecayVector = locDecayPair.second;
		for(auto locDecayCombo : locDecayVector)
		{
			auto locDecayParticles = locDecayCombo->Get_SourceParticles(true);
			locToReturnParticles.insert(locToReturnParticles.end(), locDecayParticles.begin(), locDecayParticles.end());
		}
	}

	return locToReturnParticles;
}

/*********************************************************** INLINE NAMESPACE FUNCTION DEFINITIONS ************************************************************/

inline vector<const JObject*> Get_SourceParticles(const vector<pair<Particle_t, const JObject*>>& locSourceParticles, Particle_t locPID = Unknown)
{
	//if PID is unknown, then all particles
	vector<const JObject*> locOutputParticles;
	for(const auto& locParticlePair : locSourceParticles)
	{
		if((locPID == Unknown) || (locParticlePair.first == locPID))
			locOutputParticles.push_back(locParticlePair.second);
	}
	return locOutputParticles;
}

inline vector<pair<Particle_t, const JObject*>> Get_SourceParticles_ThisVertex(const DSourceCombo* locSourceCombo)
{
	auto locSourceParticles = locSourceCombo->Get_SourceParticles(false);
	auto locFurtherDecayCombos = locSourceCombo->Get_FurtherDecayCombos();
	for(const auto& locDecayPair : locFurtherDecayCombos)
	{
		auto locDecayPID = std::get<0>(locDecayPair.first);
		if(IsDetachedVertex(locDecayPID))
			continue;
		for(auto locDecayCombo : locDecayPair.second)
		{
			auto locDecayParticles = Get_SourceParticles_ThisVertex(locDecayCombo);
			locSourceParticles.insert(locSourceParticles.end(), locDecayParticles.begin(), locDecayParticles.end());
		}
	}

	return locSourceParticles;
}

inline vector<const DSourceCombo*> Get_SourceCombos_ThisVertex(const DSourceCombo* locSourceCombo)
{
	vector<const DSourceCombo*> locVertexCombos = {locSourceCombo};
	for(const auto& locDecayPair : locSourceCombo->Get_FurtherDecayCombos())
	{
		auto locDecayPID = std::get<0>(locDecayPair.first);
		if(IsDetachedVertex(locDecayPID))
			continue;
		for(auto locDecayCombo : locDecayPair.second)
		{
			auto locDecayVertexCombos = Get_SourceCombos_ThisVertex(locDecayCombo);
			locVertexCombos.insert(locVertexCombos.end(), locDecayVertexCombos.begin(), locDecayVertexCombos.end());
		}
	}
	return locVertexCombos;
}

inline vector<pair<DSourceComboUse, vector<const DSourceCombo*>>> Get_SourceCombosAndUses_ThisVertex(const DSourceCombo* locSourceCombo)
{
	//sorted from most- to least-dependent
	vector<pair<DSourceComboUse, vector<const DSourceCombo*>>> locVertexCombosByUse;
	for(const auto& locDecayPair : locSourceCombo->Get_FurtherDecayCombos())
	{
		auto locDecayPID = std::get<0>(locDecayPair.first);
		if(IsDetachedVertex(locDecayPID))
			continue;
		locVertexCombosByUse.emplace_back(locDecayPair);
		for(auto locDecayCombo : locDecayPair.second)
		{
			auto locDecayVertexCombosByUse = Get_SourceCombosAndUses_ThisVertex(locDecayCombo);
			locVertexCombosByUse.insert(locVertexCombosByUse.end(), locDecayVertexCombosByUse.begin(), locDecayVertexCombosByUse.end());
		}
	}
	return locVertexCombosByUse;
}

inline Charge_t Get_ChargeContent(const DSourceComboInfo* locSourceComboInfo)
{
	auto locNumParticles = locSourceComboInfo->Get_NumParticles(true);

	Charge_t locCharge = d_Charged;
	auto Charge_Search = [locCharge](const pair<Particle_t, unsigned char>& locPair) -> bool
		{return Is_CorrectCharge(locPair.first, locCharge);};

	if(!std::any_of(locNumParticles.begin(), locNumParticles.end(), Charge_Search))
		return d_Neutral;

	locCharge = d_Neutral;
	if(!std::any_of(locNumParticles.begin(), locNumParticles.end(), Charge_Search))
		return d_Charged;

	return d_AllCharges;
}

inline bool Get_HasMassiveNeutrals(const DSourceComboInfo* locComboInfo) const
{
	//see if the combo info contains a massive neutral particle

	//search function
	auto Find_MassiveNeutrals = [](const pair<Particle_t, unsigned char>& locPair) -> bool
		{return (ParticleCharge(locPair.first) != 0) && (ParticleMass(locPair.first) > 0.0);};

	//do search
	auto locNumParticles = locComboInfo->Get_NumParticles(true); //true: entire chain
	return std::any_of(locNumParticles.begin(), locNumParticles.end(), Find_MassiveNeutrals);
}

const DSourceCombo* Find_Combo_AtThisStep(const DSourceCombo* locSourceCombo, const DSourceComboUse& locUseToFind, size_t locDecayInstanceIndex)
{
	for(const auto& locDecayPair : locSourceCombo->Get_FurtherDecayCombos())
	{
		const auto& locUse = locDecayPair.first;
		if(locUse == locUseToFind) //good, do stuff
			return locDecayPair.second[locDecayInstanceIndex];
		if(std::get<0>(locUse) != Unknown)
			continue; //is another step!

		//vector of combos is guaranteed to be size 1, and it's guaranteed that none of ITS further decays are unknown
		auto locComboToSearch = locDecayPair.second[0];
		for(const auto& locNestedDecayPair : locComboToSearch->Get_FurtherDecayCombos())
		{
			if(locUse == locUseToFind) //good, do stuff
				return locNestedDecayPair.second[locDecayInstanceIndex];
		}
	}
	return nullptr;
}

} //end DAnalysis namespace

#endif // DSourceCombo_h
