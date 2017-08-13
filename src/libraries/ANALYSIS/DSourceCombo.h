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
#include "DResettable.h"

#include "PID/DNeutralShower.h"

using namespace std;
using namespace jana;

namespace DAnalysis
{

/****************************************************** DEFINE LAMBDAS, USING STATEMENTS, DECLARE FUNCTIONS *******************************************************/

//forward declarations
class DSourceComboInfo;
class DSourceCombo;

//DSourceComboUse is what the combo is USED for (the decay of Particle_t (if Unknown then is just a grouping)
//signed char: vertex-z bin of the final state (combo contents)
//bool: true/false if has/doesn't-have missing decay product (is always false if decay pid == Unknown)
//last pid: Pid to exclude from calculating the invariant mass of the decay (is always Unknown if decay pid == Unknown) (e.g. if rescattering Lambda, p -> p, p, pi- then this would be Proton)
using DSourceComboUse = tuple<Particle_t, signed char, const DSourceComboInfo*, bool, Particle_t>; //e.g. Pi0, zbin, -> 2g, false, Unknown
using DSourceCombosByUse_Small = vector<pair<DSourceComboUse, vector<const DSourceCombo*>>>;

//DECLARE NAMESPACE-SCOPE FUNCTIONS
vector<const JObject*> Get_SourceParticles(const vector<pair<Particle_t, const JObject*>>& locSourceParticles, Particle_t locPID = Unknown);
vector<pair<Particle_t, const JObject*>> Get_SourceParticles_ThisVertex(const DSourceCombo* locSourceCombo, Charge_t locCharge = d_AllCharges);
vector<const DSourceCombo*> Get_SourceCombos_ThisVertex(const DSourceCombo* locSourceCombo);
vector<pair<DSourceComboUse, vector<const DSourceCombo*>>> Get_SourceCombosAndUses_ThisVertex(const DSourceCombo* locSourceCombo);
Charge_t Get_ChargeContent(const DSourceComboInfo* locSourceComboInfo);
bool Get_HasMassiveNeutrals(const DSourceComboInfo* locComboInfo);
const JObject* Get_SourceParticle_ThisStep(const DSourceCombo* locSourceCombo, Particle_t locPID, size_t locInstance, size_t& locPIDCountSoFar);

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

		//FORWARD DECLARE COMPARISON STRUCTS
		struct DCompare_ParticlePairPIDs;
		struct DCompare_FurtherDecays;

		//CONSTRUCTOR
		DSourceComboInfo(void) = delete;
		DSourceComboInfo(const vector<pair<Particle_t, unsigned char>>& locNumParticles, const vector<pair<DSourceComboUse, unsigned char>>& locFurtherDecays = {});

		//OPERATORS
		bool operator< (const DSourceComboInfo& rhs) const;

		//GET MEMBERS
		vector<pair<Particle_t, unsigned char>> Get_NumParticles(bool locEntireChainFlag = false) const;
		vector<pair<DSourceComboUse, unsigned char>> Get_FurtherDecays(void) const{return dFurtherDecays;}

		//definitions of negative values for any particle index //in order such that operator< returns order expected for string (e.g. gp->...)
		static signed char Get_VertexZIndex_OutOfRange(void){return -3;}
		static signed char Get_VertexZIndex_ZIndependent(void){return -2;}
		static signed char Get_VertexZIndex_Unknown(void){return -1;}

	private:

		//don't have decaying PID a direct member of this combo info
		//e.g. for a 2g pair, it has no idea whether or not it came from a Pi0, an eta, through direct production, etc.
		//this way, e.g., the 2g combos can be used for ANY of those possibilities, without creating new objects
		//it is the responsibility of the containing object to know what the combos are used for: DSourceComboUse

		vector<pair<Particle_t, unsigned char>> dNumParticles;
		vector<pair<DSourceComboUse, unsigned char>> dFurtherDecays; //unsigned char: # of (e.g.) pi0s, etc.
};

inline bool operator<(const DSourceComboUse& lhs, const DSourceComboUse& rhs)
{
	//this puts mixed-charge first, then fully-neutral, then fully-charged

	//first, special case of nullptr (guard against it)
	if((std::get<2>(lhs) == nullptr) || (std::get<2>(rhs) == nullptr))
	{
		if(std::get<2>(lhs) != std::get<2>(rhs))
			return (std::get<2>(lhs) == nullptr);
		if(std::get<3>(lhs) != std::get<3>(rhs))
			return (std::get<3>(lhs) < std::get<3>(rhs));
		if(std::get<4>(lhs) != std::get<4>(rhs))
			return (std::get<4>(lhs) < std::get<4>(rhs));

		if(std::get<0>(lhs) == std::get<0>(rhs))
		{
			if(std::get<1>(lhs) == std::get<1>(rhs))
				return false;
			else
				return std::get<1>(lhs) > std::get<1>(rhs);
		}
		if(ParticleMass(std::get<0>(lhs)) == ParticleMass(std::get<0>(rhs)))
			return std::get<0>(lhs) > std::get<0>(rhs);
		return (ParticleMass(std::get<0>(lhs)) > ParticleMass(std::get<0>(rhs)));
	}

	auto locChargeContent_LHS = Get_ChargeContent(std::get<2>(lhs));
	auto locChargeContent_RHS = Get_ChargeContent(std::get<2>(rhs));
	if(locChargeContent_LHS != locChargeContent_RHS)
		return locChargeContent_LHS > locChargeContent_RHS;

	if(std::get<3>(lhs) != std::get<3>(rhs))
		return (std::get<3>(lhs) < std::get<3>(rhs));
	if(std::get<4>(lhs) != std::get<4>(rhs))
		return (std::get<4>(lhs) < std::get<4>(rhs));

	//within each of those, it puts the most-massive particles first
	if(std::get<0>(lhs) == std::get<0>(rhs))
	{
		if(std::get<1>(lhs) == std::get<1>(rhs))
			return *std::get<2>(rhs) < *std::get<2>(lhs);
		else
			return std::get<1>(lhs) > std::get<1>(rhs);
	}
	if(ParticleMass(std::get<0>(lhs)) == ParticleMass(std::get<0>(rhs)))
		return std::get<0>(lhs) > std::get<0>(rhs);
	return (ParticleMass(std::get<0>(lhs)) > ParticleMass(std::get<0>(rhs)));
}

struct DSourceComboInfo::DCompare_ParticlePairPIDs
{
	bool operator()(const pair<Particle_t, unsigned char>& lhs, const pair<Particle_t, unsigned char>& rhs) const{return lhs.first < rhs.first;} //sort
	bool operator()(const pair<Particle_t, unsigned char>& lhs, Particle_t rhs) const{return lhs.first < rhs;} //lookup
	bool operator()(Particle_t lhs, const pair<Particle_t, unsigned char>& rhs) const{return lhs < rhs.first;} //lookup
};

struct DSourceComboInfo::DCompare_FurtherDecays
{
	bool operator()(const pair<DSourceComboUse, unsigned char>& lhs, const pair<DSourceComboUse, unsigned char>& rhs) const{return DAnalysis::operator<(lhs.first, rhs.first);} //sort
	bool operator()(const pair<DSourceComboUse, unsigned char>& lhs, DSourceComboUse rhs) const{return DAnalysis::operator<(lhs.first, rhs);} //lookup
	bool operator()(DSourceComboUse lhs, const pair<DSourceComboUse, unsigned char>& rhs) const{return DAnalysis::operator<(lhs, rhs.first);} //lookup
};

class DSourceCombo : public DResettable
{
	public:

		//FORWARD DECLARE COMPARISON STRUCT
		struct DCompare_FurtherDecays;

		//CONSTRUCTOR
		DSourceCombo(void) = default;
		DSourceCombo(const vector<pair<Particle_t, const JObject*>>& locSourceParticles, const DSourceCombosByUse_Small& locFurtherDecayCombos, bool locIsZIndependent = false);

		//SET MEMBERS
		void Set_Members(const vector<pair<Particle_t, const JObject*>>& locSourceParticles, const DSourceCombosByUse_Small& locFurtherDecayCombos, bool locIsZIndependent = false);
		void Reset(void);
		void Release(void){Reset();};

		//GET MEMBERS
		vector<pair<Particle_t, const JObject*>> Get_SourceParticles(bool locEntireChainFlag = false, Charge_t locCharge = d_AllCharges) const;
		DSourceCombosByUse_Small Get_FurtherDecayCombos(void) const{return dFurtherDecayCombos;}
		bool Get_IsComboingZIndependent(void) const{return dIsComboingZIndependent;}

	private:

		//FYI, if use PID is unknown: Vector is guaranteed to be size 1, and none of ITS further decays can have an unknown use PID
		//Note: Either:
		//1) The entire content is one PID
		//2) There is at most one particle of each PID
		//Otherwise, groupings of 2+ (e.g. 2pi-) will be in the further-decay section (X -> 2pi+)

		//particles & decays
		vector<pair<Particle_t, const JObject*>> dSourceParticles; //original DNeutralShower or DChargedTrack
		DSourceCombosByUse_Small dFurtherDecayCombos; //vector is e.g. size 3 if 3 pi0s needed

		//everything is z-dependent for massive neutrals.
		//however the momentum is SO z-dependent, that we can't cut on it until the end when we have the final vertex, AFTER comboing
		//a mere z-bin is not enough.
		//So, as far as COMBOING is concerned, massive neutrals are Z-INDEPENDENT
		bool dIsComboingZIndependent = false; //is false for BCAL photons
};

struct DSourceCombo::DCompare_FurtherDecays
{
	bool operator()(const pair<DSourceComboUse, vector<const DSourceCombo*>>& lhs, const pair<DSourceComboUse, vector<const DSourceCombo*>>& rhs) const{return DAnalysis::operator<(lhs.first, rhs.first);} //sort
	bool operator()(const pair<DSourceComboUse, vector<const DSourceCombo*>>& lhs, DSourceComboUse rhs) const{return DAnalysis::operator<(lhs.first, rhs);} //lookup
	bool operator()(DSourceComboUse lhs, const pair<DSourceComboUse, vector<const DSourceCombo*>>& rhs) const{return DAnalysis::operator<(lhs, rhs.first);} //lookup
};

struct DSourceComboChecker_ReusedParticle
{
	//returns true if a particle WAS reused
	bool operator()(const DSourceCombo* locCombo) const
	{
		auto locParticles = DAnalysis::Get_SourceParticles(locCombo->Get_SourceParticles(true), Unknown); //all pids
		std::sort(locParticles.begin(), locParticles.end());
		auto locUniqueIterator = std::unique(locParticles.begin(), locParticles.end());
		return (locUniqueIterator != locParticles.end());
	}
};

/*********************************************************** INLINE MEMBER FUNCTION DEFINITIONS ************************************************************/

inline DSourceComboInfo::DSourceComboInfo(const vector<pair<Particle_t, unsigned char>>& locNumParticles, const vector<pair<DSourceComboUse, unsigned char>>& locFurtherDecays) :
		dNumParticles(locNumParticles), dFurtherDecays(locFurtherDecays)
{
	std::sort(dNumParticles.begin(), dNumParticles.end(), DCompare_ParticlePairPIDs());
	std::sort(dFurtherDecays.begin(), dFurtherDecays.end(), DCompare_FurtherDecays());
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
				(*locIteratorPair.first).second += locParticlePair.second; //it exists: increase it
			else //doesn't exist. insert it
				locToReturnNumParticles.insert(locIteratorPair.first, locParticlePair);
		}
	}

	return locToReturnNumParticles;
}

inline DSourceCombo::DSourceCombo(const vector<pair<Particle_t, const JObject*>>& locSourceParticles, const DSourceCombosByUse_Small& locFurtherDecayCombos, bool locIsZIndependent) :
		dSourceParticles(locSourceParticles), dFurtherDecayCombos(locFurtherDecayCombos), dIsComboingZIndependent(locIsZIndependent) {}

inline void DSourceCombo::Reset(void)
{
	dSourceParticles.clear();
	dFurtherDecayCombos.clear();
	dIsComboingZIndependent = false;
}

inline void DSourceCombo::Set_Members(const vector<pair<Particle_t, const JObject*>>& locSourceParticles, const DSourceCombosByUse_Small& locFurtherDecayCombos, bool locIsZIndependent)
{
	dIsComboingZIndependent = locIsZIndependent;
	dSourceParticles = locSourceParticles;
	dFurtherDecayCombos = locFurtherDecayCombos;
	std::sort(dFurtherDecayCombos.begin(), dFurtherDecayCombos.end(), DCompare_FurtherDecays());
}

inline vector<pair<Particle_t, const JObject*>> DSourceCombo::Get_SourceParticles(bool locEntireChainFlag, Charge_t locCharge) const
{
	if(!locEntireChainFlag || dFurtherDecayCombos.empty())
		return dSourceParticles;

	vector<pair<Particle_t, const JObject*>> locToReturnParticles = dSourceParticles;

	auto Charge_Checker = [&locCharge](const pair<Particle_t, const JObject*>& locPair) -> bool {return !Is_CorrectCharge(locPair.first, locCharge);};
	locToReturnParticles.erase(std::remove_if(locToReturnParticles.begin(), locToReturnParticles.end(), Charge_Checker), locToReturnParticles.end());

	for(const auto& locDecayPair : dFurtherDecayCombos)
	{
		const auto& locDecayVector = locDecayPair.second;
		for(const auto& locDecayCombo : locDecayVector)
		{
			auto locDecayParticles = locDecayCombo->Get_SourceParticles(true, locCharge);
			locToReturnParticles.insert(locToReturnParticles.end(), locDecayParticles.begin(), locDecayParticles.end());
		}
	}

	return locToReturnParticles;
}

/*********************************************************** INLINE NAMESPACE FUNCTION DEFINITIONS ************************************************************/

void Print_SourceComboUse(const DSourceComboUse& locComboUse, unsigned char locNumTabs = 0, bool locIgnoreTabs = false);
inline void Print_SourceComboInfo(const DSourceComboInfo* locComboInfo, unsigned char locNumTabs = 0)
{
	if(locComboInfo == nullptr)
		return;

	auto locNumParticles = locComboInfo->Get_NumParticles(false);
	for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
	cout << "Particles: ";
	for(auto& locParticlePair : locNumParticles)
		cout << int(locParticlePair.second) << " " << ParticleType(locParticlePair.first) << ", ";
	cout << endl;

	auto locFurtherDecays = locComboInfo->Get_FurtherDecays();
	for(auto& locDecayPair : locFurtherDecays)
	{
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << int(locDecayPair.second) << " of ";
		DAnalysis::Print_SourceComboUse(locDecayPair.first, locNumTabs, true);
	}
}

inline void Print_SourceComboUse(const DSourceComboUse& locComboUse, unsigned char locNumTabs, bool locIgnoreTabs)
{
	if(!locIgnoreTabs)
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
	cout << "Use (decay pid, z-bin, combo info, has-missing-decay-product, mass-cut-pid-to-exclude): ";
	cout << ParticleType(std::get<0>(locComboUse)) << ", " << int(std::get<1>(locComboUse)) << ", " << std::get<2>(locComboUse) << ", " << std::get<3>(locComboUse) << ", " << std::get<4>(locComboUse) << ":" << endl;
	DAnalysis::Print_SourceComboInfo(std::get<2>(locComboUse), locNumTabs + 1);
}

inline void Print_SourceCombo(const DSourceCombo* locCombo, unsigned char locNumTabs = 0)
{
	if(locCombo == nullptr)
		return;

	for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
	cout << "pointer: " << locCombo << ", Z-independent?: " << locCombo->Get_IsComboingZIndependent() << ", Particles: ";
	auto locSourceParticles = locCombo->Get_SourceParticles();
	for(auto& locParticlePair : locSourceParticles)
		cout << ParticleType(locParticlePair.first) << " " << locParticlePair.second << ", ";
	cout << endl;

	DSourceCombosByUse_Small locFurtherDecayCombos = locCombo->Get_FurtherDecayCombos();
	for(auto& locDecayPair : locFurtherDecayCombos)
	{
		for(decltype(locNumTabs) locTabNum = 0; locTabNum < locNumTabs; ++locTabNum) cout << "\t";
		cout << locDecayPair.second.size() << " of ";
		DAnalysis::Print_SourceComboUse(locDecayPair.first, locNumTabs, true);
		for(auto& locCombo : locDecayPair.second)
			DAnalysis::Print_SourceCombo(locCombo, locNumTabs + 1);
	}
}

inline vector<const JObject*> Get_SourceParticles(const vector<pair<Particle_t, const JObject*>>& locSourceParticles, Particle_t locPID)
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

inline vector<const JObject*> Get_SourceParticles(const vector<pair<Particle_t, const JObject*>>& locSourceParticles, int locCharge)
{
	//ignores the charge magnitude and sign: only considers if ==/!= 0
	vector<const JObject*> locOutputParticles;
	for(const auto& locParticlePair : locSourceParticles)
	{
		auto locParticleCharge = ParticleCharge(locParticlePair.first);
		if((locParticleCharge == locCharge) && (locParticleCharge == 0))
			locOutputParticles.push_back(locParticlePair.second);
		else if(locParticleCharge*locCharge != 0)
			locOutputParticles.push_back(locParticlePair.second);
	}
	return locOutputParticles;
}

inline vector<pair<Particle_t, const JObject*>> Get_SourceParticles_ThisVertex(const DSourceCombo* locSourceCombo, Charge_t locCharge)
{
	auto locSourceParticles = locSourceCombo->Get_SourceParticles(false, locCharge);
	auto locFurtherDecayCombos = locSourceCombo->Get_FurtherDecayCombos();
	for(const auto& locDecayPair : locFurtherDecayCombos)
	{
		auto locDecayPID = std::get<0>(locDecayPair.first);
		if(IsDetachedVertex(locDecayPID))
			continue;
		for(const auto& locDecayCombo : locDecayPair.second)
		{
			auto locDecayParticles = Get_SourceParticles_ThisVertex(locDecayCombo, locCharge);
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
		for(const auto& locDecayCombo : locDecayPair.second)
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
		for(const auto& locDecayCombo : locDecayPair.second)
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
	auto Charge_Search = [&locCharge](const pair<Particle_t, unsigned char>& locPair) -> bool
		{return Is_CorrectCharge(locPair.first, locCharge);};

	if(!std::any_of(locNumParticles.begin(), locNumParticles.end(), Charge_Search))
		return d_Neutral;

	locCharge = d_Neutral;
	if(!std::any_of(locNumParticles.begin(), locNumParticles.end(), Charge_Search))
		return d_Charged;

	return d_AllCharges;
}

inline bool Get_HasMassiveNeutrals(const DSourceComboInfo* locComboInfo)
{
	//see if the combo info contains a massive neutral particle

	//search function
	auto Find_MassiveNeutrals = [](const pair<Particle_t, unsigned char>& locPair) -> bool
		{return ((ParticleCharge(locPair.first) == 0) && (ParticleMass(locPair.first) > 0.0));};

	//do search
	auto locNumParticles = locComboInfo->Get_NumParticles(true); //true: entire chain
	return std::any_of(locNumParticles.begin(), locNumParticles.end(), Find_MassiveNeutrals);
}

inline const JObject* Get_SourceParticle_ThisStep(const DSourceCombo* locSourceCombo, Particle_t locPID, size_t locInstance, size_t& locPIDCountSoFar)
{
	auto ParticleFinder = [&locPID, &locInstance, &locPIDCountSoFar](pair<Particle_t, const JObject*>& locPair) -> bool
		{return ((locPair.first != locPID) ? false : ((++locPIDCountSoFar == locInstance) ? true : false));};

	auto locParticles = locSourceCombo->Get_SourceParticles();
	auto locIterator = std::find_if(locParticles.begin(), locParticles.end(), ParticleFinder);
	if(locIterator != locParticles.end())
		return locIterator->second;

	for(const auto& locDecayPair : locSourceCombo->Get_FurtherDecayCombos())
	{
		if(std::get<0>(locDecayPair.first) != Unknown)
			continue; //a new step!
		for(const auto& locDecayCombo : locDecayPair.second)
		{
			auto locParticle = Get_SourceParticle_ThisStep(locDecayCombo, locPID, locInstance, locPIDCountSoFar);
			if(locParticle != nullptr)
				return locParticle;
		}
	}

	return nullptr;
}

inline bool Check_AreDuplicateCombos(const DSourceCombo* lhs, const DSourceCombo* rhs)
{
	//Assumes inputs are from the same source (excluding z)
	auto locParticles_lhs = lhs->Get_SourceParticles(false);
	auto locParticles_rhs = rhs->Get_SourceParticles(false);
	if(locParticles_lhs.size() != locParticles_rhs.size())
		return false;
	if(!std::is_permutation(locParticles_lhs.begin(), locParticles_lhs.end(), locParticles_rhs.begin()))
		return false;

	auto locDecayCombos_lhs = lhs->Get_FurtherDecayCombos();
	auto locDecayCombos_rhs = rhs->Get_FurtherDecayCombos();
	if(locDecayCombos_lhs.size() != locDecayCombos_rhs.size())
		return false;

	for(auto& locDecayPair : locDecayCombos_lhs)
	{
		auto locIteratorPair = std::equal_range(locDecayCombos_rhs.begin(), locDecayCombos_rhs.end(), locDecayPair.first, DSourceCombo::DCompare_FurtherDecays());
		if(locIteratorPair.first == locIteratorPair.second)
			return false; //careful, compares z's!!

		auto& locDecayCombos_lhs = locDecayPair.second;
		auto& locDecayCombos_rhs = (*locIteratorPair.first).second;
		if(locDecayCombos_lhs.size() != locDecayCombos_rhs.size())
			return false;
		if(!std::is_permutation(locDecayCombos_lhs.begin(), locDecayCombos_lhs.end(), locDecayCombos_rhs.begin(), Check_AreDuplicateCombos))
			return false;
	}
	return true;
}

} //end DAnalysis namespace

#endif // DSourceCombo_h
