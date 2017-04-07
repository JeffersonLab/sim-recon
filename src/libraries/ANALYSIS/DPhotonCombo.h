#ifndef DPhotonCombo_h
#define DPhotonCombo_h

#include <map>
#include <vector>
#include <memory>
#include <algorithm>

#include "JANA/JObject.h"
#include "JANA/JEventLoop.h"

#include "particleType.h"

#include "PID/DNeutralShower.h"

using namespace std;

namespace DAnalysis
{
/*
MEMORY, TIME CONSIDERATIONS:

//DPhotonCombo
//consider, instead of storing DNeutralShower pointer, store index to main array (unsigned char!!)

*/
/****************************************************** DEFINE LAMBDAS, USING STATEMENTS *******************************************************/

//forward declarations
class DPhotonComboInfo;
class DPhotonCombo;
class DPhotonComboer;

//DPhotonComboUse is what the combo is USED for (the decay of Particle_t (if Unknown then is just a grouping)
using DPhotonComboUse = pair<Particle_t, const DPhotonComboInfo*>; //e.g. Pi0, -> 2g
using DPhotonCombosByUse_Small = map<DPhotonComboUse, vector<const DPhotonCombo*>>; //for use WITHIN a DPhotonCombo (e.g. vector size 3 of pi0 decays)

//Compare_PhotonComboUses
auto Compare_PhotonComboUses = [](const DPhotonComboUse& lhs, const DPhotonComboUse& rhs) -> bool
{
	//this puts the least-massive particles first
	if(lhs.first == rhs.first)
		return *lhs.second < *rhs.second;
	if(ParticleMass(lhs.first) == ParticleMass(rhs.first))
		return lhs.first < rhs.first;
	return (ParticleMass(lhs.first) < ParticleMass(rhs.first));
};

//consider removing DPhotonComboInfo::dNumPhotons and instead storing in below with key = Gamma, nullptr

/************************************************************** DEFINE CLASSES ***************************************************************/

//In theory, for safety, dynamically-allocated objects should be stored in a shared_ptr
//However, the combos take a TON of memory, and a shared_ptr<T*> takes 3x the memory of a regular T* pointer
//Plus, only the DPhotonComboer class can create these objects, and it always registers them with itself.

//The info objects will exist for the life of the program, and so don't need to be recycled to a resource pool.
//The combo objects will be recycled every event into a resource pool.

//THE MOST NUMBER OF PARTICLES OF A GIVEN TYPE IS 255 (# stored in unsigned char)
class DPhotonComboInfo
{
	public:

		//FRIEND CLASS
		friend class DPhotonComboer; //so that can call the constructor

		//CONSTRUCTORS AND OPERATORS
		DPhotonComboInfo(void) = delete;
		bool operator< (const DPhotonComboInfo& rhs) const;

		//GET MEMBERS
		unsigned char Get_NumPhotons(void) const{return dNumPhotons;}
		vector<pair<DPhotonComboUse, unsigned char>> Get_FurtherDecays(void) const{return dFurtherDecays;}

	private:

		DPhotonComboInfo(unsigned char locNumPhotons, const vector<pair<DPhotonComboUse, unsigned char>>& locFurtherDecays = {});

		//don't have decaying PID a direct member of this combo info
		//e.g. for a 2g pair, it has no idea whether or not it came from a Pi0, an eta, through direct production, etc.
		//this way, e.g., the 2g combos can be used for ANY of those possibilities, without creating new objects
		//it is the responsibility of the containing object to know what the combos are used for: DPhotonComboUse

		unsigned char dNumPhotons = 0;
		//this will be sorted with Compare_PhotonComboUses
		vector<pair<DPhotonComboUse, unsigned char>> dFurtherDecays; //unsigned char: # of (e.g.) pi0s, etc.
};

class DPhotonCombo
{
	public:

		//FRIEND CLASS
		friend class DPhotonComboer; //so that can call the constructor

		DPhotonCombo(void) = delete;

		//GET MEMBERS
		vector<const DNeutralShower*> Get_PhotonShowers(bool locEntireChainFlag = false) const;
		DPhotonCombosByUse_Small Get_FurtherDecayCombos(void) const{return dFurtherDecayCombos;}

	private:
		DPhotonCombo(const vector<const DNeutralShower*>& locPhotonShowers, const DPhotonCombosByUse_Small& locFurtherDecayCombos = {});

//vector<pair<Particle_t, JObject*>>
		vector<const DNeutralShower*> dPhotonShowers;

		//DPhotonCombosByUse_Small vector: e.g. size 3 if 3 pi0s needed
		DPhotonCombosByUse_Small dFurtherDecayCombos;
};

/*********************************************************** INLINE MEMBER FUNCTION DEFINITIONS ************************************************************/

inline DPhotonComboInfo::DPhotonComboInfo(unsigned char locNumPhotons, const vector<pair<DPhotonComboUse, unsigned char>>& locFurtherDecays) :
		dNumPhotons(locNumPhotons), dFurtherDecays(locFurtherDecays)
{
	std::sort(dFurtherDecays.begin(), dFurtherDecays.end(), Compare_PhotonComboUses);
}

inline DPhotonCombo::DPhotonCombo(const vector<const DNeutralShower*>& locPhotonShowers, const DPhotonCombosByUse_Small& locFurtherDecayCombos) :
		dPhotonShowers(locPhotonShowers), dFurtherDecayCombos(locFurtherDecayCombos) {}

inline bool DPhotonComboInfo::operator< (const DPhotonComboInfo& rhs) const
{
	if(dNumPhotons != rhs.dNumPhotons)
		return dNumPhotons < rhs.dNumPhotons;

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

inline vector<const DNeutralShower*> DPhotonCombo::Get_PhotonShowers(bool locEntireChainFlag) const
{
	if(!locEntireChainFlag || dFurtherDecayCombos.empty())
		return dPhotonShowers;

	vector<const DNeutralShower*> locToReturnShowers(dPhotonShowers);
	for(const auto& locDecayPair : dFurtherDecayCombos)
	{
		const auto& locDecayVector = *(locDecayPair.second);
		for(auto locDecayCombo : locDecayVector)
		{
			auto locDecayPhotons = locDecayCombo->Get_PhotonShowers(true);
			locToReturnShowers.insert(locToReturnShowers.end(), locDecayPhotons.begin(), locDecayPhotons.end());
		}
	}

	return locToReturnShowers;
}

} //end DAnalysis namespace

#endif // DPhotonCombo_h
