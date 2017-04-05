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

/****************************************************** DEFINE LAMBDAS, USING STATEMENTS *******************************************************/

//DPhotonComboUse
//what the combo is USED for (the decay of Particle_t (if Unknown then is just a grouping)
class DPhotonComboInfo;
class DPhotonCombo;
using DPhotonComboUse = pair<Particle_t, shared_ptr<const DPhotonComboInfo>>;

//DEFINE LAMBDAS
auto Compare_PhotonComboUses = [](const DPhotonComboUse& lhs, const DPhotonComboUse& rhs) -> bool
{
	if(lhs.first == rhs.first)
		return *lhs.second < *rhs.second;
	if(ParticleMass(lhs.first) == ParticleMass(rhs.first))
		return lhs.first < rhs.first; //dunno how this would be possible ...
	return (ParticleMass(lhs.first) < ParticleMass(rhs.first));
};

//DPhotonComboUseMap
using DPhotonComboUseMap = map<DPhotonComboUse, size_t, decltype(Compare_PhotonComboUses)>; //size_t: # of (e.g.) pi0s, etc.
//MAYBE INSTEAD:
//key pair: pi0 -> 2g, 3pi0s
//value: use representing Unknown -> 3pi0 (if circular uses *this) //UGHHHH
//using DPhotonComboUseMap = map<pair<DPhotonComboUse, size_t>, DPhotonComboUse, decltype(Compare_PhotonComboUses)>; //size_t: # of (e.g.) pi0s, etc.

/************************************************************** DEFINE CLASSES ***************************************************************/

class DPhotonComboInfo
{
	public:

		//CONSTRUCTORS AND OPERATORS
		DPhotonComboInfo(void) = delete;
		DPhotonComboInfo(size_t locNumPhotons, const DPhotonComboUseMap& locFurtherDecays);
		bool operator< (const DPhotonComboInfo& rhs) const;

		//GET MEMBERS
		size_t Get_NumPhotons(void) const{return dNumPhotons;}
		DPhotonComboUseMap Get_FurtherDecays(void) const{return dFurtherDecays;}

	private:

		//don't have decaying PID a direct member of this combo info
		//e.g. for a 2g pair, it has no idea whether or not it came from a Pi0, an eta, through direct production, etc.
		//this way, e.g., the 2g combos can be used for ANY of those possibilities, without creating new objects
		//it is the responsibility of the containing object to know what the combos are used for: DPhotonComboUse

		size_t dNumPhotons = 0;
		DPhotonComboUseMap dFurtherDecays;
};

class DPhotonCombo
{
	public:
		DPhotonCombo(void) = delete;

	private:
		shared_ptr<const DPhotonComboInfo> dPhotonComboInfo = nullptr;
		vector<const DNeutralShower*> dPhotonShowers;

		//vector: e.g. size 3 if 3 pi0s needed
		map<DPhotonComboUse, vector<shared_ptr<DPhotonCombo>>> dFurtherDecays;
};

/*********************************************************** INLINE MEMBER FUNCTION DEFINITIONS ************************************************************/

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

} //end DAnalysis namespace

#endif // DPhotonCombo_h
