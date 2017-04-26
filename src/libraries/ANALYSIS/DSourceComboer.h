#ifndef DSourceComboer_h
#define DSourceComboer_h

#include <map>
#include <set>
#include <vector>
#include <memory>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <stack>

#include "TF1.h"

#include "JANA/JObject.h"
#include "JANA/JEventLoop.h"

#include "particleType.h"
#include "DANA/DApplication.h"
#include "HDGEOMETRY/DGeometry.h"

#include "PID/DNeutralShower.h"
#include "PID/DKinematicData.h"
#include "PID/DEventRFBunch.h"

#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DSourceCombo.h"
#include "ANALYSIS/DReactionStepVertexInfo.h"
#include "ANALYSIS/DReactionVertexInfo.h"
#include "ANALYSIS/DResourcePool.h"

#include "ANALYSIS/DSourceComboVertexer.h"
#include "ANALYSIS/DSourceComboP4Handler.h"
#include "ANALYSIS/DSourceComboTimeHandler.h"

using namespace std;
using namespace jana;

namespace DAnalysis
{

//BIG TO DO'S:
//fill tracks by PID
//once RF bunch is chosen, redo mass cuts involving massive neutrals
//once vertex position fully defined, place mass cuts on massive neutrals
//compute vertices using beam energy (missing mass)
//finish porting from DVertexCreator

//ANY TIME:
//Cut combo ahead of time if not enough tracks/showers
//fill in calc inv mass functions

//AT THE END:
//VERY CAREFULLY recycle resources
//CONSIDER VECTOR INSTEAD OF MAP FOR DSourceCombosByUse_Small
//store unsigned char instead of pointer???

//MUST BEWARE DUPLICATE COMBOS
//let's say a combo of charged tracks has 2 valid RF bunches
//and we need to combo 2 pi0s with them
//and the shower timing cuts are loose enough that all 4 showers satisfy both RF bunches
//if we combo the 2 rf bunches separately: WE HAVE DUPLICATE COMBOS
//and doing the duplicate check AFTER the fact takes FOREVER
//therefore, we must take the neutral showers for the 2 rfs, COMBINE THEM, and then COMBO AS A UNIT

/****************************************************** DEFINE LAMBDAS, USING STATEMENTS *******************************************************/

struct DCompare_SourceComboInfos{
	bool operator()(const DSourceComboInfo* lhs, const DSourceComboInfo* rhs) const{return *lhs < *rhs;}
};

/********************************************************** DEFINE USING STATEMENTS ***********************************************************/

//DEFINE USING STATEMENTS
using DCombosByBeamBunch = unordered_map<vector<int>, vector<const DSourceCombo*>>;
using DSourceCombosByBeamBunchByUse = unordered_map<DSourceComboUse, DCombosByBeamBunch>;
using DComboIteratorsByBeamBunch = unordered_map<vector<int>, vector<const DSourceCombo*>::const_iterator>; //vector<int>: RF bunches (empty for all)
using DSourceCombosByUse_Large = map<DSourceComboUse, vector<const DSourceCombo*>*>;

/************************************************************** DEFINE CLASSES ***************************************************************/

enum ComboingStage_t
{
	d_ChargedStage,
	d_MixedStage_ZIndependent,
	d_MixedStage
};

class DSourceComboer : public JObject
{
	public:

		DSourceComboer(void) = delete;
		DSourceComboer(JEventLoop* locEventLoop);
		DSourceComboer::~DSourceComboer(void);

		Charge_t Get_ChargeContent(const DSourceComboInfo* locSourceComboInfo) const{return dComboInfoChargeContent.find(locSourceComboInfo)->second;}
		const DSourceCombo* Get_VertexPrimaryCombo(const DSourceCombo* locReactionChargedCombo, const DReactionStepVertexInfo* locStepVertexInfo);
		DSourceComboUse Get_SourceComboUse(const DReactionStepVertexInfo* locStepVertexInfo) const{return dSourceComboUseReactionMap.find(locStepVertexInfo)->second;};

	private:

		/********************************************************** DECLARE MEMBER FUNCTIONS ***********************************************************/

		//SETUP
		void Reset_NewEvent(JEventLoop* locEventLoop);
		void Setup_NeutralShowers(JEventLoop* locEventLoop);

		//CREATE PHOTON COMBO INFOS & USES
		void Create_SourceComboInfos(const DReactionVertexInfo* locReactionVertexInfo);
		DSourceComboUse Create_ZDependentSourceComboUses(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionChargedCombo);
		DSourceComboUse Build_NewZDependentUse(const DReaction* locReaction, size_t locStepIndex, signed char locVertexZBin, const DSourceComboUse& locOrigUse, const unordered_map<size_t, DSourceComboUse>& locCreatedUseMap);

		//CREATE PHOTON COMBO INFOS & USES: UTILITY METHODS
		map<Particle_t, unsigned char> Build_ParticleMap(const DReaction* locReaction, size_t locStepIndex, Charge_t locCharge) const;
		pair<bool, map<DSourceComboUse, unsigned char>> Get_FinalStateDecayingComboUses(const DReaction* locReaction, size_t locStepIndex, const map<size_t, DSourceComboUse>& locStepComboUseMap) const;
		DSourceComboUse Make_ComboUse(Particle_t locInitPID, const map<Particle_t, unsigned char>& locNumParticles, const map<DSourceComboUse, unsigned char>& locFurtherDecays);
		const DSourceComboInfo* MakeOrGet_SourceComboInfo(const vector<pair<Particle_t, unsigned char>>& locNumParticles, const vector<pair<DSourceComboUse, unsigned char>>& locFurtherDecays);
		const DSourceComboInfo* GetOrMake_SourceComboInfo(const vector<pair<Particle_t, unsigned char>>& locNumParticles, const vector<pair<DSourceComboUse, unsigned char>>& locFurtherDecays);

		//CREATE COMBOS - GENERAL METHODS
		void Create_SourceCombos(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding);
		void Create_SourceCombos_Unknown(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding);
		//bool Do_CommonComboingTasks(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locAllBut1ComboUse, const DSourceComboUse& locSourceComboUseToAdd, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding, bool locComboingVertically);

		//COMBO VERTICALLY METHODS
		//Note that vertical comboing always takes place at the same vertex-z
		void Combo_Vertically_AllDecays(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding);
		void Combo_Vertically_NDecays(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locNMinus1ComboUse, const DSourceComboUse& locSourceComboDecayUse, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding);
		void Combo_Vertically_AllParticles(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage);
		void Combo_Vertically_NParticles(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locNMinus1ComboUse, ComboingStage_t locComboingStage);

		//COMBO HORIZONTALLY METHODS
		void Combo_Horizontally_All(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding);
		void Create_Combo_OneParticle(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage);
		void Combo_Horizontally_AddCombo(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locAllBut1ComboUse, const DSourceComboUse& locSourceComboUseToAdd, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding);
		void Combo_Horizontally_AddParticle(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locAllBut1ComboUse, Particle_t locPID, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding);

		//BUILD/RETRIEVE RESUME-AT ITERATORS
		void Build_ParticleIterators(const vector<int>& locBeamBunches, const vector<const JObject*>& locParticles);
		void Build_ComboIterators(const vector<int>& locBeamBunches, const vector<const DSourceCombo*>& locCombos, ComboingStage_t locComboingStage, signed char locVertexZBin);
		vector<const JObject*>::const_iterator Get_ResumeAtIterator_Particles(const DSourceCombo* locSourceCombo, const vector<int>& locBeamBunches) const;
		vector<const DSourceCombo*>::const_iterator Get_ResumeAtIterator_Combos(const DSourceCombo* locSourceCombo, const vector<int>& locBeamBunches, ComboingStage_t locComboingStage, signed char locVertexZBin) const;

		//GET POTENTIAL PARTICLES & COMBOS FOR COMBOING
		const vector<const DSourceCombo*>& Get_CombosForComboing(const DSourceComboUse& locComboUse, ComboingStage_t locComboingStage, const vector<int>& locBeamBunches, const DSourceCombo* locChargedCombo_WithPrevious);
		const vector<const DSourceCombo*>& Get_CombosByBeamBunch(DCombosByBeamBunch& locCombosByBunch, const vector<int>& locBeamBunches, ComboingStage_t locComboingStage, signed char locVertexZBin);

		//REGISTER VALID RF BUNCHES
		void Register_ValidRFBunches(const DSourceComboUse& locSourceComboUse, const DSourceCombo* locSourceCombo, const vector<int>& locRFBunches, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_WithNow);

		//PARTICLE UTILITY FUNCTIONS
		const vector<const JObject*>& Get_ParticlesForComboing(Particle_t locPID, ComboingStage_t locComboingStage, const vector<int>& locBeamBunches = {}, signed char locVertexZBin = 0);
		const vector<const JObject*>& Get_ShowersByBeamBunch(const vector<int>& locBeamBunches, DPhotonShowersByBeamBunch& locShowersByBunch);
		shared_ptr<const DKinematicData> Create_KinematicData(const DNeutralShower* locNeutralShower, const DVector3& locVertex) const;
		bool Get_IsZIndependent(const JObject* locObject) const;

		//COMBO UTILITY FUNCTIONS
		DSourceCombosByUse_Large& Get_CombosSoFar(ComboingStage_t locComboingStage, Charge_t locChargeContent_SearchForUse, const DSourceCombo* locChargedCombo = nullptr);
		DSourceCombosByBeamBunchByUse& Get_SourceCombosByBeamBunchByUse(Charge_t locChargeContent_SearchForUse, const DSourceCombo* locChargedCombo = nullptr);
		void Copy_ZIndependentMixedResults(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_WithNow);
		const DSourceCombo* Get_StepSourceCombo(const DReaction* locReaction, size_t locDesiredStepIndex, const DSourceCombo* locSourceCombo_Current, size_t locCurrentStepIndex = 0);
		const DSourceCombo* Get_ChargedCombo_WithNow(const DSourceCombo* locChargedCombo_Presiding) const;
		const DSourceCombo* Get_Presiding_ChargedCombo(const DSourceCombo* locChargedCombo_Presiding, const DSourceComboUse& locNextComboUse, ComboingStage_t locComboingStage, size_t locInstance) const;

		//VERTEX-Z BINNING UTILITY FUNCTIONS
		size_t Get_PhotonVertexZBin(double locVertexZ) const;
		double Get_PhotonVertexZBinCenter(signed char locVertexZBin) const;

		/************************************************************** DEFINE MEMBERS ***************************************************************/

		uint64_t dEventNumber = 0; //re-setup on new events
		string dShowerSelectionTag = "PreSelect";

		//EXPERIMENT INFORMATION
		DVector3 dTargetCenter;

		//VERTEX-DEPENDENT PHOTON INFORMATION
		//For every 10cm in vertex-z, calculate the photon p4 & time for placing mass & delta-t cuts
		//The z-range extends from the upstream end of the target - 5cm to the downstream end + 15cm
		//so for a 30-cm-long target, it's a range of 50cm: 5bins, evaluated at the center of each bin
		float dPhotonVertexZBinWidth;
		float dPhotonVertexZRangeLow;
		size_t dNumPhotonVertexZBins;

		//CHARGED TRACKS
		vector<const DChargedTrack*> dChargedTracks;
		unordered_map<Particle_t, vector<const JObject*>> dTracksByPID;

		//NEUTRAL SHOWERS
		unordered_map<signed char, DPhotonShowersByBeamBunch> dShowersByBeamBunchByZBin; //char: zbin

		//SOURCE COMBO INFOS: CREATED ONCE DURING DSOURCECOMBOER OBJECT CONSTRUCTION
			//with some exceptions (specific vertex-z, etc.)
		//want to make sure we only have one of each type: suggests using a set
		//however, after the first few events, almost all of these have already been created: vector has faster lookup time
		//therefore, use the set when creating the objects during construction, but then move the results into the vector and keep it sorted
		set<const DSourceComboInfo*, DCompare_SourceComboInfos> dSourceComboInfoSet;
		vector<const DSourceComboInfo*> dSourceComboInfos;
		unordered_map<const DSourceComboInfo*, Charge_t> dComboInfoChargeContent;
		unordered_set<const DSourceComboInfo*> dComboInfosWithMassiveNeutrals;
		//is this necessary??
		unordered_map<const DReactionVertexInfo*, DSourceComboUse> dSourceComboUseReactionMap_Primary;
		//the rest
		unordered_map<const DReactionStepVertexInfo*, DSourceComboUse> dSourceComboUseReactionMap; //primary combo info (nullptr if none)
		//combo use -> step
		unordered_map<pair<const DReactionStepVertexInfo*, DSourceComboUse>, size_t> dSourceComboInfoStepMap; //size_t: step index
		//i need to go from step -> combo use
		unordered_map<const DReaction*, map<size_t, DSourceComboUse>> dSourceComboUseReactionStepMap; //primary combo info (nullptr if none)
		//with specific vertex-z's
		unordered_map<pair<const DReactionVertexInfo*, vector<signed char>>, DSourceComboUse> dSourceComboUseVertexZMap;
		unordered_map<DSourceComboUse, DSourceComboUse> dZDependentUseToIndependentMap; //from z-dependent -> z-independent

		//SOURCE COMBOS //vector: z-bin //if attempted and all failed, DSourceCombosByUse_Large vector will be empty
		size_t dInitialComboVectorCapacity = 100;
		DSourceCombosByUse_Large dSourceCombosByUse;
		unordered_map<const DSourceCombo*, DSourceCombosByUse_Large> dMixedCombosByUseByChargedCombo; //key: charged combo //value: contains mixed & neutral combos //neutral: key is nullptr
		//also, sort by which beam bunches they are valid for: that way when comboing, we can retrieve only the combos that can possibly match the input RF bunches
		unordered_map<const DSourceCombo*, DSourceCombosByBeamBunchByUse> dSourceCombosByBeamBunchByUse; //key: charged combo //value: contains mixed & neutral combos: key is nullptr
		unordered_map<pair<const DSourceCombo*, const DReactionStepVertexInfo*>, const DSourceCombo*> dVertexPrimaryChargedComboMap; //first combo: reaction primary charged combo

		//RESUME SEARCH ITERATORS
		//e.g. if a DSourceCombo is -> 2pi0, and we want to use it as a basis for building a combo of 3pi0s,
		//then this iterator points to the first pi0 in the DSourceCombosByUse_Large vector that we want to test
		//that way we save a lot of time, since we don't have to look for it again
		//they are useful when comboing VERTICALLY, but cannot be used when comboing HORIZONTALLY
			//e.g. when comboing a pi0 (with photons = A, D) with a single photon, the photon could be B, C, or E+: no single spot to resume at
		unordered_map<pair<const JObject*, vector<int>>, vector<const JObject*>::const_iterator> dResumeSearchAfterIterators_Particles; //vector<int>: RF bunches (empty for all)
		unordered_map<pair<const DSourceCombo*, signed char>, DComboIteratorsByBeamBunch> dResumeSearchAfterIterators_Combos; //char: zbin

		//RESUME-SEARCH-AFTER OBJECTS
		//The below resume-at vectors are only useful when comboing vertically, so PID-specific versions are not needed
		unordered_map<const DSourceCombo*, const DSourceCombo*> dResumeSearchAfterMap_Combos; //key: Combo containing N combos of the type in the value //value: the last of those N combos
		unordered_map<const DSourceCombo*, const JObject*> dResumeSearchAfterMap_Particles; //key: Combo containing N particles of the type in the value //value: the last of those N particles

		//VALID RF BUNCHES BY COMBO
		unordered_map<const DSourceCombo*, vector<int>> dValidRFBunches_ByCombo;

		//RESOURCE POOLS
		DResourcePool<DSourceCombo> dResourcePool_SourceCombo;
		DResourcePool<vector<const DSourceCombo>> dResourcePool_SourceComboVector;

		//HANDLERS AND VERTEXERS
		DSourceComboVertexer* dSourceComboVertexer;
		DSourceComboP4Handler* dSourceComboP4Handler;
		DSourceComboTimeHandler* dSourceComboTimeHandler;
};

/*********************************************************** INLINE MEMBER FUNCTION DEFINITIONS ************************************************************/

inline const DSourceCombo* DSourceComboer::Get_ChargedCombo_WithNow(const DSourceCombo* locChargedCombo_Presiding) const
{
	if(locChargedCombo_Presiding == nullptr)
		return nullptr;

	for(const auto& locFurtherDecayPair : locChargedCombo_Presiding->Get_FurtherDecayCombos())
	{
		if(dComboInfoChargeContent[std::get<2>(locFurtherDecayPair.first)] == d_Charged)
			return locFurtherDecayPair.second[0]; //guaranteed to be size = 1
	}
	return nullptr; //uh oh ...
}


inline size_t DSourceComboer::Get_PhotonVertexZBin(double locVertexZ) const
{
	//given some vertex-z, what bin am I in?
	int locPhotonVertexZBin = int((locVertexZ - dPhotonVertexZRangeLow)/dPhotonVertexZBinWidth);
	if(locPhotonVertexZBin < 0)
		return 0;
	else if(locPhotonVertexZBin >= dNumPhotonVertexZBins)
		return dNumPhotonVertexZBins - 1;
	return locPhotonVertexZBin;
}

inline double DSourceComboer::Get_PhotonVertexZBinCenter(signed char locVertexZBin) const
{
	return dPhotonVertexZRangeLow + (double(locVertexZBin) + 0.5)*dPhotonVertexZBinWidth;
}

/*****
 * COMBOING PHOTONS AND RF BUNCHES
 *
 * So, this is tricky.
 * Start out by allowing ALL beam bunches, regardless of what the charged tracks want.
 * Then, as each photon is chosen, reduce the set of possible photons to choose next: only those that agree on at least one RF bunch
 * As combos are made, the valid RF bunches are saved along with the combo
 * That way, as combos are combined with other combos/particles, we make sure that only valid possibilities are chosen.
 *
 * We can't start with those only valid for the charged tracks because:
 * When we generate combos for a given info, we want to generate ALL combos at once.
 * E.g. some charged tracks may want pi0s with beam bunch = 1, but another group might want pi0s with bunch 1 OR 2.
 * Dealing with the overlap is a nightmare.  This avoids the problem entirely.
 *
 * BEWARE: Massive-neutral-particle momentum depends on the RF bunch. So a cut on the invariant mass with a neutron is effectively a cut on the RF bunches
 * Suppose: Sigma+ -> pi+ n
 * You first generate combos for -> pi+ n, and save them for the use X -> pi+, n
 * We then re-use the combos for the use Sigma+ -> pi+ n
 * But then a cut on the Sigma+ mass reduces the #valid RF bunches. So now we need a new combo!
 * We could decouple the RF bunches from the combo: e.g. save in map from combo_use -> rf bunches
 * However, this would result in many duplicate entries: e.g. X -> 2g, pi0 -> 2g, eta -> 2g, etc.
 * Users choosing final-state neutrons or KLongs is pretty rare compared to everything else: we are better off just creating new combos
 *
 * BEWARE: Massive-neutral-particle momentum depends on the RF bunch. So a cut on the invariant mass with a neutron is effectively a cut on the RF bunches.
 * So we can't actually vote on RF bunches until we choose our massive-neutral particles!!!
 */

inline void DSourceComboer::Build_ParticleIterators(const vector<int>& locBeamBunches, const vector<const JObject*>& locParticles)
{
	for(vector<const JObject*>::const_iterator locIterator = locParticles.begin(); locIterator != locParticles.end(); ++locIterator)
		dResumeSearchAfterIterators_Particles.emplace(std::make_pair(*locIterator, locBeamBunches), locIterator);
}

inline void DSourceComboer::Build_ComboIterators(const vector<int>& locBeamBunches, const vector<const DSourceCombo*>& locCombos, ComboingStage_t locComboingStage, signed char locVertexZBin)
{
	for(vector<const DSourceCombo*>::const_iterator locIterator = locCombos.begin(); locIterator != locCombos.end(); ++locIterator)
		dResumeSearchAfterIterators_Combos[std::make_pair(*locIterator, locVertexZBin)].emplace(locBeamBunches, locIterator);
}

inline vector<const JObject*>::const_iterator DSourceComboer::Get_ResumeAtIterator_Particles(const DSourceCombo* locSourceCombo, const vector<int>& locBeamBunches) const
{
	const auto& locPreviousObject = dResumeSearchAfterMap_Particles[locSourceCombo];
	return std::next(dResumeSearchAfterIterators_Particles[std::make_pair(locPreviousObject, locBeamBunches)]);
}

inline vector<const DSourceCombo*>::const_iterator DSourceComboer::Get_ResumeAtIterator_Combos(const DSourceCombo* locSourceCombo, const vector<int>& locBeamBunches, ComboingStage_t locComboingStage, signed char locVertexZBin) const
{
	const auto& locPreviousCombo = dResumeSearchAfterMap_Combos[locSourceCombo];
	return std::next(dResumeSearchAfterIterators_Combos[std::make_pair(locPreviousCombo, locVertexZBin)][locBeamBunches]);
}

inline bool DSourceComboer::Get_IsZIndependent(const JObject* locObject) const
{
	//make this virtual!
	auto locNeutralShower = dynamic_cast<const DNeutralShower*>(locObject);
	if(locNeutralShower == nullptr)
		return true;
	return (locNeutralShower->dDetectorSystem == SYS_FCAL);
}


inline DSourceCombosByUse_Large& DSourceComboer::Get_CombosSoFar(ComboingStage_t locComboingStage, Charge_t locChargeContent_SearchForUse, const DSourceCombo* locChargedCombo)
{
	//THE INPUT locChargedCombo MUST BE:
	//If reading from: Whatever charged combo you PREVIOUSLY comboed horizontally with to make the combos you're trying to get
	//If saving to (you are making a mixed): Whatever charged combo you are about to combo horizontally with to make this new, mixed combo

	//NOTE: If on mixed stage, it is NOT valid to get fully-charged combos from here! In fact, what you want is probably the input combo!
	if(locComboingStage == d_ChargedStage)
		return dSourceCombosByUse;
	else if(locChargeContent_SearchForUse == d_Neutral)
		return dMixedCombosByUseByChargedCombo[nullptr];
	return dMixedCombosByUseByChargedCombo[locChargedCombo];
}

inline DSourceCombosByBeamBunchByUse& DSourceComboer::Get_SourceCombosByBeamBunchByUse(Charge_t locChargeContent_SearchForUse, const DSourceCombo* locChargedCombo)
{
//shouldn't ever be called if charged
	//THE INPUT locChargedCombo MUST BE:
	//If reading from: Whatever charged combo you PREVIOUSLY comboed horizontally with to make the combos you're trying to get
	//If saving to (you are making a mixed): Whatever charged combo you just comboed horizontally with to make this new, mixed combo
	if(locChargeContent_SearchForUse == d_Neutral)
		return dSourceCombosByBeamBunchByUse[nullptr];
	return dSourceCombosByBeamBunchByUse[locChargedCombo];
}

inline DSourceComboer::~DSourceComboer(void)
{
//DELETE MORE THINGS!!
	for(auto locComboInfo : dSourceComboInfos)
		delete locComboInfo;
}


/*********************************************************** INLINE NAMESPACE-SCOPE FUNCTION DEFINITIONS ************************************************************/

} //end DAnalysis namespace

#endif // DSourceComboer_h
