#ifndef DSourceComboer_h
#define DSourceComboer_h

#include <map>
#include <set>
#include <vector>
#include <memory>
#include <unordered_map>
#include <cmath>
#include <algorithm>

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

using namespace std;
using namespace jana;

namespace DAnalysis
{

//BIG TO DO'S:
//Figure out how combos will be requested, and then combined together
//Don't place mass cuts on massive neutrals until vertex position fully defined
//once RF bunch is chosen, redo mass cuts involving massive neutrals
//fill in calc inv mass functions
//don't pass vert-z-bin through everything: get from combo use instead

//ANY TIME:
//Fill dShowerRFBunches_FCAL & dShowerRFBunches_Both
//compute and save Get_ChargeContent() ahead of time (and just read it)
//Cut combo ahead of time if not enough tracks/showers
//Switch DSourceCombosByUse back to Large & Small versions

//AT THE END:
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

/************************************************************** DEFINE CLASSES ***************************************************************/

enum ComboingStage_t
{
	d_ChargedStage,
	d_MixedStage_FCAL,
	d_MixedStage
};

class DSourceComboer : public JObject
{
	public:

		DSourceComboer(void) = delete;
		DSourceComboer(JEventLoop* locEventLoop);
		DSourceComboer::~DSourceComboer(void);

	private:

		/********************************************************** DEFINE USING STATEMENTS ***********************************************************/

		//DEFINE USING STATEMENTS
		using DPhotonShowersByBeamBunch = unordered_map<vector<int>, vector<const JObject*>>; //int: beam bunch n-shifts from nominal
		using DCombosByBeamBunch = unordered_map<vector<int>, vector<const DSourceCombo*>>;
		using DSourceCombosByBeamBunchByUse = unordered_map<DSourceComboUse, DCombosByBeamBunch>;
		using DComboIteratorsByBeamBunch = unordered_map<vector<int>, vector<const DSourceCombo*>::const_iterator>; //vector<int>: RF bunches (empty for all)

		/********************************************************** DECLARE MEMBER FUNCTIONS ***********************************************************/

		void Define_LooseCuts(void);
		void Reset_NewEvent(JEventLoop* locEventLoop);
		void Setup_NeutralShowers(JEventLoop* locEventLoop);

		//CREATE PHOTON COMBO INFOS
		void Create_SourceComboInfos(const DReactionVertexInfo* locReactionVertexInfo);

		//CREATE PHOTON COMBO INFOS: UTILITY METHODS
		map<Particle_t, unsigned char> Build_ParticleMap(const DReaction* locReaction, size_t locStepIndex, Charge_t locCharge) const;
		pair<bool, map<DSourceComboUse, unsigned char>> Get_FinalStateDecayingComboUses(const DReaction* locReaction, size_t locStepIndex, const map<size_t, DSourceComboUse>& locStepComboUseMap) const;
		DSourceComboUse Make_ComboUse(Particle_t locInitPID, const map<Particle_t, unsigned char>& locNumParticles, const map<DSourceComboUse, unsigned char>& locFurtherDecays);
		const DSourceComboInfo* MakeOrGet_SourceComboInfo(const vector<pair<Particle_t, unsigned char>>& locNumParticles, const vector<pair<DSourceComboUse, unsigned char>>& locFurtherDecays);
		const DSourceComboInfo* GetOrMake_SourceComboInfo(const vector<pair<Particle_t, unsigned char>>& locNumParticles, const vector<pair<DSourceComboUse, unsigned char>>& locFurtherDecays);

		//TIMING METHODS
		void Calc_PhotonBeamBunchShifts(const DNeutralShower* locNeutralShower, shared_ptr<const DKinematicData>& locKinematicData, double locRFTime,
				DPhotonShowersByBeamBunch& locShowersByBeamBunch) const;
		double Calc_MaxDeltaTError(const DNeutralShower* locNeutralShower, const shared_ptr<const DKinematicData>& locKinematicData) const;

		//CREATE COMBOS - GENERAL METHODS
		void Create_SourceCombos(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, signed char locVertexZBin, const DSourceCombo* locChargedCombo_WithNow);
		void Create_SourceCombos(const DSourceComboInfo* locSourceComboInfo, ComboingStage_t locComboingStage, signed char locVertexZBin, const DSourceCombo* locChargedCombo_WithNow);

		//COMBO VERTICALLY METHODS
		//Note that vertical comboing always takes place at the same vertex-z
		void Combo_Vertically_AllDecays(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_WithNow);
		void Combo_Vertically_NDecays(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locNMinus1ComboUse, const DSourceComboUse& locSourceComboDecayUse, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_WithNow);
		void Combo_Vertically_AllParticles(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage);
		void Combo_Vertically_NParticles(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locNMinus1ComboUse, ComboingStage_t locComboingStage);

		//COMBO HORIZONTALLY METHODS
		void Combo_Horizontally_All(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_WithNow);
		void Create_Combo_OneParticle(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage);
		void Combo_Horizontally_AddCombo(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locAllBut1ComboUse, const DSourceComboUse& locSourceComboUseToAdd, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_WithNow);
		void Combo_Horizontally_AddParticle(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locAllBut1ComboUse, Particle_t locPID, ComboingStage_t locComboingStage);

		//BUILD/RETRIEVE RESUME-AT ITERATORS
		void Build_ParticleIterators(const vector<int>& locBeamBunches, const vector<const JObject*>& locParticles);
		vector<const JObject*>::const_iterator Get_ResumeAtIterator_Particles(const DSourceCombo* locSourceCombo, const vector<int>& locBeamBunches) const;
		void Build_ComboIterators(const vector<int>& locBeamBunches, const vector<const DSourceCombo*>& locCombos, ComboingStage_t locComboingStage, signed char locVertexZBin);
		vector<const DSourceCombo*>::const_iterator Get_ResumeAtIterator_Combos(const DSourceCombo* locSourceCombo, const vector<int>& locBeamBunches, ComboingStage_t locComboingStage, signed char locVertexZBin) const;

		//GET POTENTIAL PARTICLES & COMBOS FOR COMBOING
		const vector<const JObject*>& Get_ParticlesForComboing(Particle_t locPID, ComboingStage_t locComboingStage, const vector<int>& locBeamBunches = {}, signed char locVertexZBin = 0);
		vector<const JObject*>* Get_ShowersByBeamBunch(const vector<int>& locBeamBunches, DPhotonShowersByBeamBunch& locShowersByBunch);
		const vector<const DSourceCombo*>& Get_CombosForComboing(const DSourceComboUse& locComboUse, ComboingStage_t locComboingStage, const vector<int>& locBeamBunches, const DSourceCombo* locChargedCombo_WithPrevious);
		const vector<const DSourceCombo*>& Get_CombosByBeamBunch(DCombosByBeamBunch& locCombosByBunch, const vector<int>& locBeamBunches, ComboingStage_t locComboingStage, signed char locVertexZBin);

		//GET/DETERMINE/REGISTER VALID RF BUNCHES
		const vector<int>& Get_ValidRFBunches(const DSourceComboUse& locSourceComboUse, const DSourceCombo* locSourceCombo) const;
		void Register_ValidRFBunches(const DSourceComboUse& locSourceComboUse, const DSourceCombo* locSourceCombo, const vector<int>& locRFBunches, ComboingStage_t locComboingStage, bool locHasMassiveNeutrals, const DSourceCombo* locChargedCombo_WithNow);
		vector<int> Get_CommonRFBunches(const vector<int>& locRFBunches1, const vector<int>& locRFBunches2) const;

		//PARTICLE UTILITY FUNCTIONS
		shared_ptr<const DKinematicData> Create_KinematicData(const DNeutralShower* locNeutralShower, const DVector3& locVertex) const;
		bool Get_IsFCALOnly(const JObject* locObject) const;

		//COMBO UTILITY FUNCTIONS
		DSourceCombosByUse& Get_CombosSoFar(ComboingStage_t locComboingStage, Charge_t locChargeContent_SearchForUse, const DSourceCombo* locChargedCombo = nullptr);
		DSourceCombosByBeamBunchByUse& Get_SourceCombosByBeamBunchByUse(Charge_t locChargeContent_SearchForUse, const DSourceCombo* locChargedCombo = nullptr);
		void Copy_FCALOnlyResults(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_WithNow);
		bool Get_HasMassiveNeutrals(const DSourceComboInfo* locSourceComboInfo) const;

		//VERTEX/RF UTILITY FUNCTIONS
		size_t Get_PhotonVertexZBin(double locVertexZ) const;
		double Get_PhotonVertexZBinCenter(signed char locVertexZBin) const;
		int Calc_RFBunchShift(double locTimeToStep, double locTimeToStepTo) const; //returns integer shift

		//MASS UTILITY FUNCTIONS
		bool Cut_InvariantMass(const DSourceCombo* locSourceCombo, Particle_t locDecayPID) const;
		double Calc_InvariantMass(const DSourceCombo* locSourceCombo) const;
		vector<int> Cut_InvariantMass_HasMassiveNeutral(const DSourceCombo* locSourceCombo, Particle_t locDecayPID, vector<int> locValidRFBunches) const;
		double Calc_InvariantMass_HasMassiveNeutral(const DSourceCombo* locSourceCombo, int locRFBunch) const;

		//definitions of negative values for any particle index //in order such that operator< returns order expected for string (e.g. gp->...)
		static signed char Get_VertexZIndex_FCAL(void){return -2;}
		static signed char Get_VertexZIndex_Unknown(void){return -1;}

		/************************************************************** DEFINE MEMBERS ***************************************************************/

		uint64_t dEventNumber = 0; //re-setup on new events
		string dShowerSelectionTag = "PreSelect";

		//EXPERIMENT INFORMATION
		DVector3 dTargetCenter;
		double dTargetLength = 30.0;
		double dBeamBunchPeriod = 1000.0/249.5;

		//VERTEX-DEPENDENT PHOTON INFORMATION
		//For every 10cm in vertex-z, calculate the photon p4 & time for placing mass & delta-t cuts
		//The z-range extends from the upstream end of the target - 5cm to the downstream end + 15cm
		//so for a 30-cm-long target, it's a range of 50cm: 5bins, evaluated at the center of each bin
		float dPhotonVertexZBinWidth = 10.0;
		float dPhotonVertexZRangeLow = 45.0;
		size_t dNumPhotonVertexZBins = 5;

		//due to detached vertices
		double dMaxDecayTimeOffset = 2.0;

		//CHARGED TRACKS
		vector<const DChargedTrack*> dChargedTracks;
		unordered_map<Particle_t, vector<const JObject*>> dTracksByPID;

		//NEUTRAL SHOWER DATA
		unordered_map<const JObject*, vector<int>> dShowerRFBunches_FCAL; //VECTOR MUST BE SORTED!!
		vector<unordered_map<const JObject*, vector<int>>> dShowerRFBunches_Both; //vector: vertex-z bins
		unordered_map<const DNeutralShower*, shared_ptr<const DKinematicData>> dFCALKinematics; //FCAL shower data at center of target
		vector<unordered_map<const DNeutralShower*, shared_ptr<const DKinematicData>>> dBCALKinematics; //BCAL shower data in vertex-z bins

		//SHOWERS SORTED BY RF BUNCH
		const DEventRFBunch* dInitialEventRFBunch;
		DPhotonShowersByBeamBunch dFCALPhotonShowersByBeamBunch;
		vector<DPhotonShowersByBeamBunch> dPhotonShowersByBeamBunch; //vector: vertex-z bins //FCAL + BCAL

		//CUTS
		unordered_map<DetectorSystem_t, TF1*> dPhotonTimeCutMap; //function of shower energy (p)
		unordered_map<Particle_t, pair<double, double>> dInvariantMassCuts;

		//SOURCE COMBO INFOS: CREATED ONCE DURING DSourceComboER OBJECT CONSTRUCTION
		//want to make sure we only have one of each type: suggests using a set
		//however, after the first few events, almost all of these have already been created: vector has faster lookup time
		//therefore, use the set when creating the objects during construction, but then move the results into the vector and keep it sorted
		set<const DSourceComboInfo*, DCompare_SourceComboInfos> dSourceComboInfoSet;
		vector<const DSourceComboInfo*> dSourceComboInfos;
		//is this necessary??
		unordered_map<const DReactionVertexInfo*, DSourceComboUse> dSourceComboUseReactionMap_Primary;
		//the rest
		unordered_map<shared_ptr<const DReactionStepVertexInfo>, DSourceComboUse> dSourceComboUseReactionMap; //primary combo info (nullptr if none)
		unordered_map<pair<shared_ptr<const DReactionStepVertexInfo>, DSourceComboUse>, size_t> dSourceComboInfoStepMap; //size_t: step index
		//i need to go from step -> combo use
		unordered_map<const DReaction*, map<size_t, DSourceComboUse>> dSourceComboUseReactionStepMap; //primary combo info (nullptr if none)

		//SOURCE COMBOS //vector: z-bin //if attempted and all failed, DSourceCombosByUse_Small vector will be empty
		size_t dInitialComboVectorCapacity = 100;
		DSourceCombosByUse dSourceCombosByUse;
		unordered_map<const DSourceCombo*, DSourceCombosByUse> dMixedCombosByUseByChargedCombo; //key: charged combo //value: contains mixed & neutral combos //neutral: key is nullptr
		//also, sort by which beam bunches they are valid for: that way when comboing, we can retrieve only the combos that can possibly match the input RF bunches
		unordered_map<const DSourceCombo*, DSourceCombosByBeamBunchByUse> dSourceCombosByBeamBunchByUse; //key: charged combo //value: contains mixed & neutral combos: key is nullptr

		//RESUME SEARCH ITERATORS
		//e.g. if a DSourceCombo is -> 2pi0, and we want to use it as a basis for building a combo of 3pi0s,
		//then this iterator points to the first pi0 in the DSourceCombosByUse_Small vector that we want to test
		//that way we save a lot of time, since we don't have to look for it again
		//they are useful when comboing VERTICALLY, but cannot be used when comboing HORIZONTALLY
			//e.g. when comboing a pi0 (with photons = A, D) with a single photon, the photon could be B, C, or E+: no single spot to resume at
		unordered_map<pair<const JObject*, vector<int>>, vector<const JObject*>::const_iterator> dResumeSearchAfterIterators_Particles; //vector<int>: RF bunches (empty for all)
		unordered_map<pair<const DSourceCombo*, signed char>, DComboIteratorsByBeamBunch> dResumeSearchAfterIterators_Combos; //char: zbin

		//RESUME-SEARCH-AFTER OBJECTS
		//The below resume-at vectors are only useful when comboing vertically, so PID-specific versions are not needed
		unordered_map<const DSourceCombo*, const DSourceCombo*> dResumeSearchAfterMap_Combos; //key: Combo containing N combos of the type in the value //value: the last of those N combos
		unordered_map<const DSourceCombo*, const JObject*> dResumeSearchAfterMap_Particles; //key: Combo containing N particles of the type in the value //value: the last of those N particles

		//VALID RF BUNCHES BY COMBO AND USE
		unordered_map<const DSourceCombo*, vector<int>> dValidRFBunches_ByCombo;
		unordered_map<pair<DSourceComboUse, const DSourceCombo*>, vector<int>> dValidRFBunches_ByUse; //only when massive neutral particles are somewhere in the combo

		//TOTAL FINAL STATE FOUR-MOMENTUM
		unordered_map<const DSourceCombo*, DLorentzVector> dFinalStateP4ByCombo;
		unordered_map<int, unordered_map<const DSourceCombo*, DLorentzVector>> dFinalStateP4ByCombo_HasMassiveNeutrals; //int: RF bunch
};

/*********************************************************** INLINE MEMBER FUNCTION DEFINITIONS ************************************************************/
/*
//prefer some kind of map rather than doing this over and over again
//ugh: no, would be huge
inline const DSourceCombo* DSourceComboer::Get_StepPrimaryVertexCombo_Charged(const DSourceCombo* locReactionCombo, const DReactionVertexInfo* locReactionVertexInfo, const DReactionStepVertexInfo* locStepVertexInfo)
{
	//if it's the production vertex, just return the input
	if(locStepVertexInfo->Get_ProductionVertexFlag())
		return locReactionCombo;

	//collect all combos of the desired use that are in the reaction combo
	auto locSourceComboUse = dSourceComboUseReactionMap_Charged[locStepVertexInfo];
	auto locCombos = Get_AllCombosByUse(locReactionCombo, locSourceComboUse);

	//if there is only one we can just return it.  but if there are more, we must return the correct one
	//e.g. if the primary step vertex is the N'th occurrence of the use, we must return the N'th combo
	if(locCombos.size() == 1)
		return locCombos[0];

	//get all step vertices, sort by step order
	auto locStepVertexInfos = locReactionVertexInfo->Get_StepVertexInfos();
	auto Sort_StepVertices = [](const shared_ptr<const DReactionStepVertexInfo>& lhs, const shared_ptr<const DReactionStepVertexInfo>& rhs) -> bool
			{return lhs->Get_StepIndices().front() < rhs->Get_StepIndices().front();};
	std::sort(locStepVertexInfos.begin(), locStepVertexInfos.end(), Sort_StepVertices);

	//now, find what occurrence of the use
	auto locComboIndex = size_t(0);
	for(auto locLoopStepVertexInfo : locStepVertexInfos)
	{
		if(locLoopStepVertexInfo == locStepVertexInfo)
			break;
		if(dSourceComboUseReactionMap_Charged[locLoopStepVertexInfo] == locSourceComboUse)
			++locComboIndex;
	}

	//return the corresponding combo
	return locCombos[locComboIndex];
}

inline const DSourceCombo* DSourceComboer::Get_StepCombo(const DSourceCombo* locReactionCombo, const DReaction* locReaction, size_t locStepIndex)
{
	auto locSourceComboUse = dSourceComboUseReactionMap_Charged[locStepVertexInfo];
	return Get_
}
*/

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

inline bool DSourceComboer::Get_IsFCALOnly(const JObject* locObject) const
{
	auto locNeutralShower = dynamic_cast<const DNeutralShower*>(locObject);
	if(locNeutralShower == nullptr)
		return false;
	return (locNeutralShower->dDetectorSystem == SYS_FCAL);
}

inline bool DSourceComboer::Get_HasMassiveNeutrals(const DSourceComboInfo* locComboInfo) const
{
	//see if the combo info contains a massive neutral particle

	//search function
	auto Find_MassiveNeutrals = [](const pair<Particle_t, unsigned char>& locPair) -> bool
		{return (ParticleCharge(locPair.first) != 0) && (ParticleMass(locPair.first) > 0.0);};

	//do search
	auto locNumParticles = locComboInfo->Get_NumParticles(true); //true: entire chain
	return std::any_of(locNumParticles.begin(), locNumParticles.end(), Find_MassiveNeutrals);
}

inline vector<int> DSourceComboer::Get_CommonRFBunches(const vector<int>& locRFBunches1, const vector<int>& locRFBunches2) const
{
	//check to see if one of the input vectors is empty //empty means "no idea": all possible bunches are valid
	if(locRFBunches1.empty())
		return locRFBunches2;
	else if(locRFBunches2.empty())
		return locRFBunches1;

	vector<int> locCommonRFBunches = {}; //if charged or massive neutrals, ignore (they don't choose at this stage)
	locCommonRFBunches.reserve(locRFBunches1.size() + locRFBunches2.size());
	std::set_intersection(locRFBunches1.begin(), locRFBunches1.end(), locRFBunches2.begin(), locRFBunches2.end(), std::back_inserter(locCommonRFBunches));
	return locCommonRFBunches;
}

//Mixed results are saved in: (where the keys are the charged contents of the mixed-use step)

inline DSourceCombosByUse& DSourceComboer::Get_CombosSoFar(ComboingStage_t locComboingStage, Charge_t locChargeContent_SearchForUse, const DSourceCombo* locChargedCombo)
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

inline int DSourceComboer::Calc_RFBunchShift(double locTimeToStep, double locTimeToStepTo) const
{
	double locDeltaT = locTimeToStepTo - locTimeToStep;
	return (locDeltaT > 0.0) ? int(locDeltaT/dBeamBunchPeriod + 0.5) : int(locDeltaT/dBeamBunchPeriod - 0.5);
}

inline shared_ptr<const DKinematicData> DSourceComboer::Create_KinematicData(const DNeutralShower* locNeutralShower, const DVector3& locVertex) const
{
	DVector3 locPath = locNeutralShower->dSpacetimeVertex.Vect() - locVertex;
	double locPathLength = locPath.Mag();
	double locVertexTime = locNeutralShower->dSpacetimeVertex.T() - locPathLength/29.9792458;
	DVector3 locMomentum = locNeutralShower->dEnergy*locPath.Unit();
	return std::make_shared<const DKinematicData>(Gamma, locMomentum, locVertex, locVertexTime);
}

inline bool DSourceComboer::Cut_InvariantMass(const DSourceCombo* locSourceCombo, Particle_t locDecayPID) const
{
	//Don't call if it contains massive neutrals! Call the other cut function instead!!
	auto locCutIterator = dInvariantMassCuts.find(locDecayPID);
	if(locCutIterator == dInvariantMassCuts.end())
		return true; //no cut to place!!
	auto& locMassCuts = locCutIterator->second;

	auto locInvariantMass = Calc_InvariantMass(locSourceCombo);
	return ((locInvariantMass >= locMassCuts.first) && (locInvariantMass <= locMassCuts.second));
}

inline vector<int> DSourceComboer::Cut_InvariantMass_HasMassiveNeutral(const DSourceCombo* locSourceCombo, Particle_t locDecayPID, vector<int> locValidRFBunches) const
{
	//cuts on possible RF bunches for the massive neutrals
	//if no possible rf bunch yields a massive-neutral-momentum that passes the invariant mass cut, returns an empty vector
	auto locCutIterator = dInvariantMassCuts.find(locDecayPID);
	if(locCutIterator == dInvariantMassCuts.end())
		return locSourceCombo; //no cut to place!!
	auto& locMassCuts = locCutIterator->second;

	//function for calculating and cutting the invariant mass for each rf bunch
	auto CalcAndCut_InvariantMass = [&locSourceCombo, &locMassCuts](int locRFBunch) -> bool
	{
		auto locInvariantMass = Calc_InvariantMass_HasMassiveNeutral(locSourceCombo, locRFBunch);
		return ((locInvariantMass >= locMassCuts.first) && (locInvariantMass <= locMassCuts.second));
	};

	//apply the function
	locValidRFBunches.erase(std::remove_if(locValidRFBunches.begin(), locValidRFBunches.end(), CalcAndCut_InvariantMass), locValidRFBunches.end());
	return locValidRFBunches;
}

inline double DSourceComboer::Calc_InvariantMass(const DSourceCombo* locSourceCombo, signed char locVertexZBin) const
{
	DLorentzVector locTotalP4;
	dFinalStateP4ByCombo.emplace(locSourceCombo, locTotalP4);
	return locTotalP4.M();
}

inline DLorentzVector DSourceComboer::Calc_DecayingP4_ChargedOnly(const DSourceCombo* locSourceCombo) const
{
	//see if it's already been calculated
	auto locCalcIterator = dFinalStateP4ByCombo.find(locSourceCombo);
	if(locCalcIterator != dFinalStateP4ByCombo.end())
		return *locCalcIterator;

	//the input DSourceCombo CANNOT contain any neutrals (or else this will crash!!)
	auto locSourceParticles = locSourceCombo->Get_SourceParticles(false); //false: NOT the whole chain

	DLorentzVector locTotalP4;
	for(auto locParticlePair : locSourceParticles)
	{
		auto locChargedTrack = static_cast<const DChargedTrack*>(locParticlePair.second);
		locTotalP4 += locChargedTrack->Get_Hypothesis(locParticlePair.first)->lorentzMomentum();
	}

	//now loop over decays
	auto locFurtherDecayCombos = locSourceCombo->Get_FurtherDecayCombos();
	for(const auto& locCombosByUsePair : locFurtherDecayCombos)
	{
		for(const auto& locCombo : locCombosByUsePair.second)
		{
			//this check is duplicated at the top!
			auto locIterator = dFinalStateP4ByCombo.find(locCombo);
			if(locIterator != dFinalStateP4ByCombo.end())
				locTotalP4 += locIterator->second;
			else
				locTotalP4 += Calc_DecayingP4_ChargedOnly(locCombo);
		}
	}

	dFinalStateP4ByCombo.emplace(locSourceCombo, locTotalP4);
	return locTotalP4;
}

inline double DSourceComboer::Calc_InvariantMass_HasMassiveNeutral(const DSourceCombo* locSourceCombo, int locRFBunch, signed char locVertexZBin) const
{
	auto locSourceParticles = locSourceCombo->Get_SourceParticles(false); //false: NOT the whole chain
	//vertex-z bin may be different for decay products! (detached vertex)
	//save/retrieve masses by combo instead

	DLorentzVector locTotalP4;
	for(auto locParticlePair : locSourceParticles)
	{
		auto locPID = locParticlePair.first;
		if(ParticleCharge(locPID) != 0)
		{
			auto locChargedTrack = static_cast<const DChargedTrack*>(locParticlePair.second);
			locTotalP4 += locChargedTrack->Get_Hypothesis(locPID)->lorentzMomentum();
			continue;
		}

		auto locNeutralShower = static_cast<const DNeutralShower*>(locParticlePair.second);
		if(locPID == Gamma)
		{
			auto& locKinematicsMap = (locNeutralShower->dDetectorSystem == SYS_FCAL) ? dFCALKinematics : dBCALKinematics[locVertexZBin];
			auto& locKinematicData = locKinematicsMap[locNeutralShower];
			locTotalP4 += locKinematicData->lorentzMomentum();
			continue;
		}

//handle massive neutral case!!
	}

	//now loop over decays
	auto locFurtherDecayCombos = locSourceCombo->Get_FurtherDecayCombos();
	for(const auto& locCombosByUsePair : locFurtherDecayCombos)
	{
		for(const auto& locCombo : locCombosByUsePair.second)
		{
			auto locIterator = dFinalStateP4ByCombo.find(locCombo);
			if(locIterator != dFinalStateP4ByCombo.end())
				locTotalP4 += locIterator->second;
			else
				locTotalP4 += dFinalStateP4ByCombo_HasMassiveNeutrals[locRFBunch][locCombo];
			//else call Calc_InvariantMass_HasMassiveNeutral!!
		}
	}

	dFinalStateP4ByCombo_HasMassiveNeutrals[locRFBunch].emplace(locSourceCombo, locTotalP4);
	return locTotalP4.M();
}



/*********************************************************** INLINE NAMESPACE-SCOPE FUNCTION DEFINITIONS ************************************************************/

} //end DAnalysis namespace

#endif // DSourceComboer_h
