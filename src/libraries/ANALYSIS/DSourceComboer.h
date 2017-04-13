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
//Don't place mass cuts on massive neutrals until vertex position defined
//once RF bunch is chosen, redo mass cuts involving massive neutrals
//fill in calc inv mass functions

//ANY TIME:
//Fill dShowerRFBunches_FCAL & dShowerRFBunches_Both
//Cut combo ahead of time if not enough tracks/showers

//AT THE END:
//store unsigned char instead of pointer???

/****************************************************** DEFINE LAMBDAS, USING STATEMENTS *******************************************************/

using DSourceCombosByUse_Large = unordered_map<DSourceComboUse, vector<const DSourceCombo*>*>; //for use collecting ALL of the combos in the event
auto Compare_SourceComboInfos = [](const DSourceComboInfo* lhs, const DSourceComboInfo* rhs) -> bool{return *lhs < *rhs;};

struct DCompare_SourceComboInfos{
	bool operator()(const DSourceComboInfo* lhs, const DSourceComboInfo* rhs) const{return *lhs < *rhs;}
};

/************************************************************** DEFINE CLASSES ***************************************************************/

enum ComboingStage_t
{
	d_ChargedOnly,
	d_FCALShowersOnly,
	d_AllShowers,
	d_Final
};

class DSourceComboer : public JObject
{
	public:

		DSourceComboer(void) = delete;
		DSourceComboer(JEventLoop* locEventLoop);
		DSourceComboer::~DSourceComboer(void);

		const vector<const DSourceCombo*>* Request_VertexCombos(JEventLoop* locEventLoop, const DReactionStepVertexInfo* locReactionStepVertexInfo);

		const vector<const DSourceCombo*>* Request_NeutralCombos(JEventLoop* locEventLoop, const DReactionStepVertexInfo* locReactionStepVertexInfo, DVector3 locVertex);

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

		//Create photon combo infos
		void Create_SourceComboInfos(const DReactionVertexInfo* locReactionVertexInfo);
		void Create_SourceComboInfos_Vertices(const DReactionVertexInfo* locReactionVertexInfo);
		void Create_SourceComboInfos_Neutrals(const shared_ptr<const DReactionStepVertexInfo>& locReactionStepVertexInfo);
		pair<bool, map<DSourceComboUse, unsigned char>> Get_FinalStateDecayingComboUses(const DReaction* locReaction, size_t locStepIndex, const map<size_t, DSourceComboUse>& locStepComboUseMap) const;
		const DSourceComboInfo* MakeOrGet_SourceComboInfo(const vector<pair<Particle_t, unsigned char>>& locNumParticles, const vector<pair<DSourceComboUse, unsigned char>>& locFurtherDecays);
		const DSourceComboInfo* GetOrMake_SourceComboInfo(const vector<pair<Particle_t, unsigned char>>& locNumParticles, const vector<pair<DSourceComboUse, unsigned char>>& locFurtherDecays);

		//TIMING METHODS
		void Calc_PhotonBeamBunchShifts(const DNeutralShower* locNeutralShower, shared_ptr<const DKinematicData>& locKinematicData, double locRFTime,
				DPhotonShowersByBeamBunch& locShowersByBeamBunch) const;
		double Calc_MaxDeltaTError(const DNeutralShower* locNeutralShower, const shared_ptr<const DKinematicData>& locKinematicData) const;

		//CREATE COMBOS - GENERAL METHODS
		void Create_SourceCombos(const DSourceComboUse& locSourceComboUse, ComboingStage_t locComboingStage, size_t locVertexZBin);
		void Create_SourceCombos(const DSourceComboInfo* locSourceComboInfo, ComboingStage_t locComboingStage, size_t locVertexZBin);

		//COMBO VERTICALLY METHODS
		void Combo_Vertically_AllDecays(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, size_t locVertexZBin);
		void Combo_Vertically_NDecays(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locNMinus1ComboUse, const DSourceComboUse& locSourceComboDecayUse, ComboingStage_t locComboingStage, size_t locVertexZBin);
		void Combo_Vertically_AllParticles(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, size_t locVertexZBin);
		void Combo_Vertically_NParticles(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locNMinus1ComboUse, ComboingStage_t locComboingStage, size_t locVertexZBin);

		//COMBO HORIZONTALLY METHODS
		void Combo_Horizontally_All(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, size_t locVertexZBin);
		void Create_Combo_OneParticle(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, size_t locVertexZBin);
		void Combo_Horizontally_AddCombo(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locAllBut1ComboUse, const DSourceComboUse& locSourceComboUseToAdd, ComboingStage_t locComboingStage, size_t locVertexZBin);
		void Combo_Horizontally_AddParticle(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locAllBut1ComboUse, Particle_t locPID, ComboingStage_t locComboingStage, size_t locVertexZBin);

		//BUILD/RETRIEVE RESUME-AT ITERATORS
		void Build_ParticleIterators(const vector<int>& locBeamBunches, const vector<const JObject*>& locParticles);
		vector<const JObject*>::const_iterator Get_ResumeAtIterator_Particles(const DSourceCombo* locSourceCombo, const vector<int>& locBeamBunches) const;
		void Build_ComboIterators(const vector<int>& locBeamBunches, const vector<const DSourceCombo*>& locCombos, ComboingStage_t locComboingStage, size_t locVertexZBin);
		vector<const DSourceCombo*>::const_iterator Get_ResumeAtIterator_Combos(const DSourceCombo* locSourceCombo, const vector<int>& locBeamBunches, ComboingStage_t locComboingStage, size_t locVertexZBin) const;

		//GET POTENTIAL PARTICLES & COMBOS FOR COMBOING
		const vector<const JObject*>& Get_ParticlesForComboing(Particle_t locPID, ComboingStage_t locComboingStage, const vector<int>& locBeamBunches = {}, size_t locVertexZBin = 0);
		vector<const JObject*>* Get_ShowersByBeamBunch(const vector<int>& locBeamBunches, DPhotonShowersByBeamBunch& locShowersByBunch);
		const vector<const DSourceCombo*>& Get_CombosForComboing(const DSourceComboUse& locComboUse, ComboingStage_t locComboingStage, const vector<int>& locBeamBunches, size_t locVertexZBin);
		const vector<const DSourceCombo*>& Get_CombosByBeamBunch(DCombosByBeamBunch& locCombosByBunch, const vector<int>& locBeamBunches, ComboingStage_t locComboingStage, size_t locVertexZBin);

		//GET/DETERMINE/REGISTER VALID RF BUNCHES
		const vector<int>& Get_ValidRFBunches(const DSourceComboUse& locSourceComboUse, const DSourceCombo* locSourceCombo) const;
		void Register_ValidRFBunches(const DSourceComboUse& locSourceComboUse, const DSourceCombo* locSourceCombo, const vector<int>& locRFBunches, ComboingStage_t locComboingStage, size_t locVertexZBin, bool locHasMassiveNeutrals);
		vector<int> Get_CommonRFBunches(const vector<int>& locRFBunches1, const vector<int>& locRFBunches2) const;

		//PARTICLE UTILITY FUNCTIONS
		shared_ptr<const DKinematicData> Create_KinematicData(const DNeutralShower* locNeutralShower, const DVector3& locVertex) const;
		bool Get_IsFCALOnly(const JObject* locObject) const;

		//COMBO UTILITY FUNCTIONS
		DSourceCombosByUse_Large& Get_CombosSoFar(ComboingStage_t locComboingStage, size_t locVertexZBin);
		void Copy_FCALOnlyResults(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, size_t locVertexZBin);
		bool Get_HasMassiveNeutrals(const DSourceComboInfo* locSourceComboInfo) const;

		//VERTEX/RF UTILITY FUNCTIONS
		size_t Get_PhotonVertexZBin(double locVertexZ) const;
		double Get_PhotonVertexZBinCenter(size_t locVertexZBin) const;
		int Calc_RFBunchShift(double locTimeToStep, double locTimeToStepTo) const; //returns integer shift

		//MASS UTILITY FUNCTIONS
		bool Cut_InvariantMass(const DSourceCombo* locSourceCombo, Particle_t locDecayPID) const;
		double Calc_InvariantMass(const DSourceCombo* locSourceCombo) const;
		vector<int> Cut_InvariantMass_HasMassiveNeutral(const DSourceCombo* locSourceCombo, Particle_t locDecayPID, vector<int> locValidRFBunches) const;
		double Calc_InvariantMass_HasMassiveNeutral(const DSourceCombo* locSourceCombo, int locRFBunch) const;

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
		unordered_map<shared_ptr<const DReactionStepVertexInfo>, DSourceComboUse> dSourceComboUseReactionMap_Photons; //primary combo info (nullptr if none)
		unordered_map<shared_ptr<const DReactionStepVertexInfo>, DSourceComboUse> dSourceComboUseReactionMap_Charged; //primary combo info (nullptr if none)
		unordered_map<pair<shared_ptr<const DReactionStepVertexInfo>, DSourceComboUse>, size_t> dSourceComboInfoStepMap_Photons; //size_t: step index
		unordered_map<pair<shared_ptr<const DReactionStepVertexInfo>, DSourceComboUse>, size_t> dSourceComboInfoStepMap_Charged; //size_t: step index

		//SOURCE COMBOS //vector: z-bin //if attempted and all failed, DSourceCombosByUse_Large vector will be empty
		size_t dInitialComboVectorCapacity = 100;
		DSourceCombosByBeamBunchByUse dSourceCombos_Charged;
		DSourceCombosByUse_Large dSourceCombos_FCAL;
		vector<DSourceCombosByUse_Large> dSourceCombos_Both;
		//also, sort by which beam bunches they are valid for: that way when comboing, we can retrieve only the combos that can possibly match the input RF bunches
		DSourceCombosByBeamBunchByUse dSourceCombosByBeamBunch_FCAL;
		vector<DSourceCombosByBeamBunchByUse> dSourceCombosByBeamBunch_Both; //vector: vertex-z bins //FCAL + BCAL

		//RESUME SEARCH ITERATORS
		//e.g. if a DSourceCombo is -> 2pi0, and we want to use it as a basis for building a combo of 3pi0s,
		//then this iterator points to the first pi0 in the DSourceCombosByUse_Large vector that we want to test
		//that way we save a lot of time, since we don't have to look for it again
		//they are useful when comboing VERTICALLY, but cannot be used when comboing HORIZONTALLY
			//e.g. when comboing a pi0 (with photons = A, D) with a single photon, the photon could be B, C, or E+: no single spot to resume at
		unordered_map<pair<const JObject*, vector<int>>, vector<const JObject*>::const_iterator> dResumeSearchAfterIterators_Particles; //vector<int>: RF bunches (empty for all)
		unordered_map<const DSourceCombo*, DComboIteratorsByBeamBunch> dResumeSearchAfterIterators_Combos;
		vector<unordered_map<const DSourceCombo*, DComboIteratorsByBeamBunch>> dResumeSearchAfterIterators_AllShowerCombos; //first vector: vertex-z bins

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

inline void DSourceComboer::Build_ComboIterators(const vector<int>& locBeamBunches, const vector<const DSourceCombo*>& locCombos, ComboingStage_t locComboingStage, size_t locVertexZBin)
{
	auto& locIteratorComboMap = (locComboingStage != d_AllShowers) ? dResumeSearchAfterIterators_Combos : dResumeSearchAfterIterators_AllShowerCombos[locVertexZBin];
	for(vector<const DSourceCombo*>::const_iterator locIterator = locCombos.begin(); locIterator != locCombos.end(); ++locIterator)
		locIteratorComboMap[*locIterator].emplace(locBeamBunches, locIterator);
}

inline vector<const JObject*>::const_iterator DSourceComboer::Get_ResumeAtIterator_Particles(const DSourceCombo* locSourceCombo, const vector<int>& locBeamBunches) const
{
	const auto& locPreviousObject = dResumeSearchAfterMap_Particles[locSourceCombo];
	return std::next(dResumeSearchAfterIterators_Particles[std::make_pair(locPreviousObject, locBeamBunches)]);
}

inline vector<const DSourceCombo*>::const_iterator DSourceComboer::Get_ResumeAtIterator_Combos(const DSourceCombo* locSourceCombo, const vector<int>& locBeamBunches, ComboingStage_t locComboingStage, size_t locVertexZBin) const
{
	const auto& locPreviousCombo = dResumeSearchAfterMap_Combos[locSourceCombo];
	if(locComboingStage == d_AllShowers)
		return std::next(dResumeSearchAfterIterators_AllShowerCombos[locVertexZBin][locPreviousCombo][locBeamBunches]);
	else
		return std::next(dResumeSearchAfterIterators_Combos[locPreviousCombo][locBeamBunches]);
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





inline DSourceCombosByUse_Large& DSourceComboer::Get_CombosSoFar(ComboingStage_t locComboingStage, size_t locVertexZBin)
{
	if(locComboingStage == d_ChargedOnly)
		return dSourceCombos_Charged;
	else if(locComboingStage == d_FCALShowersOnly)
		return dSourceCombos_FCAL;
	else if(locComboingStage == d_AllShowers)
		return dSourceCombos_Both[locVertexZBin];
	else
		return {};

	//CAREFUL: FINAL STAGE!!!
	//Need to save separately for: Passing charged PID cuts at production vertex


	//Assume: only one reaction, no neutrons, no detached vertices, vertex position findable
	//For a reaction: for each step/vertex (only 1 here), do charged comboing
	//Then, find the vertex and select the RF bunch(es) (if any)
		//Store results as map from use to combo to pair of vertex/rf_bunches (if cut, rf will be empty)
	//for each rf of each use/combo, find photons at production vertex (all?), ignoring charged tracks
		//do this by requesting the combo results for the primary vertex of the reaction
		//store results as map of vertex-z-bin to rf to use to combo
	//then,

	//MUST BEWARE DUPLICATE COMBOS
	//let's say a combo of charged tracks has 2 valid RF bunches
	//and we need to combo 2 pi0s with them
	//and the shower timing cuts are loose enough that all 4 showers satisfy both RF bunches
	//if we combo the 2 rf bunches separately: WE HAVE DUPLICATE COMBOS
	//and doing the duplicate check AFTER the fact takes FOREVER
	//therefore, we must take the neutral showers for the 2 rfs, COMBINE THEM, and then COMBO AS A UNIT

	//once we choose a shower that has less valid rf bunches than the total, we can reduce the size of the potential shower pool
		//we must be careful about the resume iterators!! they are saved when you are looking at (e.g.) 2 rf bunches, but on lookup may be only 1!!
			//or when we generate the "needed" vector we can apply a cut with the iterator
		//or, we can manually cut the showers by rf
//Where do I register which RF bunches are available as I'm comboing? For the saved combos?
	//you may create the combo while restricted to only one RF bunch, but it may work several others too: won't know at creation time
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

inline double DSourceComboer::Get_PhotonVertexZBinCenter(size_t locVertexZBin) const
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

inline double DSourceComboer::Calc_InvariantMass(const DSourceCombo* locSourceCombo, size_t locVertexZBin) const
{
	DLorentzVector locTotalP4;
	dFinalStateP4ByCombo.emplace(locSourceCombo, locTotalP4);
	return locTotalP4.M();
}

inline double DSourceComboer::Calc_InvariantMass_HasMassiveNeutral(const DSourceCombo* locSourceCombo, int locRFBunch, size_t locVertexZBin) const
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
		}
	}

	dFinalStateP4ByCombo_HasMassiveNeutrals[locRFBunch].emplace(locSourceCombo, locTotalP4);
	return locTotalP4.M();
}



/*********************************************************** INLINE NAMESPACE-SCOPE FUNCTION DEFINITIONS ************************************************************/

} //end DAnalysis namespace

#endif // DSourceComboer_h
