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
//Extend to photon
//Extend to charged
//Extend to BCAL + FCAL
/****************************************************** DEFINE LAMBDAS, USING STATEMENTS *******************************************************/

using DSourceCombosByUse_Large = unordered_map<DSourceComboUse, vector<const DSourceCombo*>*>; //for use collecting ALL of the combos in the event
using DSourceCombosByBeamBunch = unordered_map<int, const vector<const DSourceCombo*>*>; //int: RF bunch shift from primary
auto Compare_SourceComboInfos = [](const DSourceComboInfo* lhs, const DSourceComboInfo* rhs) -> bool{return *lhs < *rhs;};

struct DCompare_SourceComboInfos{
	bool operator()(const DSourceComboInfo* lhs, const DSourceComboInfo* rhs) const{return *lhs < *rhs;}
};

/*******
 * Call order:
 * Charged, primary vertex
 * Neutral, primary vertex
 * Charged, detached vertex
 * Neutral, detached vertex
 */

/************************************************************** DEFINE CLASSES ***************************************************************/

class DSourceComboer : public JObject
{
	public:

		DSourceComboer(void) = delete;
		DSourceComboer(JEventLoop* locEventLoop);
		DSourceComboer::~DSourceComboer(void);

		DSourceCombosByBeamBunch Request_SourceCombos(JEventLoop* locEventLoop, const DReactionStepVertexInfo* locReactionStepVertexInfo,
				DVector3 locVertex, const set<int>& locBeamBunchesToDo); //if set is empty, return for all possible bunches

	private:

		/********************************************************** DEFINE USING STATEMENTS ***********************************************************/

		//DEFINE USING STATEMENTS
		using DPhotonShowersByBeamBunch = unordered_map<int, vector<const DNeutralShower*>>; //int: beam bunch n-shifts from nominal
		using DSourceCombosByUseByBeamBunch = unordered_map<int, DSourceCombosByUse_Large>; //int: beam bunch n-shifts from nominal

		/********************************************************** DECLARE MEMBER FUNCTIONS ***********************************************************/

		void Define_LooseCuts(void);
		void Reset_NewEvent(JEventLoop* locEventLoop);
		void Setup_NeutralShowers(JEventLoop* locEventLoop);

		//Create photon combo infos
		void Create_SourceComboInfos(const DReactionVertexInfo* locReactionVertexInfo);
		void Create_SourceComboInfos_Vertices(const DReactionVertexInfo* locReactionVertexInfo);
		void Create_SourceComboInfos_Photons(const shared_ptr<const DReactionStepVertexInfo>& locReactionStepVertexInfo);
		pair<bool, map<DSourceComboUse, unsigned char>> Get_FinalStateDecayingComboUses(const DReaction* locReaction, size_t locStepIndex, const map<size_t, DSourceComboUse>& locStepComboUseMap) const;
		const DSourceComboInfo* MakeOrGet_SourceComboInfo(const vector<pair<Particle_t, unsigned char>>& locNumParticles, const vector<pair<DSourceComboUse, unsigned char>>& locFurtherDecays);
		const DSourceComboInfo* GetOrMake_SourceComboInfo(const vector<pair<Particle_t, unsigned char>>& locNumParticles, const vector<pair<DSourceComboUse, unsigned char>>& locFurtherDecays);

		//TIMING METHODS
		void Calc_PhotonBeamBunchShifts(const DNeutralShower* locNeutralShower, shared_ptr<const DKinematicData>& locKinematicData, double locRFTime,
				DPhotonShowersByBeamBunch& locShowersByBeamBunch) const;
		double Calc_MaxDeltaTError(const DNeutralShower* locNeutralShower, const shared_ptr<const DKinematicData>& locKinematicData) const;

		//CREATE COMBOS - GENERAL METHODS
		DSourceCombosByBeamBunch Create_SourceCombos(const DSourceComboUse& locSourceComboUse, const DPhotonShowersByBeamBunch& locShowersByBeamBunch,
				set<int> locBeamBunchesToDo, DSourceCombosByUseByBeamBunch& locSourceCombosByUseByBeamBunch);
		vector<const DSourceCombo*>* Create_SourceCombos(const DSourceComboUse& locSourceComboUse, const vector<const DNeutralShower*>& locShowers,
				DSourceCombosByUse_Large& locSourceCombosByUseSoFar);
		vector<const DSourceCombo*>* Create_SourceCombos(const DSourceComboInfo* locSourceComboInfo, const vector<const DNeutralShower*>& locShowers,
				DSourceCombosByUse_Large& locSourceCombosByUseSoFar);

		//COMBO VERTICALLY METHODS
		void Combo_Vertically_AllDecays(const DSourceComboUse& locComboUseToCreate, const vector<const DNeutralShower*>& locShowers, DSourceCombosByUse_Large& locSourceCombosByUseSoFar);
		void Combo_Vertically_NDecays(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locNMinus1ComboUse, const DSourceComboUse& locSourceComboDecayUse, DSourceCombosByUse_Large& locSourceCombosByUseSoFar);
		void Combo_Vertically_AllParticles(const DSourceComboUse& locComboUseToCreate, const vector<const DNeutralShower*>& locShowers, DSourceCombosByUse_Large& locSourceCombosByUseSoFar);
		void Combo_Vertically_NParticles(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locNMinus1ComboUse, const vector<const DNeutralShower*>& locShowers, DSourceCombosByUse_Large& locSourceCombosByUseSoFar);

		//COMBO HORIZONTALLY METHODS
		void Combo_Horizontally_All(const DSourceComboUse& locComboUseToCreate, const vector<const DNeutralShower*>& locShowers, DSourceCombosByUse_Large& locSourceCombosByUseSoFar);
		void Combo_Horizontally_AddCombo(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locAllBut1ComboUse, const DSourceComboUse& locSourceComboUseToAdd, DSourceCombosByUse_Large& locSourceCombosByUseSoFar);
		void Combo_Horizontally_AddParticle(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locAllBut1ComboUse, Particle_t locPID, const vector<const DNeutralShower*>& locShowers, DSourceCombosByUse_Large& locSourceCombosByUseSoFar);

		//UTILITY FUNCTIONS
		size_t Get_PhotonVertexZBin(double locVertexZ) const;
		double Get_PhotonVertexZBinCenter(size_t locVertexZBin) const;
		shared_ptr<const DKinematicData> Create_KinematicData(const DNeutralShower* locNeutralShower, const DVector3& locVertex) const;
		int Calc_RFBunchShift(double locTimeToStep, double locTimeToStepTo) const; //returns integer shift

		vector<const DSourceCombo*>* Cut_InvariantMass(const vector<const DSourceCombo*>* locSourceCombos, Particle_t locDecayPID);
		double Calc_InvariantMass(const DSourceCombo* locSourceCombo) const;

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

		//NEUTRAL SHOWER DATA
		vector<const DNeutralShower*> dFCALShowers;
		vector<const DNeutralShower*> dBCALShowers;
		unordered_map<const DNeutralShower*, shared_ptr<const DKinematicData>> dFCALKinematics; //FCAL shower data at center of target
		vector<unordered_map<const DNeutralShower*, shared_ptr<const DKinematicData>>> dBCALKinematics; //BCAL shower data in vertex-z bins

		//SHOWERS SORTED BY RF BUNCH
		const DEventRFBunch* dInitialEventRFBunch;
		DPhotonShowersByBeamBunch dFCALShowersByBeamBunch;
		vector<DPhotonShowersByBeamBunch> dBCALShowersByBeamBunch; //vector: vertex-z bins

		//CUTS
		unordered_map<DetectorSystem_t, TF1*> dPhotonTimeCutMap; //function of shower energy (p)
		unordered_map<Particle_t, pair<double, double>> dInvariantMassCuts;

		//PHOTON COMBO INFOS: CREATED ONCE DURING DSourceComboER OBJECT CONSTRUCTION
		//want to make sure we only have one of each type: suggests using a set
		//however, after the first few events, almost all of these have already been created: vector has faster lookup time
		//therefore, use the set when creating the objects during construction, but then move the results into the vector and keep it sorted
		set<const DSourceComboInfo*, DCompare_SourceComboInfos> dSourceComboInfoSet;
		vector<const DSourceComboInfo*> dSourceComboInfos;
		unordered_map<shared_ptr<const DReactionStepVertexInfo>, DSourceComboUse> dSourceComboUseReactionMap; //primary combo info (nullptr if none)
		unordered_map<pair<shared_ptr<const DReactionStepVertexInfo>, DSourceComboUse>, size_t> dSourceComboInfoStepMap; //size_t: step index

		//PHOTON COMBOS //vector: z-bin //if attempted and all failed, DSourceCombosByBeamBunch vector will be empty
		size_t dInitialComboVectorCapacity = 100;
		DSourceCombosByUse_Large dSourceCombos_Charged;
		DSourceCombosByUseByBeamBunch dSourceCombos_FCAL;
		vector<DSourceCombosByUseByBeamBunch> dSourceCombos_Both;

		//e.g. if a DSourceCombo is -> 2pi0, and we want to use it as a basis for building a combo of 3pi0s,
		//then this iterator points to the first pi0 in the DSourceCombosByUseByBeamBunch vector that we want to test
		//that way we save a lot of time, since we don't have to look for it again
		unordered_map<const DSourceCombo*, vector<const DSourceCombo*>::const_iterator> dResumeSearchIteratorMap_Combos;
		unordered_map<const DSourceCombo*, vector<const DNeutralShower*>::const_iterator> dResumeSearchIteratorMap_Showers;

};

/*********************************************************** INLINE MEMBER FUNCTION DEFINITIONS ************************************************************/

inline DSourceComboer::~DSourceComboer(void)
{
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

inline vector<const DSourceCombo*>* DSourceComboer::Cut_InvariantMass(const vector<const DSourceCombo*>* locSourceCombos, Particle_t locDecayPID)
{
	auto locCutIterator = dInvariantMassCuts.find(locDecayPID);
	if(locCutIterator == dInvariantMassCuts.end())
		return locSourceCombos; //no cut to place!!

	auto locCutSourceCombos = new vector<const DSourceCombo*>();
	locCutSourceCombos->reserve(locSourceCombos->size());
	auto& locMassCuts = locCutIterator->second;
	auto locMassCutter = [&locDecayPID, &dInvariantMassCuts](const DSourceCombo* locSourceCombo) -> bool
	{
		auto locInvariantMass = Calc_InvariantMass(locSourceCombo);
		return ((locInvariantMass >= locMassCuts.first) && (locInvariantMass <= locMassCuts.second));
	};
	std::copy_if(locSourceCombos->begin(), locSourceCombos->end(), std::back_inserter(*locCutSourceCombos), locMassCutter);
	return locCutSourceCombos;
}

inline double DSourceComboer::Calc_InvariantMass(const DSourceCombo* locSourceCombo) const
{
	return 0.0;
}

/*********************************************************** INLINE NAMESPACE-SCOPE FUNCTION DEFINITIONS ************************************************************/

} //end DAnalysis namespace

#endif // DSourceComboer_h
