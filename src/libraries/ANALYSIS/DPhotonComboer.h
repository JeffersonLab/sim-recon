#ifndef DPhotonComboer_h
#define DPhotonComboer_h

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
#include "ANALYSIS/DPhotonCombo.h"
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

using DPhotonCombosByUse_Large = unordered_map<DPhotonComboUse, vector<const DPhotonCombo*>*>; //for use collecting ALL of the combos in the event
using DPhotonCombosByBeamBunch = unordered_map<int, const vector<const DPhotonCombo*>*>; //int: RF bunch shift from primary
auto Compare_PhotonComboInfos = [](const DPhotonComboInfo* lhs, const DPhotonComboInfo* rhs) -> bool{return *lhs < *rhs;};

struct DCompare_PhotonComboInfos{
	bool operator()(const DPhotonComboInfo* lhs, const DPhotonComboInfo* rhs) const{return *lhs < *rhs;}
};

/************************************************************** DEFINE CLASSES ***************************************************************/

class DPhotonComboer : public JObject
{
	public:

		DPhotonComboer(void) = delete;
		DPhotonComboer(JEventLoop* locEventLoop);
		DPhotonComboer::~DPhotonComboer(void);

		DPhotonCombosByBeamBunch Request_PhotonCombos(JEventLoop* locEventLoop, const DReactionStepVertexInfo* locReactionStepVertexInfo,
				DVector3 locVertex, const set<int>& locBeamBunchesToDo); //if set is empty, return for all possible bunches

	private:

		/********************************************************** DEFINE USING STATEMENTS ***********************************************************/

		//DEFINE USING STATEMENTS
		using DPhotonShowersByBeamBunch = unordered_map<int, vector<const DNeutralShower*>>; //int: beam bunch n-shifts from nominal
		using DPhotonCombosByUseByBeamBunch = unordered_map<int, DPhotonCombosByUse_Large>; //int: beam bunch n-shifts from nominal

		/********************************************************** DECLARE MEMBER FUNCTIONS ***********************************************************/

		void Reset_NewEvent(JEventLoop* locEventLoop);
		void Setup_NeutralShowers(JEventLoop* locEventLoop);

		//Create photon combo infos
		void Create_PhotonComboInfos(const DReactionVertexInfo* locReactionVertexInfo);
		map<size_t, DPhotonComboUse> Create_PhotonComboInfos(const shared_ptr<const DReactionStepVertexInfo>& locReactionStepVertexInfo);
		const DPhotonComboInfo* MakeOrGet_PhotonComboInfo(size_t locNumPhotons, const vector<pair<DPhotonComboUse, unsigned char>>& locFurtherDecays);
		const DPhotonComboInfo* GetOrMake_PhotonComboInfo(size_t locNumPhotons, const vector<pair<DPhotonComboUse, unsigned char>>& locFurtherDecays);

		//TIMING METHODS
		void Calc_PhotonBeamBunchShifts(const DNeutralShower* locNeutralShower, shared_ptr<const DKinematicData>& locKinematicData, double locRFTime,
				DPhotonShowersByBeamBunch& locShowersByBeamBunch) const;
		double Calc_MaxDeltaTError(const DNeutralShower* locNeutralShower, const shared_ptr<const DKinematicData>& locKinematicData) const;

		//CREATE COMBOS METHODS
		DPhotonCombosByBeamBunch Create_PhotonCombos(const DPhotonComboUse& locPhotonComboUse, const DPhotonShowersByBeamBunch& locShowersByBeamBunch,
				set<int> locBeamBunchesToDo, DPhotonCombosByUseByBeamBunch& locPhotonCombosByUseByBeamBunch);
		vector<const DPhotonCombo*>* Create_PhotonCombos(const DPhotonComboUse& locPhotonComboUse, const vector<const DNeutralShower*>& locShowers,
				DPhotonCombosByUse_Large& locPhotonCombosByUseSoFar);
		vector<const DPhotonCombo*>* Create_PhotonCombos(const DPhotonComboInfo* locPhotonComboInfo, const vector<const DNeutralShower*>& locShowers,
				DPhotonCombosByUse_Large& locPhotonCombosByUseSoFar);

		void Combo_Vertically(const DPhotonComboUse& locNeededGroupingUse, const DPhotonComboUse& locNMinus1ComboUse, const DPhotonComboUse& locPhotonComboDecayUse, DPhotonCombosByUse_Large& locPhotonCombosByUseSoFar);
		void Combo_Horizontally(const DPhotonComboUse& locNeededGroupingUse, const DPhotonComboUse& locAllBut1ComboUse, const DPhotonComboUse& locPhotonComboDecayUse, DPhotonCombosByUse_Large& locPhotonCombosByUseSoFar);

		//UTILITY FUNCTIONS
		size_t Get_PhotonVertexZBin(double locVertexZ) const;
		double Get_PhotonVertexZBinCenter(size_t locVertexZBin) const;
		shared_ptr<const DKinematicData> Create_KinematicData(const DNeutralShower* locNeutralShower, const DVector3& locVertex) const;
		int Calc_RFBunchShift(double locTimeToStep, double locTimeToStepTo) const; //returns integer shift

		vector<const DPhotonCombo*>* Cut_InvariantMass(const vector<const DPhotonCombo*>* locPhotonCombos, Particle_t locDecayPID);
		double Calc_InvariantMass(const DPhotonCombo* locPhotonCombo) const;

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

		//PHOTON COMBO INFOS: CREATED ONCE DURING DPHOTONCOMBOER OBJECT CONSTRUCTION
		//want to make sure we only have one of each type: suggests using a set
		//however, after the first few events, almost all of these have already been created: vector has faster lookup time
		//therefore, use the set when creating the objects during construction, but then move the results into the vector and keep it sorted
		set<const DPhotonComboInfo*, DCompare_PhotonComboInfos> dPhotonComboInfoSet;
		vector<const DPhotonComboInfo*> dPhotonComboInfos;
		unordered_map<shared_ptr<const DReactionStepVertexInfo>, DPhotonComboUse> dPhotonComboUseReactionMap; //primary combo info (nullptr if none)
		unordered_map<pair<shared_ptr<const DReactionStepVertexInfo>, DPhotonComboUse>, size_t> dPhotonComboInfoStepMap; //size_t: step index

		//PHOTON COMBOS //vector: z-bin //if attempted and all failed, DPhotonCombosByBeamBunch vector will be empty
		size_t dInitialComboVectorCapacity = 100;
		DPhotonCombosByUseByBeamBunch dPhotonCombos_FCAL;
		vector<DPhotonCombosByUseByBeamBunch> dPhotonCombos_BCAL;
		vector<DPhotonCombosByUseByBeamBunch> dPhotonCombos_Both;

		//e.g. if a DPhotonCombo is -> 2pi0, and we want to use it as a basis for building a combo of 3pi0s,
		//then this iterator points to the first pi0 in the DPhotonCombosByUseByBeamBunch vector that we want to test
		//that way we save a lot of time, since we don't have to look for it again
		unordered_map<const DPhotonCombo*, vector<const DPhotonCombo*>::const_iterator> dResumeSearchIteratorMap;
};

/*********************************************************** INLINE MEMBER FUNCTION DEFINITIONS ************************************************************/

inline DPhotonComboer::~DPhotonComboer(void)
{
	for(auto locComboInfo : dPhotonComboInfos)
		delete locComboInfo;
}

inline size_t DPhotonComboer::Get_PhotonVertexZBin(double locVertexZ) const
{
	//given some vertex-z, what bin am I in?
	int locPhotonVertexZBin = int((locVertexZ - dPhotonVertexZRangeLow)/dPhotonVertexZBinWidth);
	if(locPhotonVertexZBin < 0)
		return 0;
	else if(locPhotonVertexZBin >= dNumPhotonVertexZBins)
		return dNumPhotonVertexZBins - 1;
	return locPhotonVertexZBin;
}

inline double DPhotonComboer::Get_PhotonVertexZBinCenter(size_t locVertexZBin) const
{
	return dPhotonVertexZRangeLow + (double(locVertexZBin) + 0.5)*dPhotonVertexZBinWidth;
}

inline int DPhotonComboer::Calc_RFBunchShift(double locTimeToStep, double locTimeToStepTo) const
{
	double locDeltaT = locTimeToStepTo - locTimeToStep;
	return (locDeltaT > 0.0) ? int(locDeltaT/dBeamBunchPeriod + 0.5) : int(locDeltaT/dBeamBunchPeriod - 0.5);
}

inline shared_ptr<const DKinematicData> DPhotonComboer::Create_KinematicData(const DNeutralShower* locNeutralShower, const DVector3& locVertex) const
{
	DVector3 locPath = locNeutralShower->dSpacetimeVertex.Vect() - locVertex;
	double locPathLength = locPath.Mag();
	double locVertexTime = locNeutralShower->dSpacetimeVertex.T() - locPathLength/29.9792458;
	DVector3 locMomentum = locNeutralShower->dEnergy*locPath.Unit();
	return std::make_shared<const DKinematicData>(Gamma, locMomentum, locVertex, locVertexTime);
}

inline vector<const DPhotonCombo*>* DPhotonComboer::Cut_InvariantMass(const vector<const DPhotonCombo*>* locPhotonCombos, Particle_t locDecayPID)
{
	auto locCutIterator = dInvariantMassCuts.find(locDecayPID);
	if(locCutIterator == dInvariantMassCuts.end())
		return locPhotonCombos; //no cut to place!!

	auto locCutPhotonCombos = new vector<const DPhotonCombo*>();
	locCutPhotonCombos->reserve(locPhotonCombos->size());
	auto& locMassCuts = locCutIterator->second;
	auto locMassCutter = [&locDecayPID, &dInvariantMassCuts](const DPhotonCombo* locPhotonCombo) -> bool
	{
		auto locInvariantMass = Calc_InvariantMass(locPhotonCombo);
		return ((locInvariantMass >= locMassCuts.first) && (locInvariantMass <= locMassCuts.second));
	};
	std::copy_if(locPhotonCombos->begin(), locPhotonCombos->end(), std::back_inserter(*locCutPhotonCombos), locMassCutter);
	return locCutPhotonCombos;
}

inline double DPhotonComboer::Calc_InvariantMass(const DPhotonCombo* locPhotonCombo) const
{
	return 0.0;
}

/*********************************************************** INLINE NAMESPACE-SCOPE FUNCTION DEFINITIONS ************************************************************/

} //end DAnalysis namespace

#endif // DPhotonComboer_h
