#ifndef DSourceComboTimeHandler_h
#define DSourceComboTimeHandler_h

#include <unordered_map>
#include <vector>

#include "particleType.h"
#include "DLorentzVector.h"
#include "PID/DKinematicData.h"
#include "PID/DNeutralShower.h"
#include "PID/DEventRFBunch.h"
#include "ANALYSIS/DSourceCombo.h"
#include "ANALYSIS/DSourceComboVertexer.h"

using namespace std;
using namespace jana;

namespace DAnalysis
{

using DPhotonKinematicsByZBin = unordered_map<signed char, unordered_map<const DNeutralShower*, shared_ptr<const DKinematicData>>>; //char: z-bin
using DPhotonShowersByBeamBunch = unordered_map<vector<int>, vector<const JObject*>>; //int: beam bunch n-shifts from nominal

class DSourceComboTimeHandler
{
	public:
		DSourceComboTimeHandler(void) = delete;
		DSourceComboTimeHandler(JEventLoop* locEventLoop, const DSourceComboer* locSourceComboer, const DSourceComboVertexer* locSourceComboVertexer);

		//SETUP
		void Reset(void);
		void Setup_NeutralShowers(const vector<const DNeutralShower*>& locNeutralShowers, const DEventRFBunch* locInitialEventRFBunch);

		//GET SETUP RESULTS
		DPhotonKinematicsByZBin Get_PhotonKinematics(void) const{return dPhotonKinematics;}
		unordered_map<signed char, DPhotonShowersByBeamBunch> Get_ShowersByBeamBunchByZBin(void) const{return dShowersByBeamBunchByZBin;}

		//GET VALID/COMMON RF BUNCHES
		//Note that combining with an empty vector does NOT remove all bunches!!
		//It is assumed that an empty vector means "unknown," and that you just want to use the other set instead
		vector<int> Get_ValidRFBunches(const JObject* locObject, signed char locVertexZBin) const{return dShowerRFBunches[locVertexZBin][locObject];}
		vector<int> Get_CommonRFBunches(const vector<int>& locRFBunches1, const JObject* locObject, signed char locVertexZBin) const;
		vector<int> Get_CommonRFBunches(const vector<int>& locRFBunches1, const vector<int>& locRFBunches2) const;

		//UTILITY FUNCTIONS
		int Calc_RFBunchShift(double locTimeToStep, double locTimeToStepTo) const;

		//SELECT RF BUNCHES
		vector<int> Select_RFBunches_Charged(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locChargedCombo) const;

	private:

		shared_ptr<const DKinematicData> Create_KinematicData_Photon(const DNeutralShower* locNeutralShower, const DVector3& locVertex) const;

		void Calc_PhotonBeamBunchShifts(const DNeutralShower* locNeutralShower, shared_ptr<const DKinematicData>& locKinematicData, double locRFTime, signed char locZBin);
		vector<int> Calc_BeamBunchShifts(double locVertexTime, double locRFTime, double locDeltaTCut, bool locIncludeDecayTimeOffset) const;

		double Calc_MaxDeltaTError(const DNeutralShower* locNeutralShower, const shared_ptr<const DKinematicData>& locKinematicData) const;

		//HANDLERS
		const DSourceComboer* dSourceComboer;
		const DSourceComboVertexer* dSourceComboVertexer;
		const DAnalysisUtilities* dAnalysisUtilities;

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

		//SHOWERS SORTED BY RF BUNCH
		const DEventRFBunch* dInitialEventRFBunch = nullptr;

		//zbins: showers are stored based on what zbins they are VALID for
		//that means that the FCAL showers will appear in EVERY zbin: duplicate entries
		//even though this takes less memory, it's faster than merging vectors each time they are requested
		//Note that the RF-shift vectors are sorted
		DPhotonKinematicsByZBin dPhotonKinematics;
		//shower -> RFs
		unordered_map<signed char, unordered_map<const JObject*, vector<int>>> dShowerRFBunches; //char: z-bin
		//RF -> showers
		unordered_map<signed char, DPhotonShowersByBeamBunch> dShowersByBeamBunchByZBin; //char: zbin

		unordered_map<pair<const DKinematicData*, vector<const DKinematicData*>>, DLorentzVector> dChargedParticlePOCAToVertexX4;

		unordered_map<const DSourceCombo*, vector<int>> dChargedComboRFBunches; //empty vector = FAILED cuts

		//CUTS
		unordered_map<Particle_t, unordered_map<DetectorSystem_t, TF1*>> dPIDTimingCuts; //function of shower energy (p)
};

inline void DSourceComboTimeHandler::Reset(void)
{
	dInitialEventRFBunch = nullptr;
	dPhotonKinematics.clear(); //can probably delete instead of shared_ptr!!
	dShowerRFBunches.clear();
	dChargedComboRFBunches.clear();
	dShowersByBeamBunchByZBin.clear();
	dChargedParticlePOCAToVertexX4.clear();
}

inline shared_ptr<const DKinematicData> DSourceComboTimeHandler::Create_KinematicData_Photon(const DNeutralShower* locNeutralShower, const DVector3& locVertex) const
{
	DVector3 locPath = locNeutralShower->dSpacetimeVertex.Vect() - locVertex;
	double locPathLength = locPath.Mag();
	double locVertexTime = locNeutralShower->dSpacetimeVertex.T() - locPathLength/29.9792458;
	DVector3 locMomentum = locNeutralShower->dEnergy*locPath.Unit();
	return std::make_shared<const DKinematicData>(Gamma, locMomentum, locVertex, locVertexTime);
}

inline vector<int> DSourceComboTimeHandler::Get_ValidRFBunches(const JObject* locObject, signed char locVertexZBin) const
{
	const auto& locBunchesByObject = dShowerRFBunches.find(locVertexZBin)->second;
	auto locIterator = locBunchesByObject.find(locObject);
	if(locIterator == locBunchesByObject.end())
		return {};
	return locIterator->second;
}

inline int DSourceComboTimeHandler::Calc_RFBunchShift(double locTimeToStep, double locTimeToStepTo) const
{
	double locDeltaT = locTimeToStepTo - locTimeToStep;
	return (locDeltaT > 0.0) ? int(locDeltaT/dBeamBunchPeriod + 0.5) : int(locDeltaT/dBeamBunchPeriod - 0.5);
}

inline vector<int> DSourceComboTimeHandler::Get_CommonRFBunches(const vector<int>& locRFBunches1, const JObject* locObject, signed char locVertexZBin) const
{
	return Get_CommonRFBunches(locRFBunches1, Get_ValidRFBunches(locObject, locVertexZBin));
}

inline vector<int> DSourceComboTimeHandler::Get_CommonRFBunches(const vector<int>& locRFBunches1, const vector<int>& locRFBunches2) const
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


}

#endif // DSourceComboTimeHandler_h
