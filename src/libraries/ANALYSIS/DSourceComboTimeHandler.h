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
		DSourceComboTimeHandler(JEventLoop* locEventLoop, DSourceComboer* locSourceComboer, const DSourceComboVertexer* locSourceComboVertexer);

		//SETUP
		void Reset(void);
		void Setup_NeutralShowers(const vector<const DNeutralShower*>& locNeutralShowers, const DEventRFBunch* locInitialEventRFBunch);
		void Set_BeamParticles(const vector<const DBeamPhoton*>& locBeamParticles);

		//GET SETUP RESULTS
		DPhotonKinematicsByZBin Get_PhotonKinematics(void) const{return dPhotonKinematics;}
		unordered_map<signed char, DPhotonShowersByBeamBunch> Get_ShowersByBeamBunchByZBin(void) const{return dShowersByBeamBunchByZBin;}
		vector<const DKinematicData*> Get_BeamParticlesByRFBunch(int locRFBunch) const;

		//GET VALID/COMMON RF BUNCHES
		//Note that combining with an empty vector does NOT remove all bunches!!
		//It is assumed that an empty vector means "unknown," and that you just want to use the other set instead
		vector<int> Get_ValidRFBunches(const JObject* locObject, signed char locVertexZBin) const{return dShowerRFBunches[locVertexZBin][locObject];}
		vector<int> Get_CommonRFBunches(const vector<int>& locRFBunches1, const JObject* locObject, signed char locVertexZBin) const;
		vector<int> Get_CommonRFBunches(const vector<int>& locRFBunches1, const vector<int>& locRFBunches2) const;

		//UTILITY FUNCTIONS
		int Calc_RFBunchShift(double locTimeToStepTo) const{return Calc_RFBunchShift(dInitialEventRFBunch->dTime, locTimeToStepTo);}
		int Calc_RFBunchShift(double locTimeToStep, double locTimeToStepTo) const;
		double Calc_PropagatedRFTime(double locPrimaryVertexZ, int locNumRFBunchShifts, double locDetachedVertexTimeOffset) const;

		//SELECT RF BUNCHES
		bool Select_RFBunches_Charged(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locChargedCombo, vector<int>& locValidRFBunches);
		bool Select_RFBunches_PhotonVertices(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionChargedCombo, const DSourceCombo* locReactionFullCombo, vector<int>& locValidRFBunches);
		int Select_RFBunch_Full(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, const DSourceCombo* locReactionChargedCombo, const vector<int>& locRFBunches);
		bool Cut_Timing_MissingMassVertices(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionChargedCombo, const DSourceCombo* locReactionFullCombo, int locRFBunch);

	private:

		shared_ptr<const DKinematicData> Create_KinematicData_Photon(const DNeutralShower* locNeutralShower, const DVector3& locVertex) const;

		void Calc_PhotonBeamBunchShifts(const DNeutralShower* locNeutralShower, shared_ptr<const DKinematicData>& locKinematicData, double locRFTime, signed char locZBin);
		vector<int> Calc_BeamBunchShifts(double locVertexTime, double locRFTime, double locDeltaTCut, bool locIncludeDecayTimeOffset) const;

		double Calc_MaxDeltaTError(const DNeutralShower* locNeutralShower, const shared_ptr<const DKinematicData>& locKinematicData) const;

		vector<int> Get_RFBunches_ChargedTrack(const DChargedTrack* locChargedTrack, Particle_t locPID, bool locIsProductionVertex, const DSourceCombo* locVertexPrimaryCombo, DVector3 locVertex, double locTimeOffset, double locPropagatedRFTime);

		double Calc_RFDeltaTChiSq(const DNeutralShower* locNeutralShower, int locNumShifts, const TVector3& locVertex, double locPropagatedRFTime) const;
		double Calc_RFDeltaTChiSq(const DChargedTrackHypothesis* locHypothesis, int locNumShifts, double locVertexTime, double locPropagatedRFTime) const;

		//HANDLERS AND UTILITIES
		DSourceComboer* dSourceComboer;
		const DSourceComboVertexer* dSourceComboVertexer;
		const DAnalysisUtilities* dAnalysisUtilities;

		//EXPERIMENT INFORMATION
		DVector3 dTargetCenter;
		double dTargetLength = 30.0;
		double dBeamBunchPeriod = 1000.0/249.5;
		bool dUseSigmaForRFSelectionFlag = false;

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

		unordered_map<const DSourceCombo*, vector<int>> dChargedComboRFBunches; //empty vector = FAILED cuts //combo: charged
		unordered_map<const DSourceCombo*, vector<int>> dPhotonVertexRFBunches; //combo: full
		unordered_map<const DSourceCombo*, int> dFullComboRFBunches; //combo: full
		unordered_map<const DSourceCombo*, bool> dFullComboTimeCutResults; //combo: full

		unordered_map<int, vector<const DKinematicData*>> dBeamParticlesByRFBunch;

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

inline void DSourceComboTimeHandler::Set_BeamParticles(const vector<const DBeamPhoton*>& locBeamParticles)
{
	//ASSUME that the tagger timing is good enough to determine which bunch it belongs too!
	for(auto locBeamParticle : locBeamParticles)
	{
		auto locRFBunch = Calc_RFBunchShift(dInitialEventRFBunch->dTime, locBeamParticle->time());
		dBeamParticlesByRFBunch[locRFBunch].push_back(locBeamParticle);
	}
}

inline vector<const DKinematicData*> DSourceComboTimeHandler::Get_BeamParticlesByRFBunch(int locRFBunch) const
{
	auto locIterator = dBeamParticlesByRFBunch.find(locRFBunch);
	if(locIterator != dBeamParticlesByRFBunch.end())
		return locIterator->second;
	return {};
}

inline double DSourceComboTimeHandler::Calc_PropagatedRFTime(double locPrimaryVertexZ, int locNumRFBunchShifts, double locDetachedVertexTimeOffset) const
{
	//propagate rf time to vertex and add time offset (faster to just do it here rather than for each particle)
	return dInitialEventRFBunch->dTime + (locPrimaryVertexZ - dTargetCenter.Z())/SPEED_OF_LIGHT + locNumRFBunchShifts*dBeamBunchPeriod + locDetachedVertexTimeOffset;
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

inline double DSourceComboTimeHandler::Calc_RFDeltaTChiSq(const DNeutralShower* locNeutralShower, int locNumShifts, const TVector3& locVertex, double locPropagatedRFTime) const
{
	//calc vertex time, get delta-t cut
	double locPathLength = (locNeutralShower->dSpacetimeVertex.Vect() - locVertex).Mag();
	double locVertexTime = locNeutralShower->dSpacetimeVertex.T() - locPathLength/29.9792458 - locTimeOffset;
	double locVertexTimeVariance = dUseSigmaForRFSelectionFlag ? locNeutralShower->dCovarianceMatrix(4, 4) : 1.0;

	double locDeltaT = locVertexTime - locPropagatedRFTime;
	return locDeltaT*locDeltaT/locVertexTimeVariance;
}

inline double DSourceComboTimeHandler::Calc_RFDeltaTChiSq(const DChargedTrackHypothesis* locHypothesis, int locNumShifts, double locVertexTime, double locPropagatedRFTime) const
{
	//calc vertex time, get delta-t cut
	double locVertexTimeVariance = dUseSigmaForRFSelectionFlag ? locHypothesis->errorMatrix(6, 6) : 1.0;
	double locDeltaT = locVertexTime - locPropagatedRFTime;
	return locDeltaT*locDeltaT/locVertexTimeVariance;
}

}

#endif // DSourceComboTimeHandler_h
