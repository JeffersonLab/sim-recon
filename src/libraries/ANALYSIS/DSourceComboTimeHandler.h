#ifndef DSourceComboTimeHandler_h
#define DSourceComboTimeHandler_h

#include <unordered_map>
#include <vector>
#include <map>
#include <memory>

#include "TF1.h"
#include "TH2I.h"
#include "TROOT.h"
#include "TFile.h"
#include "TDirectoryFile.h"

#include "particleType.h"
#include "DLorentzVector.h"
#include "PID/DKinematicData.h"
#include "PID/DNeutralShower.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DEventRFBunch.h"
#include "PID/DDetectorMatches.h"
#include "ANALYSIS/DSourceCombo.h"
#include "ANALYSIS/DReactionVertexInfo.h"
#include "ANALYSIS/DAnalysisUtilities.h"

using namespace std;
using namespace jana;

namespace DAnalysis
{

using DPhotonKinematicsByZBin = unordered_map<signed char, unordered_map<const DNeutralShower*, shared_ptr<const DKinematicData>>>; //char: z-bin
using DPhotonShowersByBeamBunch = map<vector<int>, vector<const JObject*>>; //int: beam bunch n-shifts from nominal
class DSourceComboer;
class DSourceComboVertexer;

class DSourceComboTimeHandler
{
	public:
		DSourceComboTimeHandler(void) = delete;
		DSourceComboTimeHandler(JEventLoop* locEventLoop, DSourceComboer* locSourceComboer, const DSourceComboVertexer* locSourceComboVertexer);
		~DSourceComboTimeHandler(void){Fill_Histograms();}

		//SETUP
		void Reset(void);
		void Set_DebugLevel(int locDebugLevel){dDebugLevel = locDebugLevel;}
		void Setup(const vector<const DNeutralShower*>& locNeutralShowers, const DEventRFBunch* locInitialEventRFBunch, const DDetectorMatches* locDetectorMatches);
		void Set_BeamParticles(const vector<const DBeamPhoton*>& locBeamParticles);

		//GET SETUP RESULTS
		const DEventRFBunch* Get_InitialEventRFBunch(void) const{return dInitialEventRFBunch;}
		DPhotonKinematicsByZBin Get_PhotonKinematics(void) const{return dPhotonKinematics;}
		unordered_map<signed char, DPhotonShowersByBeamBunch> Get_ShowersByBeamBunchByZBin(void) const{return dShowersByBeamBunchByZBin;}
		vector<const DKinematicData*> Get_BeamParticlesByRFBunch(int locRFBunch, unsigned int locPlusMinusBunchRange) const;

		//GET CUTS
		bool Get_TimeCut(Particle_t locPID, DetectorSystem_t locSystem, TF1* locTimeCut_ns) const;

		//GET VALID/COMMON RF BUNCHES
		//Note that combining with an empty vector does NOT remove all bunches!!
		//It is assumed that an empty vector means "unknown," and that you just want to use the other set instead
		vector<int> Get_ValidRFBunches(const JObject* locObject, signed char locVertexZBin) const;
		vector<int> Get_CommonRFBunches(const vector<int>& locRFBunches1, const JObject* locObject, signed char locVertexZBin) const;
		vector<int> Get_CommonRFBunches(const vector<int>& locRFBunches1, const vector<int>& locRFBunches2) const;

		//UTILITY FUNCTIONS
		int Calc_RFBunchShift(double locTimeToStepTo) const{return Calc_RFBunchShift(dInitialEventRFBunch->dTime, locTimeToStepTo);}
		int Calc_RFBunchShift(double locTimeToStep, double locTimeToStepTo) const;
		double Calc_RFTime(int locNumRFBunchShifts) const;
		double Calc_PropagatedRFTime(double locPrimaryVertexZ, int locNumRFBunchShifts, double locDetachedVertexTimeOffset) const;
		double Get_BeamBunchPeriod(void) const{return dBeamBunchPeriod;}

		//SELECT RF BUNCHES
		bool Select_RFBunches_Charged(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locChargedCombo, vector<int>& locValidRFBunches);
		bool Select_RFBunches_PhotonVertices(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, vector<int>& locValidRFBunches);
		bool Select_RFBunches_AllVerticesUnknown(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, Charge_t locCharge, vector<int>& locValidRFBunches);
		int Select_RFBunch_Full(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, const vector<int>& locRFBunches);
		bool Cut_Timing_MissingMassVertices(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, const DKinematicData* locBeamParticle, int locRFBunch);

		//GET POCA TO VERTEX
		DLorentzVector Get_ChargedPOCAToVertexX4(const DChargedTrackHypothesis* locHypothesis, bool locIsProductionVertex, const DSourceCombo* locReactionFullCombo, const DSourceCombo* locVertexPrimaryCombo, const DKinematicData* locBeamPhoton, DVector3 locVertex);

		//VERTEX-Z BINNING UTILITY FUNCTIONS
		size_t Get_NumVertexZBins(void) const{return dNumPhotonVertexZBins;}
		signed char Get_PhotonVertexZBin(double locVertexZ) const;
		double Get_PhotonVertexZBinCenter(signed char locVertexZBin) const;
		size_t Get_VertexZBin_TargetCenter(void) const{return Get_PhotonVertexZBin(dTargetCenter.Z());}

		void Vote_OldMethod(const DSourceCombo* locReactionFullCombo, vector<int>& locValidRFBunches);

	private:

		pair<DVector3, double> Calc_Photon_Kinematics(const DNeutralShower* locNeutralShower, const DVector3& locVertex) const;
		shared_ptr<const DKinematicData> Create_KinematicData_Photon(const DNeutralShower* locNeutralShower, const DVector3& locVertex) const;

		void Calc_PhotonBeamBunchShifts(const DNeutralShower* locNeutralShower, shared_ptr<const DKinematicData>& locKinematicData, double locRFTime, signed char locZBin);
		vector<int> Calc_BeamBunchShifts(double locVertexTime, double locOrigRFBunchPropagatedTime, double locDeltaTCut, bool locIncludeDecayTimeOffset, Particle_t locPID, DetectorSystem_t locSystem, double locP);

		double Calc_MaxDeltaTError(const DNeutralShower* locNeutralShower, const shared_ptr<const DKinematicData>& locKinematicData) const;

		bool Get_RFBunches_ChargedTrack(const DChargedTrackHypothesis* locHypothesis, bool locIsProductionVertex, const DSourceCombo* locVertexPrimaryCombo, DVector3 locVertex, double locPropagatedRFTime, bool locOnlyTrackFlag, vector<int>& locRFBunches);
		TF1* Get_TimeCutFunction(Particle_t locPID, DetectorSystem_t locSystem) const;

		bool Compute_RFChiSqs_UnknownVertices(const DSourceCombo* locReactionFullCombo, Charge_t locCharge, const vector<int>& locRFBunches, unordered_map<int, double>& locChiSqByRFBunch, map<int, map<Particle_t, map<DetectorSystem_t, vector<pair<float, float>>>>>& locRFDeltaTsForHisting);
		bool Cut_PhotonPID(const DNeutralShower* locNeutralShower, const DVector3& locVertex, double locPropagatedRFTime, bool locTargetCenterFlag);
		bool Cut_TrackPID(const DChargedTrackHypothesis* locHypothesis, bool locIsProductionVertex, const DSourceCombo* locFullReactionCombo, const DSourceCombo* locVertexPrimaryCombo, const DKinematicData* locBeamPhoton, DVector3 locVertex, double locPropagatedRFTime);

		pair<double, double> Calc_RFDeltaTChiSq(const DNeutralShower* locNeutralShower, const TVector3& locVertex, double locPropagatedRFTime) const;
		pair<double, double> Calc_RFDeltaTChiSq(const DChargedTrackHypothesis* locHypothesis, double locVertexTime, double locPropagatedRFTime) const;

		void Fill_Histograms(void);


		//HANDLERS AND UTILITIES
		DSourceComboer* dSourceComboer;
		const DSourceComboVertexer* dSourceComboVertexer;
		const DAnalysisUtilities* dAnalysisUtilities;
		int dDebugLevel = 0;

		//EXPERIMENT INFORMATION
		DVector3 dTargetCenter;
		double dTargetLength = 30.0;
		double dBeamBunchPeriod = 1000.0/249.5;
		bool dUseSigmaForRFSelectionFlag = false;

		//VERTEX-DEPENDENT PHOTON INFORMATION
		//For every 10cm in vertex-z, calculate the photon p4 & time for placing mass & delta-t cuts
		//The z-range extends from the upstream end of the target - 5cm to the downstream end + 25cm
		//so for a 30-cm-long target, it's a range of 50cm: 5bins, evaluated at the center of each bin
		//Make sure that the center of the target is the center of a zbin!!!
		float dPhotonVertexZBinWidth = 10.0;
		float dPhotonVertexZRangeLow = 45.0;
		size_t dNumPhotonVertexZBins = 6;
		//due to detached vertices
		double dMaxDecayTimeOffset = 0.0; //changed in constructor if needed
		double dDecayTimeUncertainty = 0.0; //due to uncertain path length

		//SHOWERS SORTED BY RF BUNCH
		const DEventRFBunch* dInitialEventRFBunch = nullptr;
		const DDetectorMatches* dDetectorMatches = nullptr;

		//zbins: showers are stored based on what zbins they are VALID for
		//that means that the FCAL showers will appear in EVERY zbin: duplicate entries
		//even though this takes less memory, it's faster than merging vectors each time they are requested
		//Note that the RF-shift vectors are sorted
		DPhotonKinematicsByZBin dPhotonKinematics;
		//shower -> RFs
		unordered_map<signed char, unordered_map<const JObject*, vector<int>>> dShowerRFBunches; //char: z-bin
		//RF -> showers
		unordered_map<signed char, DPhotonShowersByBeamBunch> dShowersByBeamBunchByZBin; //char: zbin

		map<pair<const DKinematicData*, vector<const DKinematicData*>>, DLorentzVector> dChargedParticlePOCAToVertexX4; //pair: charged hypo, then vertex particles

		unordered_map<const DSourceCombo*, pair<bool, vector<int>>> dChargedComboRFBunches; //bool: passed/failed cuts (can pass with empty vector if no timing info) //combo: charged
		unordered_map<const DSourceCombo*, pair<bool, vector<int>>> dPhotonVertexRFBunches; //bool: passed/failed cuts (can pass with empty vector if no timing info) //combo: full
		unordered_map<const DSourceCombo*, int> dFullComboRFBunches; //combo: full

		unordered_map<int, vector<const DKinematicData*>> dBeamParticlesByRFBunch;

		//CUTS
		//Unknown: initial RF selection for photons (at beginning of event, prior to vertex) //can be separate cut function
		map<Particle_t, map<DetectorSystem_t, TF1*>> dPIDTimingCuts; //function of p

		//HISTOGRAMS
		//Unknown: initial RF selection for photons (at beginning of event, prior to vertex) //can be separate cut function
		TH2* dHist_BeamRFDeltaTVsBeamE;
		map<Particle_t, map<DetectorSystem_t, TH2*>> dHistMap_RFDeltaTVsP_BestRF; //PID Unknown: photons prior to vertex selection
		map<Particle_t, map<DetectorSystem_t, TH2*>> dHistMap_RFDeltaTVsP_AllRFs; //PID Unknown: photons prior to vertex selection
		vector<pair<float, float>> dBeamRFDeltaTs;
		map<Particle_t, map<DetectorSystem_t, vector<pair<float, float>>>> dSelectedRFDeltaTs; //first float is p, 2nd is delta-t //PID Unknown: photons prior to vertex selection
		map<Particle_t, map<DetectorSystem_t, vector<pair<float, float>>>> dAllRFDeltaTs; //first float is p, 2nd is delta-t //PID Unknown: photons prior to vertex selection
};

inline void DSourceComboTimeHandler::Reset(void)
{
	Fill_Histograms();

	dInitialEventRFBunch = nullptr;

	dPhotonKinematics.clear(); //can probably delete instead of shared_ptr!!
	for(auto& locZPair : dShowerRFBunches)
		locZPair.second.clear();
	for(auto& locZPair : dShowersByBeamBunchByZBin)
		locZPair.second.clear();
	dChargedParticlePOCAToVertexX4.clear();
	dBeamParticlesByRFBunch.clear();

	dChargedComboRFBunches.clear();
	dPhotonVertexRFBunches.clear();
	dFullComboRFBunches.clear();
}

inline signed char DSourceComboTimeHandler::Get_PhotonVertexZBin(double locVertexZ) const
{
	//given some vertex-z, what bin am I in?
	if(locVertexZ < dPhotonVertexZRangeLow)
		return DSourceComboInfo::Get_VertexZIndex_OutOfRange();
	int locPhotonVertexZBin = int((locVertexZ - dPhotonVertexZRangeLow)/dPhotonVertexZBinWidth);
	return ((locPhotonVertexZBin >= int(dNumPhotonVertexZBins)) ? DSourceComboInfo::Get_VertexZIndex_OutOfRange() : locPhotonVertexZBin);
}

inline double DSourceComboTimeHandler::Get_PhotonVertexZBinCenter(signed char locVertexZBin) const
{
	return dPhotonVertexZRangeLow + (double(locVertexZBin) + 0.5)*dPhotonVertexZBinWidth;
}

inline void DSourceComboTimeHandler::Set_BeamParticles(const vector<const DBeamPhoton*>& locBeamParticles)
{
	for(const auto& locBeamParticle : locBeamParticles)
	{
		auto locRFBunch = Calc_RFBunchShift(dInitialEventRFBunch->dTime, locBeamParticle->time());
		dBeamParticlesByRFBunch[locRFBunch].push_back(locBeamParticle);
		dBeamRFDeltaTs.emplace_back(locBeamParticle->energy(), locBeamParticle->time() - dInitialEventRFBunch->dTime);
	}
}

inline bool DSourceComboTimeHandler::Get_TimeCut(Particle_t locPID, DetectorSystem_t locSystem, TF1* locTimeCut_ns) const
{
	auto locIterator = dPIDTimingCuts.find(locPID);
	if(locIterator == dPIDTimingCuts.end())
		return false;

	auto locSystemMap = locIterator->second;
	auto locSystemIterator = locSystemMap.find(locSystem);
	if(locSystemIterator == locSystemMap.end())
		return false;

	locTimeCut_ns = locSystemIterator->second;
	return true;
}

inline vector<const DKinematicData*> DSourceComboTimeHandler::Get_BeamParticlesByRFBunch(int locCenterRFBunch, unsigned int locPlusMinusBunchRange) const
{
	vector<const DKinematicData*> locBeamParticles;
	for(auto locRFBunch = locCenterRFBunch - int(locPlusMinusBunchRange); locRFBunch <= locCenterRFBunch + int(locPlusMinusBunchRange); ++locRFBunch)
	{
		auto locIterator = dBeamParticlesByRFBunch.find(locRFBunch);
		if(locIterator != dBeamParticlesByRFBunch.end())
			locBeamParticles.insert(locBeamParticles.end(), locIterator->second.begin(), locIterator->second.end());
	}
	return locBeamParticles;
}

inline double DSourceComboTimeHandler::Calc_RFTime(int locNumRFBunchShifts) const
{
	return dInitialEventRFBunch->dTime + locNumRFBunchShifts*dBeamBunchPeriod;
}

inline double DSourceComboTimeHandler::Calc_PropagatedRFTime(double locPrimaryVertexZ, int locNumRFBunchShifts, double locDetachedVertexTimeOffset) const
{
	//propagate rf time to vertex and add time offset (faster to just do it here rather than for each particle)
	return Calc_RFTime(locNumRFBunchShifts) + (locPrimaryVertexZ - dTargetCenter.Z())/SPEED_OF_LIGHT + locDetachedVertexTimeOffset;
}

inline shared_ptr<const DKinematicData> DSourceComboTimeHandler::Create_KinematicData_Photon(const DNeutralShower* locNeutralShower, const DVector3& locVertex) const
{
	auto locKinematicsPair = Calc_Photon_Kinematics(locNeutralShower, locVertex);
	return std::make_shared<const DKinematicData>(Gamma, locKinematicsPair.first, locVertex, locKinematicsPair.second);
}

inline pair<DVector3, double> DSourceComboTimeHandler::Calc_Photon_Kinematics(const DNeutralShower* locNeutralShower, const DVector3& locVertex) const
{
	//returns momentum, vertex time
	auto locPath = locNeutralShower->dSpacetimeVertex.Vect() - locVertex;
	auto locPathLength = locPath.Mag();
	auto locVertexTime = locNeutralShower->dSpacetimeVertex.T() - locPathLength/29.9792458;
	auto locMomentum = locNeutralShower->dEnergy*locPath.Unit();
	return std::make_pair(locMomentum, locVertexTime);
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

inline pair<double, double> DSourceComboTimeHandler::Calc_RFDeltaTChiSq(const DNeutralShower* locNeutralShower, const TVector3& locVertex, double locPropagatedRFTime) const
{
	//calc vertex time, get delta-t cut
	auto locPathLength = (locNeutralShower->dSpacetimeVertex.Vect() - locVertex).Mag();
	auto locVertexTime = locNeutralShower->dSpacetimeVertex.T() - locPathLength/29.9792458;
	auto locVertexTimeVariance = dUseSigmaForRFSelectionFlag ? (*(locNeutralShower->dCovarianceMatrix))(4, 4) : 1.0;

	auto locDeltaT = locVertexTime - locPropagatedRFTime;
	auto locChiSq = locDeltaT*locDeltaT/locVertexTimeVariance;
	if(dDebugLevel >= 5)
		cout << "neutral Calc_RFDeltaTChiSq(): vertex time, rf time, delta-t, chisq = " << locVertexTime << ", " << locPropagatedRFTime << ", " << locDeltaT << ", " << locChiSq << endl;

	return std::make_pair(locDeltaT, locChiSq);
}

inline pair<double, double> DSourceComboTimeHandler::Calc_RFDeltaTChiSq(const DChargedTrackHypothesis* locHypothesis, double locVertexTime, double locPropagatedRFTime) const
{
	//calc vertex time, get delta-t cut
	auto locErrorMatrix = locHypothesis->errorMatrix();
	auto locVertexTimeVariance = (dUseSigmaForRFSelectionFlag && (locErrorMatrix != nullptr)) ? (*locErrorMatrix)(6, 6) : 1.0;
	auto locDeltaT = locVertexTime - locPropagatedRFTime;
	auto locChiSq = locDeltaT*locDeltaT/locVertexTimeVariance;
	if(dDebugLevel >= 5)
		cout << "charged Calc_RFDeltaTChiSq(): vertex time, rf time, delta-t, chisq = " << locVertexTime << ", " << locPropagatedRFTime << ", " << locDeltaT << ", " << locChiSq << endl;

	return std::make_pair(locDeltaT, locChiSq);
}

inline TF1* DSourceComboTimeHandler::Get_TimeCutFunction(Particle_t locPID, DetectorSystem_t locSystem) const
{
	auto locPIDIterator = dPIDTimingCuts.find(locPID);
	if(locPIDIterator == dPIDTimingCuts.end())
		return nullptr;

	auto& locSystemMap = locPIDIterator->second;
	auto locSystemIterator = locSystemMap.find(locSystem);
	if(locSystemIterator == locSystemMap.end())
		return nullptr;

	return locSystemIterator->second;
}

}

#endif // DSourceComboTimeHandler_h
