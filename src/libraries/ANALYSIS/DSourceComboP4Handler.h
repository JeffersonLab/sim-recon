#ifndef DSourceComboP4Handler_h
#define DSourceComboP4Handler_h

#include <unordered_map>
#include <vector>
#include <map>
#include <set>

#include "TF1.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TROOT.h"
#include "TFile.h"
#include "TDirectoryFile.h"

#include "particleType.h"
#include "DLorentzVector.h"
#include "PID/DKinematicData.h"
#include "PID/DNeutralShower.h"
#include "PID/DChargedTrack.h"
#include "ANALYSIS/DSourceCombo.h"
#include "ANALYSIS/DSourceComboTimeHandler.h"

using namespace std;
using namespace jana;

namespace DAnalysis
{

class DSourceComboVertexer;
class DSourceComboer;

class DSourceComboP4Handler
{
	public:

		//CONSTRUCTORS
		DSourceComboP4Handler(void) = delete;
		DSourceComboP4Handler(DSourceComboer* locSourceComboer, bool locCreateHistsFlag = true);
		~DSourceComboP4Handler(void){Fill_Histograms();}

		//SET HANDLERS
		void Set_SourceComboVertexer(const DSourceComboVertexer* locSourceComboVertexer){dSourceComboVertexer = locSourceComboVertexer;}
		void Set_SourceComboTimeHandler(const DSourceComboTimeHandler* locSourceComboTimeHandler){dSourceComboTimeHandler = locSourceComboTimeHandler;}

		//SETUP FOR EACH EVENT
		void Reset(void);
		void Set_PhotonKinematics(const DPhotonKinematicsByZBin& locPhotonKinematics){dPhotonKinematics = locPhotonKinematics;}

		//GET CUTS
		bool Get_InvariantMassCuts(Particle_t locPID, pair<float, float>& locMinMaxCuts_GeV) const;
		bool Get_MissingMassSquaredCuts(Particle_t locPID, pair<TF1*, TF1*>& locMinMaxCuts_GeVSq) const;

		//GET/CALC PARTICLE P4
		DLorentzVector Get_P4_NotMassiveNeutral(Particle_t locPID, const JObject* locObject, const DVector3& locVertex, bool locAccuratePhotonsFlag) const;
		DLorentzVector Calc_MassiveNeutralP4(const DNeutralShower* locNeutralShower, Particle_t locPID, const DVector3& locVertex, double locVertexTime) const;

		//GET/CALC COMBO P4
		DLorentzVector Calc_P4_NoMassiveNeutrals(const DSourceCombo* locVertexCombo, const DVector3& locVertex, signed char locVertexZBin, const DKinematicData* locBeamParticle, bool locAccuratePhotonsFlag);
		DLorentzVector Calc_P4_SourceParticles(const DSourceCombo* locVertexCombo, const DVector3& locVertex, double locRFVertexTime, bool locAccuratePhotonsFlag);
		bool Calc_P4_Decay(bool locIsProductionVertex, const DSourceCombo* locReactionFullCombo, const DSourceComboUse& locDecayUse, const DSourceCombo* locDecayCombo, DVector3 locVertex, int locRFBunch, double locRFVertexTime, DLorentzVector& locDecayP4, const DKinematicData* locBeamParticle, bool locAccuratePhotonsFlag);
		bool Calc_P4_HasMassiveNeutrals(bool locIsProductionVertex, const DSourceCombo* locReactionFullCombo, const DSourceCombo* locVertexCombo, DVector3 locVertex, int locRFBunch, double locRFVertexTime, const DSourceComboUse& locToExcludeUse, DLorentzVector& locTotalP4, const DKinematicData* locBeamParticle, bool locAccuratePhotonsFlag);
		DLorentzVector Get_P4_HasMassiveNeutrals(bool locIsProductionVertex, const DSourceCombo* locReactionFullCombo, const DSourceCombo* locVertexCombo, int locRFBunch) const; //ASSUMES ALREADY CALCULATED!!

		//CUT
		//use this method when the combo DOES NOT contain massive neutrals
		bool Cut_InvariantMass_NoMassiveNeutrals(const DSourceCombo* locVertexCombo, Particle_t locDecayPID, const DVector3& locVertex, signed char locVertexZBin, bool locAccuratePhotonsFlag);
		bool Cut_InvariantMass_HasMassiveNeutral_OrPhotonVertex(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, vector<int>& locValidRFBunches);
		bool Cut_InvariantMass_HasMassiveNeutral(bool locIsProductionVertex, const DSourceCombo* locReactionFullCombo, const DSourceCombo* locVertexCombo, Particle_t locDecayPID, double locPrimaryVertexZ, const DVector3& locVertex, double locTimeOffset, vector<int>& locValidRFBunches, const DKinematicData* locBeamParticle, bool locAccuratePhotonsFlag);
		bool Cut_InvariantMass_MissingMassVertex(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, const DKinematicData* locBeamParticle, int locRFBunch);
		bool Cut_InvariantMass_AccuratePhotonKinematics(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, const DKinematicData* locBeamParticle, int locRFBunch);
		bool Cut_MissingMass(const DReaction* locReaction, const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, const DKinematicData* locBeamParticle, int locRFBunch);

	private:
		DLorentzVector Get_P4(Particle_t locPID, const JObject* locObject, signed char locVertexZBin, int locRFBunch);
		bool Get_InvariantMassCut(const DSourceCombo* locSourceCombo, Particle_t locDecayPID, bool locAccuratePhotonsFlag, pair<float, float>& locMinMaxMassCuts_GeV) const;
		void Fill_Histograms(void);

		int dDebugLevel = 0;

		//HANDLERS
		DSourceComboer* dSourceComboer; //for quickly determining whether a combo has a massive neutral or not, and to facilitate access to combo vertices
		const DSourceComboVertexer* dSourceComboVertexer = nullptr; //for getting vertex positions & time offsets for massive neutral p4 calculations
		const DSourceComboTimeHandler* dSourceComboTimeHandler = nullptr; //for getting the propagated RF time for massive neutral p4 calculations

		//NEUTRAL SHOWER DATA
		DPhotonKinematicsByZBin dPhotonKinematics; //FCAL shower data at center of target, BCAL in vertex-z bins

		//TOTAL FINAL STATE FOUR-MOMENTUM
		map<pair<const DSourceCombo*, signed char>, DLorentzVector> dFinalStateP4ByCombo; //signed char: vertex-z bin
		map<tuple<bool, const DSourceCombo*, const DSourceCombo*, int>, DLorentzVector> dFinalStateP4ByCombo_HasMassiveNeutrals; //int: RF bunch //bool: is prod vertex //first combo: reaction full

		//CUTS
		double dMaxMassiveNeutralBeta = 0.99999;
		double d2PhotonInvariantMassCutError = 0.02;
		map<Particle_t, pair<double, double>> dInvariantMassCuts;
		map<Particle_t, pair<TF1*, TF1*>> dMissingMassSquaredCuts; //cuts are function of beam energy //For none missing, Particle_t = unknown

		//HISTOGRAMS
		TH2* dHist_NoneMissing_MissingPzVsBeamEnergy_PreMissMassSqCut;
		vector<pair<float, float>> dMissingPzVsBeamEnergy_PreMissMassSqCut;
		TH2* dHist_NoneMissing_MissingPzVsBeamEnergy_PostMissMassSqCut;
		vector<pair<float, float>> dMissingPzVsBeamEnergy_PostMissMassSqCut;
		TH2* dHist_NoneMissing_MissingPtVsMissingPz_PreMissMassSqCut;
		vector<pair<float, float>> dMissingPtVsMissingPz_PreMissMassSqCut;
		TH2* dHist_NoneMissing_MissingPtVsMissingPz_PostMissMassSqCut;
		vector<pair<float, float>> dMissingPtVsMissingPz_PostMissMassSqCut;

		set<tuple<Particle_t, const DSourceCombo*>> dInvariantMassFilledSet; //only filled once per combo, even if used for different vertex-z bins!!!!
		set<tuple<Particle_t, bool, const DSourceCombo*, const DSourceCombo*, int>> dInvariantMassFilledSet_MassiveNeutral; //int: RF bunch //bool: is prod vertex //first combo: reaction full
		map<DetectorSystem_t, TH1*> dHistMap_2GammaMass; //SYS_BCAL/FCAL if both in bcal/fcal, SYS_NULL if one each
		map<DetectorSystem_t, vector<float>> d2GammaInvariantMasses;
		map<Particle_t, TH1*> dHistMap_InvariantMass;
		map<Particle_t, TH2*> dHistMap_MissingMassSquaredVsBeamEnergy; //none missing = Unknown pid
		map<Particle_t, vector<float>> dInvariantMasses;
		map<Particle_t, vector<pair<float, float>>> dMissingMassPairs;
};

inline void DSourceComboP4Handler::Reset(void)
{
	Fill_Histograms();

	dInvariantMassFilledSet.clear();
	dInvariantMassFilledSet_MassiveNeutral.clear();
	dPhotonKinematics.clear();
	dFinalStateP4ByCombo.clear();
	dFinalStateP4ByCombo_HasMassiveNeutrals.clear();
}

inline bool DSourceComboP4Handler::Get_InvariantMassCuts(Particle_t locPID, pair<float, float>& locMinMaxCuts_GeV) const
{
	auto locIterator = dInvariantMassCuts.find(locPID);
	if(locIterator == dInvariantMassCuts.end())
		return false;
	locMinMaxCuts_GeV = locIterator->second;
	return true;
}

inline bool DSourceComboP4Handler::Get_MissingMassSquaredCuts(Particle_t locPID, pair<TF1*, TF1*>& locMinMaxCuts_GeVSq) const
{
	auto locIterator = dMissingMassSquaredCuts.find(locPID);
	if(locIterator == dMissingMassSquaredCuts.end())
		return false;
	locMinMaxCuts_GeVSq = locIterator->second;
	return true;
}

inline DLorentzVector DSourceComboP4Handler::Get_P4_HasMassiveNeutrals(bool locIsProductionVertex, const DSourceCombo* locReactionFullCombo, const DSourceCombo* locVertexCombo, int locRFBunch) const
{
	auto locTuple = std::make_tuple(locIsProductionVertex, locReactionFullCombo, locVertexCombo, locRFBunch);
	return dFinalStateP4ByCombo_HasMassiveNeutrals.find(locTuple)->second;
}

}

#endif // DSourceComboP4Handler_h
