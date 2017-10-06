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
		void Set_DebugLevel(int locDebugLevel){dDebugLevel = locDebugLevel;}

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
		//locRFVertexTime: the RF time propagated to the vertex, through any decaying particles if necessary
		DLorentzVector Calc_P4_NoMassiveNeutrals(const DSourceCombo* locReactionCombo, const DSourceCombo* locVertexCombo, const DVector3& locVertex, signed char locVertexZBin, const DKinematicData* locBeamParticle, const DSourceComboUse& locToExcludeUse, size_t locInstanceToExclude, bool locAccuratePhotonsFlag);
		DLorentzVector Calc_P4_SourceParticles(const DSourceCombo* locVertexCombo, const DVector3& locVertex, double locRFVertexTime, bool locAccuratePhotonsFlag);
		bool Calc_P4_Decay(bool locIsProductionVertex, bool locIsPrimaryProductionVertex, const DSourceCombo* locReactionFullCombo, const DSourceComboUse& locDecayUse, const DSourceCombo* locDecayCombo, DVector3 locVertex, double locTimeOffset, int locRFBunch, double locRFVertexTime, DLorentzVector& locDecayP4, const DKinematicData* locBeamParticle, const DSourceComboUse& locToExcludeUse, size_t locInstanceToExclude, bool locAccuratePhotonsFlag);
		bool Calc_P4_HasMassiveNeutrals(bool locIsProductionVertex, bool locIsPrimaryProductionVertex, const DSourceCombo* locReactionFullCombo, const DSourceCombo* locCurrentCombo, DVector3 locVertex, double locTimeOffset, int locRFBunch, double locRFVertexTime, const DSourceComboUse& locToExcludeUse, size_t locInstanceToExclude, DLorentzVector& locTotalP4, const DKinematicData* locBeamParticle, bool locAccuratePhotonsFlag);

		//CUT
		//use this method when the combo DOES NOT contain massive neutrals
		bool Cut_InvariantMass_NoMassiveNeutrals(const DSourceCombo* locVertexCombo, Particle_t locDecayPID, Particle_t locTargetPIDToSubtract, const DVector3& locVertex, signed char locVertexZBin, bool locAccuratePhotonsFlag);
		bool Cut_InvariantMass_HasMassiveNeutral_OrPhotonVertex(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, vector<int>& locValidRFBunches);
		bool Cut_InvariantMass_HasMassiveNeutral(bool locIsProductionVertex, bool locIsPrimaryProductionVertex, const DSourceCombo* locReactionFullCombo, const DSourceCombo* locVertexCombo, Particle_t locDecayPID, Particle_t locTargetPIDToSubtract, double locPrimaryVertexZ, const DVector3& locVertex, double locTimeOffset, vector<int>& locValidRFBunches, const DKinematicData* locBeamParticle, bool locAccuratePhotonsFlag);
		bool Cut_InvariantMass_MissingMassVertex(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, const DKinematicData* locBeamParticle, int locRFBunch);
		bool Cut_InvariantMass_AccuratePhotonKinematics(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, const DKinematicData* locBeamParticle, int locRFBunch);
		bool Cut_MissingMassSquared(const DReaction* locReaction, const DReactionVertexInfo* locReactionVertexInfo, const DSourceComboUse& locReactionFullComboUse, const DSourceCombo* locReactionFullCombo, const DKinematicData* locBeamParticle, int locRFBunch);

	private:
		DLorentzVector Get_P4(Particle_t locPID, const JObject* locObject, signed char locVertexZBin, int locRFBunch);
		bool Get_InvariantMassCut(const DSourceCombo* locSourceCombo, Particle_t locDecayPID, bool locAccuratePhotonsFlag, pair<float, float>& locMinMaxMassCuts_GeV) const;
		bool Cut_MissingMassSquared(const DReaction* locReaction, const DReactionVertexInfo* locReactionVertexInfo, const DSourceComboUse& locReactionFullComboUse, const DSourceCombo* locFullCombo, Particle_t locMissingPID, int locStepIndex, int locDecayStepIndex, const DLorentzVector& locInitialStateP4, int locRFBunch, const DKinematicData* locBeamParticle, DLorentzVector& locMissingP4);
		void Fill_Histograms(void);

		int dDebugLevel = 0;
		bool dPrintCutFlag = false;

		//Get/Set Cuts
		void Define_DefaultCuts(void);
		void Get_CommandLineCuts_MM2(void);
		void Get_CommandLineCuts_IM(void);
		void Get_CommandLineCuts_MissingEnergy(void);
		void Create_CutFunctions(void);

		//HANDLERS
		DSourceComboer* dSourceComboer; //for quickly determining whether a combo has a massive neutral or not, and to facilitate access to combo vertices
		const DSourceComboVertexer* dSourceComboVertexer = nullptr; //for getting vertex positions & time offsets for massive neutral p4 calculations
		const DSourceComboTimeHandler* dSourceComboTimeHandler = nullptr; //for getting the propagated RF time for massive neutral p4 calculations

		//NEUTRAL SHOWER DATA
		DPhotonKinematicsByZBin dPhotonKinematics; //FCAL shower data at center of target, BCAL in vertex-z bins

		//TOTAL FINAL STATE FOUR-MOMENTUM
		map<pair<const DSourceCombo*, signed char>, DLorentzVector> dFinalStateP4ByCombo; //signed char: vertex-z bin
		//int: RF bunch //bool: is prod vertex //first combo: reaction full //kindata: beam //use: use to exclude //size_t: instance to exclude
		map<tuple<bool, const DSourceCombo*, const DSourceCombo*, int, const DKinematicData*, DSourceComboUse, size_t>, DLorentzVector> dFinalStateP4ByCombo_HasMassiveNeutrals;

		//CUT DEFAULTS
		string dDefaultMissingMassSquaredCutFunctionString = "[0]";
		map<Particle_t, pair<string, string>> dMissingMassSquaredCuts_TF1FunctionStrings; //pair: low bound, high bound
		map<Particle_t, pair<vector<double>, vector<double>>> dMissingMassSquaredCuts_TF1Params; //pair: low bound, high bound
		pair<string, string> dMissingEnergyCuts_TF1FunctionStrings; //pair: low bound, high bound
		pair<vector<double>, vector<double>> dMissingEnergyCuts_TF1Params; //pair: low bound, high bound

		//CUTS
		double dMaxMassiveNeutralBeta = 0.99999;
		double d2PhotonInvariantMassCutError = 0.02; //see derivation at top of .cc file
		map<Particle_t, pair<double, double>> dInvariantMassCuts;
		map<Particle_t, pair<TF1*, TF1*>> dMissingMassSquaredCuts; //cuts are function of beam energy //For none missing, Particle_t = unknown
		pair<TF1*, TF1*> dMissingECuts; //for no-missing-particle only

		//HISTOGRAMS
		TH2* dHist_NoneMissing_MissingEVsBeamEnergy_PreMissMassSqCut;
		vector<pair<float, float>> dMissingEVsBeamEnergy_PreMissMassSqCut;
		TH2* dHist_NoneMissing_MissingEVsBeamEnergy_PostMissMassSqCut;
		vector<pair<float, float>> dMissingEVsBeamEnergy_PostMissMassSqCut;
		TH2* dHist_NoneMissing_MissingPtVsMissingE_PreMissMassSqCut;
		vector<pair<float, float>> dMissingPtVsMissingE_PreMissMassSqCut;
		TH2* dHist_NoneMissing_MissingPtVsMissingE_PostMissMassSqCut;
		vector<pair<float, float>> dMissingPtVsMissingE_PostMissMassSqCut;

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

}

#endif // DSourceComboP4Handler_h
