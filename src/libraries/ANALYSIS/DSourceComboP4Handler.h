#ifndef DSourceComboP4Handler_h
#define DSourceComboP4Handler_h

#include <unordered_map>
#include <vector>
#include <map>

#include "TF1.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TROOT.h"
#include "TFile.h"
#include "TDirectoryFile.h"

#include "JANA/JEventLoop.h"

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
		DSourceComboP4Handler(JEventLoop* locEventLoop, DSourceComboer* locSourceComboer);
		~DSourceComboP4Handler(void){Fill_Histograms();}

		//SET HANDLERS
		void Set_SourceComboVertexer(const DSourceComboVertexer* locSourceComboVertexer){dSourceComboVertexer = locSourceComboVertexer;}
		void Set_SourceComboTimeHandler(const DSourceComboTimeHandler* locSourceComboTimeHandler){dSourceComboTimeHandler = locSourceComboTimeHandler;}

		//SETUP FOR EACH EVENT
		void Reset(void);
		void Set_PhotonKinematics(const DPhotonKinematicsByZBin& locPhotonKinematics){dPhotonKinematics = locPhotonKinematics;}

		//GET/CALC PARTICLE P4
		DLorentzVector Get_P4_NotMassiveNeutral(Particle_t locPID, const JObject* locObject, signed char locVertexZBin) const;
		DLorentzVector Calc_MassiveNeutralP4(const DNeutralShower* locNeutralShower, Particle_t locPID, const DVector3& locVertex, double locVertexTime) const;

		//GET/CALC COMBO P4
		DLorentzVector Calc_P4_NoMassiveNeutrals(const DSourceCombo* locVertexCombo, signed char locVertexZBin);
		DLorentzVector Calc_P4_SourceParticles(const DSourceCombo* locVertexCombo, signed char locVertexZBin, const DVector3& locVertex, double locRFVertexTime);
		bool Calc_P4_Decay(bool locIsProductionVertex, const DSourceCombo* locReactionFullCombo, const DSourceCombo* locVertexCombo, const DSourceComboUse& locDecayUse, const DSourceCombo* locDecayCombo, DVector3 locVertex, int locRFBunch, double locRFVertexTime, DLorentzVector& locDecayP4);
		bool Calc_P4_HasMassiveNeutrals(bool locIsProductionVertex, const DSourceCombo* locReactionFullCombo, const DSourceCombo* locVertexCombo, DVector3 locVertex, int locRFBunch, double locRFVertexTime, const DSourceComboUse& locToExcludeUse, DLorentzVector& locTotalP4);
		DLorentzVector Get_P4_HasMassiveNeutrals(bool locIsProductionVertex, const DSourceCombo* locReactionFullCombo, const DSourceCombo* locVertexCombo, int locRFBunch) const; //ASSUMES ALREADY CALCULATED!!
		bool Cut_MissingMass(const DReaction* locReaction, const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, const DKinematicData* locBeamParticle, int locRFBunch);

		//CUT
		//use this method when the combo DOES NOT contain massive neutrals
		bool Cut_InvariantMass_NoMassiveNeutrals(const DSourceCombo* locVertexCombo, Particle_t locDecayPID, signed char locVertexZBin);
		bool Cut_InvariantMass_HasMassiveNeutral_OrPhotonVertex(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, vector<int>& locValidRFBunches);
		bool Cut_InvariantMass_HasMassiveNeutral(bool locIsProductionVertex, const DSourceCombo* locReactionFullCombo, const DSourceCombo* locVertexCombo, Particle_t locDecayPID, double locPrimaryVertexZ, const DVector3& locVertex, double locTimeOffset, vector<int>& locValidRFBunches);
		bool Cut_InvariantMass_MissingMassVertex(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, const DKinematicData* locBeamParticle, int locRFBunch);

	private:
		DLorentzVector Get_P4(Particle_t locPID, const JObject* locObject, signed char locVertexZBin, int locRFBunch);
		void Fill_Histograms(void);

		//HANDLERS
		DSourceComboer* dSourceComboer; //for quickly determining whether a combo has a massive neutral or not, and to facilitate access to combo vertices
		const DSourceComboVertexer* dSourceComboVertexer = nullptr; //for getting vertex positions & time offsets for massive neutral p4 calculations
		const DSourceComboTimeHandler* dSourceComboTimeHandler = nullptr; //for getting the propagated RF time for massive neutral p4 calculations
		size_t dDebugLevel = 0;

		//NEUTRAL SHOWER DATA
		DPhotonKinematicsByZBin dPhotonKinematics; //FCAL shower data at center of target, BCAL in vertex-z bins

		//TOTAL FINAL STATE FOUR-MOMENTUM
		map<pair<const DSourceCombo*, signed char>, DLorentzVector> dFinalStateP4ByCombo; //signed char: vertex-z bin
		map<tuple<bool, const DSourceCombo*, const DSourceCombo*, int>, DLorentzVector> dFinalStateP4ByCombo_HasMassiveNeutrals; //int: RF bunch //bool: is prod vertex //first combo: reaction full

		//CUTS
		double dMaxMassiveNeutralBeta = 0.99999;
		map<Particle_t, pair<double, double>> dInvariantMassCuts;
		map<Particle_t, pair<TF1*, TF1*>> dMissingMassSquaredCuts; //cuts are function of beam energy //For none missing, Particle_t = unknown

		//HISTOGRAMS
		map<Particle_t, TH1*> dHistMap_InvariantMass;
		map<Particle_t, TH2*> dHistMap_MissingMassSquaredVsBeamEnergy; //none missing = Unknown pid
		map<Particle_t, vector<float>> dInvariantMasses;
		map<Particle_t, vector<pair<float, float>>> dMissingMassPairs;
};

inline void DSourceComboP4Handler::Reset(void)
{
	Fill_Histograms();
	for(auto& locPIDPair : dInvariantMasses)
		locPIDPair.second.clear();
	for(auto& locPIDPair : dMissingMassPairs)
		locPIDPair.second.clear();

	dPhotonKinematics.clear();
	dFinalStateP4ByCombo.clear();
	dFinalStateP4ByCombo_HasMassiveNeutrals.clear();
}

inline bool DSourceComboP4Handler::Cut_InvariantMass_NoMassiveNeutrals(const DSourceCombo* locVertexCombo, Particle_t locDecayPID, signed char locVertexZBin)
{
	//Don't call if it contains massive neutrals! Call the other cut function instead!!
	auto locCutIterator = dInvariantMassCuts.find(locDecayPID);
	if(locCutIterator == dInvariantMassCuts.end())
		return true; //no cut to place!!
	auto& locMassCuts = locCutIterator->second;

	auto locInvariantMass = Calc_P4_NoMassiveNeutrals(locVertexCombo, locVertexZBin).M();

	//save and cut
	dInvariantMasses[locDecayPID].emplace_back(locInvariantMass);
	return ((locInvariantMass >= locMassCuts.first) && (locInvariantMass <= locMassCuts.second));
}

inline DLorentzVector DSourceComboP4Handler::Get_P4_HasMassiveNeutrals(bool locIsProductionVertex, const DSourceCombo* locReactionFullCombo, const DSourceCombo* locVertexCombo, int locRFBunch) const
{
	auto locTuple = std::make_tuple(locIsProductionVertex, locReactionFullCombo, locVertexCombo, locRFBunch);
	return dFinalStateP4ByCombo_HasMassiveNeutrals.find(locTuple)->second;
}

}

#endif // DSourceComboP4Handler_h
