#ifndef DSourceComboP4Handler_h
#define DSourceComboP4Handler_h

#include <unordered_map>
#include <vector>

#include "JANA/JEventLoop.h"

#include "particleType.h"
#include "DLorentzVector.h"
#include "PID/DKinematicData.h"
#include "PID/DNeutralShower.h"
#include "ANALYSIS/DSourceCombo.h"

using namespace std;
using namespace jana;

namespace DAnalysis
{

class DSourceComboVertexer;
class DSourceComboTimeHandler;

class DSourceComboP4Handler
{
	public:

		//CONSTRUCTORS
		DSourceComboP4Handler(void) = delete;
		DSourceComboP4Handler(JEventLoop* locEventLoop, DSourceComboer* locSourceComboer);

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
		DLorentzVector Calc_P4_NoMassiveNeutrals(const DSourceCombo* locSourceCombo, signed char locVertexZBin);
		DLorentzVector Calc_P4_SourceParticles(const DSourceCombo* locSourceCombo, signed char locVertexZBin, const DVector3& locVertex, double locRFVertexTime);
		bool Calc_P4_Decay(bool locIsProductionVertex, const DSourceCombo* locReactionChargedCombo, const DSourceCombo* locSourceCombo, const DSourceComboUse& locDecayUse, const DSourceCombo* locDecayCombo, DVector3 locVertex, int locRFBunch, double locRFVertexTime, DLorentzVector& locDecayP4);
		bool Calc_P4_HasMassiveNeutrals(bool locIsProductionVertex, const DSourceCombo* locReactionChargedCombo, const DSourceCombo* locSourceCombo, DVector3 locVertex, int locRFBunch, double locRFVertexTime, const DSourceComboUse& locToExcludeUse, DLorentzVector& locTotalP4);
		DLorentzVector Get_P4_HasMassiveNeutrals(bool locIsProductionVertex, const DSourceCombo* locReactionChargedCombo, const DSourceCombo* locSourceCombo, int locRFBunch) const; //ASSUMES ALREADY CALCULATED!!

		//CUT
		//use this method when the combo DOES NOT contain massive neutrals
		bool Cut_InvariantMass_NoMassiveNeutrals(const DSourceCombo* locSourceCombo, Particle_t locDecayPID, signed char locVertexZBin);
		bool Cut_InvariantMass_HasMassiveNeutral_OrPhotonVertex(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionChargedCombo, vector<int>& locValidRFBunches);
		bool Cut_InvariantMass_HasMassiveNeutral(bool locIsProductionVertex, const DSourceCombo* locReactionChargedCombo, const DSourceCombo* locSourceCombo, Particle_t locDecayPID, double locPrimaryVertexZ, const DVector3& locVertex, double locTimeOffset, vector<int>& locValidRFBunches);
		bool Cut_InvariantMass_MissingMassVertex(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, const DKinematicData* locBeamParticle, int locRFBunch);

	private:
		DLorentzVector Get_P4(Particle_t locPID, const JObject* locObject, signed char locVertexZBin, int locRFBunch);

		//HANDLERS
		DSourceComboer* dSourceComboer; //for quickly determining whether a combo has a massive neutral or not, and to facilitate access to combo vertices
		const DSourceComboVertexer* dSourceComboVertexer = nullptr; //for getting vertex positions & time offsets for massive neutral p4 calculations
		const DSourceComboTimeHandler* dSourceComboTimeHandler = nullptr; //for getting the propagated RF time for massive neutral p4 calculations

		//NEUTRAL SHOWER DATA
		DPhotonKinematicsByZBin dPhotonKinematics; //FCAL shower data at center of target, BCAL in vertex-z bins

		//TOTAL FINAL STATE FOUR-MOMENTUM
		unordered_map<pair<const DSourceCombo*, signed char>, DLorentzVector> dFinalStateP4ByCombo; //signed char: vertex-z bin
		unordered_map<tuple<bool, const DSourceCombo*, const DSourceCombo*, int>, DLorentzVector> dFinalStateP4ByCombo_HasMassiveNeutrals; //int: RF bunch //bool: is prod vertex //first combo: reaction full

		//CUTS
		double dMaxMassiveNeutralBeta = 0.99999;
		unordered_map<Particle_t, pair<double, double>> dInvariantMassCuts;
};

inline void DSourceComboP4Handler::Reset(void)
{
	dPhotonKinematics.clear();
	dFinalStateP4ByCombo.clear();
	dFinalStateP4ByCombo_HasMassiveNeutrals.clear();
}

inline bool DSourceComboP4Handler::Cut_InvariantMass_NoMassiveNeutrals(const DSourceCombo* locSourceCombo, Particle_t locDecayPID, signed char locVertexZBin)
{
	//Don't call if it contains massive neutrals! Call the other cut function instead!!
	auto locCutIterator = dInvariantMassCuts.find(locDecayPID);
	if(locCutIterator == dInvariantMassCuts.end())
		return true; //no cut to place!!
	auto& locMassCuts = locCutIterator->second;

	auto locInvariantMass = Calc_P4_NoMassiveNeutrals(locSourceCombo, locVertexZBin).M();
	return ((locInvariantMass >= locMassCuts.first) && (locInvariantMass <= locMassCuts.second));
}

inline DLorentzVector DSourceComboP4Handler::Get_P4_HasMassiveNeutrals(bool locIsProductionVertex, const DSourceCombo* locReactionChargedCombo, const DSourceCombo* locSourceCombo, int locRFBunch) const
{
	auto locTuple = std::make_tuple(locIsProductionVertex, locReactionChargedCombo, locSourceCombo, locRFBunch);
	return dFinalStateP4ByCombo_HasMassiveNeutrals.find(locTuple)->second;
}

}

#endif // DSourceComboP4Handler_h
