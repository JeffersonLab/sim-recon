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
		DSourceComboP4Handler(JEventLoop* locEventLoop, const DSourceComboer* locSourceComboer);

		//SET HANDLERS
		void Set_SourceComboVertexer(const DSourceComboVertexer* locSourceComboVertexer){dSourceComboVertexer = locSourceComboVertexer;}
		void Set_SourceComboTimeHandler(const DSourceComboTimeHandler* locSourceComboTimeHandler){dSourceComboTimeHandler = locSourceComboTimeHandler;}

		//SETUP FOR EACH EVENT
		void Reset(void);
		void Set_PhotonKinematics(const DPhotonKinematicsByZBin& locPhotonKinematics){dPhotonKinematics = locPhotonKinematics;}


		DLorentzVector Get_P4(const DSourceCombo* locSourceCombo, signed char locVertexZBin = DSourceComboInfo::Get_VertexZIndex_Unknown());

		//CUT
		//use this method when the combo DOES NOT contain massive neutrals
		bool Cut_InvariantMass(const DSourceCombo* locSourceCombo, Particle_t locDecayPID, signed char locVertexZBin) const;
		//use this method when the combo contains massive neutrals
		vector<int> Cut_InvariantMass(const DSourceCombo* locSourceCombo, Particle_t locDecayPID, vector<int> locValidRFBunches, const DVector3& locVertex, double locNoShiftRFVertexTime) const;

		DLorentzVector Calc_P4(const DSourceCombo* locSourceCombo, signed char locVertexZBin) const;
		DLorentzVector Calc_P4(const DSourceCombo* locSourceCombo, DVector3 locVertex, int locRFBunch) const;
		DLorentzVector Calc_MassiveNeutralP4(const DNeutralShower* locNeutralShower, Particle_t locPID, const DVector3& locVertex, double locVertexTime) const;

	private:
		DLorentzVector Get_P4(Particle_t locPID, const JObject* locObject, signed char locVertexZBin, int locRFBunch);

		//HANDLERS
		const DSourceComboer* dSourceComboer; //for quickly determining whether a combo has a massive neutral or not, and to facilitate access to combo vertices
		const DSourceComboVertexer* dSourceComboVertexer; //for getting vertex positions & time offsets for massive neutral p4 calculations
		const DSourceComboTimeHandler* dSourceComboTimeHandler; //for getting the propagated RF time for massive neutral p4 calculations

		//NEUTRAL SHOWER DATA
		DPhotonKinematicsByZBin dPhotonKinematics; //FCAL shower data at center of target, BCAL in vertex-z bins

		//TOTAL FINAL STATE FOUR-MOMENTUM
		unordered_map<pair<const DSourceCombo*, signed char>, DLorentzVector> dFinalStateP4ByCombo; //signed char: vertex-z bin
		unordered_map<int, unordered_map<const DSourceCombo*, DLorentzVector>> dFinalStateP4ByCombo_HasMassiveNeutrals; //int: RF bunch

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

inline DLorentzVector DSourceComboP4Handler::Get_P4(const DSourceCombo* locSourceCombo, signed char locVertexZBin)
{
	auto locIterator = dFinalStateP4ByCombo.find(std::make_pair(locSourceCombo, locVertexZBin));
	if(locIterator != dFinalStateP4ByCombo.end())
		return *locIterator;

	return Calc_P4(locSourceCombo, locVertexZBin);
}

inline bool DSourceComboP4Handler::Cut_InvariantMass(const DSourceCombo* locSourceCombo, Particle_t locDecayPID, signed char locVertexZBin) const
{
	//Don't call if it contains massive neutrals! Call the other cut function instead!!
	auto locCutIterator = dInvariantMassCuts.find(locDecayPID);
	if(locCutIterator == dInvariantMassCuts.end())
		return true; //no cut to place!!
	auto& locMassCuts = locCutIterator->second;

	auto locInvariantMass = Calc_P4(locSourceCombo, locVertexZBin).M();
	return ((locInvariantMass >= locMassCuts.first) && (locInvariantMass <= locMassCuts.second));
}

}

#endif // DSourceComboP4Handler_h
