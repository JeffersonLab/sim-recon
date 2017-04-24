#ifndef DSourceComboP4Handler_h
#define DSourceComboP4Handler_h

#include <unordered_map>
#include <vector>

#include "particleType.h"
#include "DLorentzVector.h"
#include "PID/DKinematicData.h"
#include "ANALYSIS/DSourceCombo.h"

using namespace std;
using namespace jana;

namespace DAnalysis
{

//OK, there are several times we might call this:
//Charged only
//With neutrals, but vertex not chosen yet
//With neutrals, vertex chosen, but RF bunch not chosen yet (choosing it for massive neutrals!)
//We will not bother calling/cutting a final time for photons with a finalized vertex
class DSourceComboP4Handler
{
	public:
		DSourceComboP4Handler(void);

		void Reset(void);

		DLorentzVector Get_P4(const DSourceCombo* locSourceCombo, signed char locVertexZBin = DSourceComboInfo::Get_VertexZIndex_Unknown());

		//CUT
		//use this method when the combo DOES NOT contain massive neutrals
		bool Cut_InvariantMass(const DSourceCombo* locSourceCombo, Particle_t locDecayPID, signed char locVertexZBin) const;
		//use this method when the combo contains massive neutrals
		vector<int> Cut_InvariantMass(const DSourceCombo* locSourceCombo, Particle_t locDecayPID, DVector3 locVertex, vector<int> locValidRFBunches) const;

		DLorentzVector Calc_P4(const DSourceCombo* locSourceCombo, signed char locVertexZBin) const;
		DLorentzVector Calc_P4(const DSourceCombo* locSourceCombo, DVector3 locVertex, int locRFBunch) const;

	private:
		DLorentzVector Get_P4(Particle_t locPID, const JObject* locObject);

		//TOTAL FINAL STATE FOUR-MOMENTUM
		unordered_map<pair<const DSourceCombo*, signed char>, DLorentzVector> dFinalStateP4ByCombo; //signed char: vertex-z bin
		unordered_map<int, unordered_map<const DSourceCombo*, DLorentzVector>> dFinalStateP4ByCombo_HasMassiveNeutrals; //int: RF bunch

		//CUTS
		unordered_map<Particle_t, pair<double, double>> dInvariantMassCuts;
};

inline DSourceComboP4Handler::DSourceComboP4Handler(void)
{
	//INVARIANT MASS CUTS: MESONS
	dInvariantMassCuts.emplace(Pi0, std::make_pair(0.08, 0.19));
	dInvariantMassCuts.emplace(KShort, std::make_pair(0.3, 0.7));
	dInvariantMassCuts.emplace(Eta, std::make_pair(0.3, 0.8));
	dInvariantMassCuts.emplace(omega, std::make_pair(0.4, 1.2));
	dInvariantMassCuts.emplace(EtaPrime, std::make_pair(0.6, 1.3));
	dInvariantMassCuts.emplace(phiMeson, std::make_pair(0.8, 1.2));
	dInvariantMassCuts.emplace(Jpsi, std::make_pair(2.7, 3.5));

	//INVARIANT MASS CUTS: BARYONS
	dInvariantMassCuts.emplace(Lambda, std::make_pair(1.0, 1.2));
	dInvariantMassCuts.emplace(Sigma0, std::make_pair(1.1, 1.3));
	dInvariantMassCuts.emplace(SigmaPlus, dInvariantMassCuts[Sigma0]);
	dInvariantMassCuts.emplace(SigmaMinus, dInvariantMassCuts[Sigma0]);
	dInvariantMassCuts.emplace(XiMinus, std::make_pair(1.1, 1.5));
	dInvariantMassCuts.emplace(Xi0, dInvariantMassCuts[XiMinus]);
}

inline void DSourceComboP4Handler::Reset(void)
{
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

inline DLorentzVector DSourceComboP4Handler::Get_P4(Particle_t locPID, const JObject* locObject)
{
	if(ParticleCharge(locPID) != 0)
	{
		auto locChargedTrack = static_cast<const DChargedTrack*>(locObject);
		return locChargedTrack->Get_Hypothesis(locPID)->lorentzMomentum();
	}

	auto locNeutralShower = static_cast<const DNeutralShower*>(locObject);
	if(locPID == Gamma)
	{
		auto& locKinematicsMap = (locNeutralShower->dDetectorSystem == SYS_FCAL) ? dFCALKinematics : dBCALKinematics[locVertexZBin];
		auto& locKinematicData = locKinematicsMap[locNeutralShower];
		return locKinematicData->lorentzMomentum();
	}

	//HANDLE MASSIVE NEUTRAL CASE!!!
}

inline bool DSourceComboP4Handler::Cut_InvariantMass(const DSourceCombo* locSourceCombo, Particle_t locDecayPID) const
{
	//Don't call if it contains massive neutrals! Call the other cut function instead!!
	auto locCutIterator = dInvariantMassCuts.find(locDecayPID);
	if(locCutIterator == dInvariantMassCuts.end())
		return true; //no cut to place!!
	auto& locMassCuts = locCutIterator->second;

	auto locInvariantMass = Calc_P4(locSourceCombo).M();
	return ((locInvariantMass >= locMassCuts.first) && (locInvariantMass <= locMassCuts.second));
}

inline vector<int> DSourceComboP4Handler::Cut_InvariantMass(const DSourceCombo* locSourceCombo, Particle_t locDecayPID, vector<int> locValidRFBunches) const
{
	//cuts on possible RF bunches for the massive neutrals
	//if no possible rf bunch yields a massive-neutral-momentum that passes the invariant mass cut, returns an empty vector
	auto locCutIterator = dInvariantMassCuts.find(locDecayPID);
	if(locCutIterator == dInvariantMassCuts.end())
		return locSourceCombo; //no cut to place!!
	auto& locMassCuts = locCutIterator->second;

	//function for calculating and cutting the invariant mass for each rf bunch
	auto CalcAndCut_InvariantMass = [&locSourceCombo, &locMassCuts](int locRFBunch) -> bool
	{
		auto locInvariantMass = Calc_P4(locSourceCombo, locRFBunch).M();
		return ((locInvariantMass >= locMassCuts.first) && (locInvariantMass <= locMassCuts.second));
	};

	//apply the function
	locValidRFBunches.erase(std::remove_if(locValidRFBunches.begin(), locValidRFBunches.end(), CalcAndCut_InvariantMass), locValidRFBunches.end());
	return locValidRFBunches;
}

inline DLorentzVector DSourceComboP4Handler::Calc_P4(const DSourceCombo* locSourceCombo, signed char locVertexZBin) const
{
	DLorentzVector locTotalP4;

	//loop over particles
	auto locSourceParticles = locSourceCombo->Get_SourceParticles(false); //false: NOT the whole chain
	for(auto locParticlePair : locSourceParticles)
		locTotalP4 += Get_P4(locParticlePair.first, locParticlePair.second);

	//loop over decays
	auto locFurtherDecayCombos = locSourceCombo->Get_FurtherDecayCombos();
	for(const auto& locCombosByUsePair : locFurtherDecayCombos)
	{
		for(const auto& locCombo : locCombosByUsePair.second)
			locTotalP4 += Get_P4(locCombo);
	}

	//save the results and return
	dFinalStateP4ByCombo.emplace(std::make_pair(locSourceCombo, locVertexZBin), locTotalP4);
	return locTotalP4;
}


inline DLorentzVector DSourceComboP4Handler::Calc_P4(const DSourceCombo* locSourceCombo, DVector3 locVertex, int locRFBunch) const
{
	auto locSourceParticles = locSourceCombo->Get_SourceParticles(false); //false: NOT the whole chain
	//vertex-z bin may be different for decay products! (detached vertex)
	//save/retrieve masses by combo instead

	DLorentzVector locTotalP4;
	for(auto locParticlePair : locSourceParticles)
	{

//handle massive neutral case!!
	}

	//now loop over decays
	auto locFurtherDecayCombos = locSourceCombo->Get_FurtherDecayCombos();
	for(const auto& locCombosByUsePair : locFurtherDecayCombos)
	{
		for(const auto& locCombo : locCombosByUsePair.second)
		{
			auto locIterator = dFinalStateP4ByCombo.find(locCombo);
			if(locIterator != dFinalStateP4ByCombo.end())
				locTotalP4 += locIterator->second;
			else
				locTotalP4 += dFinalStateP4ByCombo_HasMassiveNeutrals[locRFBunch][locCombo];
			//else call Calc_InvariantMass_HasMassiveNeutral!!
		}
	}

	dFinalStateP4ByCombo_HasMassiveNeutrals[locRFBunch].emplace(locSourceCombo, locTotalP4);
	return locTotalP4;
}

}

#endif // DSourceComboP4Handler_h
