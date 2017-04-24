#include "ANALYSIS/DSourceComboP4Handler.h"

/************************************************ 2-PHOTON MASS ERROR DERIVATION ***********************************************
*
* For this error estimate, consider the decay pi0 -> 2g  (or eta -> 2g).
* The equation for the invariant mass squared of the 2 photons is:
* m^2 = (E1 + E2)^2 - (px1 + px2)^2 - (py1 + py2)^2 - (pz1 + pz2)^2
*
* The difference in mass squared due to the error is: (using "_g" prefix for guess and "_t" prefix for "true")
* delta(m^2) = m_g^2 - m^2_t
* delta(m^2) = (px1 + px2)^2_t + (py1 + py2)^2_t + (pz1 + pz2)^2_t - (px1 + px2)^2_g - (py1 + py2)^2_g - (pz1 + pz2)^2_g
*
* Simplifying, we will choose the 2 photons to have opposite phi's at 0 & 180 degrees
* Thus, py1 = py2 = 0 for guess & true momentum.
* delta(m^2) = (px1 + px2)^2_t + (pz1 + pz2)^2_t - (px1 + px2)^2_g - (pz1 + pz2)^2_g
* Also, px1_unit = -px2_unit, but we'll use this a little later.
*
* Now, assume that both photons have the same E & theta: Whatever the worst-case is, it will be the same for both photons.
* Since from before px1_unit = -px2_unit, Now px1_t = -px2_t and px1_g = -px2_g. Also, pz1_t = pz2_t = pz_t, & pz1_g = pz2_g = pz_g
* delta(m^2) = 4*pz_t^2 - 4*pz_g^2
*
* Now for photons, pz = E*dz/d, where d is the distance between the shower and the vertex, and dz is the z-component of d.
* Plugging this in gives:
* delta(m^2) = 4*E^2*[(dz_t/d_t)^2 - (dz_g/d_g)^2]
* Using dz_g/d_g = cos(theta_g) gives:
* delta(m^2) = 4*E^2*[(dz_t/d_t)^2 - cos^2(theta_g)]
*
* Now, m_g^2 = (E1 + E2)^2 - (px1 + px2)^2_g - (py1 + py2)^2_g - (pz1 + pz2)^2_g
* However, we've defined our photon guesses pretty narrowly: py = 0, px1 = -px2, pz1 = pz2, E1 = E2
* Thus, m_g^2 = (2E)^2 - (2pz_g)^2
* Plugging in pz_g = E*dz_g/d_g yields:
* Thus, m_g^2 = 4*E^2*[1 - (dz_g/d_g)^2]
* Using dz_g/d_g = cos(theta_g)
* Thus, m_g^2 = 4*E^2*sin^2(theta_g)
* Rearranging gives: 4*E^2 = m_g^2/sin^2(theta_g)
*
* Plugging the above into delta(m^2)
* delta(m^2) = m_g^2*[(dz_t/d_t)^2 - cos^2(theta_g)]/sin^2(theta_g)
*
* But, we want delta_m, not delta(m^2)
* delta_m = m_g - m_t, m_t = sqrt(m_g^2 - delta(m^2))
* delta_m = m_g - sqrt(m_g^2 - delta(m^2))
* delta_m = m_g - sqrt(m_g^2 - m_g^2*[(dz_t/d_t)^2 - cos^2(theta_g)]/sin^2(theta_g))
*
* Rearrange and cancel terms
* delta_m = m_g - m_g*sqrt(1 - [(dz_t/d_t)^2 - cos^2(theta_g)]/sin^2(theta_g))
* delta_m = m_g - m_g*sqrt([sin^2(theta_g) - (dz_t/d_t)^2 + cos^2(theta_g)]/sin^2(theta_g))
* delta_m = m_g - m_g*sqrt[1 - (dz_t/d_t)^2]/sin(theta_g)
* delta_m = m_g - m_g*sqrt[(d_t^2 - dz_t^2)/d_t^2]/sin(theta_g)
*
* Note that d_t^2 - dz_t^2 = dx^2 (since dy is zero)
* delta_m = m_g - m_g*sqrt[dx^2/d_t^2]/sin(theta_g)
* delta_m = m_g - m_g*dx/(d_t*sin(theta_g))
*
* Getting the true dz_t & d_t in terms of some z_error:
* d_t = sqrt(dx^2 + dz_t^2), dz_t = dz_g + z_error
* delta_m = m_g - m_g*dx/(sin(theta_g)*sqrt(dx^2 + (dz_g + z_error)^2))
* And dz_g = dx/tan(theta_g)
* delta_m = m_g - m_g*dx/(sin(theta_g)*sqrt(dx^2 + (dx/tan(theta_g) + z_error)^2))
* delta_m = m_g - m_g/(sin(theta_g)*sqrt(1 + (1/tan(theta_g) + z_error/dx)^2))
*
* For BCAL, dx = 65.
* For the FCAL, dx = dz*tan(theta), dz = 650 - z_error (approx 650):
* delta_m = m_g - m_g/(sin(theta_g)*sqrt(1 + (1 + z_error/(650 - z_error))^2/tan^2(theta_g)))
* delta_m = m_g - m_g/(cos(theta_g)*sqrt(tan^2(theta_g) + (1 + z_error/(650 - z_error))^2))
*
* delta_m Is larger at higher m_g, max at 45 degrees (and is thus small for FCAL)
* In fact, for the FCAL, delta_m is ~25 MeV for the eta mass when the z_error is 30cm (max: center of target + detached vertex)
* Therefore, if the center of the target is used, the error is negligible compared to the width of the mass cut.
*
* For the BCAL:
* With m_g near eta mass, z_error = 15: delta_m_max = ~60 MeV
* With m_g near eta mass, z_error = 10: delta_m_max = ~40 MeV
* With m_g near eta mass, z_error = 5: delta_m_max = ~20 MeV
* With m_g near eta mass, z_error = 3: delta_m_max = ~15 MeV
* With m_g near eta mass, z_error = 2: delta_m_max = ~9 MeV
* Instead of the above, you can of course plot the delta_m for real data, and get something similar
* So, choose a z_error of 5: compute at center of 10-cm-wide bins.
*
* OK, but what about eta -> 3pi0?  Or Lambda -> pi0, n?
* eta -> 3pi0: Max error for a given pi0 is ~small, not bad when combined with 3 others: it's fine as long as cut is wide.
* pi0, n: Neutron is likelys slow, similar to charged tracks: Error is much larger, cannot combo massive neutrals without exact vertex position
* Well, OK, we can COMBO them, but we can't place mass cuts.
*
*******************************************************************************************************************************/

namespace DAnalysis
{

DSourceComboP4Handler::DSourceComboP4Handler(JEventLoop* locEventLoop)
{
	//GET THE GEOMETRY
	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry* locGeometry = locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());

	//TARGET INFORMATION
	double locTargetCenterZ = 65.0;
	locGeometry->GetTargetZ(locTargetCenterZ);
	dTargetCenter.SetXYZ(0.0, 0.0, locTargetCenterZ);
	locGeometry->GetTargetLength(dTargetLength);

	//INITIALIZE PHOTON VERTEX-Z EVALUATION BINNING
	dPhotonVertexZBinWidth = 10.0;
	dPhotonVertexZRangeLow = dTargetCenter.Z() - dTargetLength/2.0 - 5.0;
	dNumPhotonVertexZBins = round((dTargetLength + 20.0)/dPhotonVertexZBinWidth);

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

void DSourceComboP4Handler::Setup_NeutralShowers(const vector<const DNeutralShower*>& locNeutralShowers)
{
	//Precompute kinematics for the neutral showers, before comboing
	//Even if it turns out some of this isn't technically needed,
	//it's still faster than doing a check to see if this has been done or not for every single photon-combo request

	//CALCULATE KINEMATICS
	//FCAL: at target center
	for(auto& locShower : locNeutralShowers)
	{
		if(locNeutralShower->dDetectorSystem == SYS_FCAL)
		{
			dPhotonKinematics[DSourceComboInfo::Get_VertexZIndex_FCAL()].emplace(locShower, Create_KinematicData(locShower, dTargetCenter));
			continue;
		}
		//BCAL: in vertex-z bins
		for(size_t loc_i = 0; loc_i < dNumPhotonVertexZBins; ++loc_i)
		{
			DVector3 locBinCenter(0.0, 0.0, Get_PhotonVertexZBinCenter(loc_i));
			dPhotonKinematics[loc_i].emplace(locShower, Create_KinematicData(locShower, locBinCenter));
		}
	}
}

DLorentzVector DSourceComboP4Handler::Get_P4(Particle_t locPID, const JObject* locObject, signed char locVertexZBin)
{
	if(ParticleCharge(locPID) != 0)
	{
		auto locChargedTrack = static_cast<const DChargedTrack*>(locObject);
		return locChargedTrack->Get_Hypothesis(locPID)->lorentzMomentum();
	}

	auto locNeutralShower = static_cast<const DNeutralShower*>(locObject);
	if(locPID == Gamma)
	{
		if(locNeutralShower->dDetectorSystem == SYS_FCAL)
			locVertexZBin = DSourceComboInfo::Get_VertexZIndex_FCAL();
		auto& locKinematicData = dPhotonKinematics[locVertexZBin][locNeutralShower];
		return locKinematicData->lorentzMomentum();
	}

	//HANDLE MASSIVE NEUTRAL CASE!!!
}

vector<int> DSourceComboP4Handler::Cut_InvariantMass(const DSourceCombo* locSourceCombo, Particle_t locDecayPID, vector<int> locValidRFBunches) const
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

DLorentzVector DSourceComboP4Handler::Calc_P4(const DSourceCombo* locSourceCombo, signed char locVertexZBin) const
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

DLorentzVector DSourceComboP4Handler::Calc_P4(const DSourceCombo* locSourceCombo, DVector3 locVertex, int locRFBunch) const
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
