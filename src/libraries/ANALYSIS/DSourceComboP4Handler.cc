#include "ANALYSIS/DSourceComboP4Handler.h"


/*************************************************** PHOTON P4 RECONSTRUCTION **************************************************
*
* The exact photon momentum is unknown until its production vertex is known.
* However, that vertex is combo-dependent. We'd like to make cuts on the pi0 mass globally in advance, rather than on a combo-by-combo basis.
* This would be a huge savings in time and memory.
*
* The momentum of the hypotheses calculated in DNeutralParticle is based on the DVertex (class) position.
* This vertex is determined from all of the "good" charged tracks in the event, typically by doing a kinematic fit.
*
* However, there a several potential problems with using this vertex:
* 1) There may have been extra (junk, accidental) tracks reconstructed in the event. These will throw off the vertex position.
*    And, if there are only 2 tracks, it can throw it off considerably.
* 2) The track position & momentum are different for each hypothesis, and it's not clear in advance which hypothesis should be used.
* 3) The photons may have come from a detached vertex, which can be 10+ cm downstream of the main vertex.
*
* So while the DVertex is OK for someone looking for a quick estimate, it should not be used in an actual analysis.
*
* Now, as we said before, the true photon momentum is combo-dependent, and we want to do loose mass cuts in a combo-independent way.
* So, we can compute all of the p4's at a specific vertex position (e.g. center of the target), rather than separately for each combo.
* But, how much of an impact will a given error in the vertex position have on the calculated 2-photon invariant mass?
*
* The calculation below determines that the maximum 2-photon-mass error occurs when both photons are at 45 degrees near the eta peak (less near pi0).
* Specifically: z_err = 5cm yields a mass error of ~20 MeV, 10cm -> ~40 MeV, 15cm -> ~60 MeV, etc.
*
* So, what is the maximum delta_m we can tolerate with our loose cuts?
* The idea is that the loose cuts should be wide enough around the signal peak to provide enough statistics for sideband subtraction.
* So, no matter what we choose, it won't affect the signal peak.  But we also don't want to affect the sidebands too much.
* E.g. pi0 mass peak is from ~110 -> ~160 MeV, loose cut is 50 -> 220 MeV at the moment
* Therefore, you probably want to keep the maximum delta_m at around the 20 MeV level.
* This means that the max z_error should be about 5 cm
*
* Therefore, every 10 cm, from 5cm upstream of the target to 15 cm downstream of the target (detached vertex) (5 bins):
* Compute p4s at the centers of these vertex-z bins and do loose mass cuts
*
*******************************************************************************************************************************/

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

DSourceComboP4Handler::DSourceComboP4Handler(JEventLoop* locEventLoop, DSourceComboer* locSourceComboer) : dSourceComboer(locSourceComboer)
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

DLorentzVector DSourceComboP4Handler::Get_P4_NonMassiveNeutral(Particle_t locPID, const JObject* locObject, signed char locVertexZBin) const
{
	if(ParticleCharge(locPID) != 0)
	{
		auto locChargedTrack = static_cast<const DChargedTrack*>(locObject);
		return locChargedTrack->Get_Hypothesis(locPID)->lorentzMomentum();
	}

	//assume is photon!
	auto locNeutralShower = static_cast<const DNeutralShower*>(locObject);
	if(locNeutralShower->dDetectorSystem == SYS_FCAL)
		locVertexZBin = DSourceComboInfo::Get_VertexZIndex_FCAL();
	auto& locKinematicData = dPhotonKinematics.find(locVertexZBin)->second.find(locNeutralShower)->second;
	return locKinematicData->lorentzMomentum();
}

DLorentzVector DSourceComboP4Handler::Calc_MassiveNeutralP4(const DNeutralShower* locNeutralShower, Particle_t locPID, const DVector3& locVertex, double locRFVertexTime) const
{
	//locRFVertexTime: the RF time propagated to the vertex, through any decaying particles if necessary
	auto locMass = ParticleMass(locPID);
	auto locPath = locNeutralShower->dSpacetimeVertex.Vect() - locVertex;
	double locDeltaT = locNeutralShower->dSpacetimeVertex.T() - locRFVertexTime;

	double locBeta = locPath.Mag()/(locDeltaT*29.9792458);
	if(locBeta >= 1.0)
		locBeta = dMaxMassiveNeutralBeta;
	if(locBeta < 0.0)
		locBeta = 0.0;

	auto locGamma = 1.0/sqrt(1.0 - locBeta*locBeta);
	auto locPMag = locGamma*locBeta*locMass;
	locPath.SetMag(locPMag); //is now the momentum!

	auto locEnergy = sqrt(locPMag*locPMag + locMass*locMass);
	return DLorentzVector(locPath, locEnergy);
}


DLorentzVector DSourceComboP4Handler::Calc_P4_NoMassiveNeutrals(const DSourceCombo* locSourceCombo, signed char locVertexZBin)
{
	auto locIterator = dFinalStateP4ByCombo.find(std::make_pair(locSourceCombo, locVertexZBin));
	if(locIterator != dFinalStateP4ByCombo.end())
		return *locIterator;

	DLorentzVector locTotalP4;

	//loop over particles
	//vertex-z bin may be different for decay products! (detached vertex)
	//save/retrieve masses by combo instead
	auto locSourceParticles = locSourceCombo->Get_SourceParticles(false); //false: NOT the whole chain
	for(auto locParticlePair : locSourceParticles)
		locTotalP4 += Get_P4_NonMassiveNeutral(locParticlePair.first, locParticlePair.second, locVertexZBin);

	//loop over decays
	auto locFurtherDecayCombos = locSourceCombo->Get_FurtherDecayCombos();
	for(const auto& locCombosByUsePair : locFurtherDecayCombos)
	{
		auto locDecayVertexZBin = std::get<1>(locCombosByUsePair.first);
		for(const auto& locCombo : locCombosByUsePair.second)
			locTotalP4 += Calc_P4_NoMassiveNeutrals(locCombo, locDecayVertexZBin);
	}

	//save the results and return
	dFinalStateP4ByCombo.emplace(std::make_pair(locSourceCombo, locVertexZBin), locTotalP4);
	return locTotalP4;
}

DLorentzVector DSourceComboP4Handler::Calc_P4_SourceParticles(const DSourceCombo* locSourceCombo, signed char locVertexZBin, const DVector3& locVertex, double locRFVertexTime)
{
	DLorentzVector locTotalP4(0.0, 0.0, 0.0, 0.0);
	auto locSourceParticles = locSourceCombo->Get_SourceParticles(false); //false: NOT the whole chain
	for(auto locParticlePair : locSourceParticles)
	{
		auto locPID = locParticlePair.first;
		if((ParticleCharge(locPID) == 0) && (ParticleMass(locPID) > 0.0))
		{
			auto locNeutralShower = static_cast<const DNeutralShower*>(locParticlePair.second);
			locTotalP4 += Calc_MassiveNeutralP4(locNeutralShower, locPID, locVertex, locRFVertexTime);
		}
		else
			locTotalP4 += Get_P4_NonMassiveNeutral(locParticlePair.first, locParticlePair.second, locVertexZBin);
	}
	return locTotalP4;
}

bool DSourceComboP4Handler::Calc_P4_Decay(bool locIsProductionVertex, const DSourceCombo* locReactionChargedCombo, const DSourceCombo* locSourceCombo, const DSourceComboUse& locDecayUse, const DSourceCombo* locDecayCombo, signed char locVertexZBin, DVector3 locVertex, int locRFBunch, double locRFVertexTime, DLorentzVector& locDecayP4)
{
	auto locDecayPID = std::get<0>(locDecayUse);
	auto locDecayVertexZBin = std::get<1>(locDecayUse);
	auto locHasMassiveNeutrals = dSourceComboer->Get_HasMassiveNeutrals(std::get<2>(locDecayUse));

	if(!IsDetachedVertex(locDecayPID))
	{
		if(!locHasMassiveNeutrals)
			locDecayP4 = Calc_P4_NoMassiveNeutrals(locSourceCombo, locDecayVertexZBin);
		else if(!Calc_P4_HasMassiveNeutrals(locIsProductionVertex, locReactionChargedCombo, locDecayCombo, locVertexZBin, locVertex, locRFBunch, locRFVertexTime, DSourceComboUse(Unknown, 0, nullptr), locDecayP4))
			return false;
		return true;
	}

	//detached vertex!
	if(!locHasMassiveNeutrals)
		locDecayP4 = Calc_P4_NoMassiveNeutrals(locSourceCombo, locDecayVertexZBin);
	else
	{
		//p4 better have already been calculated: if not, it means it was uncalcable (e.g. vertex unknown)
		auto locDecayP4LookupTuple = std::make_tuple(locIsProductionVertex, locReactionChargedCombo, locDecayCombo, locRFBunch);
		auto locDecayIterator = dFinalStateP4ByCombo_HasMassiveNeutrals.find(locDecayP4LookupTuple);
		if(locDecayIterator == dFinalStateP4ByCombo_HasMassiveNeutrals.end())
			return false; //failed!
		locDecayP4 = *locDecayIterator;
	}

	return true;
}

bool DSourceComboP4Handler::Calc_P4_HasMassiveNeutrals(bool locIsProductionVertex, const DSourceCombo* locReactionChargedCombo, const DSourceCombo* locSourceCombo, signed char locVertexZBin, DVector3 locVertex, int locRFBunch, double locRFVertexTime, const DSourceComboUse& locToExcludeUse, DLorentzVector& locTotalP4)
{
	auto locP4LookupTuple = std::make_tuple(locIsProductionVertex, locReactionChargedCombo, locSourceCombo, locRFBunch);
	auto locIterator = dFinalStateP4ByCombo_HasMassiveNeutrals.find(locP4LookupTuple);
	if(locIterator != dFinalStateP4ByCombo_HasMassiveNeutrals.end())
		return *locIterator;

	//final state particles
	locTotalP4 += Calc_P4_SourceParticles(locSourceCombo, locVertexZBin, locVertex, locRFVertexTime);

	//now loop over decays
	auto locFurtherDecayCombos = locSourceCombo->Get_FurtherDecayCombos();
	for(const auto& locCombosByUsePair : locFurtherDecayCombos)
	{
		if(locCombosByUsePair.first == locToExcludeUse)
			continue; //e.g. when calc'ing missing p4

		for(const auto& locDecayCombo : locCombosByUsePair.second)
		{
			DLorentzVector locDecayP4;
			if(!Calc_P4_Decay(locIsProductionVertex, locReactionChargedCombo, locSourceCombo, locCombosByUsePair.first, locDecayCombo, locVertexZBin, locVertex, locRFBunch, locRFVertexTime, locDecayP4))
				return false;
			locTotalP4 += locDecayP4;
		}
	}

	dFinalStateP4ByCombo_HasMassiveNeutrals.emplace(locP4LookupTuple, locTotalP4);
	return true;
}

bool DSourceComboP4Handler::Cut_InvariantMass_HasMassiveNeutral(bool locIsProductionVertex, const DSourceCombo* locReactionChargedCombo, const DSourceCombo* locSourceCombo, signed char locVertexZBin, Particle_t locDecayPID, double locPrimaryVertexZ, const DVector3& locVertex, double locTimeOffset, vector<int>& locValidRFBunches)
{
	//cuts on possible RF bunches for the massive neutrals
	//if no possible rf bunch yields a massive-neutral-momentum that passes the invariant mass cut, returns an empty vector
	auto locCutIterator = dInvariantMassCuts.find(locDecayPID);
	if(locCutIterator == dInvariantMassCuts.end())
		return true; //no cut to place!!
	auto& locMassCuts = locCutIterator->second;

	//function for calculating and cutting the invariant mass for each rf bunch
	auto CalcAndCut_InvariantMass = [&](int locRFBunch) -> bool
	{
		auto locRFVertexTime = dSourceComboTimeHandler->Calc_PropagatedRFTime(locPrimaryVertexZ, locRFBunch, locTimeOffset);

		DLorentzVector locTotalP4(0.0, 0.0, 0.0, 0.0);
		if(!Calc_P4_HasMassiveNeutrals(locIsProductionVertex, locReactionChargedCombo, locSourceCombo, locVertexZBin, locVertex, locRFBunch, locRFVertexTime, DSourceComboUse(Unknown, 0, nullptr), locTotalP4))
			return true; //can't cut it yet!
		return ((locTotalP4.M() >= locMassCuts.first) && (locTotalP4.M() <= locMassCuts.second));
	};

	//apply the function
	locValidRFBunches.erase(std::remove_if(locValidRFBunches.begin(), locValidRFBunches.end(), CalcAndCut_InvariantMass), locValidRFBunches.end());
	return !locValidRFBunches.empty();
}

bool DSourceComboP4Handler::Cut_InvariantMass_HasMassiveNeutral(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, const DSourceCombo* locReactionChargedCombo, vector<int>& locValidRFBunches)
{
	auto locPrimaryComboUse = dSourceComboer->Get_PrimaryComboUse(locReactionVertexInfo);
	auto locPrimaryVertexZ = dSourceComboVertexer->Get_PrimaryVertex(locReactionVertexInfo, locReactionChargedCombo).Z();

	//loop over vertices in reverse step order //dependency for calculating invariant mass
	auto locStepVertexInfos = DAnalysis::Get_StepVertexInfos_ReverseOrderByStep(locReactionVertexInfo);
	for(auto locStepVertexInfo : locReactionVertexInfo->Get_StepVertexInfos())
	{
		if(locStepVertexInfo->Get_DanglingVertexFlag())
			continue; //unknown position!

		//if this vertex doesn't contain a massive neutral, don't need to do anything special:
		//all p4 cuts have already been done, and any p4 calcs not yet performed will be done if they are needed later
		auto locVertexComboUse = dSourceComboer->Get_SourceComboUse(locStepVertexInfo);
		if(!dSourceComboer->Get_HasMassiveNeutrals(std::get<2>(locVertexComboUse)))
			continue;
		auto locVertexZBin = std::get<1>(locVertexComboUse);

		//see if vertex position is determinable
		auto locIsProductionVertex = locStepVertexInfo->Get_ProductionVertexFlag();
		auto locVertexPrimaryChargedCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionChargedCombo, locStepVertexInfo);
		if(!dSourceComboVertexer->Get_VertexDeterminableWithCharged(locIsProductionVertex, locVertexPrimaryChargedCombo))
			continue; //vertex position indeterminate at this stage: don't include these particles

		//get vertex position and time offset
		auto locVertex = dSourceComboVertexer->Get_Vertex(locIsProductionVertex, locVertexPrimaryChargedCombo);
		auto locTimeOffset = dSourceComboVertexer->Get_TimeOffset(locIsProductionVertex, locReactionChargedCombo, locVertexPrimaryChargedCombo);

		//get all combos at this vertex, with their uses, in reverse dependency order
		auto locVertexPrimaryFullCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionFullCombo, locStepVertexInfo);
		vector<pair<DSourceComboUse, vector<const DSourceCombo*>>> locSourceCombosAndUses_ThisVertex = DAnalysis::Get_SourceCombosAndUses_ThisVertex(locVertexPrimaryFullCombo);
		locSourceCombosAndUses_ThisVertex.emplace(locSourceCombosAndUses_ThisVertex.begin(), locVertexComboUse, locVertexPrimaryFullCombo); //not ideal ...

		//loop over combos & uses at this vertex in dependency order (in reverse)
		for(auto locIterator = locSourceCombosAndUses_ThisVertex.rbegin(); locIterator != locSourceCombosAndUses_ThisVertex.rend(); ++locIterator)
		{
			auto locDecayComboUse = locIterator->first;
			auto locDecayPID = std::get<0>(locDecayComboUse);
			if(locDecayPID == Unknown)
				continue; //no mass cut to place!
			if(!dSourceComboer->Get_HasMassiveNeutrals(std::get<2>(locDecayComboUse)))
				continue; //already cut!

			//loop over combos
			for(const auto& locDecayCombo : locIterator->second)
			{
				if(!Cut_InvariantMass_HasMassiveNeutral(locIsProductionVertex, locReactionChargedCombo, locDecayCombo, locVertexZBin, locDecayPID, locPrimaryVertexZ, locVertex, locTimeOffset, locValidRFBunches))
					return false; //failed mass cut for all possible rf bunches!
			}
		}
	}

	return true;
}

}

#endif // DSourceComboP4Handler_h
