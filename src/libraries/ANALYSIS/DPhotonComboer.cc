#include "ANALYSIS/DPhotonComboer.h"


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

	/**************************************************** PHOTON-RF DELTA-T CUTS ***************************************************
	*
	* We would also like to place photon timing cuts in advance as well.  However, these also depend on the vertex position.
	*
	* The error on the delta-t (delta_delta_t) is derived in a separate section below. The result is:
	* All: delta_delta_t = dr*[1/sin(theta) - sqrt(1 + (1/tan(theta) - z_error/dr)^2)]/c - z_error/c
	* BCAL: delta_delta_t = 65*[1/sin(theta) - sqrt(1 + (1/tan(theta) - z_error/65)^2)]/c - z_error/c
	* FCAL: delta_delta_t = (~650 - z_error)*[1/cos(theta) - sqrt(tan^2(theta) + (1 - z_error/(~650 - z_error))^2)]/c - z_error/c
	*
	* For the FCAL, at a z_error of 30cm (center of target + detached-z), the maximum error in delta_delta_t is 22ps
	* At a z_error of 5cm, delta_delta_t is 3.5ps
	*
	* For the BCAL (dr = 65cm), at a z_error of 30cm (center of target + detached-z), the maximum error in delta_delta_t is 1.8ns (at 140 degrees)
	* For a z_error of 5cm, delta_delta_t is 300ps
	*
	* So, since we are evaluating the photons at z-vertex steps of 10cm anyway (for invariant mass), do the same for the photons.
	* Then, when performing delta-t cuts, increase the width of the cut by the max time error.
	* This error is strongly theta-dependent, so we can calculate the max error on a shower-by-shower basis using the above equations.
	*
	*******************************************************************************************************************************/

	/*********************************************** PHOTON-RF DELTA-T CUT DERIVATION **********************************************
	*
	* The error on the photon time difference (delta_delta_t) can be found as:
	* delta_t = t_prop_shower - t_beam
	* delta_delta_t = delta_t_actual - delta_t_guess
	* delta_delta_t = (t_prop_shower_actual - t_beam_actual) - (t_prop_shower_guess - t_beam_guess)
	* delta_delta_t = (t_prop_shower_actual - t_prop_shower_guess) - (t_beam_actual - t_beam_guess)
	*
	* t_prop_shower = t_shower - path_shower/c
	* delta_delta_t = ((t_shower - path_shower_actual/c) - (t_shower - path_shower_guess/c)) - (t_beam_actual - t_beam_guess)
	* delta_delta_t = (path_shower_guess - path_shower_actual)/c - (t_beam_actual - t_beam_guess)
	*
	* t_beam = t_RF_targcenter + (vertz - targz)/c
	* delta_delta_t = (path_shower_guess - path_shower_actual)/c - ((t_RF_targcenter + (vertz_actual - targz)/c) - (t_RF_targcenter + (vertz_guess - targz)/c))
	* delta_delta_t = [path_shower_guess - path_shower_actual - vertz_actual + vertz_guess]/c
	*
	* define z_error = vertz_actual - vertz_guess
	* delta_delta_t = [path_shower_guess - path_shower_actual - z_error]/c
	*
	* path_shower_guess = shower_pos.Perp()/sin(theta)
	* path_shower_actual = sqrt(shower_pos.Perp()^2 + path_z_shower_actual^2)
	*
	* path_z_shower_actual = shower_z - vertz_actual, vertz_actual = z_error + vertz_guess
	* path_z_shower_actual = shower_z - (z_error + vertz_guess)
	* path_z_shower_guess = shower_z - vertz_guess
	* path_z_shower_actual = path_z_shower_guess - z_error
	* So, path_shower_actual = sqrt(shower_pos.Perp()^2 + (path_z_shower_guess - z_error)^2)
	*
	* path_z_shower_guess = shower_pos.Perp()/tan(theta)
	* So, path_shower_actual = sqrt(shower_pos.Perp()^2 + (shower_pos.Perp()/tan(theta) - z_error)^2)
	*
	* Using dr for shower_pos.Perp() and reducing gives:
	* path_shower_actual = dr*sqrt(1 + (1/tan(theta) - z_error/dr)^2)
	* and path_shower_guess = dr/sin(theta)
	*
	* So: delta_delta_t = dr*[1/sin(theta) - sqrt(1 + (1/tan(theta) - z_error/dr)^2)]/c - z_error/c
	*
	* For the FCAL, dr = dz*tan(theta), dz = 650 - z_error (approx 650)
	* So: delta_delta_t = (650 - z_error)*[1/cos(theta) - sqrt(tan^2(theta) + (1 - z_error/(650 - z_error))^2)]/c - z_error/c
	*
	* For the FCAL, the delta_delta_t is at most 23ps when the z_error is 30cm (center of target + detached vertex) (12 degrees)
	* Therefore, the z_error is irrelevant: Just choose the center of the target and increase the width of the delta_t cut as needed.
	*
	* However, for the BCAL the time error is rather large (largest at large angles (e.g. 140 degrees))
	* Even with a z_error of 5cm, the delta_delta_t is 300ps.  Therefore use 10-cm-wide vertex-z bins and increase the delta_t cut width as needed.
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
	*******************************************************************************************************************************/

void DVertexCreator::Register_ProductionPhotons(const DReaction* locReaction, const vector<int>& locStepIndices)
{
	for(auto& locStep : locReactionSteps)
	{
		//if any decay products are charged or missing, continue
		Particle_t locMissingPID;
		if(!locStep->Get_MissingPID(locMissingPID))
			continue;

		deque<Particle_t> locPIDs;
		locStep->Get_NonMissingFinalChargedPIDs(locPIDs);
		if(!locPIDs.empty())
			continue;

		locStep->Get_NonMissingFinalNeutralPIDs(locPIDs)
		DParticleDecay locPair(locStep->Get_InitialParticleID(), );
	}
	using DParticleDecay = pair<Particle_t, map<Particle_t, int> >; //if first PID is unknown then is direct #photons (not a decay)
	auto dProductionVertexNeutralsSet = set<pair<DParticleDecay, int> >(DParticleDecayComparer);
}
