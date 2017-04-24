#include "ANALYSIS/DSourceComboer.h"

namespace DAnalysis
{

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
* We would also like to place photon timing cuts in advance as well (these narrow down the possible RF bunch).
* However, these also depend on the vertex position.
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
* Because we are comboing as channel-independent as possible, when doing a photon combo we don't know if it's at a detached vertex or not.
* So, we must also include a max timing offset when doing time cuts (to select RF bunches) due for detached vertices.
* This time offset is also derived below
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
* Now, for the max time offset due to detached vertices:
* Worst case: Xi- -> Lambda -> pi0, n: ~30cm, say avg beta = 1/3 (beta can't be too small or ~30cm distance not likely!)
* max_time_offset = delta_t_neutral - delta_t_rf
* delta_t_neutral = delta_x/(beta*c) = 30/(30*1/3) = 3 ns
* delta_t_rf = delta_x/(beta*c) = 30/(1*30) = 1ns
* max_time_offset = 3ns - 1ns = 2ns
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


/*************************************************** CHARGED TRACK TIMING CUTS *************************************************
*
* Charged time cuts are dependent on combo vertex, especially for low-theta tracks.
* Wherever the combo vertex is, the track won't pass through it until after the kinfit, so do "final" PID cuts at the end
*
* Once we have a vertex for the combo, compute POCA to the vertex and do pre-kinfit time cuts there.
* This will be pretty accurate, except for slow-single-track channels at low-theta, for which there's nothing you can do anyway.
*
* We don't want to wait until we have a combo vertex to do some PID timing cuts, but unfortunately we have to.
* Since computing neutral timing every 10cm, we can try to do the same for charged tracks as well, but this doesn't work.
*
* The maximum error associated with this is:
*
* delta_t = t_prop_track - t_beam
* delta_delta_t = delta_t_actual - delta_t_guess
* delta_delta_t = (t_prop_track_actual - t_beam_actual) - (t_prop_track_guess - t_beam_guess)
* delta_delta_t = (t_prop_track_actual - t_prop_track_guess) - (t_beam_actual - t_beam_guess)
*
* t_prop_track = t_track - path_track/(beta_track*c)
* delta_delta_t = ((t_track - path_track_actual/(beta_track*c)) - (t_track - path_track_guess/(beta_track*c))) - (t_beam_actual - t_beam_guess)
* delta_delta_t = (path_track_guess - path_track_actual)/(beta_track*c) - (t_beam_actual - t_beam_guess)
*
* t_beam = t_RF_targcenter + (vertz - targz)/c
* delta_delta_t = (path_track_guess - path_track_actual)/(beta_track*c) - ((t_RF_targcenter + (vertz_actual - targz)/c) - (t_RF_targcenter + (vertz_guess - targz)/c))
* delta_delta_t = (path_track_guess - path_track_actual)/(beta_track*c) + (vertz_guess - vertz_actual)/c
*
* define z_error = vertz_actual - vertz_guess
* delta_delta_t = (path_track_guess - path_track_actual)/(beta_track*c) - z_error/c
*
* From here, assume track is straight over the distance z_error/2:
*
* FCAL:
* path_track_guess = path_z_guess/cos(theta)
* path_track_actual = path_z_actual/cos(theta), path_z_actual = path_z_guess - z_error
* path_track_guess - path_track_actual = path_z_guess/cos(theta) - (path_z_guess - z_error)/cos(theta) = z_error/cos(theta)
* delta_delta_t = z_error/(cos(theta)*beta_track*c) - z_error/c
* delta_delta_t = (z_error/c) * [1/(cos(theta)*beta_track) - 1]
*
* BCAL:
* path_track_guess = path_r/sin(theta)
* path_track_actual = sqrt(path_z_actual*path_z_actual + path_r*path_r)
* path_z_actual = path_z_guess - z_error
* path_z_guess = path_r/tan(theta)
* path_z_actual = path_r/tan(theta) - z_error
* path_track_actual = sqrt((path_r/tan(theta) - z_error)^2 + path_r*path_r)
* path_track_actual = path_r*sqrt((1/tan(theta) - z_error/path_r)^2 + 1)
* delta_delta_t = path_r*(1/sin(theta) - sqrt((1/tan(theta) - z_error/path_r)^2 + 1))/(beta_track*c) - z_error/c
*
* These errors are too large:
* For slow tracks the errors are huge, and for fast tracks the errors are small.
* However, for fast tracks the timing isn't good enough to tell one PID from another anyway.
* So this does not gain much.
*
* Instead, charged track timing cuts cannot be placed until the vertex position is found.
*
*******************************************************************************************************************************/

/********************************************************************* CONSTRUCTOR **********************************************************************/

DSourceComboer::DSourceComboer(JEventLoop* locEventLoop)
{
	//GET THE GEOMETRY
	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry* locGeometry = locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());

	//BEAM BUNCH PERIOD
	vector<double> locBeamPeriodVector;
	locEventLoop->GetCalib("PHOTON_BEAM/RF/beam_period", locBeamPeriodVector);
	dBeamBunchPeriod = locBeamPeriodVector[0];

	//TARGET INFORMATION
	double locTargetCenterZ = 65.0;
	locGeometry->GetTargetZ(locTargetCenterZ);
	dTargetCenter.SetXYZ(0.0, 0.0, locTargetCenterZ);
	locGeometry->GetTargetLength(dTargetLength);

	//Get preselect tag
	gPARMS->SetDefaultParameter("COMBO:SHOWER_SELECT_TAG", dShowerSelectionTag);

	dBCALKinematics.resize(dNumPhotonVertexZBins);

	//DEFINE LOOSE CUTS
	Define_LooseCuts();

	//CREATE DSourceComboINFO'S
	vector<const DReactionVertexInfo*> locVertexInfos;
	locEventLoop->Get(locVertexInfos);
	for(auto locVertexInfo : locVertexInfos)
		Create_SourceComboInfos(locVertexInfo);

	//TRANSFER INFOS FROM SET TO VECTOR
	dSourceComboInfos.reserve(dSourceComboInfoSet.size());
	std::copy(dSourceComboInfoSet.begin(), dSourceComboInfoSet.end(), std::back_inserter(dSourceComboInfos));
	dSourceComboInfoSet.clear(); //free up the memory
}

void DSourceComboer::Define_LooseCuts(void)
{
	//Create photon/RF delta-t cuts
	dPhotonTimeCutMap.emplace(SYS_BCAL, new TF1("df_BCALShowerRFTimeComboingCut", "[0]", 0.0, 12.0));
	dPhotonTimeCutMap[SYS_BCAL]->SetParameter(0, 2.0);
	dPhotonTimeCutMap.emplace(SYS_FCAL, new TF1("df_FCALShowerRFTimeComboingCut", "[0]", 0.0, 12.0));
	dPhotonTimeCutMap[SYS_FCAL]->SetParameter(0, 2.0);

	// Timing Cuts: Photon
	dPIDTimingCuts[Gamma][SYS_BCAL] = 3.0;
	dPIDTimingCuts[Gamma][SYS_FCAL] = 2.5;

	// Timing Cuts: Leptons
	dPIDTimingCuts[Electron][SYS_BCAL] = 1.0;
	dPIDTimingCuts[Electron][SYS_FCAL] = 2.5;
	dPIDTimingCuts[Electron][SYS_TOF] = 2.5;
	dPIDTimingCuts[Positron] = dPIDTimingCuts[Electron];
	dPIDTimingCuts[MuonMinus] = dPIDTimingCuts[Electron];
	dPIDTimingCuts[MuonPlus] = dPIDTimingCuts[Electron];

	// Timing Cuts: Mesons
	dPIDTimingCuts[PiPlus][SYS_BCAL] = 2.0;
	dPIDTimingCuts[PiPlus][SYS_FCAL] = 2.5;
	dPIDTimingCuts[PiPlus][SYS_TOF] = 2.5;

	dPIDTimingCuts[PiMinus][SYS_BCAL] = 2.0;
	dPIDTimingCuts[PiMinus][SYS_FCAL] = 2.5;
	dPIDTimingCuts[PiMinus][SYS_TOF] = 2.5;

	dPIDTimingCuts[KPlus][SYS_BCAL] = 0.75;
	dPIDTimingCuts[KPlus][SYS_FCAL] = 2.5;
	dPIDTimingCuts[KPlus][SYS_TOF] = 2.0;
	dPIDTimingCuts[KMinus] = dPIDTimingCuts[KPlus];

	// Timing Cuts: Baryons
	dPIDTimingCuts[Proton][SYS_BCAL] = 2.5;
	dPIDTimingCuts[Proton][SYS_FCAL] = 2.5;
	dPIDTimingCuts[Proton][SYS_TOF] = 2.5;

	dPIDTimingCuts[AntiProton] = dPIDTimingCuts[Proton];
}

/******************************************************************* CREATE DSOURCOMBOINFO'S ********************************************************************/

void DSourceComboer::Create_SourceComboInfos(const DReactionVertexInfo* locReactionVertexInfo)
{
	//FULL combo use: Segregate each step into (up to) 3 combos: a fully charged, a fully neutral, and a mixed
	//that way we will combo each separately before combining them horizontally: maximum re-use, especially of time-intensive neutral comboing

	//We will register what steps these combos are created for
	map<size_t, DSourceComboUse> locStepComboUseMap; //size_t = step index

	//loop over steps in reverse order
	auto locReaction = locReactionVertexInfo->Get_Reaction();
	auto locReactionSteps = locReaction->Get_ReactionSteps();
	for(auto locStepIterator = locReactionSteps.rbegin(); locStepIterator != locReactionSteps.rend(); ++locStepIterator)
	{
		auto locStep = *locStepIterator;
		auto locStepIndex = locReaction->Get_NumReactionSteps() - std::distance(locReactionSteps.rbegin(), locStepIterator) - 1;

		//create combo uses for all charged, all neutral, then for any mixed decays
		map<Particle_t, unsigned char> locChargedParticleMap = Build_ParticleMap(locReaction, locStepIndex, d_Charged);
		map<Particle_t, unsigned char> locNeutralParticleMap = Build_ParticleMap(locReaction, locStepIndex, d_Neutral);

		//get combo infos for final-state decaying particles //if not present, ignore parent
		auto locFinalStateDecayingComboUsesPair = Get_FinalStateDecayingComboUses(locReaction, locStepIndex, locStepComboUseMap);
		auto locIncludeParentFlag = locFinalStateDecayingComboUsesPair.first;
		auto& locFurtherDecays = locFinalStateDecayingComboUsesPair.second;

		//split up further-decays into all-charged, all-neutral, and mixed
		map<DSourceComboUse, unsigned char> locFurtherDecays_Charged, locFurtherDecays_Neutral, locFurtherDecays_Mixed;
		for(const auto& locDecayPair : locFurtherDecays)
		{
			auto locChargeContent = dComboInfoChargeContent[std::get<2>(locDecayPair.first)];
			if(locChargeContent == d_Charged)
				locFurtherDecays_Charged.emplace(locDecayPair);
			else if(locChargeContent == d_Neutral)
				locFurtherDecays_Neutral.emplace(locDecayPair);
			else
				locFurtherDecays_Mixed.emplace(locDecayPair);
		}

		//determine whether to include the decay itself in the comboing (or just the products)
		//only include if can make an invariant mass cut (what it's used for here)
		//we will still group these separately from the rest of the particles
		if((locStepIndex != 0) || !DAnalysis::Get_IsFirstStepBeam(locReaction)) //decay
		{
			//ignore parent if products include missing particles
			if(DAnalysis::Check_IfMissingDecayProduct(locReaction, locStepIndex))
				locIncludeParentFlag = false;
		}
		else //direct production
			locIncludeParentFlag = false;

		//create combo uses for each case
		Particle_t locInitPID = locIncludeParentFlag ? locStep->Get_InitialPID() : Unknown;
		bool locNoChargedFlag = (locChargedParticleMap.empty() && locFurtherDecays_Charged.empty());
		bool locNoNeutralFlag = (locNeutralParticleMap.empty() && locFurtherDecays_Neutral.empty());

		DSourceComboUse locPrimaryComboUse(Unknown, DSourceComboInfo::Get_VertexZIndex_Unknown(), nullptr);
		if(locNoChargedFlag && locNoNeutralFlag) //only mixed
			locPrimaryComboUse = Make_ComboUse(locInitPID, {}, locFurtherDecays_Mixed);
		else if(locNoNeutralFlag && locFurtherDecays_Mixed.empty()) //only charged
			locPrimaryComboUse = Make_ComboUse(locInitPID, locChargedParticleMap, locFurtherDecays_Charged);
		else if(locNoChargedFlag && locFurtherDecays_Mixed.empty()) //only neutral
			locPrimaryComboUse = Make_ComboUse(locInitPID, locNeutralParticleMap, locFurtherDecays_Neutral);
		else //some combination
		{
			auto locFurtherDecays_All = locFurtherDecays_Mixed;
			//create a combo for each charged group, with init pid = unknown
			if(!locNoChargedFlag)
			{
				auto locComboUse_Charged = Make_ComboUse(Unknown, locChargedParticleMap, locFurtherDecays_Charged);
				locFurtherDecays_All.emplace(locComboUse_Charged, 1);
			}
			if(!locNoNeutralFlag)
			{
				auto locComboUse_Neutral = Make_ComboUse(Unknown, locNeutralParticleMap, locFurtherDecays_Neutral);
				locFurtherDecays_All.emplace(locComboUse_Neutral, 1);
			}

			locPrimaryComboUse = Make_ComboUse(locInitPID, {}, locFurtherDecays_All);
		}

		locStepComboUseMap.emplace(locStepIndex, locPrimaryComboUse);
	}

	//Register the results!!
	for(const auto& locStepVertexInfo : locReactionVertexInfo->Get_StepVertexInfos())
		dSourceComboUseReactionMap.emplace(locStepVertexInfo, locStepComboUseMap[locStepVertexInfo->Get_StepIndices().front()]);
	for(const auto& locUseStepPair : locStepComboUseMap)
		dSourceComboInfoStepMap.emplace(std::make_pair(locReactionVertexInfo->Get_StepVertexInfo(locUseStepPair.first), locUseStepPair.second), locUseStepPair.first);
	dSourceComboUseReactionStepMap.emplace(locReaction, locStepComboUseMap);

	dSourceComboUseReactionMap_Primary[locReactionVertexInfo] = locStepComboUseMap[0];
}

DSourceComboUse DSourceComboer::Create_ZDependentSourceComboUses(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionChargedCombo)
{
	//this creates new uses, with the specific vertex-z bins needed
	//note that the use can have a different structure from the charged!! (although not likely)
	//E.g. if something crazy like 2 KShorts -> 3pi, each at a different vertex-z bin, then they will no longer be grouped together vertically (separate uses: horizontally instead)

	//see if they've already been created.  if so, just return it.
	auto locVertexZBins = dSourceComboVertexer.Get_VertexZBins(locReactionVertexInfo, locReactionChargedCombo);
	auto locCreationPair = std::make_pair(locReactionVertexInfo, locVertexZBins);
	auto locUseIterator = dSourceComboUseVertexZMap.find(locCreationPair);
	if(locUseIterator != dSourceComboUseVertexZMap.end())
		return locUseIterator->second; //already created! we are done

	auto locReaction = locReactionVertexInfo->Get_Reaction();
	auto locStepVertexInfos = locReactionVertexInfo->Get_StepVertexInfos();

	//sort vertex infos in reverse-step order
	auto Comparator_ReverseOrderByStep = [](const DReactionStepVertexInfo* lhs, const DReactionStepVertexInfo* rhs) -> bool
		{return lhs->Get_StepIndices().front() > rhs->Get_StepIndices().front();} // >: reverse order
	std::sort(locStepVertexInfos.begin(), locStepVertexInfos.end(), Comparator_ReverseOrderByStep);

	//loop over vertex infos in reverse-step order
	unordered_map<size_t, DSourceComboUse> locCreatedUseMap; //size_t: step index
	for(auto locStepVertexInfo : locStepVertexInfos)
	{
		//handle indeterminate vertex case
		if(!dSourceComboVertexer.Get_VertexDeterminableWithCharged(locStepVertexInfo))
		{
			//no different than before: vertex bin remains unknown
			//however, it's possible that decay products down the road have a determined vertex-z
			//so, we can't just automatically reuse the infos/uses
		}

		//for this vertex, get the vertex z bin
		auto locVertexZBin = dSourceComboVertexer.Get_VertexZBin(locReactionChargedCombo, locStepVertexInfo);

		//loop over the steps at this vertex z bin, in reverse order
		auto locStepIndices = locStepVertexInfo->Get_StepIndices();
		for(auto locStepIterator = locStepIndices.rbegin(); locStepIterator != locStepIndices.rend(); ++locStepIterator)
		{
			auto locStepIndex = *locStepIterator;
			auto locStepOrigUse = dSourceComboUseReactionStepMap[locReaction][locStepIndex];

			//build new use for the further decays, setting the vertex z-bins
			auto locNewComboUse = Build_NewZDependentUse(locReaction, locStepIndex, locVertexZBin, locStepOrigUse, locCreatedUseMap);
			locCreatedUseMap.emplace(locStepIndex, locNewComboUse);
		}
	}

	dSourceComboUseVertexZMap.emplace(locCreationPair, locCreatedUseMap[0]);
	return locCreatedUseMap[0];
}

DSourceComboUse DSourceComboer::Build_NewZDependentUse(const DReaction* locReaction, size_t locStepIndex, signed char locVertexZBin, const DSourceComboUse& locOrigUse, const unordered_map<size_t, DSourceComboUse>& locCreatedUseMap)
{
	//each step can be broken up into combo infos with a depth of 2 (grouping charges separately)
	auto locStep = locReaction->Get_ReactionStep(locStepIndex);
	auto locOrigInfo = std::get<2>(locOrigUse);
	if(dComboInfoChargeContent[locOrigInfo] == d_Charged)
	{
		dZDependentUseToIndependentMap.emplace(locOrigUse, locOrigUse);
		return locOrigUse; //no need to change!: no neutrals anyway
	}

	map<DSourceComboUse, unsigned char> locNewFurtherDecays;
	auto locOrigFurtherDecays = locOrigInfo->Get_FurtherDecays();
	for(const auto& locDecayPair : locOrigFurtherDecays)
	{
		const auto& locOrigDecayUse = locDecayPair.first;
		auto locDecayPID = std::get<0>(locOrigDecayUse);
		if(locDecayPID != Unknown)
		{
			//these decays are represented by other steps, and have already been saved
			for(unsigned char locInstance = 1; locInstance <= locDecayPair.second; ++locInstance)
			{
				auto locParticleIndex = DAnalysis::Get_ParticleIndex(locStep, locDecayPID, locInstance);
				auto locDecayStepIndex = DAnalysis::Get_DecayStepIndex(locReaction, locStepIndex, locParticleIndex);
				const auto& locSavedDecayUse = locCreatedUseMap.find(locDecayStepIndex)->second; //is same as locOrigDecayUse, except different zbins along chain

				//save the use for this decay
				auto locUseIterator = locNewFurtherDecays.find(locSavedDecayUse);
				if(locUseIterator != locNewFurtherDecays.end())
					++(locUseIterator->second);
				else
					locNewFurtherDecays.emplace(locSavedDecayUse, 1);
			}
		}
		else //is unknown (and guaranteed to be size 1 since has unknown parent)
		{
			//must dig down, but only one level: their decays must terminate at new steps (or end)
			auto locNewComboUse = Build_NewZDependentUse(locReaction, locStepIndex, locVertexZBin, locOrigDecayUse, locCreatedUseMap);
			//save the use for this decay
			auto locUseIterator = locNewFurtherDecays.find(locNewComboUse);
			if(locUseIterator != locNewFurtherDecays.end())
				++(locUseIterator->second);
			else
				locNewFurtherDecays.emplace(locNewComboUse, 1);
		}
	}

	//build and save new info, use, and return
	auto locNewComboInfo = locNewFurtherDecays.empty() ? locOrigInfo : GetOrMake_SourceComboInfo(locOrigInfo->Get_NumParticles(), locNewFurtherDecays);
	DSourceComboUse locNewComboUse(std::get<0>(locOrigUse), locVertexZBin, locNewComboInfo);
	dZDependentUseToIndependentMap.emplace(locNewComboUse, locOrigUse);
	return locNewComboUse;
}

pair<bool, map<DSourceComboUse, unsigned char>> DSourceComboer::Get_FinalStateDecayingComboUses(const DReaction* locReaction, size_t locStepIndex, const map<size_t, DSourceComboUse>& locStepComboUseMap) const
{
	//get combo infos for final-state decaying particles //if one is not present, ignore parent
	auto locIncludeParentFlag = true; //unless changed below
	map<DSourceComboUse, unsigned char> locFurtherDecays;
	auto locStep = locReaction->Get_ReactionStep(locStepIndex);
	for(size_t loc_i = 0; loc_i < locStep->Get_NumFinalPIDs(); ++loc_i)
	{
		int locDecayStepIndex = DAnalysis::Get_DecayStepIndex(locReaction, locStepIndex, loc_i);
		if(locDecayStepIndex < 0)
			continue;
		auto locUseIterator = locStepComboUseMap.find(size_t(locDecayStepIndex));
		if(locUseIterator == locStepComboUseMap.end())
			locIncludeParentFlag = false;
		else
		{
			//save decay
			auto& locSourceComboUse = locUseIterator->second;
			auto locDecayIterator = locFurtherDecays.find(locSourceComboUse);
			if(locDecayIterator == locFurtherDecays.end())
				locFurtherDecays.emplace(locSourceComboUse, 1);
			else
				++(locDecayIterator->second);
		}
	}

	return std::make_pair(locIncludeParentFlag, locFurtherDecays);
}

map<Particle_t, unsigned char> DSourceComboer::Build_ParticleMap(const DReaction* locReaction, size_t locStepIndex, Charge_t locCharge) const
{
	//build map of charged particles
	map<Particle_t, unsigned char> locNumParticles;
	auto locParticles = locReaction->Get_FinalPIDs(locStepIndex, false, false, locCharge, true); //no missing or decaying, include duplicates
	for(const auto& locPID : locParticles)
	{
		auto locPIDIterator = locNumParticles.find(locPID);
		if(locPIDIterator != locNumParticles.end())
			++(locPIDIterator->second);
		else
			locNumParticles.emplace(locPID, 1);
	}

	return locNumParticles;
}

DSourceComboUse DSourceComboer::Make_ComboUse(Particle_t locInitPID, const map<Particle_t, unsigned char>& locNumParticles, const map<DSourceComboUse, unsigned char>& locFurtherDecays)
{
	//convert locFurtherDecays map to a vector
	vector<pair<DSourceComboUse, unsigned char>> locDecayVector;
	locDecayVector.reserve(locFurtherDecays.size());
	std::copy(locFurtherDecays.begin(), locFurtherDecays.end(), std::back_inserter(locDecayVector));

	//convert locNumParticles map to a vector
	vector<pair<Particle_t, unsigned char>> locParticleVector;
	locParticleVector.reserve(locNumParticles.size());
	std::copy(locNumParticles.begin(), locNumParticles.end(), std::back_inserter(locParticleVector));

	//make or get the combo info
	auto locComboInfo = MakeOrGet_SourceComboInfo(locParticleVector, locDecayVector);
	return DSourceComboUse(locInitPID, DSourceComboInfo::Get_VertexZIndex_Unknown(), locComboInfo);
}

const DSourceComboInfo* DSourceComboer::MakeOrGet_SourceComboInfo(const vector<pair<Particle_t, unsigned char>>& locNumParticles, const vector<pair<DSourceComboUse, unsigned char>>& locFurtherDecays)
{
	//to be called (indirectly) by constructor: during the stage when primarily making
	//create the object on the stack
	DSourceComboInfo locSearchForInfo(locNumParticles, locFurtherDecays);

	//then search through the set to retrieve the pointer to the corresponding object if it already exists
	auto locInfoIterator = dSourceComboInfoSet.find(&locSearchForInfo);
	if(locInfoIterator != dSourceComboInfoSet.end())
		return *locInfoIterator; //it exists: return it

	//doesn't exist, make it and insert it into the sorted vector in the correct spot
	auto locComboInfo = new DSourceComboInfo(locNumParticles, locFurtherDecays);
	dSourceComboInfoSet.insert(locComboInfo);
	dComboInfoChargeContent.emplace(locComboInfo, DAnalysis::Get_ChargeContent(locComboInfo));
	if(DAnalysis::Get_HasMassiveNeutrals(locComboInfo))
		dComboInfosWithMassiveNeutrals.insert(locComboInfo);
	return locComboInfo;
}

const DSourceComboInfo* DSourceComboer::GetOrMake_SourceComboInfo(const vector<pair<Particle_t, unsigned char>>& locNumParticles, const vector<pair<DSourceComboUse, unsigned char>>& locFurtherDecays)
{
	//to be called when making combos: during the stage when primarily getting
	//create the object on the stack
	DSourceComboInfo locSearchForInfo(locNumParticles, locFurtherDecays);

	//then search through the vector to retrieve the pointer to the corresponding object if it already exists
	auto locIteratorPair = std::equal_range(dSourceComboInfos.begin(), dSourceComboInfos.end(), locSearchForInfo, DCompare_SourceComboInfos());
	if(locIteratorPair.first != locIteratorPair.second)
		return *(locIteratorPair.first); //it exists: return it

	//doesn't exist, make it and insert it into the sorted vector in the correct spot
	auto locComboInfo = new DSourceComboInfo(locNumParticles, locFurtherDecays);
	dSourceComboInfos.emplace(locIteratorPair.first, locComboInfo);
	dComboInfoChargeContent.emplace(locComboInfo, DAnalysis::Get_ChargeContent(locComboInfo));
	if(DAnalysis::Get_HasMassiveNeutrals(locComboInfo))
		dComboInfosWithMassiveNeutrals.insert(locComboInfo);
	return locComboInfo;
}

/********************************************************************** SETUP FOR NEW EVENT ***********************************************************************/

void DSourceComboer::Reset_NewEvent(JEventLoop* locEventLoop)
{
	//check if it's actually a new event
	auto locEventNumber = locEventLoop->GetJEvent().GetEventNumber();
	if(locEventNumber == dEventNumber)
		return; //nope
	dEventNumber = locEventNumber;

	//CLEAR OLD DATA
	dBCALShowers.clear();
	dFCALShowers.clear();
	dFCALKinematics.clear();
	for(auto& locMap : dBCALKinematics)
		locMap.clear();

	//SETUP NEUTRAL SHOWERS
	Setup_NeutralShowers(locEventLoop);

//CLEAR MORE THINGS!!
	//RECYCLE COMBO & VECTOR POINTERS
	//be careful! don't recycle combos with a use pid != unknown, because they are just copies! not unique pointers!

	dSourceComboP4Handler.Reset();
	dSourceComboVertexer.Reset();
}

void DSourceComboer::Setup_NeutralShowers(JEventLoop* locEventLoop)
{
	//Precompute a few things for the neutral showers, before comboing
	//Even if it turns out some of this isn't technically needed,
	//it's still faster than doing a check to see if this has been done or not for every single photon-combo request

	//GET RF BUNCH
	locEventLoop->GetSingle(dInitialEventRFBunch);

	//GET NEUTRAL SHOWERS
	vector<const DNeutralShower*> locNeutralShowers;
	locEventLoop->Get(locNeutralShowers);

	//ARRANGE NEUTRAL SHOWERS
	for(auto& locShower : locNeutralShowers)
	{
		auto& locContainer = (locShower->dDetectorSystem == SYS_BCAL) ? dBCALShowers : dFCALShowers;
		locContainer.push_back(locShower);
	}

	//CALCULATE KINEMATICS
	//FCAL: at target center
	for(auto& locShower : dFCALShowers)
		dFCALKinematics.emplace(locShower, Create_KinematicData(locShower, dTargetCenter));
	//BCAL: in vertex-z bins
	for(size_t loc_i = 0; loc_i < dNumPhotonVertexZBins; ++loc_i)
	{
		DVector3 locBinCenter(0.0, 0.0, Get_PhotonVertexZBinCenter(loc_i));
		for(auto& locShower : dBCALShowers)
			dBCALKinematics[loc_i].emplace(locShower, Create_KinematicData(locShower, locBinCenter));
	}

	//DETERMINE WHICH RF BUNCHES ARE VALID
	//FCAL: at target center
	for(auto& locShower : dFCALShowers)
		Calc_PhotonBeamBunchShifts(locShower, dFCALKinematics[locShower], dInitialEventRFBunch->dTime, dFCALPhotonShowersByBeamBunch);
	//BCAL + FCAL: in vertex-z bins
	for(size_t loc_i = 0; loc_i < dNumPhotonVertexZBins; ++loc_i)
	{
		const auto& locAllPhotonsByBeamBunch = dPhotonShowersByBeamBunch[loc_i];

		//propagate RF time to vertex position
		double locPropagatedRFTime = dInitialEventRFBunch->dTime + (Get_PhotonVertexZBinCenter(loc_i) - dTargetCenter.Z())/29.9792458;
		for(auto& locShower : dBCALShowers)
			Calc_PhotonBeamBunchShifts(locShower, dBCALKinematics[loc_i][locShower], locPropagatedRFTime, locAllPhotonsByBeamBunch);

		//insert the previously-done FCAL photons
		for(auto locBeamBunchPair : dFCALPhotonShowersByBeamBunch)
		{
			const auto& locShowers = locBeamBunchPair.second;
			auto locBothIterator = locAllPhotonsByBeamBunch.find(locBeamBunchPair.first);
			if(locBothIterator != locAllPhotonsByBeamBunch.end())
			{
				auto locPhotonVector = locBothIterator->second;
				locPhotonVector.insert(locPhotonVector.end(), locShowers.begin(), locShowers.end());
			}
			else
				locAllPhotonsByBeamBunch.emplace(locBeamBunchPair);
		}
	}
}

void DSourceComboer::Calc_PhotonBeamBunchShifts(const DNeutralShower* locNeutralShower, shared_ptr<const DKinematicData>& locKinematicData,
		double locRFTime, DPhotonShowersByBeamBunch& locShowersByBeamBunch) const
{
	//get delta-t cut
	DetectorSystem_t locSystem = locNeutralShower->dDetectorSystem;
	double locDeltaTCut = dPhotonTimeCutMap[locSystem]->Eval(locNeutralShower->dEnergy) + Calc_MaxDeltaTError(locNeutralShower, locKinematicData);

	//prepare for loop over possible #-RF-shifts
	double locVertexTime = locKinematicData->time();
	int locOrigNumShifts = Calc_RFBunchShift(locRFTime, locVertexTime); //get best shift

	//start with best-shift, then loop up until fails cut
	int locNumShifts = locOrigNumShifts;
	double locDeltaT = locVertexTime - (locRFTime + locNumShifts*dBeamBunchPeriod);
	while(fabs(locDeltaT) < locDeltaTCut)
	{
		locShowersByBeamBunch[locNumShifts].push_back(locNeutralShower);
		++locNumShifts;
		locDeltaT = locVertexTime - (locRFTime + locNumShifts*dBeamBunchPeriod);
	}

	//now loop down in n-shifts
	locNumShifts = locOrigNumShifts - 1;
	locDeltaT = locVertexTime - (locRFTime + locNumShifts*dBeamBunchPeriod);
	while(fabs(locDeltaT) < locDeltaTCut)
	{
		locShowersByBeamBunch[locNumShifts].push_back(locNeutralShower);
		--locNumShifts;
		locDeltaT = locVertexTime - (locRFTime + locNumShifts*dBeamBunchPeriod);
	}

	//this means we may need to accept EARLIER RF bunches
	//continue down-shift loop, this time including time offset
	locDeltaT = locVertexTime - (locRFTime + locNumShifts*dBeamBunchPeriod + dMaxDecayTimeOffset); //+dMaxTimeOffset: takes longer for "RF" time to get there (due to slow decaying particle)
	while(fabs(locDeltaT) < locDeltaTCut)
	{
		locShowersByBeamBunch[locNumShifts].push_back(locNeutralShower);
		--locNumShifts;
		locDeltaT = locVertexTime - (locRFTime + locNumShifts*dBeamBunchPeriod + dMaxDecayTimeOffset);
	}
}

double DSourceComboer::Calc_MaxDeltaTError(const DNeutralShower* locNeutralShower, const shared_ptr<const DKinematicData>& locKinematicData) const
{
	double locTheta = locKinematicData->momentum().Theta();
	if(locNeutralShower->dDetectorSystem == SYS_BCAL)
	{
		float& locZError = dPhotonVertexZBinWidth/2.0; //evaluated at center of bin
		double locR = locNeutralShower->dSpacetimeVertex.Vect().Perp();
		double locPathError = locR*(1.0/sin(locTheta) - sqrt(1.0 + pow(1.0/tan(locTheta) - locZError/locR, 2.0))) - locZError;
		return locPathError/SPEED_OF_LIGHT;
	}

	//FCAL
	double locDeltaZ = locNeutralShower->dSpacetimeVertex.Z() - dTargetCenter.Z();
	double locMaxZError = dTargetLength/2.0 + 15.0; //center of target + detached vertex
	//delta_delta_t = (650 - z_error)*[1/cos(theta) - sqrt(tan^2(theta) + (1 - z_error/(650 - z_error))^2)]/c - z_error/c
	double locPathErrorTerm = 1.0/cos(locTheta) - sqrt(pow(tan(locTheta), 2.0) + pow(1.0 - locMaxZError/(locDeltaZ - locMaxZError), 2.0));
	double locPathError = (locDeltaZ - locMaxZError)*locPathErrorTerm - locMaxZError;
	return locPathError/SPEED_OF_LIGHT;
}

/******************************************************************* CREATE DSOURCOMBOINFO'S ********************************************************************/
/******************************************************************* CREATE DSOURCOMBOINFO'S ********************************************************************/

void DSourceComboer::Build_ParticleCombos(JEventLoop* locEventLoop, const DReactionVertexInfo* locReactionVertexInfo)
{
	Reset_NewEvent(locEventLoop); //does nothing if not actually a new event

	//This builds the combos and creates DParticleCombo & DParticleComboSteps (doing whatever is necessary)


	//What to do about unknown vertices???:
		//Will know later (e.g. post beam selection): Use all possible showers, postpone BCAL shower PID & mass cuts (FCAL still fine)
		//Will never know: Use previous position for






	//FINAL DISCUSSION:

	//charged stage: charged only, no neutrals in infos

	//when on mixed stage (existing charged + neutrals, comboing into fully-neutral & mixed):
	//loop over charged combos: calc vertices, then build/convert FULL combo use with given vertex z-bins
	//then, just build the whole combo all it once, almost as before. however, some things are different
		//get charged particles to combo: choice is REDUCED to those from that vertex in the input combo
		//get charged combos to combo: if sub-combo is fully-charged, choice is REDUCED to be the input charged combo contents (almost always size 1)
			//thus we don't use ANY of the saved charged combos any more
			//and when we retrieve mixed combos for further comboing, they are specific (temporary) to this charged combo
				//Mixed results are saved in: unordered_map<mixed_use, unordered_map<charged_combo, vector<mixed_combo>>> (where the keys are the charged contents of the mixed-use step)
				//So that way we can re-use between channels
				//But how to RETRIEVE from here?, we need to get the charged combo from the given use //tricky, but we can do it
		//we do these because we don't want to rebuild the charged combos from scratch: wastes time, choices are restricted by vertex-z, we don't want to recompute vertex-z, we don't want dupe combos

	//combo the mixed stage in two stages:
	//FCAL showers only: z-bin any
	//All showers
		//here, they are comboed with uses having a specific vertex-z set
		//fully-neutral combos saved-to/retrieved-from charged-independent area for re-use (use already contains specific vertex-z bin)
		//first grab fcal combo results from fcal-only use area (or mixed area), first setting the z-bin to -1

	//Massive neutrals: Just combo with the rest of the neutrals
		//Just hold off on both mass cuts and timing cuts

	//get step vertex infos (sorted in dependency order)
	auto locStepVertexInfos = locReactionVertexInfo->Get_StepVertexInfos();
	auto locPrimaryStepVertexInfo = locReactionVertexInfo->Get_StepVertexInfo(0);
	auto locPrimaryComboUse = dSourceComboUseReactionMap[locPrimaryStepVertexInfo];


	//handle special case of no charged tracks
	if(dComboInfoChargeContent[std::get<2>(locPrimaryComboUse)] == d_Neutral)
	{
		//do stuff
		return;
	}

	//Build vertex combos (returns those for the primary vertex, others are stored)
	Create_SourceCombos(locPrimaryComboUse, d_ChargedStage, nullptr);
	const auto& locChargedCombos = *(Get_CombosSoFar(d_ChargedStage, d_Charged, nullptr)[locPrimaryComboUse]);

	//loop over primary vertex combos //each contains decay combos except when dangling
	for(auto locChargedCombo : locChargedCombos)
	{
		//Calc all the vertex positions and time offsets for the vertices for these combos (where possible without beam energy)
		dSourceComboVertexer.Calc_VertexTimeOffsets(locReactionVertexInfo, locChargedCombo);

		//For the charged tracks, apply timing cuts to determine which RF bunches are possible
//handle unknown vertex case!!
		vector<int> locPossibleBeamBunches_Charged;

		//deal with special case of FULLY charged
		auto locChargeContent = dComboInfoChargeContent[std::get<2>(locPrimaryComboUse)];
		if(locChargeContent == d_Charged)
		{
		}

		//Create full source-particle combos (including neutrals): First using only FCAL showers, then using all showers
		Create_SourceCombos(locPrimaryComboUse, d_MixedStage_ZIndependent, locChargedCombo);
		auto locZDependentComboUse = Create_ZDependentSourceComboUses(locReactionVertexInfo, locChargedCombo);
		Create_SourceCombos(locZDependentComboUse, d_MixedStage, locChargedCombo);

		//Then, get the full combos, but only those that satisfy the charged RF bunches
		const auto& locPrimaryFullCombos = Get_CombosForComboing(locZDependentComboUse, d_MixedStage, locPossibleBeamBunches_Charged, locChargedCombo);

	}


}

/**************************************************************** BUILD SOURCE COMBOS - GENERAL *****************************************************************/


/****************************************************** COMBOING STRATEGY ******************************************************
*
* Combos are not created in advance.  They are created on-demand, when needed.
*
* The BCAL photons are evaluated in different vertex-z bins for calculating their kinematics (momentum & timing).
* This is because their kinematics have a strong dependence on vertex-z, while the FCAL showers do not (see above derivations).
* Whereas the FCAL photons have only a small dependence, so their kinematics are regardless of vertex-z.
*
* The key to this being efficient (besides splitting the BCAL photons into vertex-z bins and placing timing cuts) is combo re-use.
* For example, suppose a channel needs 3 pi0s.
* First this will build all combos for 1 pi0, then all combos for 2 pi0s, then 3.  Placing mass cuts along the way.
* The results after each of these steps is saved.  That way, if someone then requests 2 pi0s, we merely have to return the results from the previous work.
* Also, if someone later requests 4pi0s, then we just take the 3pi0 results and expand them by 1 pi0.
*
* Ultimately, this results in a clusterfuck of recursive calls.
* Also, because of how the combo-info classes are structured (decaying PID NOT a member), you have be extremely careful not to get into an infinite loop.
* So, modify this code at your own peril. Just try not to take the rest of the collaboration down with you.
*
* Now, technically, when we construct combos for a (e.g.) pi0, we are saving 2 different results:
*    The combos of 2 photons, and which of those combos survive the pi0 mass cut.
* That way, if later someone wants to build an eta, all we have to do is take 2-photon combos and place eta mass cuts.
*
* Note that combos are constructed separately for different beam bunches.
* This is because photons only survive their timing cuts for certain beam bunches.
* Comboing only within a given beam bunch reduces the #photons we need to combo, and is thus faster.
*
* When comboing, first all of the FCAL showers alone are used to build the requested combos.
* Then, the BCAL showers surviving the timing cuts within the input vertex-z bin are used to build the requested combos.
* Finally, combos are created using a mix of these BCAL & FCAL showers.
* The results from this comboing is saved for all cases, that way they can be easily retrieved and combined as needed for similar requests.
*
* Note, if we combo vertically (e.g. 3pi0s, 2pi+'s, etc.), they are created with a use that may be a subset of the old.
* Then, when we combo them horizontally, they are promoted out of the vertical combo, at the same level as everything else in the new horizontal combo.
* This reduces the depth-complexity of the combos.
*
*******************************************************************************************************************************/

/*
 * suppose reaction is 0) g, p -> omega, p
 *                     1)         omega -> 3pi
 *                     2)                   pi0 -> 2g
 *
 * The purpose of passing through the charged combo:
 * 1) To retrieve the correct charged combo when comboing it to neutrals to create mixed
 * 2) To save the mixed comboing results in a way that they can be reused
 *
 * Note that, for a given step, the particles are grouped together as:
 * Decay PID -> All_Charged, All_Neutral, any mixed decays
 * Where All_Charged, All_Neutral are separate uses containing the entirety of the contents for that charge (X -> charged, X -> neutral)
 * Whereas any mixed decays are on the same level as the X -> charged & X -> neutral "decays"
 *
 * It will have uses like:
 * 0) X -> omega, A (mixed, further decays)		//presiding = 0, withnow = A
 *    A) X -> p (charged)						//both = nullptr
 * 1) omega -> B, C (mixed, further decays)		//presiding = 1, withnow = B
 *    B) X -> pi+, pi- (charged)				//both = nullptr
 *    C) X -> pi0 (neutral)						//both = nullptr
 * 2) pi0 -> 2g									//both = nullptr
 *
 */

void DSourceComboer::Create_SourceCombos(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding)
{
	//if on mixed stage, it is impossible for this function to be called with a fully-charged use (already exists!!)
	const auto& locDecayPID = std::get<0>(locComboUseToCreate);
	const auto& locVertexZBin = std::get<1>(locComboUseToCreate);
	const auto& locSourceComboInfo = std::get<2>(locComboUseToCreate);

	//we will create these combos for an "Unknown" decay (i.e. no decay, just a direct grouping)
	//then, when we return from this function, we can cut on the invariant mass of the system for any decay we might need it for
	DSourceComboUse locUnknownComboUse(Unknown, locVertexZBin, locSourceComboInfo);
	Create_SourceCombos_Unknown(locUnknownComboUse, locComboingStage, locChargedCombo_Presiding);

	//if all we want is a direct grouping (unknown), then the combos have already been made: return
	if(locDecayPID == Unknown)
		return;

	//Get combos so far
	auto locChargedCombo_WithNow = Get_ChargedCombo_WithNow(locChargedCombo_Presiding);
	auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, dComboInfoChargeContent[locSourceComboInfo], locChargedCombo_WithNow);
	auto locInfoChargeContent = dComboInfoChargeContent[locSourceComboInfo];

	//get the combos that we just created
	auto locSourceCombos = locSourceCombosByUseSoFar[locUnknownComboUse];

	if((locComboingStage == d_ChargedStage) && (locInfoChargeContent != d_Charged))
	{
		//don't cut yet! we don't have the neutrals! just copy results and return
		locSourceCombosByUseSoFar.emplace(locComboUseToCreate, locSourceCombos);
		return;
	}

	//get combos by beam bunch
	auto* locSourceCombosByBeamBunchByUse = (locComboingStage != d_ChargedStage) ? &(Get_SourceCombosByBeamBunchByUse(locInfoChargeContent, locChargedCombo_WithNow)) : nullptr;
	auto* locCombosByBeamBunch = (locComboingStage != d_ChargedStage) ? &((*locSourceCombosByBeamBunchByUse)[locComboUseToCreate]) : nullptr;

	//initialize vector for storing results
	locSourceCombosByUseSoFar.emplace(locComboUseToCreate, dResourcePool_SourceComboVector.Get_Resource());
	locSourceCombosByUseSoFar[locComboUseToCreate]->reserve(dInitialComboVectorCapacity);

	//place an invariant mass cut & save the results, but handle differently for massive neutrals
	if(dComboInfosWithMassiveNeutrals.find(locSourceComboInfo) == dComboInfosWithMassiveNeutrals.end())
	{
		for(auto locSourceCombo : *locSourceCombos)
		{
			if(!dSourceComboP4Handler.Cut_InvariantMass(locSourceCombo, locDecayPID, locVertexZBin))
				continue;

			//save the results
			locSourceCombosByUseSoFar[locComboUseToCreate]->push_back(locSourceCombo);
			if(locComboingStage == d_ChargedStage)
				continue;

			//register beam bunches
			const auto& locBeamBunches = Get_ValidRFBunches(locUnknownComboUse, locSourceCombo);
			for(const auto& locBeamBunch : locBeamBunches)
				(*locCombosByBeamBunch)[{locBeamBunch}].push_back(locSourceCombo);
			if(locBeamBunches.empty())
				(*locCombosByBeamBunch)[locBeamBunches].push_back(locSourceCombo);
		}
	}
	else //has massive neutrals
	{
		for(auto locSourceCombo : *locSourceCombos)
		{
			auto locBeamBunches = Get_ValidRFBunches(locUnknownComboUse, locSourceCombo);
			//if locBeamBunches is empty, don't cut: contains only massive neutrals, no determination of rf bunch yet
			if(!locBeamBunches.empty())
			{
				locBeamBunches = dSourceComboP4Handler.Cut_InvariantMass(locSourceCombo, locDecayPID, locBeamBunches);
				if(locBeamBunches.empty())
					continue; //failed the cut
			}

			//save the results
			locSourceCombosByUseSoFar[locComboUseToCreate]->push_back(locSourceCombo);

			//register beam bunches
			dValidRFBunches_ByUse.emplace(std::make_pair(locComboUseToCreate, locSourceCombo), locBeamBunches);
			for(const auto& locBeamBunch : locBeamBunches)
				(*locCombosByBeamBunch)[{locBeamBunch}].push_back(locSourceCombo);
			if(locBeamBunches.empty())
				(*locCombosByBeamBunch)[locBeamBunches].push_back(locSourceCombo);
		}
	}
}

void DSourceComboer::Create_SourceCombos_Unknown(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding)
{

	//First combo VERTICALLY, and then HORIZONTALLY
	//What does this mean?
	//Vertically: Make combos of size N of each PID needed (e.g. 3 pi0s)
	//Horizontally: Make combos of different PIDs (e.g. 2pi0, pi+, pi-, p)

	//Why start with vertical comboing?
	//because the thing that takes the most time is when someone decides to analyze (e.g.) 2pi0, 3pi0, then 3pi0 eta, 3pi0 something else, 4pi0, etc.
	//we want to make the Npi0 combos as needed, then reuse the Npi0s when making combos of other types
	//thus we want to build vertically (pi0s together, then etas together), and THEN horizontally (combine pi0s & etas, etc)
	//plus, when building vertically, it's easier to keep track of things since the PID / decay-parent is the same

	//Build all possible combos for all NEEDED GROUPINGS for each of the FURTHER DECAYS (if not done already)
	//this becomes a series of recursive calls
	//e.g. if need 3 pi0s, call for 2pi0s, which calls for 1pi0, which calls for 2g
		//then do the actual pi0 groupings on the return

	Combo_Vertically_AllDecays(locComboUseToCreate, locComboingStage, locChargedCombo_Presiding);
	if((locComboingStage == d_ChargedStage) || (dComboInfoChargeContent[std::get<2>(locComboUseToCreate)] == d_Neutral))
		Combo_Vertically_AllParticles(locComboUseToCreate, locComboingStage); //no such thing as a "mixed" particle

	//OK, now build horizontally!! //group particles with different PIDs
	Combo_Horizontally_All(locComboUseToCreate, locComboingStage, locChargedCombo_Presiding);
}
/*
bool DSourceComboer::Do_CommonComboingTasks(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locAllBut1ComboUse, const DSourceComboUse& locSourceComboUseToAdd, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_WithNow, bool locComboingVertically)
{
	//Returns true if comboing is done, false if there's more work to do
	auto& locSourceCombosByUseToSaveTo = Get_CombosSoFar(locComboingStage, dComboInfoChargeContent[locComboUseToCreate], locChargedCombo_WithNow);

	if(!locComboingVertically)
	{
		//this check will never be needed if comboing vertically (can't mix combo-charges)
		auto locChargeContent = dComboInfoChargeContent[std::get<2>(locSourceComboUseToAdd)];
		if((locComboingStage == d_ChargedStage) && (locChargeContent == d_Neutral))
		{
			//can't add neutrals, so we are already done! just copy the results to the new vector
			auto& locSourceCombosByUseSoFar = Get_CombosSoFar(d_ChargedStage, d_Charged, nullptr);
			auto locCombos_AllBut1 = locSourceCombosByUseSoFar[locAllBut1ComboUse];
			locSourceCombosByUseToSaveTo.emplace(locComboUseToCreate, locCombos_AllBut1);
			return true;
		}
	}

	//if on the all-showers stage, first copy over ALL fcal-only results
	if(locComboingStage == d_MixedStage)
		Copy_FCALOnlyResults(locComboUseToCreate, locComboingStage, locChargedCombo_WithNow);
	else //initialize vector for storing results
	{
		locSourceCombosByUseToSaveTo.emplace(locComboUseToCreate, dResourcePool_SourceComboVector.Get_Resource());
		locSourceCombosByUseToSaveTo[locComboUseToCreate]->reserve(dInitialComboVectorCapacity);
	}

	return false;
}
*/
/************************************************************** BUILD PHOTON COMBOS - VERTICALLY ****************************************************************/

void DSourceComboer::Combo_Vertically_AllDecays(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding)
{
	auto locChargedCombo_WithNow = Get_ChargedCombo_WithNow(locChargedCombo_Presiding);

	//Get combos so far
	auto locComboInfo = std::get<2>(locComboUseToCreate);
	auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, dComboInfoChargeContent[locComboInfo], locChargedCombo_WithNow);

	//get combo use contents
	auto locVertexZBin = std::get<1>(locComboUseToCreate);
	auto locNumParticlesNeeded = locComboInfo->Get_NumParticles();
	auto locFurtherDecays = locComboInfo->Get_FurtherDecays();

	//for each further decay map entry (e.g. pi0, 3), this is a collection of the uses representing those groupings //e.g. Unknown -> 3pi0
	for(const auto& locFurtherDecayPair : locFurtherDecays)
	{
		auto& locSourceComboDecayUse = locFurtherDecayPair.first; //e.g. pi0, -> 2g
		auto& locNumDecaysNeeded = locFurtherDecayPair.second; //N of the above decay (e.g. pi0s)
		auto locDecayChargeContent = dComboInfoChargeContent[std::get<2>(locSourceComboDecayUse)];

		if((locComboingStage == d_ChargedStage) && (locDecayChargeContent == d_Neutral))
			continue; //skip for now!!

		if(locNumDecaysNeeded == 1)
		{
			//if on a mixed stage, and the to-build combo info is fully charged, skip it: it's already been done
			if((locComboingStage != d_ChargedStage) && (locDecayChargeContent == d_Charged))
				continue;

			//final dependency: no more "further decays" after this
			//e.g. the input locSourceComboInfo here has no photons and 1 further decay use (pi0)
			//so, just build the pi0 combos directly
			if(locSourceCombosByUseSoFar.find(locSourceComboDecayUse) != locSourceCombosByUseSoFar.end()) //if not done already!
			{
				//must dive down to get the next charged combo
				//building for the first time: the first one (later ones will be grabbed when building these combos vertically (in Combo_Vertically_NDecays))
				auto locChargedCombo_NextPresiding = Get_Presiding_ChargedCombo(locChargedCombo_Presiding, locSourceComboDecayUse, locComboingStage, 1);

				//must return to top-level combo function, as this may have any structure
				Create_SourceCombos(locSourceComboDecayUse, locComboingStage, locChargedCombo_NextPresiding);
			}
			continue;
		}

		//OK, so we need a grouping of N > 1 decays (e.g. pi0s)
		//so, let's create a use of Unknown -> N pi0s (e.g.)
		//if we can just utilize the use from the input combo-info, then we will. if not, we'll make a new one
		auto locNeededGroupingUse = locComboUseToCreate;
		if((locFurtherDecays.size() > 1) || !locNumParticlesNeeded.empty()) //if true: can't use the input
		{
			auto locGroupingComboInfo = GetOrMake_SourceComboInfo({}, {locSourceComboDecayUse, locNumDecaysNeeded}); // -> N pi0s (e.g.)
			locNeededGroupingUse = (Unknown, locVertexZBin, locGroupingComboInfo); // Unknown -> Npi0s (e.g.)
		}

		// Now, see whether the combos for this grouping have already been done
		if(locSourceCombosByUseSoFar.find(locNeededGroupingUse) != locSourceCombosByUseSoFar.end())
			continue; //it's already done!!

		//it's not already done.  darn it.
		//build an info and a use for a direct grouping of N - 1 decays //e.g. 2pi0s
		auto locNMinus1ComboUse = locSourceComboDecayUse; //initialize (is valid if #needed == 2, otherwise will create it)
		if(locNumDecaysNeeded > 2)
		{
			auto locNMinus1Info = GetOrMake_SourceComboInfo({}, {locSourceComboDecayUse, locNumDecaysNeeded - 1}); // 0 detected particles, N - 1 pi0s (e.g.)
			locNMinus1ComboUse(Unknown, locVertexZBin, locNMinus1Info); // Unknown -> N - 1 pi0s (e.g.)
		}

		// Now, see whether the combos for the direct N - 1 grouping have already been done.  If not, create them
		if(locSourceCombosByUseSoFar.find(locNMinus1ComboUse) != locSourceCombosByUseSoFar.end())
			Combo_Vertically_AllDecays(locNMinus1ComboUse, locComboingStage, locChargedCombo_WithNow); //no need to go to top-level combo function since just N - 1: can re-call this one

		//Finally, we can actually DO the grouping, between the N - 1 combos and the one-off combos
		Combo_Vertically_NDecays(locNeededGroupingUse, locNMinus1ComboUse, locSourceComboDecayUse, locComboingStage, locChargedCombo_WithNow);
	}
}

void DSourceComboer::Combo_Vertically_NDecays(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locNMinus1ComboUse, const DSourceComboUse& locSourceComboDecayUse, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding)
{
	auto locVertexZBin = std::get<1>(locComboUseToCreate);

	//Get combos so far
	auto locChargeContent = dComboInfoChargeContent[std::get<2>(locComboUseToCreate)];
	auto locChargedCombo_WithNow = Get_ChargedCombo_WithNow(locChargedCombo_Presiding);
	auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, locChargeContent, locChargedCombo_WithNow);

	//e.g. we are grouping 1 pi0 with N - 1 pi0s to make a combo of N pi0s
	//so, let's get the combos for (e.g.) 1 pi0 and for N - 1 pi0s
	const auto& locCombos_NMinus1 = *locSourceCombosByUseSoFar[locNMinus1ComboUse]; //Combos are a vector of (e.g.): -> N - 1 pi0s
	if(locCombos_NMinus1.empty())
		return; //bail!

	//if on the all-showers stage, first copy over ALL fcal-only results
	if(locComboingStage == d_MixedStage)
		Copy_ZIndependentMixedResults(locComboUseToCreate, locComboingStage, locChargedCombo_WithNow);
	else //initialize vector for storing results
	{
		locSourceCombosByUseSoFar.emplace(locComboUseToCreate, dResourcePool_SourceComboVector.Get_Resource());
		locSourceCombosByUseSoFar[locComboUseToCreate]->reserve(dInitialComboVectorCapacity);
	}

	//if comboing N mixed combos (locComboUseToCreate) (which are thus all used in the same step), do this:
	//locChargedCombo_WithNow corresponds to N mixed combos
	const DSourceCombo* locChargedCombo_WithPrevious = nullptr; //corresponds to ONE mixed combo: the next one //will set it below
	auto locInstance = locCombos_NMinus1.front()->Get_FurtherDecayCombos()[locSourceComboDecayUse].size() + 1;
	const DSourceCombo* locChargedCombo_WithPrevious = Get_ChargedCombo_WithNow(Get_Presiding_ChargedCombo(locChargedCombo_Presiding, locSourceComboDecayUse, locComboingStage, locInstance));

	//now, for each combo of N - 1 (e.g.) pi0s, see which of the single-decay combos are a valid grouping
	//valid grouping:
		//TEST 1: If (e.g.) pi0s have names "A", "B", "C", don't include the grouping "ABA", and don't include "ACB" if we already have "ABC"
		//TEST 2: Also, don't re-use a shower we've already used (e.g. if A & C each contain the same photon, don't group them together)
		//Technically, if we pass Test 2 we automatically pass Test 1.
		//However, validating for Test 1 is much faster, as discussed below.
	for(auto locCombo_NMinus1 : locCombos_NMinus1)
	{
		//loop over potential combos to add to the group, creating a new combo for each valid (non-duplicate) grouping
		//however, we don't have to loop over all of the combos!!

		//first of all, get the potential combos that satisfy the RF bunches for the N - 1 combo
		const auto& locValidRFBunches_NMinus1 = Get_ValidRFBunches(locNMinus1ComboUse, locCombo_NMinus1);
		const auto& locDecayCombos_1 = Get_CombosForComboing(locSourceComboDecayUse, locComboingStage, locValidRFBunches_NMinus1, locChargedCombo_WithPrevious);

		//now, note that all of the combos are stored in the order in which they were created (e.g. A, B, C, D)
		//so (e.g.), groupings of 2 will be created and saved in the order: AB, AC, AD, BC, BD, CD
		//above, on the B-loop, we start the search at "C," not at A, because this was already tested on an earlier pass
		//therefore, start the search one AFTER the LAST (e.g. -> 2 photon) combo of the N - 1 group
		//this will guarantee we pass "TEST 1" without ever checking

		//actually, we already saved the iterator to the first (e.g.) pi0 to test when we saved the N - 1 combo, so just retrieve it
		auto locComboSearchIterator = Get_ResumeAtIterator_Combos(locCombo_NMinus1, locValidRFBunches_NMinus1, locComboingStage, locVertexZBin);
		if(locComboSearchIterator == std::end(locDecayCombos_1))
			continue; //e.g. this combo is "AD" and there are only 4 reconstructed combos (ABCD): no potential matches! move on to the next N - 1 combo

		//before we loop, first get all of the showers used to make the N - 1 grouping, and sort it so that we can quickly search it
		auto locUsedParticles_NMinus1 = DAnalysis::Get_SourceParticles(locCombo_NMinus1->Get_SourceParticles(true)); //true: entire chain
		std::sort(locUsedParticles_NMinus1.begin(), locUsedParticles_NMinus1.end()); //must sort, because when retrieving entire chain is unsorted

		//this function will do our "TEST 2"
		auto Search_Duplicates = [&locUsedParticles_NMinus1](const JObject* locParticle) -> bool
				{return std::binary_search(locUsedParticles_NMinus1.begin(), locUsedParticles_NMinus1.end(), locParticle);};

		auto locIsZIndependent_NMinus1 = locCombo_NMinus1->Get_IsZIndependent();

		//now loop over the potential combos
		for(; locComboSearchIterator != locDecayCombos_1.end(); ++locComboSearchIterator)
		{
			const auto locDecayCombo_1 = *locComboSearchIterator;

			//If on all-showers stage, and combo is fcal-only, don't save (combo already created!!)
			auto locIsZIndependent = locIsZIndependent_NMinus1 && locDecayCombo_1->Get_IsZIndependent();
			if((locComboingStage == d_MixedStage) && locIsZIndependent)
				continue; //this combo has already been created (assuming it was valid): during the FCAL-only stage

			//conduct "TEST 2" search: search the N - 1 shower vector to see if any of the showers in this combo are duplicated
			auto locUsedParticles_1 = DAnalysis::Get_SourceParticles(locDecayCombo_1->Get_SourceParticles(true)); //true: entire chain
			if(std::any_of(locUsedParticles_1.begin(), locUsedParticles_1.end(), Search_Duplicates))
				continue; //at least one photon was a duplicate, this combo won't work

			//no duplicates: this combo is unique.  build a new combo!

			//See which RF bunches match up //guaranteed to be at least one, due to selection in Get_ParticlesForComboing() function
			const auto& locValidRFBunches_1 = Get_ValidRFBunches(locSourceComboDecayUse, locDecayCombo_1);
			auto locValidRFBunches = Get_CommonRFBunches(locValidRFBunches_NMinus1, locValidRFBunches_1);

			//take the vector of N - 1 (e.g. -> 2g) combos and add the new one
			auto locAllDecayCombos = locCombo_NMinus1->Get_FurtherDecayCombos()[locSourceComboDecayUse];
			locAllDecayCombos.push_back(locDecayCombo_1);

			//then create the new combo
			auto locFurtherDecayCombos(locSourceComboDecayUse, locAllDecayCombos); //arguments (e.g.): (pi0, -> 2g), N combos of: -> 2g
			auto locCombo = dResourcePool_SourceCombo.Get_Resource();
			locCombo->Set_Members({}, locFurtherDecayCombos, locIsZIndependent); // 1 combo of N (e.g.) pi0s

			//save it! //in creation order!
			locSourceCombosByUseSoFar[locComboUseToCreate]->push_back(locCombo);
			Register_ValidRFBunches(locComboUseToCreate, locCombo, locValidRFBunches, locComboingStage, locChargedCombo_WithNow);

			//finally, in case we add more (e.g.) pi0s later (N + 1), save the last pi0
			//so that we will start the search for the next (e.g.) pi0 in the location after the last one
			dResumeSearchAfterMap_Combos[locCombo] = locDecayCombo_1;
		}
	}
}

void DSourceComboer::Combo_Vertically_AllParticles(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage)
{
	//get combo use contents
	auto locVertexZBin = std::get<1>(locComboUseToCreate);
	auto locNumParticlesNeeded = std::get<2>(locComboUseToCreate)->Get_NumParticles();
	auto locFurtherDecays = std::get<2>(locComboUseToCreate)->Get_FurtherDecays();

	//Get combos so far //guaranteed not to be mixed
	auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, d_Neutral); //if not neutral then is on charged stage: argument doesn't matter

	//for each further decay map entry (e.g. pi0, 3), this is a collection of the uses representing those groupings //e.g. Unknown -> 3pi0
	for(const auto& locParticlePair : locNumParticlesNeeded)
	{
		//get PID information
		auto& locPID = locParticlePair.first; //e.g. pi0, -> 2g
		auto& locNumPIDNeeded = locParticlePair.second; //N of the above decay (e.g. pi0s)

		if(locNumPIDNeeded == 1)
			continue; //nothing to do vertically; we will combo this horizontally later

		if((locComboingStage == d_ChargedStage) && (ParticleCharge(locPID) == 0))
			continue; //skip for now!!

		//OK, so we need a grouping of N > 1 particles with the same PID (e.g. g's)
		//so, let's create a use of Unknown -> N g's (e.g.)
		//if we can just utilize the use from the input combo-info, then we will. if not, we'll make a new one
		DSourceComboUse locNeededGroupingUse = locComboUseToCreate;
		if((locNumParticlesNeeded.size() > 1) || !locFurtherDecays.empty()) //if true: can't use the input
		{
			auto locGroupingComboInfo = GetOrMake_SourceComboInfo({locPID, locNumPIDNeeded}, {}); // -> N g's (e.g.)
			locNeededGroupingUse = (Unknown, locVertexZBin, locGroupingComboInfo); // Unknown -> N g's (e.g.)
		}

		//See whether the combos for this grouping have already been done
		if(locSourceCombosByUseSoFar.find(locNeededGroupingUse) != locSourceCombosByUseSoFar.end())
			continue; //it's already done!!

		//it's not already done.  darn it.
		//if it's a direct combo of 2 particles, just make it and continue
		if(locNumPIDNeeded == 2)
		{
			Combo_Vertically_NParticles(locNeededGroupingUse, DSourceComboUse(), locComboingStage);
			continue;
		}

		//build an info and a use for a direct grouping of N - 1 particles //e.g. 3 g's
		auto locNMinus1Info = GetOrMake_SourceComboInfo({locPID, locNumPIDNeeded - 1}, {}); // N - 1 g's (e.g.), no decaying particles
		DSourceComboUse locNMinus1ComboUse(Unknown, locVertexZBin, locNMinus1Info); // Unknown -> N - 1 g's (e.g.)

		// Now, see whether the combos for the direct N - 1 grouping have already been done.  If not, create them
		if(locSourceCombosByUseSoFar.find(locNMinus1ComboUse) != locSourceCombosByUseSoFar.end())
			Combo_Vertically_AllParticles(locNMinus1ComboUse, locComboingStage); //no need to go to top-level combo function since just N - 1: can re-call this one

		//Finally, we can actually DO the grouping, between the N - 1 particles and one more particle
		Combo_Vertically_NParticles(locNeededGroupingUse, locNMinus1ComboUse, locComboingStage);
	}
}

void DSourceComboer::Combo_Vertically_NParticles(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locNMinus1ComboUse, ComboingStage_t locComboingStage)
{
	auto locVertexZBin = std::get<1>(locComboUseToCreate);

	//either: combining two particles with the same PID to create a new combo, or combining a combo of N particles (with the same PID) with one more particle
	auto locComboInfo = std::get<2>(locComboUseToCreate);
	auto locParticlePair = locComboInfo->Get_NumParticles().back(); //is guaranteed to be size 1
	auto locPID = locParticlePair.first;
	auto locNumParticles = locParticlePair.second;

	//Get combos so far //guaranteed not to be mixed
	auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, d_Neutral); //if not neutral then is on charged stage: argument doesn't matter

	//if on the all-showers stage, first copy over ALL fcal-only results
	if(locComboingStage == d_MixedStage)
		Copy_ZIndependentMixedResults(locComboUseToCreate, locComboingStage, nullptr);
	else //initialize vector for storing results
	{
		locSourceCombosByUseSoFar.emplace(locComboUseToCreate, dResourcePool_SourceComboVector.Get_Resource());
		locSourceCombosByUseSoFar[locComboUseToCreate]->reserve(dInitialComboVectorCapacity);
	}

	//For checking RF bunches (ignored if on charged stage)
	const unordered_map<const JObject*, vector<int>>& locShowerRFMap = (locComboingStage == d_MixedStage_ZIndependent) ? dShowerRFBunches_FCAL : dShowerRFBunches_Both[locVertexZBin];

	if(locNumParticles == 2)
	{
		//Get particles for comboing
		const auto& locParticles = Get_ParticlesForComboing(locPID, locComboingStage, {}, locVertexZBin);

		auto locLastIteratorToCheck = std::prev(locParticles.end());
		for(auto locFirstIterator = locParticles.begin(); locFirstIterator != locLastIteratorToCheck; ++locFirstIterator)
		{
			const auto locRFBunches_First = (locPID == Gamma) ? &(locShowerRFMap[*locFirstIterator]) : nullptr;
			for(auto locSecondIterator = std::next(locFirstIterator); locSecondIterator != locParticles.end(); ++locSecondIterator)
			{
				auto locIsZIndependent = (locComboingStage == d_MixedStage_ZIndependent) || (Get_IsZIndependent(*locFirstIterator) && Get_IsZIndependent(*locSecondIterator));
				if((locComboingStage == d_MixedStage) && locIsZIndependent)
					continue; //this combo has already been created (assuming it was valid): during the FCAL-only stage

				//See which RF bunches match up, if any
				vector<int> locValidRFBunches = {};
				if(locPID == Gamma) //if charged or massive neutrals, ignore (they don't choose at this stage)
				{
					const auto& locRFBunches_Second = locShowerRFMap[*locSecondIterator];
					locValidRFBunches = Get_CommonRFBunches(*locRFBunches_First, locRFBunches_Second);
					if(locValidRFBunches.empty())
						continue;
				}

				auto locCombo = dResourcePool_SourceCombo.Get_Resource();
				locCombo->Set_Members({std::make_pair(locPID, *locFirstIterator), std::make_pair(locPID, *locSecondIterator)}, {}, locIsZIndependent);
				locSourceCombosByUseSoFar[locComboUseToCreate]->push_back(locCombo); //save it //in creation order

				Register_ValidRFBunches(locComboUseToCreate, locCombo, locValidRFBunches, locComboingStage, nullptr);

				//in case we add more particles with the same PID later (N + 1), save last object with this PID
				//so that we will start the search for the next particle one spot after it
				dResumeSearchAfterMap_Particles[locCombo] = *locSecondIterator;
			}
		}
		return;
	}

	//create combo of N same-PID-particles by adding one particle to previously-created combos of N - 1 same-PID-particles
	const auto& locCombos_NMinus1 = *locSourceCombosByUseSoFar[locNMinus1ComboUse]; //Each combo contains a vector of N - 1 same-PID-particles
	for(auto locCombo_NMinus1 : locCombos_NMinus1)
	{
		//Get particles for comboing
		const auto& locValidRFBunches_NMinus1 = Get_ValidRFBunches(locNMinus1ComboUse, locCombo_NMinus1);
		const auto& locParticles = Get_ParticlesForComboing(locPID, locComboingStage, locValidRFBunches_NMinus1, locVertexZBin);

		//retrieve where to begin the search
		auto locParticleSearchIterator = Get_ResumeAtIterator_Particles(locCombo_NMinus1, locValidRFBunches_NMinus1);
		if(locParticleSearchIterator == std::end(locParticles))
			continue; //e.g. this combo is "AD" and there are only 4 reconstructed combos (ABCD): no potential matches! move on to the next N - 1 combo

		auto locIsZIndependent_NMinus1 = locCombo_NMinus1->Get_IsZIndependent();

		for(; locParticleSearchIterator != locParticles.end(); ++locParticleSearchIterator)
		{
			auto locIsZIndependent = (locComboingStage == d_MixedStage_ZIndependent) || (locIsZIndependent_NMinus1 && Get_IsZIndependent(*locParticleSearchIterator));
			if((locComboingStage == d_MixedStage) && locIsZIndependent)
				continue; //this combo has already been created (assuming it was valid): during the FCAL-only stage

			//See which RF bunches match up //guaranteed to be at least one, due to selection in Get_ParticlesForComboing() function
			vector<int> locValidRFBunches = {}; //if charged or massive neutrals, ignore (they don't choose at this stage)
			if(locPID == Gamma)
			{
				const vector<int>& locRFBunches_Second = locShowerRFMap[*locParticleSearchIterator];
				locValidRFBunches = Get_CommonRFBunches(locValidRFBunches_NMinus1, locRFBunches_Second);
			}

			auto locComboParticlePairs = locCombo_NMinus1->Get_SourceParticles();
			locComboParticlePairs.emplace_back(locPID, *locParticleSearchIterator);
			auto locCombo = dResourcePool_SourceCombo.Get_Resource();
			locCombo->Set_Members(locComboParticlePairs, {}, locIsZIndependent);
			locSourceCombosByUseSoFar[locComboUseToCreate]->push_back(locCombo); //save it //in creation order

			Register_ValidRFBunches(locComboUseToCreate, locCombo, locValidRFBunches, locComboingStage, nullptr);

			//in case we add more particles with the same PID later (N + 1), save last object with this PID
			//so that we will start the search for the next particle one spot after it
			dResumeSearchAfterMap_Particles[locCombo] = *locParticleSearchIterator;
		}
	}
}

/************************************************************* BUILD PHOTON COMBOS - HORIZONTALLY ***************************************************************/

void DSourceComboer::Combo_Horizontally_All(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding)
{
	//get combo use contents
	auto locVertexZBin = std::get<1>(locComboUseToCreate);
	const auto& locComboInfo = std::get<2>(locComboUseToCreate);
	auto locNumParticlesNeeded = locComboInfo->Get_NumParticles();
	auto locFurtherDecays = locComboInfo->Get_FurtherDecays();

	//first handle special cases:
	if(locNumParticlesNeeded.empty() && (locFurtherDecays.size() == 1))
		return; //e.g. we just need N pi0s together: already done when comboing vertically!!
	if(locFurtherDecays.empty() && (locNumParticlesNeeded.size() == 1))
	{
		//we just need N (e.g.) photons together
		auto& locParticlePair = locNumParticlesNeeded.front();
		if(locParticlePair.second > 1)
			return; //already done when comboing vertically!!

		//not much of a combo if there's only 1, is it? //e.g. 1 charged track at a vertex
		if((locComboingStage == d_ChargedStage) && (ParticleCharge(locParticlePair.first) == 0))
			return; //skip for now!!
		Create_Combo_OneParticle(locComboUseToCreate, locComboingStage);
		return;
	}

	//see if there is another combo that already exists that is a subset of what we requested
	//e.g. if we need a charged combo, a neutral combo, and a mixed: search for:
		//charged + neutral (no mixed)
		//charged + mixed (no neutral)
		//neutral + mixed (no charged)
	//e.g. if we need 2pi0s, one omega, and 1g: search for:
		//2pi0s, one omega: if exists, just combo that with 1g
		//2pi0s, one photon: if exists, just combo with one omega
		//etc.

	//save in case need to create these
	DSourceComboUse locComboUse_SubsetToBuild(Unknown, locVertexZBin, nullptr);

	//for each further decay map entry (e.g. pi0, 3), this is a collection of the uses representing those groupings //e.g. Unknown -> 3pi0
	//decays are sorted by: mixed-charge first, then fully-neutral, then fully-charged
	//within a charge: loop from heaviest-mass to least (most likely to be missing)
	auto locChargedCombo_WithNow = Get_ChargedCombo_WithNow(locChargedCombo_Presiding);
	for(auto locDecayIterator = locFurtherDecays.begin(); locDecayIterator != locFurtherDecays.end(); ++locDecayIterator)
	{
		//build a DSourceComboUse with everything EXCEPT this set of decays, and see if it already exists
		//build the further-decays, removing this decay
		auto locFurtherDecaysToSearchFor = locFurtherDecays;
		const auto& locSourceComboUse_ThisDecay = locDecayIterator->first;
		auto locChargeContent_ThisDecay = dComboInfoChargeContent[std::get<2>(locSourceComboUse_ThisDecay)];
		locFurtherDecaysToSearchFor.erase(locFurtherDecaysToSearchFor.begin() + std::distance(locFurtherDecays.begin(), locDecayIterator));

		//build the DSourceComboUse
		auto locAllBut1ComboInfo = GetOrMake_SourceComboInfo(locNumParticlesNeeded, locFurtherDecaysToSearchFor);
		auto locAllBut1ChargeContent = dComboInfoChargeContent[locAllBut1ComboInfo];
		if((locComboingStage == d_ChargedStage) && (locAllBut1ChargeContent == d_Neutral))
			continue; //this won't be done yet!
		DSourceComboUse locAllBut1ComboUse(Unknown, locVertexZBin, locAllBut1ComboInfo); // Unknown -> everything but this decay

		if((locComboingStage != d_ChargedStage) && (locAllBut1ChargeContent == d_Charged))
		{
			//yes, it's already been done!
			//just combo the All-but-1 combos to those from this decay and return the results
			Combo_Horizontally_AddCombo(locComboUseToCreate, locAllBut1ComboUse, locSourceComboUse_ThisDecay, locComboingStage, locChargedCombo_WithNow);
			return;
		}

		//Get combos so far
		auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, dComboInfoChargeContent[locAllBut1ComboInfo], locChargedCombo_WithNow);

		// Now, see whether the combos for this grouping have already been done
		if(locSourceCombosByUseSoFar.find(locAllBut1ComboUse) == locSourceCombosByUseSoFar.end()) //if true: not yet
		{
			//if on the first one (heaviest mass), save this subset in case we need to create it (if nothing else already done)
			if(locDecayIterator == locFurtherDecays.begin())
				locComboUse_SubsetToBuild = locAllBut1ComboUse;
			continue; // try the next decay
		}

		//yes, it's already been done!
		//just combo the All-but-1 combos to those from this decay and save the results
		if((locComboingStage == d_ChargedStage) && (locChargeContent_ThisDecay == d_Neutral))
		{
			//this won't be done yet! just copy the all-but-1 as the desired combos
			auto& locAllBut1Combos = locSourceCombosByUseSoFar[locAllBut1ComboUse];
			locSourceCombosByUseSoFar.emplace(locComboUseToCreate, locAllBut1Combos);
		}
		else
			Combo_Horizontally_AddCombo(locComboUseToCreate, locAllBut1ComboUse, locSourceComboUse_ThisDecay, locComboingStage, locChargedCombo_WithNow);
		return;
	}

	//ok, none of the subsets without a decay has yet been created. let's try subsets without detected particles
	if((locComboingStage == d_ChargedStage) || (dComboInfoChargeContent[locComboInfo] == d_Neutral)) //no loose particles when mixing charged & neutral
	{
		for(auto locParticleIterator = locNumParticlesNeeded.begin(); locParticleIterator != locNumParticlesNeeded.end(); ++locParticleIterator)
		{
			//build a DSourceComboUse with everything EXCEPT this set of particles, and see if it already exists
			//combo the particle horizontally, removing this PID
			auto locNumParticlesToSearchFor = locNumParticlesNeeded;
			const auto& locParticlePair = *locParticleIterator;
			locNumParticlesToSearchFor.erase(locNumParticlesToSearchFor.begin() + std::distance(locNumParticlesNeeded.begin(), locParticleIterator));

			//build the DSourceComboUse
			auto locAllBut1ComboInfo = GetOrMake_SourceComboInfo(locNumParticlesToSearchFor, locFurtherDecays);
			if((locComboingStage == d_ChargedStage) && (dComboInfoChargeContent[locAllBut1ComboInfo] == d_Neutral))
				continue; //this won't be done yet!
			DSourceComboUse locAllBut1ComboUse(Unknown, locVertexZBin, locAllBut1ComboInfo); // Unknown -> everything but these particles

			//Get combos so far
			auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, d_Neutral); //if not neutral then is on charged stage: argument doesn't matter

			// Now, see whether the combos for this grouping have already been done
			if(locSourceCombosByUseSoFar.find(locAllBut1ComboUse) == locSourceCombosByUseSoFar.end()) //if true: not yet
			{
				//if on the first one and there's no decays, save this subset in case we need to create it (if nothing else already done)
				if((locParticleIterator == locNumParticlesNeeded.begin()) && locFurtherDecays.empty())
					locComboUse_SubsetToBuild = locAllBut1ComboUse;
				continue; // try the next PID
			}

			//yes, it's already been done!
			//just combo the All-but-1 combos to those from this particle and return the results
			if((locComboingStage == d_ChargedStage) && (ParticleCharge(locParticlePair.first) == 0))
			{
				//this won't be done yet! just copy the all-but-1 as the desired combos
				auto& locAllBut1Combos = locSourceCombosByUseSoFar[locAllBut1ComboUse];
				locSourceCombosByUseSoFar.emplace(locComboUseToCreate, locAllBut1Combos);
				return;
			}

			if(locParticlePair.second > 1)
			{
				//create a combo use for X -> N particles of this type
				auto locSourceInfo_NParticles = GetOrMake_SourceComboInfo({locParticlePair}, {});
				DSourceComboUse locSourceComboUse_NParticles(Unknown, locSourceInfo_NParticles);
				Combo_Horizontally_AddCombo(locComboUseToCreate, locAllBut1ComboUse, locSourceComboUse_NParticles, locComboingStage, locChargedCombo_WithNow);
			}
			else
				Combo_Horizontally_AddParticle(locComboUseToCreate, locAllBut1ComboUse, locParticlePair.first, locComboingStage, locChargedCombo_WithNow);
			return;
		}
	}

	//none of the possible immediate subsets have been created
	//therefore, create one of them (the one without the heaviest particle), and then do the remaining combo
	Combo_Horizontally_All(locComboUse_SubsetToBuild, locComboingStage, locChargedCombo_WithNow);

	//do the final combo!
	if(locFurtherDecays.empty())
	{
		//subset was missing a detected PID
		const auto& locParticlePair = locNumParticlesNeeded.front();
		if((locComboingStage == d_ChargedStage) && (ParticleCharge(locParticlePair.first) == 0))
		{
			//this won't be done yet! just copy the all-but-1 as the desired combos
			auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, d_Charged);
			auto& locAllBut1Combos = locSourceCombosByUseSoFar[locComboUse_SubsetToBuild];
			locSourceCombosByUseSoFar.emplace(locComboUseToCreate, locAllBut1Combos);
			return;
		}
		if(locParticlePair.second > 1)
		{
			//create a combo use for X -> N particles of this type
			auto locSourceInfo_NParticles = GetOrMake_SourceComboInfo({locParticlePair}, {});
			DSourceComboUse locSourceComboUse_NParticles(Unknown, locVertexZBin, locSourceInfo_NParticles);
			Combo_Horizontally_AddCombo(locComboUseToCreate, locComboUse_SubsetToBuild, locSourceComboUse_NParticles, locComboingStage, locChargedCombo_WithNow);
		}
		else
			Combo_Horizontally_AddParticle(locComboUseToCreate, locComboUse_SubsetToBuild, locParticlePair.first, locComboingStage, locChargedCombo_WithNow);
	}
	else //subset was missing a decay PID
	{
		auto locComboUseToAdd = locFurtherDecays.front().first;
		if((locComboingStage == d_ChargedStage) && (std::get<2>(locComboUseToAdd) == d_Neutral))
		{
			//this won't be done yet! just copy the all-but-1 as the desired combos
			auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, d_Charged);
			auto& locAllBut1Combos = locSourceCombosByUseSoFar[locComboUse_SubsetToBuild];
			locSourceCombosByUseSoFar.emplace(locComboUseToCreate, locAllBut1Combos);
		}
		else
			Combo_Horizontally_AddCombo(locComboUseToCreate, locComboUse_SubsetToBuild, locComboUseToAdd, locComboingStage, locChargedCombo_WithNow);
	}
}

void DSourceComboer::Create_Combo_OneParticle(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage)
{
	//not much of a combo if there's only 1, is it? //e.g. 1 charged track at a vertex

	//Get combos so far
	auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, d_Neutral); //if not neutral then is on charged stage: argument doesn't matter

	//get combo use contents
	auto locParticlePair = std::get<2>(locComboUseToCreate)->Get_NumParticles().front();
	auto locVertexZBin = std::get<1>(locComboUseToCreate);

	//if on the mixed stage, must be doing all neutrals: first copy over ALL fcal-only results
	if(locComboingStage == d_MixedStage)
		Copy_ZIndependentMixedResults(locComboUseToCreate, locComboingStage, nullptr);
	else //initialize vector for storing results
	{
		locSourceCombosByUseSoFar.emplace(locComboUseToCreate, dResourcePool_SourceComboVector.Get_Resource());
		locSourceCombosByUseSoFar[locComboUseToCreate]->reserve(dInitialComboVectorCapacity);
	}

	//For checking RF bunches (ignored if on charged stage)
	const unordered_map<const JObject*, vector<int>>& locShowerRFMap = (locComboingStage == d_MixedStage_ZIndependent) ? dShowerRFBunches_FCAL : dShowerRFBunches_Both[locVertexZBin];

	auto locPID = locParticlePair.first;

	//Get particles for comboing
	const auto& locParticles = Get_ParticlesForComboing(locPID, locComboingStage, {}, locVertexZBin);
	for(auto locParticle : locParticles)
	{
		auto locIsZIndependent = Get_IsZIndependent(locParticle);
		if((locComboingStage == d_MixedStage) && locIsZIndependent)
			continue; //this combo has already been created (assuming it was valid): during the FCAL-only stage

		auto locCombo = dResourcePool_SourceCombo.Get_Resource();
		locCombo->Set_Members({std::make_pair(locPID, locParticle)}, {}, locIsZIndependent);
		locSourceCombosByUseSoFar[locComboUseToCreate]->push_back(locCombo); //save it //in creation order
		if(locPID == Gamma)
			Register_ValidRFBunches(locComboUseToCreate, locCombo, locShowerRFMap.find(locParticle)->second, locComboingStage, nullptr);
		else
			Register_ValidRFBunches(locComboUseToCreate, locCombo, {}, locComboingStage, nullptr);
	}
}

void DSourceComboer::Combo_Horizontally_AddCombo(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locAllBut1ComboUse, const DSourceComboUse& locSourceComboUseToAdd, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding)
{
	//e.g. we are grouping N pi0s and M photons (> 1) with L etas (>= 1), etc. to make combos
	//so, let's get the combos for the main grouping

	//Get combos so far
	auto locChargeContent_AllBut1 = dComboInfoChargeContent[std::get<2>(locAllBut1ComboUse)];
	auto locChargedCombo_WithNow = Get_ChargedCombo_WithNow(locChargedCombo_Presiding);
	auto& locSourceCombosByUseToSaveTo = Get_CombosSoFar(locComboingStage, dComboInfoChargeContent[std::get<2>(locComboUseToCreate)], locChargedCombo_WithNow);

	bool locGetFromSoFarFlag = (locComboingStage == d_ChargedStage) || (locChargeContent_AllBut1 != d_Charged);
	auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, locChargeContent_AllBut1, locChargedCombo_WithNow);

	vector<const DSourceCombo*> locChargedComboVector = {locChargedCombo_WithNow}; //ugh
	auto locCombos_AllBut1 = locGetFromSoFarFlag ? locSourceCombosByUseSoFar[locAllBut1ComboUse] : &locChargedComboVector; //Combos are a vector of (e.g.): -> N pi0s

	auto locChargeContent = dComboInfoChargeContent[std::get<2>(locSourceComboUseToAdd)];
	if((locComboingStage == d_ChargedStage) && (locChargeContent == d_Neutral))
	{
		//can't add neutrals, so we are already done! just copy the results to the new vector
		locSourceCombosByUseToSaveTo.emplace(locComboUseToCreate, locCombos_AllBut1);
		return;
	}

	//if on the all-showers stage, first copy over ALL fcal-only results
	if(locComboingStage == d_MixedStage)
		Copy_ZIndependentMixedResults(locComboUseToCreate, locComboingStage, locChargedCombo_WithNow);
	else //initialize vector for storing results
	{
		locSourceCombosByUseToSaveTo.emplace(locComboUseToCreate, dResourcePool_SourceComboVector.Get_Resource());
		locSourceCombosByUseToSaveTo[locComboUseToCreate]->reserve(dInitialComboVectorCapacity);
	}

	auto locDecayPID_UseToAdd = std::get<0>(locSourceComboUseToAdd);
	auto locComboInfo_UseToAdd = std::get<2>(locSourceComboUseToAdd);

	//check if on mixed stage but comboing to charged
	if((locComboingStage != d_ChargedStage) && (locChargeContent == d_Charged))
	{
		//only one valid option: locChargedCombo_WithNow: create all combos immediately
		for(auto locCombo_AllBut1 : *locCombos_AllBut1)
		{
			//first of all, get the potential combos that satisfy the RF bunches for the all-but-1 combo
			const auto& locValidRFBunches = Get_ValidRFBunches(locAllBut1ComboUse, locCombo_AllBut1);

			auto locIsZIndependent = locCombo_AllBut1->Get_IsZIndependent();
			if((locComboingStage == d_MixedStage) && locIsZIndependent)
				continue; //this combo has already been created (assuming it was valid): during the FCAL-only stage

			//get contents of the all-but-1 so that we can add to them
			auto locFurtherDecayCombos_Needed = locCombo_AllBut1->Get_FurtherDecayCombos(); //the all-but-1 combo contents by use
			auto locComboParticles = locCombo_AllBut1->Get_SourceParticles();

			//if the combo-to-add is a direct grouping of measured particles (e.g. X -> 2g), promote them
			if((locDecayPID_UseToAdd == Unknown) && locComboInfo_UseToAdd->Get_FurtherDecays().empty()) //is impossible unless in add-particles phase
			{
				auto locUsedParticlePairs_ToAdd = locChargedCombo_WithNow->Get_SourceParticles(false);
				locComboParticles.insert(locComboParticles.end(), locUsedParticlePairs_ToAdd.begin(), locUsedParticlePairs_ToAdd.end());
			}
			else //add the decay-combo to the further-decays map
			{
				//first building the further-decays for it)
				auto locFurtherDecayCombos_ToAdd = locChargedCombo_WithNow->Get_FurtherDecayCombos(); //the to-add combo contents by use
				locFurtherDecayCombos_Needed.emplace(locSourceComboUseToAdd, locFurtherDecayCombos_ToAdd[locSourceComboUseToAdd]); //add to it the new PID
			}

			//create and save it! //in creation order!
			auto locCombo = dResourcePool_SourceCombo.Get_Resource();
			locCombo->Set_Members(locComboParticles, locFurtherDecayCombos_Needed, locIsZIndependent); // create combo with all PIDs
			locSourceCombosByUseToSaveTo[locComboUseToCreate]->push_back(locCombo);
			Register_ValidRFBunches(locComboUseToCreate, locCombo, locValidRFBunches, locComboingStage, locChargedCombo_WithNow);
		}
	}

	//get the previous charged combo (if needed)
	const DSourceCombo* locChargedCombo_WithPrevious = Get_ChargedCombo_WithNow(Get_Presiding_ChargedCombo(locChargedCombo_Presiding, locSourceComboUseToAdd, locComboingStage, 1));

	//now, for each combo of all-but-1-PIDs, see which of the to-add combos we can group to it
	//valid grouping: Don't re-use a shower we've already used
	for(auto locCombo_AllBut1 : *locCombos_AllBut1)
	{
		//first of all, get the potential combos that satisfy the RF bunches for the all-but-1 combo
		const auto& locValidRFBunches_AllBut1 = Get_ValidRFBunches(locAllBut1ComboUse, locCombo_AllBut1);
		const auto& locDecayCombos_ToAdd = Get_CombosForComboing(locSourceComboUseToAdd, locComboingStage, locValidRFBunches_AllBut1, locChargedCombo_WithPrevious);

		//before we loop, first get all of the showers used to make the all-but-1 grouping, and sort it so that we can quickly search it
		auto locUsedParticles_AllBut1 = DAnalysis::Get_SourceParticles(locCombo_AllBut1->Get_SourceParticles(true)); //true: entire chain
		std::sort(locUsedParticles_AllBut1.begin(), locUsedParticles_AllBut1.end()); //must sort, because when retrieving entire chain is unsorted

		//this function will do our validity test
		auto Search_Duplicates = [&locUsedParticles_AllBut1](const JObject* locParticle) -> bool
			{return std::binary_search(locUsedParticles_AllBut1.begin(), locUsedParticles_AllBut1.end(), locParticle);};

		auto locIsZIndependent_AllBut1 = locCombo_AllBut1->Get_IsZIndependent();

		//loop over potential combos to add to the group, creating a new combo for each valid (non-duplicate) grouping
		for(const auto& locDecayCombo_ToAdd : locDecayCombos_ToAdd)
		{
			auto locIsZIndependent = (locIsZIndependent_AllBut1 && locDecayCombo_ToAdd->Get_IsZIndependent());
			if((locComboingStage == d_MixedStage) && locIsZIndependent)
				continue; //this combo has already been created (assuming it was valid): during the FCAL-only stage

			//search the all-but-1 shower vector to see if any of the showers in this combo are duplicated
			auto locUsedParticles_ToAdd = DAnalysis::Get_SourceParticles(locDecayCombo_ToAdd->Get_SourceParticles(true)); //true: entire chain

			//conduct search
			if(std::any_of(locUsedParticles_ToAdd.begin(), locUsedParticles_ToAdd.end(), Search_Duplicates))
				continue; //at least one photon was a duplicate, this combo won't work

			//no duplicates: this combo is unique.  build a new combo

			//See which RF bunches match up //guaranteed to be at least one, due to selection in Get_CombosForComboing() function
			vector<int> locValidRFBunches = {}; //if charged or massive neutrals, ignore (they don't choose at this stage)
			if(locComboingStage != d_ChargedStage)
			{
				const auto& locValidRFBunches_1 = Get_ValidRFBunches(locSourceComboUseToAdd, locDecayCombo_ToAdd);
				locValidRFBunches = Get_CommonRFBunches(locValidRFBunches_AllBut1, locValidRFBunches_1);
			}

			//get contents of the all-but-1 so that we can add to them
			auto locFurtherDecayCombos_Needed = locCombo_AllBut1->Get_FurtherDecayCombos(); //the all-but-1 combo contents by use
			auto locComboParticles = locCombo_AllBut1->Get_SourceParticles();

			//if the combo-to-add is a direct grouping of measured particles (e.g. X -> 2g), promote them
			if((locDecayPID_UseToAdd == Unknown) && locComboInfo_UseToAdd->Get_FurtherDecays().empty()) //is impossible unless in add-particles phase
			{
				auto locUsedParticlePairs_ToAdd = locDecayCombo_ToAdd->Get_SourceParticles(false);
				locComboParticles.insert(locComboParticles.end(), locUsedParticlePairs_ToAdd.begin(), locUsedParticlePairs_ToAdd.end());
			}
			else //add the decay-combo to the further-decays map
			{
				//first building the further-decays for it)
				auto locFurtherDecayCombos_ToAdd = locDecayCombo_ToAdd->Get_FurtherDecayCombos(); //the to-add combo contents by use
				locFurtherDecayCombos_Needed.emplace(locSourceComboUseToAdd, locFurtherDecayCombos_ToAdd[locSourceComboUseToAdd]); //add to it the new PID
			}

			//create and save it! //in creation order!
			auto locCombo = dResourcePool_SourceCombo.Get_Resource();
			locCombo->Set_Members(locComboParticles, locFurtherDecayCombos_Needed, locIsZIndependent); // create combo with all PIDs
			locSourceCombosByUseToSaveTo[locComboUseToCreate]->push_back(locCombo);
			Register_ValidRFBunches(locComboUseToCreate, locCombo, locValidRFBunches, locComboingStage, locChargedCombo_WithNow);
		}
	}
}

void DSourceComboer::Combo_Horizontally_AddParticle(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locAllBut1ComboUse, Particle_t locPID, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding)
{
	auto locVertexZBin = std::get<1>(locComboUseToCreate);

	//Get combos so far
	auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, d_Neutral); //if not neutral then is on charged stage: argument doesn't matter

	//e.g. we are grouping a whole bunch of particles and decays with a lone particle to make new combos
	//so, let's get the combos for this initial grouping
	auto locChargedCombo_WithNow = Get_ChargedCombo_WithNow(locChargedCombo_Presiding);
	vector<const DSourceCombo*> locChargedComboVector = {locChargedCombo_WithNow}; //ugh
	auto locChargeContent_AllBut1 = dComboInfoChargeContent[std::get<2>(locAllBut1ComboUse)];
	bool locGetFromSoFarFlag = (locComboingStage == d_ChargedStage) || (locChargeContent_AllBut1 != d_Charged);
	auto locCombos_AllBut1 = locGetFromSoFarFlag ? locSourceCombosByUseSoFar[locAllBut1ComboUse] : &locChargedComboVector; //Combos are a vector of (e.g.): -> N pi0s

	if((locComboingStage == d_ChargedStage) && (ParticleCharge(locPID) == 0))
	{
		//can't add neutrals, so we are already done! just copy the results to the new vector
		locSourceCombosByUseSoFar[locComboUseToCreate] = locCombos_AllBut1;
		return;
	}

	//if on the all-showers stage, first copy over ALL fcal-only results
	if(locComboingStage == d_MixedStage)
		Copy_ZIndependentMixedResults(locComboUseToCreate, locComboingStage, nullptr);
	else //initialize vector for storing results
	{
		locSourceCombosByUseSoFar.emplace(locComboUseToCreate, dResourcePool_SourceComboVector.Get_Resource());
		locSourceCombosByUseSoFar[locComboUseToCreate]->reserve(dInitialComboVectorCapacity);
	}

	//For checking RF bunches (ignored if on charged stage)
	const unordered_map<const JObject*, vector<int>>& locShowerRFMap = (locComboingStage == d_MixedStage_ZIndependent) ? dShowerRFBunches_FCAL : dShowerRFBunches_Both[locVertexZBin];

	//loop over the combos
	for(auto locCombo_AllBut1 : *locCombos_AllBut1)
	{
		//now, for each combo of all-but-1-PIDs, see which of the particles can group to it
		//valid grouping: Don't re-use a particle we've already used

		//before we loop, first get all of the particles of the given PID used to make the all-but-1 grouping, and sort it so that we can quickly search it
		auto locUsedParticlePairs_AllBut1 = locCombo_AllBut1->Get_SourceParticles(true);
		auto locUsedParticles_AllBut1 = DAnalysis::Get_SourceParticles(locUsedParticlePairs_AllBut1, locPID); //true: entire chain
		std::sort(locUsedParticles_AllBut1.begin(), locUsedParticles_AllBut1.end()); //necessary: may be out of order due to comboing of different decays

		//also, pre-get the further decays & FCAL-only flag, as we'll need them to build new combos
		auto locFurtherDecays = locCombo_AllBut1->Get_FurtherDecayCombos(); //the all-but-1 combo contents by use
		auto locIsZIndependent_AllBut1 = locCombo_AllBut1->Get_IsZIndependent();

		//Get potential particles for comboing
		const auto& locValidRFBunches_AllBut1 = Get_ValidRFBunches(locAllBut1ComboUse, locCombo_AllBut1);
		const auto& locParticles = Get_ParticlesForComboing(locPID, locComboingStage, locValidRFBunches_AllBut1, locVertexZBin);

		//loop over potential showers to add to the group, creating a new combo for each valid (non-duplicate) grouping
		for(const auto& locParticle : locParticles)
		{
			auto locIsZIndependent = (locComboingStage == d_MixedStage_ZIndependent) || (locIsZIndependent_AllBut1 && Get_IsZIndependent(locParticle));
			if((locComboingStage == d_MixedStage) && locIsZIndependent)
				continue; //this combo has already been created (assuming it was valid): during the FCAL-only stage

			//conduct search
			if(std::binary_search(locUsedParticles_AllBut1.begin(), locUsedParticles_AllBut1.end(), locParticle))
				continue; //this shower has already been used, this combo won't work

			//See which RF bunches match up //guaranteed to be at least one, due to selection in Get_CombosForComboing() function
			vector<int> locValidRFBunches = {}; //if charged or massive neutrals, ignore (they don't choose at this stage)
			if(locComboingStage != d_ChargedStage)
			{
				const vector<int>& locRFBunches_Second = locShowerRFMap[locParticle];
				locValidRFBunches = Get_CommonRFBunches(locValidRFBunches_AllBut1, locRFBunches_Second);
			}

			//no duplicates: this combo is unique.  build a new combo
			auto locComboParticles = locUsedParticlePairs_AllBut1;
			locComboParticles.emplace_back(locPID, locParticle);
			auto locCombo = dResourcePool_SourceCombo.Get_Resource();
			locCombo->Set_Members(locComboParticles, locFurtherDecays, locIsZIndependent); // create combo with all PIDs

			//save it! //in creation order!
			locSourceCombosByUseSoFar[locComboUseToCreate]->push_back(locCombo);
			Register_ValidRFBunches(locComboUseToCreate, locCombo, locValidRFBunches, locComboingStage, nullptr);
		}
	}
}

/***************************************************************** PARTICLE UTILITY FUNCTIONS *****************************************************************/

const vector<int>& DSourceComboer::Get_ValidRFBunches(const DSourceComboUse& locSourceComboUse, const DSourceCombo* locSourceCombo) const
{
	//in general, the valid rf bunches are only combo-dependent, not use-dependent
	//however, they ARE use-dependent if the combo contains a massive-neutral particle
	//this is because the timing is needed to get the massive-neutral momentum, which is used in invariant mass cuts for the use
	auto locRFComboIterator = dValidRFBunches_ByCombo.find(locSourceCombo);
	if(locRFComboIterator != dValidRFBunches_ByCombo.end())
		return locRFComboIterator->second; //combo-dependent (no massive neutrals)

	auto locRFComboUseIterator = dValidRFBunches_ByUse.find(std::make_pair(locSourceComboUse, locSourceCombo));
	if(locRFComboUseIterator == dValidRFBunches_ByUse.end())
		return {}; //shouldn't be possible!!!

	return locRFComboUseIterator->second; //depends on massive neutrals
}

const vector<const JObject*>& DSourceComboer::Get_ParticlesForComboing(Particle_t locPID, ComboingStage_t locComboingStage, const vector<int>& locBeamBunches, signed char locVertexZBin)
{
	//find all particles that have an overlapping beam bunch with the input
	if(ParticleCharge(locPID) != 0) //charged tracks
		return dTracksByPID[locPID]; //rf bunch & vertex-z are irrelevant
	else if(locPID != Gamma) //massive neutrals
		return dPhotonShowersByBeamBunch[0][{}]; //all neutrals: cannot do PID at all, and cannot do mass cuts until a specific vertex is chosen, so vertex-z doesn't matter
	else if(locComboingStage == d_MixedStage_ZIndependent) //fcal photons
	{
		auto locGroupBunchIterator = dFCALPhotonShowersByBeamBunch.find(locBeamBunches);
		if(locGroupBunchIterator != dFCALPhotonShowersByBeamBunch.end())
			return locGroupBunchIterator->second;
		return Get_ShowersByBeamBunch(locBeamBunches, dFCALPhotonShowersByBeamBunch);
	}
	else //bcal + fcal photons
	{
		auto locGroupBunchIterator = dPhotonShowersByBeamBunch[locVertexZBin].find(locBeamBunches);
		if(locGroupBunchIterator != dPhotonShowersByBeamBunch[locVertexZBin].end())
			return locGroupBunchIterator->second;
		return Get_ShowersByBeamBunch(locBeamBunches, dPhotonShowersByBeamBunch[locVertexZBin]);
	}
}

vector<const JObject*>* DSourceComboer::Get_ShowersByBeamBunch(const vector<int>& locBeamBunches, DPhotonShowersByBeamBunch& locShowersByBunch)
{
	//find all particles that have an overlapping beam bunch with the input
	//this won't happen often (max probably tens of times each event), so we can be a little inefficient
	vector<int> locBunchesSoFar = {*locBeamBunches.begin()};
	for(auto locBunchIterator = std::next(locBeamBunches.begin()); locBunchIterator != locBeamBunches.end(); ++locBunchIterator)
	{
		const auto& locComboShowers = locShowersByBunch[locBunchesSoFar];
		const auto& locBunchShowers = locShowersByBunch[vector<int>(*locBunchIterator)];
		locBunchesSoFar.push_back(*locBunchIterator);
		if(locBunchShowers.empty())
		{
			locShowersByBunch.emplace(locBunchesSoFar, locComboShowers);
			continue;
		}

		//merge and move-emplace
		vector<const JObject*> locMergeResult;
		locMergeResult.reserve(locComboShowers.size() + locBunchShowers.size());
		std::set_union(locComboShowers.begin(), locComboShowers.end(), locBunchShowers.begin(), locBunchShowers.end(), std::back_inserter(locMergeResult));
		locShowersByBunch.emplace(locBunchesSoFar, std::move(locMergeResult));
		Build_ParticleIterators(locBeamBunches, locShowersByBunch[locBeamBunches]);
	}
	return &(locShowersByBunch[locBeamBunches]);
}

/******************************************************************* COMBO UTILITY FUNCTIONS ******************************************************************/

void DSourceComboer::Register_ValidRFBunches(const DSourceComboUse& locSourceComboUse, const DSourceCombo* locSourceCombo, const vector<int>& locRFBunches, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_WithNow)
{
	//THE INPUT locChargedCombo MUST BE:
	//Whatever charged combo you just comboed horizontally with to make this new, mixed combo

	//search and register
	auto locComboInfo = std::get<2>(locSourceComboUse);
	if(dComboInfosWithMassiveNeutrals.find(locComboInfo) != dComboInfosWithMassiveNeutrals.end())
		dValidRFBunches_ByUse.emplace(std::make_pair(locSourceComboUse, locSourceCombo), locRFBunches); //contains a massive neutral
	else
		dValidRFBunches_ByCombo.emplace(locSourceCombo, locRFBunches);

	//also, register for each individual bunch: so that we can get valid combos for some input rf bunches later
	auto locVertexZBin = std::get<1>(locSourceComboUse);
	if(locComboingStage != d_ChargedStage)
	{
		auto& locSourceCombosByBeamBunchByUse = Get_SourceCombosByBeamBunchByUse(dComboInfoChargeContent[locComboInfo], locChargedCombo_WithNow);
		auto& locCombosByBeamBunch = locSourceCombosByBeamBunchByUse[locSourceComboUse];
		for(const auto& locBeamBunch : locRFBunches)
		{
			auto& locComboVector = locCombosByBeamBunch[{locBeamBunch}];
			locComboVector.push_back(locSourceCombo);
			dResumeSearchAfterIterators_Combos[std::make_pair(locSourceCombo, locVertexZBin)].emplace({locBeamBunch}, std::prev(locComboVector.end()));
		}
	}
	if(locRFBunches.empty()) //all //don't need to save the by-beam-bunch, but still need to save the resume-after iterator
	{
		auto& locComboVector = *(Get_CombosSoFar(locComboingStage, dComboInfoChargeContent[std::get<2>(locSourceComboUse)], locChargedCombo_WithNow)[locSourceComboUse]);
		dResumeSearchAfterIterators_Combos[std::make_pair(locSourceCombo, locVertexZBin)].emplace(locRFBunches, std::prev(locComboVector.end()));
	}
}

const vector<const DSourceCombo*>& DSourceComboer::Get_CombosForComboing(const DSourceComboUse& locComboUse, ComboingStage_t locComboingStage, const vector<int>& locBeamBunches, const DSourceCombo* locChargedCombo_WithPrevious)
{
	//THE INPUT locChargedCombo MUST BE:
	//Whatever charged combo you PREVIOUSLY comboed horizontally with to make the combos you're trying to get

	//find all combos for the given use that have an overlapping beam bunch with the input
	auto locChargeContent = dComboInfoChargeContent[std::get<2>(locComboUse)];
	if(locBeamBunches.empty() || (locChargeContent == d_Charged)) //e.g. fully charged, or a combo of 2 KLongs (RF bunches not saved for massive neutrals)
		return Get_CombosSoFar(locComboingStage, locChargeContent, locChargedCombo_WithPrevious)[locComboUse];

	auto& locSourceCombosByBeamBunchByUse = Get_SourceCombosByBeamBunchByUse(locChargeContent, locChargedCombo_WithPrevious);
	auto locGroupBunchIterator = locSourceCombosByBeamBunchByUse[locComboUse].find(locBeamBunches);
	if(locGroupBunchIterator != locSourceCombosByBeamBunchByUse[locComboUse].end())
		return locGroupBunchIterator->second;
	return Get_CombosByBeamBunch(locSourceCombosByBeamBunchByUse[locComboUse], locBeamBunches, locComboingStage, std::get<1>(locComboUse));
}

const vector<const DSourceCombo*>& DSourceComboer::Get_CombosByBeamBunch(DCombosByBeamBunch& locCombosByBunch, const vector<int>& locBeamBunches, ComboingStage_t locComboingStage, signed char locVertexZBin)
{
	//find all combos for the given use that have an overlapping beam bunch with the input
	//this shouldn't be called very many times per event
	vector<int> locBunchesSoFar = {*locBeamBunches.begin()};
	for(auto locBunchIterator = std::next(locBeamBunches.begin()); locBunchIterator != locBeamBunches.end(); ++locBunchIterator)
	{
		const auto& locComboShowers = locCombosByBunch[locBunchesSoFar];
		const auto& locBunchShowers = locCombosByBunch[vector<int>(*locBunchIterator)];
		locBunchesSoFar.push_back(*locBunchIterator);
		if(locBunchShowers.empty())
		{
			locCombosByBunch.emplace(locBunchesSoFar, locComboShowers);
			continue;
		}

		//merge and move-emplace
		vector<const DSourceCombo*> locMergeResult;
		locMergeResult.reserve(locComboShowers.size() + locBunchShowers.size());
		std::set_union(locComboShowers.begin(), locComboShowers.end(), locBunchShowers.begin(), locBunchShowers.end(), std::back_inserter(locMergeResult));
		locCombosByBunch.emplace(locBunchesSoFar, std::move(locMergeResult));
		Build_ComboIterators(locBeamBunches, locCombosByBunch[locBeamBunches], locComboingStage, locVertexZBin);
	}
	return &(locCombosByBunch[locBeamBunches]);
}

void DSourceComboer::Copy_ZIndependentMixedResults(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_WithNow)
{
	//Copy the results from the FCAL-only stage through to the both stage (that way we don't have to repeat them)

	//THE INPUT locChargedCombo MUST BE:
	//Whatever charged combo you are about to combo horizontally with to make this new, mixed combo

	//Get combos so far
	auto locVertexZBin = std::get<1>(locComboUseToCreate);
	auto locComboInfo = std::get<2>(locComboUseToCreate);
	auto locChargeContent = dComboInfoChargeContent[locComboInfo];
	auto& locSourceCombosByUseSoFar = Get_CombosSoFar(locComboingStage, locChargeContent, locChargedCombo_WithNow);

	//Get the combo vectors
	auto locComboUseFCAL = std::make_tuple(std::get<0>(locComboUseToCreate), DSourceComboInfo::Get_VertexZIndex_FCAL(), locComboInfo);
	const auto& locFCALComboVector = *(locSourceCombosByUseSoFar[locComboUseFCAL]);
	auto& locBothComboVector = *(locSourceCombosByUseSoFar[locComboUseToCreate]);

	//Copy over the combos
	locBothComboVector.reserve(locFCALComboVector.size() + dInitialComboVectorCapacity);
	locBothComboVector.assign(locFCALComboVector.begin(), locFCALComboVector.end());

	//Copy over the combos-by-beam-bunch
	auto& locSourceCombosByBeamBunchByUse = Get_SourceCombosByBeamBunchByUse(locChargeContent, locChargedCombo_WithNow);

	const auto& locCombosByBeamBunch = locSourceCombosByBeamBunchByUse[locComboUseFCAL];
	for(const auto& locComboBeamBunchPair : locCombosByBeamBunch)
	{
		if(locComboBeamBunchPair.first.size() == 1) //don't copy the overlap ones: they are not complete & need to be filled on the fly
			locSourceCombosByBeamBunchByUse[locComboUseToCreate].emplace(locComboBeamBunchPair);
	}

	//Copy over the resume-after iterators
	for(vector<const DSourceCombo*>::const_iterator locComboIterator = locBothComboVector.begin(); locComboIterator != locBothComboVector.end(); ++locComboIterator)
	{
		const auto& locRFBunches = Get_ValidRFBunches(locComboUseToCreate, *locComboIterator);
		for(const auto& locBeamBunch : locRFBunches)
			dResumeSearchAfterIterators_Combos[std::make_pair(*locComboIterator, locVertexZBin)].emplace(locBeamBunch, locComboIterator);
		if(locRFBunches.empty()) //all
			dResumeSearchAfterIterators_Combos[std::make_pair(*locComboIterator, locVertexZBin)].emplace(locRFBunches, locComboIterator);
	}
}

const DSourceCombo* DSourceComboer::Get_StepSourceCombo(const DReaction* locReaction, size_t locDesiredStepIndex, const DSourceCombo* locSourceCombo_Current, size_t locCurrentStepIndex)
{
	//Get the list of steps we need to traverse //particle pair: step index, particle instance index
	vector<pair<size_t, int>> locParticleIndices = {std::make_pair(locDesiredStepIndex, DReactionStep::Get_ParticleIndex_Initial())};
	while(locParticleIndices.back().first != locCurrentStepIndex)
	{
		auto locParticlePair = DAnalysis::Get_InitialParticleDecayFromIndices(locReaction, locParticleIndices.back().first);
		auto locStep = locReaction->Get_ReactionStep(locParticlePair.first);
		auto locInstanceIndex = DAnalysis::Get_ParticleInstanceIndex(locStep, locParticlePair.second);
		locParticleIndices.emplace_back(locParticlePair.first, locInstanceIndex);
	}

	//start from back of locParticleIndices, searching
	while(true)
	{
		auto locNextStep = locParticleIndices[locParticleIndices.size() - 2].first;
		auto locInstanceToFind = locParticleIndices.back().second;
		const auto& locUseToFind = dSourceComboUseReactionStepMap[locReaction][locNextStep];
		locSourceCombo_Current = DAnalysis::Find_Combo_AtThisStep(locSourceCombo_Current, locUseToFind, locInstanceToFind);
		if(locSourceCombo_Current == nullptr)
			return nullptr; //e.g. entirely neutral step when input is charged
		if(locNextStep == locDesiredStepIndex)
			return locSourceCombo_Current;
		locParticleIndices.pop_back();
	}

	return nullptr;
}

const DSourceCombo* DSourceComboer::Get_Presiding_ChargedCombo(const DSourceCombo* locChargedCombo_Presiding, const DSourceComboUse& locNextComboUse, ComboingStage_t locComboingStage, size_t locInstance) const
{
	//locInstance starts from ONE!!
	if(locComboingStage == d_ChargedStage)
		return nullptr;
	if(locChargedCombo_Presiding == nullptr)
		return nullptr;
	if(dComboInfoChargeContent[std::get<2>(locNextComboUse)] != d_AllCharges)
		return nullptr; //not needed

	auto locFurtherDecayCombos = locChargedCombo_Presiding->Get_FurtherDecayCombos();

	auto locUseToFind = (locComboingStage == d_MixedStage_ZIndependent) ? locNextComboUse : dZDependentUseToIndependentMap.find(locNextComboUse)->second;
	auto locUseIterator = locFurtherDecayCombos.find(locUseToFind);

	//check if the use you are looking for is a temporary (e.g. vertical grouping of 2KShorts when comboing horizontally)
	if(locUseIterator == locFurtherDecayCombos.end())
		return locChargedCombo_Presiding; //temporary: the presiding is still the same!

	//get the vector of potential charged combos
	auto locNextChargedComboVector = locUseIterator->second;

	//if on z-independent, don't need to do anything fancy, just return the requested instance
	if(locComboingStage == d_MixedStage_ZIndependent)
		return locNextChargedComboVector[locInstance - 1];

	//there might be multiple combos (e.g. K0 decays), each at a different vertex-z
	//so, we must retrieve the N'th charged combo with the correct vertex-z bin
	size_t locCount = 0;
	auto locDesiredVertexZBin = std::get<1>(locNextComboUse);
	for(auto locNextPotentialCombo : locNextChargedComboVector)
	{
		auto locNextVertexZBin = dSourceComboVertexer.Get_VertexZBin(locNextPotentialCombo);
		if(locNextVertexZBin != locDesiredVertexZBin)
			continue;
		if(++locCount == locInstance)
			return locNextPotentialCombo;
	}

	return nullptr; //uh oh ...
}

} //end DAnalysis namespace
