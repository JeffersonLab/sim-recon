#include "ANALYSIS/DSourceComboTimeHandler.h"
#include "ANALYSIS/DSourceComboer.h"


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

namespace DAnalysis
{

DSourceComboTimeHandler::DSourceComboTimeHandler(JEventLoop* locEventLoop, const DSourceComboer* locSourceComboer, const DSourceComboVertexer* locSourceComboVertexer) :
		dSourceComboer(locSourceComboer), dSourceComboVertexer(locSourceComboVertexer)
{
	//UTILITIES
	locEventLoop->GetSingle(dAnalysisUtilities);

	//BEAM BUNCH PERIOD
	vector<double> locBeamPeriodVector;
	locEventLoop->GetCalib("PHOTON_BEAM/RF/beam_period", locBeamPeriodVector);
	dBeamBunchPeriod = locBeamPeriodVector[0];

	//These functions can have the same name because we are no longer adding them to the global ROOT list of functions

	// Timing Cuts: Photon
	dPIDTimingCuts[Gamma].emplace(SYS_BCAL, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[Gamma][SYS_BCAL]->SetParameter(0, 3.0);
	dPIDTimingCuts[Gamma].emplace(SYS_FCAL, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[Gamma][SYS_FCAL]->SetParameter(0, 2.5);

	// Timing Cuts: Leptons
	dPIDTimingCuts[Electron].emplace(SYS_BCAL, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[Electron][SYS_BCAL]->SetParameter(0, 1.0);
	dPIDTimingCuts[Electron].emplace(SYS_FCAL, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[Electron][SYS_FCAL]->SetParameter(0, 2.5);
	dPIDTimingCuts[Electron].emplace(SYS_TOF, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[Electron][SYS_TOF]->SetParameter(0, 2.5);

	dPIDTimingCuts.emplace(Positron, dPIDTimingCuts[Electron]);
	dPIDTimingCuts.emplace(MuonMinus, dPIDTimingCuts[Electron]);
	dPIDTimingCuts.emplace(MuonPlus, dPIDTimingCuts[Electron]);

	// Timing Cuts: Mesons
	dPIDTimingCuts[PiPlus].emplace(SYS_BCAL, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[PiPlus][SYS_BCAL]->SetParameter(0, 2.0);
	dPIDTimingCuts[PiPlus].emplace(SYS_FCAL, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[PiPlus][SYS_FCAL]->SetParameter(0, 2.5);
	dPIDTimingCuts[PiPlus].emplace(SYS_TOF, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[PiPlus][SYS_TOF]->SetParameter(0, 2.5);
	dPIDTimingCuts.emplace(PiMinus, dPIDTimingCuts[PiPlus]);

	dPIDTimingCuts[KPlus].emplace(SYS_BCAL, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[KPlus][SYS_BCAL]->SetParameter(0, 0.75);
	dPIDTimingCuts[KPlus].emplace(SYS_FCAL, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[KPlus][SYS_FCAL]->SetParameter(0, 2.5);
	dPIDTimingCuts[KPlus].emplace(SYS_TOF, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[KPlus][SYS_TOF]->SetParameter(0, 2.0);
	dPIDTimingCuts.emplace(KMinus, dPIDTimingCuts[KPlus]);

	// Timing Cuts: Baryons
	dPIDTimingCuts[Proton].emplace(SYS_BCAL, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[Proton][SYS_BCAL]->SetParameter(0, 2.5);
	dPIDTimingCuts[Proton].emplace(SYS_FCAL, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[Proton][SYS_FCAL]->SetParameter(0, 2.5);
	dPIDTimingCuts[Proton].emplace(SYS_TOF, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[Proton][SYS_TOF]->SetParameter(0, 2.5);
	dPIDTimingCuts.emplace(AntiProton, dPIDTimingCuts[Proton]);
}

void DSourceComboTimeHandler::Setup_NeutralShowers(const vector<const DNeutralShower*>& locNeutralShowers, const DEventRFBunch* locInitialEventRFBunch)
{
	//Precompute a few things for the neutral showers, before comboing
	//Even if it turns out some of this isn't technically needed,
	//it's still faster than doing a check to see if this has been done or not for every single photon-combo request

	//GET RF BUNCH
	dInitialEventRFBunch = locInitialEventRFBunch;

	//ARRANGE NEUTRAL SHOWERS
	vector<const DNeutralShower*> locBCALShowers, locFCALShowers;
	for(auto& locShower : locNeutralShowers)
	{
		auto& locContainer = (locShower->dDetectorSystem == SYS_BCAL) ? locBCALShowers : locFCALShowers;
		locContainer.push_back(locShower);
	}

	//CALCULATE KINEMATICS
	//FCAL: at target center
	auto locFCALZBin = DSourceComboInfo::Get_VertexZIndex_FCAL();
	for(auto& locShower : locFCALShowers)
		dPhotonKinematics[locFCALZBin].emplace(locShower, Create_KinematicData_Photon(locShower, dTargetCenter));

	//BCAL: in vertex-z bins
	for(size_t loc_i = 0; loc_i < dNumPhotonVertexZBins; ++loc_i)
	{
		DVector3 locBinCenter(0.0, 0.0, dSourceComboer->Get_PhotonVertexZBinCenter(loc_i));
		for(auto& locShower : locBCALShowers)
			dPhotonKinematics[loc_i].emplace(locShower, Create_KinematicData_Photon(locShower, locBinCenter));
	}

	//DETERMINE WHICH RF BUNCHES ARE VALID
	//FCAL: at target center
	for(auto& locShower : locFCALShowers)
		Calc_PhotonBeamBunchShifts(locShower, dPhotonKinematics[locFCALZBin][locShower], dInitialEventRFBunch->dTime, locFCALZBin);

	//BCAL + FCAL: in vertex-z bins
	for(size_t loc_i = 0; loc_i < dNumPhotonVertexZBins; ++loc_i)
	{
		//propagate RF time to vertex position
		double locPropagatedRFTime = dInitialEventRFBunch->dTime + (dSourceComboer->Get_PhotonVertexZBinCenter(loc_i) - dTargetCenter.Z())/SPEED_OF_LIGHT;
		for(auto& locShower : locBCALShowers)
			Calc_PhotonBeamBunchShifts(locShower, dPhotonKinematics[loc_i][locShower], locPropagatedRFTime, loc_i);

		//insert the previously-done FCAL photons
		dShowersByBeamBunchByZBin[loc_i].insert(dShowersByBeamBunchByZBin[locFCALZBin].begin(), dShowersByBeamBunchByZBin[locFCALZBin].end());
		for(auto locBeamBunchPair : dShowersByBeamBunchByZBin[locFCALZBin])
		{
			const auto& locShowers = locBeamBunchPair.second;
			auto locBothIterator = dShowersByBeamBunchByZBin[loc_i].find(locBeamBunchPair.first);
			if(locBothIterator != dShowersByBeamBunchByZBin[loc_i].end())
			{
				auto locPhotonVector = locBothIterator->second;
				locPhotonVector.insert(locPhotonVector.end(), locShowers.begin(), locShowers.end());
			}
			else
				dShowersByBeamBunchByZBin[loc_i].emplace(locBeamBunchPair);
		}
	}
}

void DSourceComboTimeHandler::Calc_PhotonBeamBunchShifts(const DNeutralShower* locNeutralShower, shared_ptr<const DKinematicData>& locKinematicData, double locRFTime, signed char locZBin)
{
	//get delta-t cut
	DetectorSystem_t locSystem = locNeutralShower->dDetectorSystem;
	double locDeltaTCut = dPIDTimingCuts[Gamma][locSystem]->Eval(locNeutralShower->dEnergy) + Calc_MaxDeltaTError(locNeutralShower, locKinematicData);

	//do loop over possible #-RF-shifts
	double locVertexTime = locKinematicData->time();
	auto locRFShifts = Calc_BeamBunchShifts(locVertexTime, locRFTime, locDeltaTCut, true);

	auto locJObject = static_cast<const JObject*>(locNeutralShower);
	dShowerRFBunches[locZBin].emplace(locJObject, locRFShifts);
	for(auto locNumShifts : locRFShifts)
		dShowersByBeamBunchByZBin[locZBin][{locNumShifts}].push_back(locJObject);
}

vector<int> DSourceComboTimeHandler::Calc_BeamBunchShifts(double locVertexTime, double locPropagatedRFTime, double locDeltaTCut, bool locIncludeDecayTimeOffset) const
{
	auto locOrigNumShifts = Calc_RFBunchShift(locPropagatedRFTime, locVertexTime); //get best shift
	vector<int> locRFShifts;

	//start with best-shift, then loop up until fails cut
	int locNumShifts = locOrigNumShifts;
	double locDeltaT = locVertexTime - (locPropagatedRFTime + locNumShifts*dBeamBunchPeriod);
	while(fabs(locDeltaT) < locDeltaTCut)
	{
		locRFShifts.push_back(locNumShifts);
		++locNumShifts;
		locDeltaT = locVertexTime - (locPropagatedRFTime + locNumShifts*dBeamBunchPeriod);
	}

	//now loop down in n-shifts
	locNumShifts = locOrigNumShifts - 1;
	locDeltaT = locVertexTime - (locPropagatedRFTime + locNumShifts*dBeamBunchPeriod);
	while(fabs(locDeltaT) < locDeltaTCut)
	{
		locRFShifts.push_back(locNumShifts);
		--locNumShifts;
		locDeltaT = locVertexTime - (locPropagatedRFTime + locNumShifts*dBeamBunchPeriod);
	}

	//due to detached vertices, we may need to accept EARLIER RF bunches
	if(locIncludeDecayTimeOffset)
	{
		//continue down-shift loop, this time including time offset
		//+dMaxTimeOffset: takes longer for "RF" time to get there (due to slow decaying particle)
		locDeltaT = locVertexTime - (locPropagatedRFTime + locNumShifts*dBeamBunchPeriod + dMaxDecayTimeOffset);
		while(fabs(locDeltaT) < locDeltaTCut)
		{
			locRFShifts.push_back(locNumShifts);
			--locNumShifts;
			locDeltaT = locVertexTime - (locPropagatedRFTime + locNumShifts*dBeamBunchPeriod + dMaxDecayTimeOffset);
		}
	}

	std::sort(locRFShifts.begin(), locRFShifts.end());
	return locRFShifts;
}

double DSourceComboTimeHandler::Calc_MaxDeltaTError(const DNeutralShower* locNeutralShower, const shared_ptr<const DKinematicData>& locKinematicData) const
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

vector<int> DSourceComboTimeHandler::Select_RFBunches_Charged(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionChargedCombo) const
{
	auto locRFIterator = dChargedComboRFBunches.find(locReactionChargedCombo);
	if(locRFIterator != dChargedComboRFBunches.end())
		return locRFIterator->second; //already computed, return results!!

	//All charged tracks vote, even those not at the primary vertex
	//loop over vertices, get all charged particles at that vertex, utilize that + time offset

	//loop over vertices
	vector<int> locValidRFBunches;
	for(auto locStepVertexInfo : locReactionVertexInfo->Get_StepVertexInfos())
	{
		auto locIsProductionVertex = locStepVertexInfo->Get_ProductionVertexFlag();
		auto locVertexPrimaryCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionChargedCombo, locStepVertexInfo);
		if(!dSourceComboVertexer->Get_VertexDeterminableWithCharged(locIsProductionVertex, locVertexPrimaryCombo))
			continue; //vertex position indeterminate at this stage: don't include these tracks

		//get combo, vertex position, and time offset from RF bunch
		auto locVertex = dSourceComboVertexer->Get_Vertex(locIsProductionVertex, locVertexPrimaryCombo);
		auto locTimeOffset = dSourceComboVertexer->Get_TimeOffset(locIsProductionVertex, locReactionChargedCombo, locVertexPrimaryCombo);
		double locPropagatedRFTime = dInitialEventRFBunch->dTime + (locVertex.Z() - dTargetCenter.Z())/SPEED_OF_LIGHT;

		//loop over charged particles
		auto locChargedParticles = DAnalysis::Get_SourceParticles_ThisVertex(locVertexPrimaryCombo);
		for(auto locParticlePair : locChargedParticles)
		{
			auto locPID = locParticlePair.first;
			auto locChargedTrack = static_cast<const DChargedTrack*>(locParticlePair.second);
			auto locHypothesis = locChargedTrack->Get_Hypothesis(locPID);

			//evaluate timing at the POCA to the vertex
			auto locX4(locHypothesis->x4());
			auto locP4(locHypothesis->lorentzMomentum());

			auto locPOCAPair = std::make_pair(locHypothesis, dSourceComboVertexer->Get_VertexParticles(locIsProductionVertex, locVertexPrimaryCombo));
			auto locPOCAIterator = dChargedParticlePOCAToVertexX4.find(locPOCAPair);
			if(locPOCAIterator != dChargedParticlePOCAToVertexX4.end())
				locX4 = locPOCAIterator->second;
			else //do this perhaps time intensive propagation
			{
				dAnalysisUtilities->Propagate_Track(locHypothesis->charge(), locVertex, locX4, locP4, nullptr);
				dChargedParticlePOCAToVertexX4.emplace(locPOCAPair, locX4); //save results so we don't have to do it again
			}

			auto locVertexTime = locX4.T() - locTimeOffset;
			auto locSystem = locHypothesis->t1_detector();
			auto locDeltaTCut = dPIDTimingCuts[locPID][locSystem]->Eval(locP4.P());

			//do loop over possible #-RF-shifts
			auto locParticleRFBunches = Calc_BeamBunchShifts(locVertexTime, locPropagatedRFTime, locDeltaTCut, false);
			if(locParticleRFBunches.empty())
				return {};
			locValidRFBunches = Get_CommonRFBunches(locValidRFBunches, locParticleRFBunches);
			if(locValidRFBunches.empty())
				return {};
		}
	}

	dChargedComboRFBunches.emplace(locReactionChargedCombo, locValidRFBunches);
	return locValidRFBunches;
}

int DSourceComboTimeHandler::Select_RFBunch_Full(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, const DSourceCombo* locReactionChargedCombo, const vector<int>& locRFBunches)
{
	//initialize chisq's
	unordered_map<int, double> locChiSqByRFBunch;
	for(auto locRFBunch : locRFBunches)
		locChiSqByRFBunch.emplace(locRFBunch, 0.0);

	//loop over vertices
	auto locPrimaryVertexZ = dSourceComboVertexer->Get_PrimaryVertex(locReactionVertexInfo, locReactionChargedCombo).Z();
	for(auto locStepVertexInfo : locReactionVertexInfo->Get_StepVertexInfos())
	{
		if(locStepVertexInfo->Get_DanglingVertexFlag())
			continue; //unknown position!

		auto locIsProductionVertex = locStepVertexInfo->Get_ProductionVertexFlag();
		auto locVertexPrimaryChargedCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionChargedCombo, locStepVertexInfo);
		auto locVertexPrimaryFullCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionFullCombo, locStepVertexInfo);
		if(!dSourceComboVertexer->Get_VertexDeterminableWithCharged(locIsProductionVertex, locVertexPrimaryChargedCombo))
			continue; //vertex position indeterminate at this stage: don't include these particles

		//get combo, vertex position, and time offset from RF bunch
		auto locVertex = dSourceComboVertexer->Get_Vertex(locIsProductionVertex, locVertexPrimaryChargedCombo);
		auto locTimeOffset = dSourceComboVertexer->Get_TimeOffset(locIsProductionVertex, locReactionChargedCombo, locVertexPrimaryChargedCombo);
		auto locParticles = DAnalysis::Get_SourceParticles_ThisVertex(locVertexPrimaryFullCombo);

		for(auto locRFBunch : locRFBunches)
		{
			//propagate rf time to vertex and add time offset (faster to just do it here rather than for each particle)
			auto locPropagatedRFTime = Calc_PropagatedRFTime(locPrimaryVertexZ, locRFBunch, locTimeOffset);

			//loop over all particles
			for(auto locParticlePair : locParticles)
			{
				auto locPID = locParticlePair.first;
				if((ParticleCharge(locPID) == 0) && (ParticleMass(locPID) > 0.0))
					continue; //ignore massive neutrals: timing defines their momentum, cannot be used

				if(ParticleCharge(locPID) == 0)
				{
					auto locNeutralShower = static_cast<const DNeutralShower*>(locParticlePair.second);
					locChiSqByRFBunch[locRFBunch] += Calc_RFDeltaTChiSq(locNeutralShower, locRFBunch, locVertex, locPropagatedRFTime);
				}
				else //charged
				{
					auto locChargedTrack = static_cast<const DChargedTrack*>(locParticlePair.second);
					auto locHypothesis = locChargedTrack->Get_Hypothesis(locPID);

					//get the timing at the POCA to the vertex (computed previously!)
					auto locPOCAPair = std::make_pair(locHypothesis, dSourceComboVertexer->Get_VertexParticles(locIsProductionVertex, locVertexPrimaryChargedCombo));
					auto locVertexTime = dChargedParticlePOCAToVertexX4.find(locPOCAPair)->second.T();
					locChiSqByRFBunch[locRFBunch] += Calc_RFDeltaTChiSq(locHypothesis, locRFBunch, locVertexTime, locPropagatedRFTime);
				}
			}
		}
	}

	//ok, total chisq's are computed, pick the one that is the best!
	auto Compare_RFChiSqs = [](const pair<int, double>& lhs, const pair<int, double>& rhs) -> bool {return lhs.second < rhs.second;};
	return std::max_element(locChiSqByRFBunch.begin(), locChiSqByRFBunch.end(), Compare_RFChiSqs)->first;
}

}

#endif // DSourceComboTimeHandler_h
