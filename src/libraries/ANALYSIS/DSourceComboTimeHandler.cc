#include "ANALYSIS/DSourceComboTimeHandler.h"
#include "ANALYSIS/DSourceComboer.h"
#include "ANALYSIS/DSourceComboVertexer.h"


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

//What about particles with unknown vertices? They shouldn't vote (unless those are the only particles), but still need to do a PID cut

//Charged tracks:
//If vertex is known, can select RF bunch and apply final (only) cut and fill hists:
//If vertex is unknown, do PID cuts at end (after beam selection) using the best vertex available: NEED A NEW FUNCTION

//Photons:
//Initially, without vertex: Hist/cut separately (as unknown)
//If vertex is known, can select RF bunch and apply final (only) cut and fill hists:
//If vertex is unknown, do PID cuts at end (after beam selection) using the best vertex available: NEED A NEW FUNCTION

namespace DAnalysis
{

DSourceComboTimeHandler::DSourceComboTimeHandler(JEventLoop* locEventLoop, DSourceComboer* locSourceComboer, const DSourceComboVertexer* locSourceComboVertexer) :
		dSourceComboer(locSourceComboer), dSourceComboVertexer(locSourceComboVertexer)
{
	gPARMS->SetDefaultParameter("COMBO:DEBUG_LEVEL", dDebugLevel);

	//UTILITIES
	locEventLoop->GetSingle(dAnalysisUtilities);

	//GET THE GEOMETRY
	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry* locGeometry = locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());

	//TARGET INFORMATION
	double locTargetCenterZ = 65.0;
	locGeometry->GetTargetZ(locTargetCenterZ);
	dTargetCenter.SetXYZ(0.0, 0.0, locTargetCenterZ);

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
	dSelectedRFDeltaTs[Gamma][SYS_BCAL].reserve(1000);
	dSelectedRFDeltaTs[Gamma][SYS_FCAL].reserve(1000);

	//Unknown: initial RF selection for photons (at beginning of event, prior to vertex) //can be separate cut function
	dPIDTimingCuts.emplace(Unknown, dPIDTimingCuts[Gamma]);
	dSelectedRFDeltaTs.emplace(Unknown, dSelectedRFDeltaTs[Gamma]);

	// Timing Cuts: Leptons
	dPIDTimingCuts[Electron].emplace(SYS_BCAL, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[Electron][SYS_BCAL]->SetParameter(0, 1.0);
	dPIDTimingCuts[Electron].emplace(SYS_TOF, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[Electron][SYS_TOF]->SetParameter(0, 2.5);
	dPIDTimingCuts[Electron].emplace(SYS_FCAL, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[Electron][SYS_FCAL]->SetParameter(0, 2.5);
	dSelectedRFDeltaTs[Electron][SYS_BCAL].reserve(1000);
	dSelectedRFDeltaTs[Electron][SYS_FCAL].reserve(1000);
	dSelectedRFDeltaTs[Electron][SYS_TOF].reserve(1000);
	dSelectedRFDeltaTs[Electron][SYS_START].reserve(1000);

	dPIDTimingCuts.emplace(Positron, dPIDTimingCuts[Electron]);
	dPIDTimingCuts.emplace(MuonMinus, dPIDTimingCuts[Electron]);
	dPIDTimingCuts.emplace(MuonPlus, dPIDTimingCuts[Electron]);
	dSelectedRFDeltaTs.emplace(Positron, dSelectedRFDeltaTs[Electron]);
	dSelectedRFDeltaTs.emplace(MuonMinus, dSelectedRFDeltaTs[Electron]);
	dSelectedRFDeltaTs.emplace(MuonPlus, dSelectedRFDeltaTs[Electron]);

	// Timing Cuts: Mesons
	dPIDTimingCuts[PiPlus].emplace(SYS_BCAL, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[PiPlus][SYS_BCAL]->SetParameter(0, 2.0);
	dPIDTimingCuts[PiPlus].emplace(SYS_TOF, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[PiPlus][SYS_TOF]->SetParameter(0, 2.5);
	dPIDTimingCuts[PiPlus].emplace(SYS_FCAL, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[PiPlus][SYS_FCAL]->SetParameter(0, 2.5);
	dPIDTimingCuts.emplace(PiMinus, dPIDTimingCuts[PiPlus]);
	dSelectedRFDeltaTs.emplace(PiPlus, dSelectedRFDeltaTs[Electron]);
	dSelectedRFDeltaTs.emplace(PiMinus, dSelectedRFDeltaTs[Electron]);

	dPIDTimingCuts[KPlus].emplace(SYS_BCAL, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[KPlus][SYS_BCAL]->SetParameter(0, 0.75);
	dPIDTimingCuts[KPlus].emplace(SYS_TOF, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[KPlus][SYS_TOF]->SetParameter(0, 2.0);
	dPIDTimingCuts[KPlus].emplace(SYS_FCAL, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[KPlus][SYS_FCAL]->SetParameter(0, 2.5);
	dPIDTimingCuts.emplace(KMinus, dPIDTimingCuts[KPlus]);
	dSelectedRFDeltaTs.emplace(KPlus, dSelectedRFDeltaTs[Electron]);
	dSelectedRFDeltaTs.emplace(KMinus, dSelectedRFDeltaTs[Electron]);

	// Timing Cuts: Baryons
	dPIDTimingCuts[Proton].emplace(SYS_BCAL, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[Proton][SYS_BCAL]->SetParameter(0, 2.5);
	dPIDTimingCuts[Proton].emplace(SYS_TOF, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[Proton][SYS_TOF]->SetParameter(0, 2.0);
	dPIDTimingCuts[Proton].emplace(SYS_FCAL, new TF1("df_TimeCut", "[0]", 0.0, 12.0));
	dPIDTimingCuts[Proton][SYS_FCAL]->SetParameter(0, 2.5);
	dPIDTimingCuts.emplace(AntiProton, dPIDTimingCuts[Proton]);
	dSelectedRFDeltaTs.emplace(Proton, dSelectedRFDeltaTs[Electron]);
	dSelectedRFDeltaTs.emplace(AntiProton, dSelectedRFDeltaTs[Electron]);
	dAllRFDeltaTs = dSelectedRFDeltaTs;

	vector<DetectorSystem_t> locTimingSystems_Charged {SYS_TOF, SYS_BCAL, SYS_FCAL, SYS_START};
	vector<DetectorSystem_t> locTimingSystems_Neutral {SYS_BCAL, SYS_FCAL};
	vector<Particle_t> locPIDs {Unknown, Gamma, Electron, Positron, MuonPlus, MuonMinus, PiPlus, PiMinus, KPlus, KMinus, Proton, AntiProton};

	//get file name
	string locOutputFileName = "hd_root.root";
	if(gPARMS->Exists("OUTPUT_FILENAME"))
		gPARMS->GetParameter("OUTPUT_FILENAME", locOutputFileName);

	//CREATE HISTOGRAMS
	japp->RootWriteLock(); //to prevent undefined behavior due to directory changes, etc.
	{
		//get and change to the base (file/global) directory
		TDirectory* locCurrentDir = gDirectory;
		TFile* locFile = (TFile*)gROOT->FindObject(locOutputFileName.c_str());
		if(locFile != NULL)
			locFile->cd("");
		else
			gDirectory->cd("/");

		string locDirName = "Independent";
		TDirectoryFile* locDirectoryFile = static_cast<TDirectoryFile*>(gDirectory->GetDirectory(locDirName.c_str()));
		if(locDirectoryFile == NULL)
			locDirectoryFile = new TDirectoryFile(locDirName.c_str(), locDirName.c_str());
		locDirectoryFile->cd();

		locDirName = "Combo_Construction";
		locDirectoryFile = static_cast<TDirectoryFile*>(gDirectory->GetDirectory(locDirName.c_str()));
		if(locDirectoryFile == NULL)
			locDirectoryFile = new TDirectoryFile(locDirName.c_str(), locDirName.c_str());
		locDirectoryFile->cd();

		for(auto locPID : locPIDs)
		{
			auto locPIDString = string((locPID != Unknown) ? ParticleType(locPID) : "Photons_PreVertex");

			locDirName = string("PID_") + locPIDString;
			locDirectoryFile = static_cast<TDirectoryFile*>(gDirectory->GetDirectory(locDirName.c_str()));
			if(locDirectoryFile == NULL)
				locDirectoryFile = new TDirectoryFile(locDirName.c_str(), locDirName.c_str());
			locDirectoryFile->cd();

			auto& locTimingSystems = (ParticleCharge(locPID) == 0) ? locTimingSystems_Neutral : locTimingSystems_Charged;
			for(auto locSystem : locTimingSystems)
			{
				auto locHistName = string("All_RFDeltaTVsP_") + string(SystemName(locSystem));
				auto locHist = gDirectory->Get(locHistName.c_str());
				if(locHist == nullptr)
				{
					auto locHistTitle = string((locPID != Unknown) ? ParticleName_ROOT(locPID) : "Photons_PreVertex");
					locHistTitle += string(" Candidates, ") + string(SystemName(locSystem)) + string(";p (GeV/c);#Deltat_{Particle - All RFs}");
					dHistMap_RFDeltaTVsP_AllRFs[locPID][locSystem] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 400, 0.0, 12.0, 1400, -7.0, 7.0);
				}
				else
					dHistMap_RFDeltaTVsP_AllRFs[locPID][locSystem] = static_cast<TH2*>(locHist);

				if(locPID == Unknown)
					continue;

				locHistName = string("Best_RFDeltaTVsP_") + string(SystemName(locSystem));
				locHist = gDirectory->Get(locHistName.c_str());
				if(locHist == nullptr)
				{
					auto locHistTitle = string((locPID != Unknown) ? ParticleName_ROOT(locPID) : "Photons_PreVertex");
					locHistTitle += string(" Candidates, ") + string(SystemName(locSystem)) + string(";p (GeV/c);#Deltat_{Particle - Best RF}");
					dHistMap_RFDeltaTVsP_BestRF[locPID][locSystem] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 400, 0.0, 12.0, 1400, -7.0, 7.0);
				}
				else
					dHistMap_RFDeltaTVsP_BestRF[locPID][locSystem] = static_cast<TH2*>(locHist);
			}
			gDirectory->cd("..");
		}

		//Beam-RF delta-t
		string locHistName = "BeamRFDeltaTVsBeamE";
		auto locHist = gDirectory->Get(locHistName.c_str());
		if(locHist == nullptr)
			dHist_BeamRFDeltaTVsBeamE = new TH2I(locHistName.c_str(), ";Beam Energy;#Deltat_{Beam - RF}", 400, 0.0, 12.0, 3600, -18.0, 18.0);
		else
			dHist_BeamRFDeltaTVsBeamE = static_cast<TH2*>(locHist);

		locCurrentDir->cd();
	}
	japp->RootUnLock(); //unlock
}

DLorentzVector DSourceComboTimeHandler::Get_ChargedParticlePOCAToVertexX4(const DChargedTrackHypothesis* locHypothesis, bool locIsProductionVertex, const DSourceCombo* locVertexPrimaryFullCombo, const DKinematicData* locBeamParticle) const
{
	auto locPOCAPair = std::make_pair(locHypothesis, dSourceComboVertexer->Get_ConstrainingParticles(locIsProductionVertex, locVertexPrimaryFullCombo, locBeamParticle));
	auto locIterator = dChargedParticlePOCAToVertexX4.find(locPOCAPair);
	return ((locIterator != dChargedParticlePOCAToVertexX4.end()) ? locIterator->second : DLorentzVector());
}

void DSourceComboTimeHandler::Setup_NeutralShowers(const vector<const DNeutralShower*>& locNeutralShowers, const DEventRFBunch* locInitialEventRFBunch)
{
	//Precompute a few things for the neutral showers, before comboing
	//Even if it turns out some of this isn't technically needed,
	//it's still faster than doing a check to see if this has been done or not for every single photon-combo request

	//GET RF BUNCH
	dInitialEventRFBunch = locInitialEventRFBunch;

	//ARRANGE NEUTRAL SHOWERS
	//also, save to unknown-z, unknown-rf (all showers)
	vector<const DNeutralShower*> locBCALShowers, locFCALShowers;
	auto locUnknownZBin = DSourceComboInfo::Get_VertexZIndex_Unknown();
	for(const auto& locShower : locNeutralShowers)
	{
		auto& locContainer = (locShower->dDetectorSystem == SYS_BCAL) ? locBCALShowers : locFCALShowers;
		locContainer.push_back(locShower);
		dShowersByBeamBunchByZBin[locUnknownZBin][{}].push_back(locShower);
	}

	//CALCULATE KINEMATICS
	//FCAL: at target center
	auto locFCALZBin = DSourceComboInfo::Get_VertexZIndex_ZIndependent();
	for(const auto& locShower : locFCALShowers)
		dPhotonKinematics[locFCALZBin].emplace(locShower, Create_KinematicData_Photon(locShower, dTargetCenter));

	//BCAL: in vertex-z bins
	for(size_t loc_i = 0; loc_i < dNumPhotonVertexZBins; ++loc_i)
	{
		DVector3 locBinCenter(0.0, 0.0, dSourceComboer->Get_PhotonVertexZBinCenter(loc_i));
		for(const auto& locShower : locBCALShowers)
			dPhotonKinematics[loc_i].emplace(locShower, Create_KinematicData_Photon(locShower, locBinCenter));
	}

	//DETERMINE WHICH RF BUNCHES ARE VALID
	//FCAL: at target center
	for(const auto& locShower : locFCALShowers)
		Calc_PhotonBeamBunchShifts(locShower, dPhotonKinematics[locFCALZBin][locShower], dInitialEventRFBunch->dTime, locFCALZBin);

	//BCAL + FCAL: in vertex-z bins
	for(size_t loc_i = 0; loc_i < dNumPhotonVertexZBins; ++loc_i)
	{
		//propagate RF time to vertex position
		double locPropagatedRFTime = dInitialEventRFBunch->dTime + (dSourceComboer->Get_PhotonVertexZBinCenter(loc_i) - dTargetCenter.Z())/SPEED_OF_LIGHT;
		for(const auto& locShower : locBCALShowers)
			Calc_PhotonBeamBunchShifts(locShower, dPhotonKinematics[loc_i][locShower], locPropagatedRFTime, loc_i);

		//insert the previously-done FCAL photons
		dShowersByBeamBunchByZBin[loc_i].insert(dShowersByBeamBunchByZBin[locFCALZBin].begin(), dShowersByBeamBunchByZBin[locFCALZBin].end());
		for(const auto& locBeamBunchPair : dShowersByBeamBunchByZBin[locFCALZBin])
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

	//remove duplicates in the z-unknown vector
	for(auto& locRFShowerPair : dShowersByBeamBunchByZBin[locUnknownZBin])
	{
		if(locRFShowerPair.first.empty())
			continue; //no chance that this has duplicates
		auto& locShowerVector = locRFShowerPair.second;
		std::sort(locShowerVector.begin(), locShowerVector.end());
		locShowerVector.erase(std::unique(locShowerVector.begin(), locShowerVector.end()), locShowerVector.end());
	}
}

void DSourceComboTimeHandler::Calc_PhotonBeamBunchShifts(const DNeutralShower* locNeutralShower, shared_ptr<const DKinematicData>& locKinematicData, double locRFTime, signed char locZBin)
{
	//get delta-t cut
	DetectorSystem_t locSystem = locNeutralShower->dDetectorSystem;
	double locDeltaTCut = dPIDTimingCuts[Gamma][locSystem]->Eval(locNeutralShower->dEnergy) + Calc_MaxDeltaTError(locNeutralShower, locKinematicData);

	//do loop over possible #-RF-shifts
	double locVertexTime = locKinematicData->time();
	auto locRFShifts = Calc_BeamBunchShifts(locVertexTime, locRFTime, locDeltaTCut, true, Unknown, locSystem, locNeutralShower->dEnergy);

	auto locJObject = static_cast<const JObject*>(locNeutralShower);
	dShowerRFBunches[locZBin].emplace(locJObject, locRFShifts);
	for(const auto& locNumShifts : locRFShifts)
	{
		dShowersByBeamBunchByZBin[locZBin][{locNumShifts}].push_back(locJObject);
		dShowersByBeamBunchByZBin[DSourceComboInfo::Get_VertexZIndex_Unknown()][{locNumShifts}].push_back(locJObject);
	}
}

vector<int> DSourceComboTimeHandler::Calc_BeamBunchShifts(double locVertexTime, double locOrigRFBunchPropagatedTime, double locDeltaTCut, bool locIncludeDecayTimeOffset, Particle_t locPID, DetectorSystem_t locSystem, double locP)
{
	if(dDebugLevel >= 10)
		cout << "DSourceComboTimeHandler::Calc_BeamBunchShifts(): PID, system, period, orig time, vertex time, orig prop rf time, delta-t cut, include offset: " << locPID << ", " << SystemName(locSystem) << ", " << dBeamBunchPeriod << ", " << dInitialEventRFBunch->dTime << ", " << locVertexTime << ", " << locOrigRFBunchPropagatedTime << ", " << locDeltaTCut << ", " << locIncludeDecayTimeOffset << endl;
	auto locOrigNumShifts = Calc_RFBunchShift(locOrigRFBunchPropagatedTime, locVertexTime); //get best shift
	vector<int> locRFShifts;

	//start with best-shift, then loop up until fails cut
	auto locNumShifts = locOrigNumShifts;
	auto locDeltaT = locVertexTime - (locOrigRFBunchPropagatedTime + locNumShifts*dBeamBunchPeriod);
	//cout << "num shifts, delta-t = " << locNumShifts << ", " << locDeltaT << endl;
	dAllRFDeltaTs[locPID][locSystem].push_back(std::make_pair(locP, locDeltaT));
	while((fabs(locDeltaT) < locDeltaTCut) || (locOrigNumShifts == locNumShifts)) //extra condition for histogramming purposes only
	{
		if(fabs(locDeltaT) < locDeltaTCut)
			locRFShifts.push_back(locNumShifts);
		++locNumShifts;
		locDeltaT = locVertexTime - (locOrigRFBunchPropagatedTime + locNumShifts*dBeamBunchPeriod);
		dAllRFDeltaTs[locPID][locSystem].push_back(std::make_pair(locP, locDeltaT));
	}

	//now loop down in n-shifts
	locNumShifts = locOrigNumShifts - 1;
	locDeltaT = locVertexTime - (locOrigRFBunchPropagatedTime + locNumShifts*dBeamBunchPeriod);
	dAllRFDeltaTs[locPID][locSystem].push_back(std::make_pair(locP, locDeltaT));
	while(fabs(locDeltaT) < locDeltaTCut)
	{
		locRFShifts.push_back(locNumShifts);
		--locNumShifts;
		locDeltaT = locVertexTime - (locOrigRFBunchPropagatedTime + locNumShifts*dBeamBunchPeriod);
		dAllRFDeltaTs[locPID][locSystem].push_back(std::make_pair(locP, locDeltaT));
	}

	//due to detached vertices, we may need to accept EARLIER RF bunches
	if(locIncludeDecayTimeOffset)
	{
		//continue down-shift loop, this time including time offset
		//+dMaxTimeOffset: takes longer for "RF" time to get there (due to slow decaying particle)
		//Note that in the delta-t histograms, the dMaxDecayTimeOffset would cause a shift in the distribution for positive delta-t's, so we won't include these in the histograms
		locDeltaT = locVertexTime - (locOrigRFBunchPropagatedTime + locNumShifts*dBeamBunchPeriod + dMaxDecayTimeOffset);
		while(fabs(locDeltaT) < locDeltaTCut)
		{
			locRFShifts.push_back(locNumShifts);
			--locNumShifts;
			locDeltaT = locVertexTime - (locOrigRFBunchPropagatedTime + locNumShifts*dBeamBunchPeriod + dMaxDecayTimeOffset);
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
		auto locZError = dPhotonVertexZBinWidth/2.0; //evaluated at center of bin
		auto locR = locNeutralShower->dSpacetimeVertex.Vect().Perp();
		auto locPathError = locR*(1.0/sin(locTheta) - sqrt(1.0 + pow(1.0/tan(locTheta) - locZError/locR, 2.0))) - locZError;
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

bool DSourceComboTimeHandler::Select_RFBunches_Charged(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionChargedCombo, vector<int>& locValidRFBunches)
{
	if(dDebugLevel >= 10)
		cout << "DSourceComboTimeHandler::Select_RFBunches_Charged()" << endl;
	auto locRFIterator = dChargedComboRFBunches.find(locReactionChargedCombo);
	if(locRFIterator != dChargedComboRFBunches.end())
	{
		locValidRFBunches = locRFIterator->second; //already computed, return results!!
		if(dDebugLevel >= 10)
			cout << "Previously done, # valid bunches = " << locValidRFBunches.size() << endl;
		return (!locValidRFBunches.empty());
	}

	//All charged tracks that are at a known vertex position vote, even those not at the primary vertex
	//loop over vertices, get all charged particles at that vertex, utilize that + time offset
	//Applying the RF time cuts is effectively applying PID timing cuts

	//loop over vertices
	locValidRFBunches.clear();
	for(const auto& locStepVertexInfo : locReactionVertexInfo->Get_StepVertexInfos())
	{
		auto locIsProductionVertex = locStepVertexInfo->Get_ProductionVertexFlag();
		auto locVertexPrimaryCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionChargedCombo, locStepVertexInfo);
		if(!dSourceComboVertexer->Get_VertexDeterminableWithCharged(locStepVertexInfo))
			continue; //vertex position indeterminate at this stage: don't include these tracks

		//get combo, vertex position, and time offset from RF bunch
		auto locVertex = dSourceComboVertexer->Get_Vertex(locIsProductionVertex, locVertexPrimaryCombo, nullptr);
		auto locTimeOffset = dSourceComboVertexer->Get_TimeOffset(locIsProductionVertex, locReactionChargedCombo, locVertexPrimaryCombo, nullptr);
		double locPropagatedRFTime = dInitialEventRFBunch->dTime + (locVertex.Z() - dTargetCenter.Z())/SPEED_OF_LIGHT;
		//cout << "vertex z, targ z, init time, prop time = " << locVertex.Z() << ", " << dTargetCenter.Z() << ", " << dInitialEventRFBunch->dTime << ", " << locPropagatedRFTime << endl;

		//loop over charged particles
		auto locChargedParticles = DAnalysis::Get_SourceParticles_ThisVertex(locVertexPrimaryCombo);
		for(const auto& locParticlePair : locChargedParticles)
		{
			auto locChargedHypo = static_cast<const DChargedTrack*>(locParticlePair.second)->Get_Hypothesis(locParticlePair.first);
			vector<int> locParticleRFBunches;
			if(!Get_RFBunches_ChargedTrack(locChargedHypo, locIsProductionVertex, locVertexPrimaryCombo, locVertex, locTimeOffset, locPropagatedRFTime, locParticleRFBunches))
			{
				if(dDebugLevel >= 10)
					cout << "pid, has no timing info = " << locChargedHypo->PID() << endl;
				continue;
			}
			if(dDebugLevel >= 10)
				cout << "pid, # valid bunches = " << locChargedHypo->PID() << ", " << locParticleRFBunches.size() << endl;

			if(locParticleRFBunches.empty())
			{
				dChargedComboRFBunches.emplace(locReactionChargedCombo, vector<int>{});
				return false;
			}

			locValidRFBunches = Get_CommonRFBunches(locValidRFBunches, locParticleRFBunches);
			if(dDebugLevel >= 10)
				cout << "#common bunches = " << locValidRFBunches.size() << endl;
			if(locValidRFBunches.empty())
			{
				dChargedComboRFBunches.emplace(locReactionChargedCombo, vector<int>{});
				return false;
			}
		}
	}

	dChargedComboRFBunches.emplace(locReactionChargedCombo, locValidRFBunches);
	if(dDebugLevel >= 10)
		cout << "# valid bunches = " << locValidRFBunches.size() << endl;
	return true;
}

bool DSourceComboTimeHandler::Select_RFBunches_PhotonVertices(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, vector<int>& locValidRFBunches)
{
	//also, do PID cuts of photons at charged vertices while we're at it
	auto locRFIterator = dPhotonVertexRFBunches.find(locReactionFullCombo);
	if(locRFIterator != dPhotonVertexRFBunches.end())
	{
		locValidRFBunches = locRFIterator->second; //already computed, return results!!
		return (!locValidRFBunches.empty());
	}

	//loop over vertices
	for(const auto& locStepVertexInfo : locReactionVertexInfo->Get_StepVertexInfos())
	{
		if(locStepVertexInfo->Get_DanglingVertexFlag())
			continue;

		auto locIsProductionVertex = locStepVertexInfo->Get_ProductionVertexFlag();
		bool locVertexDeterminableWithCharged = dSourceComboVertexer->Get_VertexDeterminableWithCharged(locStepVertexInfo);

		auto locVertexPrimaryFullCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionFullCombo, locStepVertexInfo);
		if(dSourceComboVertexer->Get_ConstrainingParticles(locIsProductionVertex, locVertexPrimaryFullCombo).empty())
			continue; //vertex not found yet

		//get combo, vertex position, and time offset from RF bunch
		auto locVertex = dSourceComboVertexer->Get_Vertex(locIsProductionVertex, locVertexPrimaryFullCombo, nullptr);
		auto locVertexZBin = dSourceComboVertexer->Get_VertexZBin(locIsProductionVertex, locVertexPrimaryFullCombo, nullptr);

		auto locTimeOffset = dSourceComboVertexer->Get_TimeOffset(locIsProductionVertex, locReactionFullCombo, locVertexPrimaryFullCombo, nullptr);
		double locPropagatedRFTime = dInitialEventRFBunch->dTime + (locVertex.Z() - dTargetCenter.Z())/SPEED_OF_LIGHT;

		//loop over particles at this vertex: BCAL photons & charged tracks get to vote (FCAL photons already voted, but faster)
		auto locSourceParticles = DAnalysis::Get_SourceParticles_ThisVertex(locVertexPrimaryFullCombo);
		for(const auto& locParticlePair : locSourceParticles)
		{
			auto locPID = locParticlePair.first;
			vector<int> locParticleRFBunches;
			if(ParticleCharge(locPID) == 0)
			{
				if(ParticleMass(locPID) > 0.0)
					continue; //massive neutral, can't vote
				auto locNeutralShower = static_cast<const DNeutralShower*>(locParticlePair.second);
				locParticleRFBunches = (locNeutralShower->dDetectorSystem == SYS_FCAL) ? locValidRFBunches : dShowerRFBunches[locVertexZBin][locParticlePair.second]; //FCAL has already voted!

				//Do PID cut
				auto PhotonCutter = [&](int locRFBunch) -> bool {return Cut_PhotonPID(locNeutralShower, locVertex, locPropagatedRFTime + locRFBunch*dBeamBunchPeriod, false);};
				locParticleRFBunches.erase(std::remove_if(locParticleRFBunches.begin(), locParticleRFBunches.end(), PhotonCutter), locParticleRFBunches.end());
			}
			else if(locVertexDeterminableWithCharged) //charged, previously cut
				continue;
			else //charged, a new vertex: do PID cuts
			{
				auto locChargedHypo = static_cast<const DChargedTrack*>(locParticlePair.second)->Get_Hypothesis(locParticlePair.first);
				if(!Get_RFBunches_ChargedTrack(locChargedHypo, locIsProductionVertex, locVertexPrimaryFullCombo, locVertex, locTimeOffset, locPropagatedRFTime, locParticleRFBunches))
				{
					if(dDebugLevel >= 10)
						cout << "pid, has no timing info = " << locChargedHypo->PID() << endl;
					continue;
				}
			}

			if(locParticleRFBunches.empty())
			{
				dPhotonVertexRFBunches.emplace(locReactionFullCombo, vector<int>{});
				return false;
			}
			if(locValidRFBunches.empty())
			{
				locValidRFBunches = locParticleRFBunches;
				continue;
			}

			//get common rf bunches
			locValidRFBunches = Get_CommonRFBunches(locValidRFBunches, locParticleRFBunches);
			if(locValidRFBunches.empty())
			{
				dPhotonVertexRFBunches.emplace(locReactionFullCombo, vector<int>{});
				return false;
			}
		}
	}

	dPhotonVertexRFBunches.emplace(locReactionFullCombo, locValidRFBunches);
	return true;
}

bool DSourceComboTimeHandler::Select_RFBunches_AllVerticesUnknown(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, Charge_t locCharge, vector<int>& locValidRFBunches)
{
	//this function is called to select RF bunches if every vertex is unknown
	auto locRFIterator = dPhotonVertexRFBunches.find(locReactionFullCombo);
	if(locRFIterator != dPhotonVertexRFBunches.end())
	{
		locValidRFBunches = locRFIterator->second; //already computed, return results!!
		return (!locValidRFBunches.empty());
	}

	locValidRFBunches.clear();

	//get and loop over all particles
	auto locSourceParticles = locReactionFullCombo->Get_SourceParticles(true, locCharge);
	auto locVertexZBin = dSourceComboer->Get_VertexZBin_TargetCenter();
	for(const auto& locParticlePair : locSourceParticles)
	{
		auto locPID = locParticlePair.first;
		vector<int> locParticleRFBunches;
		if(ParticleCharge(locPID) == 0)
		{
			if(ParticleMass(locPID) > 0.0)
				continue; //massive neutral, can't vote
			auto locNeutralShower = static_cast<const DNeutralShower*>(locParticlePair.second);
			locParticleRFBunches = (locNeutralShower->dDetectorSystem == SYS_FCAL) ? locValidRFBunches : dShowerRFBunches[locVertexZBin][locParticlePair.second]; //FCAL has already voted!

			//Do PID cut
			auto PhotonCutter = [&](int locRFBunch) -> bool {return Cut_PhotonPID(locNeutralShower, dTargetCenter, dInitialEventRFBunch->dTime + locRFBunch*dBeamBunchPeriod, true);};
			locParticleRFBunches.erase(std::remove_if(locParticleRFBunches.begin(), locParticleRFBunches.end(), PhotonCutter), locParticleRFBunches.end());
		}
		else //charged, a new vertex: do PID cuts
		{
			auto locChargedHypo = static_cast<const DChargedTrack*>(locParticlePair.second)->Get_Hypothesis(locParticlePair.first);
			auto locVertex = locChargedHypo->position();
			auto locPropagatedRFTime = dInitialEventRFBunch->dTime + (locVertex.Z() - dTargetCenter.Z())/SPEED_OF_LIGHT;
			if(!Get_RFBunches_ChargedTrack(locChargedHypo, true, nullptr, locVertex, 0.0, locPropagatedRFTime, locParticleRFBunches))
			{
				if(dDebugLevel >= 10)
					cout << "pid, has no timing info = " << locChargedHypo->PID() << endl;
				continue;
			}
		}

		if(locParticleRFBunches.empty())
		{
			dPhotonVertexRFBunches.emplace(locReactionFullCombo, vector<int>{});
			return false;
		}

		//get common rf bunches
		locValidRFBunches = Get_CommonRFBunches(locValidRFBunches, locParticleRFBunches);
		if(locValidRFBunches.empty())
		{
			dPhotonVertexRFBunches.emplace(locReactionFullCombo, vector<int>{});
			return false;
		}
	}

	dPhotonVertexRFBunches.emplace(locReactionFullCombo, locValidRFBunches);
	return true;
}

int DSourceComboTimeHandler::Select_RFBunch_Full(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, const vector<int>& locRFBunches)
{
	//The prior functions narrowed down the possible RF bunches.  This function actually selects the RF bunch.

	if(dDebugLevel >= 10)
		cout << "DSourceComboTimeHandler::Select_RFBunch_Full()" << endl;

	auto locRFIterator = dFullComboRFBunches.find(locReactionFullCombo);
	if(locRFIterator != dFullComboRFBunches.end())
	{
		if(dDebugLevel >= 10)
			cout << "bunch previously determined: " << locRFIterator->second << endl;
		return locRFIterator->second; //already computed, return results!!
	}

	//note: even if there's only 1 bunch, we still need to cut on photon timing (if any) (and histogram): proceed with the function
	if(locRFBunches.empty())
		return 0; //no information, hope that the default was correct //only detected particles are massive neutrals: e.g. g, p -> K_L, (p)

	//initialize chisq's
	unordered_map<int, double> locChiSqByRFBunch;
	for(const auto& locRFBunch : locRFBunches)
		locChiSqByRFBunch.emplace(locRFBunch, 0.0);

	//voting:
	//first use only particles at vertices that have been determined
	//then use charged particles at their POCAs to the beamline
	//then use photons at the target center
	bool locVotedFlag = false;

	//loop over vertices, using only particles at the vertices that have been determined
	auto locPrimaryVertexZ = dSourceComboVertexer->Get_PrimaryVertex(locReactionVertexInfo, locReactionFullCombo, nullptr).Z();
	if(dDebugLevel >= 10)
		cout << "primary vertex z: " << locPrimaryVertexZ << endl;
	map<int, map<Particle_t, map<DetectorSystem_t, vector<pair<float, float>>>>> locRFDeltaTsForHisting; //first float is p, 2nd is delta-t //PID Unknown: photons prior to vertex selection
	for(const auto& locStepVertexInfo : locReactionVertexInfo->Get_StepVertexInfos())
	{
		if(dDebugLevel >= 10)
			cout << "Step: " << locStepVertexInfo->Get_StepIndices().front() << endl;

		if(locStepVertexInfo->Get_DanglingVertexFlag())
			continue; //unknown position!

		auto locIsProductionVertex = locStepVertexInfo->Get_ProductionVertexFlag();
		auto locVertexPrimaryFullCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionFullCombo, locStepVertexInfo);
		if(!dSourceComboVertexer->Get_VertexDeterminableWithCharged(locStepVertexInfo) && !dSourceComboVertexer->Get_VertexDeterminableWithPhotons(locStepVertexInfo))
			continue; //vertex position indeterminate at this stage: don't include these particles

		//get combo, vertex position, and time offset from RF bunch
		auto locVertex = dSourceComboVertexer->Get_Vertex(locIsProductionVertex, locVertexPrimaryFullCombo, nullptr);
		auto locTimeOffset = dSourceComboVertexer->Get_TimeOffset(locIsProductionVertex, locReactionFullCombo, locVertexPrimaryFullCombo, nullptr);
		auto locParticles = DAnalysis::Get_SourceParticles_ThisVertex(locVertexPrimaryFullCombo);
		if(dDebugLevel >= 10)
			cout << "this vertex z, time offset: " << locVertex.Z() << ", " << locTimeOffset << endl;

		for(const auto& locRFBunch : locRFBunches)
		{
			//propagate rf time to vertex and add time offset (faster to just do it here rather than for each particle)
			auto locPropagatedRFTime = Calc_PropagatedRFTime(locPrimaryVertexZ, locRFBunch, locTimeOffset);
			if(dDebugLevel >= 10)
				cout << "bunch, propagated time: " << locRFBunch << ", " << locPropagatedRFTime << endl;

			//loop over all particles
			for(const auto& locParticlePair : locParticles)
			{
				auto locPID = locParticlePair.first;
				if((ParticleCharge(locPID) == 0) && (ParticleMass(locPID) > 0.0))
					continue; //ignore massive neutrals: timing defines their momentum, cannot be used

				locVotedFlag = true;
				if(ParticleCharge(locPID) == 0)
				{
					auto locNeutralShower = static_cast<const DNeutralShower*>(locParticlePair.second);
					auto locRFDeltaTPair = Calc_RFDeltaTChiSq(locNeutralShower, locVertex, locPropagatedRFTime);
					locChiSqByRFBunch[locRFBunch] += locRFDeltaTPair.second;
					locRFDeltaTsForHisting[locRFBunch][locPID][locNeutralShower->dDetectorSystem].emplace_back(locNeutralShower->dEnergy, locRFDeltaTPair.first);
				}
				else //charged
				{
					auto locChargedHypo = static_cast<const DChargedTrack*>(locParticlePair.second)->Get_Hypothesis(locParticlePair.first);

					//get the timing at the POCA to the vertex (computed previously!)
					auto locPOCAPair = std::make_pair(locChargedHypo, dSourceComboVertexer->Get_ConstrainingParticles(locIsProductionVertex, locVertexPrimaryFullCombo, nullptr));
					auto locVertexTime = dChargedParticlePOCAToVertexX4.find(locPOCAPair)->second.T();
					auto locRFDeltaTPair = Calc_RFDeltaTChiSq(locChargedHypo, locVertexTime, locPropagatedRFTime);
					locChiSqByRFBunch[locRFBunch] += locRFDeltaTPair.second;
					locRFDeltaTsForHisting[locRFBunch][locPID][locChargedHypo->t1_detector()].emplace_back(locChargedHypo->momentum().Mag(), locRFDeltaTPair.first);
				}
			}
			if(dDebugLevel >= 10)
				cout << "bunch, total delta-t chisq = " << locRFBunch << ", " << locChiSqByRFBunch[locRFBunch] << endl;
		}
	}

	//if not voted yet, use charged particles at their POCAs to the beamline, else use photons at the target center
	if(!locVotedFlag)
	{
		if(!Compute_RFChiSqs_UnknownVertices(locReactionFullCombo, d_Charged, locRFBunches, locChiSqByRFBunch, locRFDeltaTsForHisting))
			Compute_RFChiSqs_UnknownVertices(locReactionFullCombo, d_Neutral, locRFBunches, locChiSqByRFBunch, locRFDeltaTsForHisting);
	}

	//ok, total chisq's are computed, pick the one that is the best!
	auto Compare_RFChiSqs = [](const pair<int, double>& lhs, const pair<int, double>& rhs) -> bool {return lhs.second < rhs.second;};
	auto locRFBunch = std::max_element(locChiSqByRFBunch.begin(), locChiSqByRFBunch.end(), Compare_RFChiSqs)->first;

	if(dDebugLevel >= 10)
		cout << "chosen bunch: " << locRFBunch << endl;

	//propagate chosen bunch timing into hist staging vector
	for(auto& locPIDPair : locRFDeltaTsForHisting[locRFBunch])
	{
		for(auto& locSystemPair : locPIDPair.second)
		{
			auto& locSaveToVector = dSelectedRFDeltaTs[locPIDPair.first][locSystemPair.first];
			auto& locMoveFromVector = locSystemPair.second;
			std::move(locMoveFromVector.begin(), locMoveFromVector.end(), std::back_inserter(locSaveToVector));
		}
	}

	//save it and return
	dFullComboRFBunches.emplace(locReactionFullCombo, locRFBunch);
	return locRFBunch;
}

bool DSourceComboTimeHandler::Compute_RFChiSqs_UnknownVertices(const DSourceCombo* locReactionFullCombo, Charge_t locCharge, const vector<int>& locRFBunches, unordered_map<int, double>& locChiSqByRFBunch, map<int, map<Particle_t, map<DetectorSystem_t, vector<pair<float, float>>>>>& locRFDeltaTsForHisting)
{
	//get and loop over all particles
	auto locSourceParticles = locReactionFullCombo->Get_SourceParticles(true, locCharge);
	bool locVotedFlag = false;
	for(const auto& locRFBunch : locRFBunches)
	{
		//loop over all particles
		for(const auto& locParticlePair : locSourceParticles)
		{
			auto locPID = locParticlePair.first;
			if((ParticleCharge(locPID) == 0) && (ParticleMass(locPID) > 0.0))
				continue; //ignore massive neutrals: timing defines their momentum, cannot be used

			locVotedFlag = true;
			if(ParticleCharge(locPID) == 0)
			{
				auto locNeutralShower = static_cast<const DNeutralShower*>(locParticlePair.second);
				auto locPropagatedRFTime = Calc_PropagatedRFTime(dTargetCenter.Z(), locRFBunch, 0.0);
				auto locRFDeltaTPair = Calc_RFDeltaTChiSq(locNeutralShower, dTargetCenter, locPropagatedRFTime);
				locChiSqByRFBunch[locRFBunch] += locRFDeltaTPair.second;
				locRFDeltaTsForHisting[locRFBunch][locPID][locNeutralShower->dDetectorSystem].emplace_back(locNeutralShower->dEnergy, locRFDeltaTPair.first);
			}
			else //charged
			{
				auto locChargedHypo = static_cast<const DChargedTrack*>(locParticlePair.second)->Get_Hypothesis(locParticlePair.first);
				auto locPropagatedRFTime = Calc_PropagatedRFTime(locChargedHypo->position().Z(), locRFBunch, 0.0);
				auto locRFDeltaTPair = Calc_RFDeltaTChiSq(locChargedHypo, locChargedHypo->time(), locPropagatedRFTime);
				locChiSqByRFBunch[locRFBunch] += locRFDeltaTPair.second;
				locRFDeltaTsForHisting[locRFBunch][locPID][locChargedHypo->t1_detector()].emplace_back(locChargedHypo->momentum().Mag(), locRFDeltaTPair.first);
			}
		}
		if(dDebugLevel >= 10)
			cout << "bunch, total delta-t chisq = " << locRFBunch << ", " << locChiSqByRFBunch[locRFBunch] << endl;
	}
	return locVotedFlag;
}

bool DSourceComboTimeHandler::Cut_Timing_MissingMassVertices(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, const DKinematicData* locBeamParticle, int locRFBunch)
{
	if(dDebugLevel >= 10)
		cout << "DSourceComboTimeHandler::Cut_Timing_MissingMassVertices()" << endl;

	auto locTimeCutIterator = dFullComboTimeCutResults.find(locReactionFullCombo);
	if(locTimeCutIterator != dFullComboTimeCutResults.end())
	{
		if(dDebugLevel >= 10)
			cout << "previously computed: " << locTimeCutIterator->second << endl;
		return locTimeCutIterator->second; //already computed, return results!!
	}

	//All charged tracks vote, even those not at the primary vertex
	//loop over vertices, get all charged particles at that vertex, utilize that + time offset
	auto locPrimaryVertexZ = dSourceComboVertexer->Get_PrimaryVertex(locReactionVertexInfo, locReactionFullCombo, locBeamParticle).Z();

	//loop over vertices
	for(const auto& locStepVertexInfo : locReactionVertexInfo->Get_StepVertexInfos())
	{
		if(dDebugLevel >= 10)
			cout << "Step: " << locStepVertexInfo->Get_StepIndices().front() << ", dangling-flag = " << locStepVertexInfo->Get_DanglingVertexFlag() << endl;
		if(locStepVertexInfo->Get_DanglingVertexFlag())
			continue;

		if(dDebugLevel >= 10)
			cout << "determinable flags: " << dSourceComboVertexer->Get_VertexDeterminableWithCharged(locStepVertexInfo) << ", " << dSourceComboVertexer->Get_VertexDeterminableWithPhotons(locStepVertexInfo) << endl;
		if(dSourceComboVertexer->Get_VertexDeterminableWithCharged(locStepVertexInfo))
			continue; //these have already been cut!
		if(dSourceComboVertexer->Get_VertexDeterminableWithPhotons(locStepVertexInfo))
			continue; //these have already been cut!

		//get combo, vertex position, and time offset from RF bunch
		auto locIsProductionVertex = locStepVertexInfo->Get_ProductionVertexFlag();
		auto locVertexPrimaryFullCombo = dSourceComboer->Get_VertexPrimaryCombo(locReactionFullCombo, locStepVertexInfo);
		auto locVertex = dSourceComboVertexer->Get_Vertex(locIsProductionVertex, locVertexPrimaryFullCombo, locBeamParticle);

		auto locTimeOffset = dSourceComboVertexer->Get_TimeOffset(locIsProductionVertex, locReactionFullCombo, locVertexPrimaryFullCombo, locBeamParticle);
		auto locPropagatedRFTime = Calc_PropagatedRFTime(locPrimaryVertexZ, locRFBunch, locTimeOffset);

		//loop over particles at this vertex: BCAL photons & charged tracks will get cut (FCAL photons already voted!)
		auto locSourceParticles = DAnalysis::Get_SourceParticles_ThisVertex(locVertexPrimaryFullCombo);
		for(const auto& locParticlePair : locSourceParticles)
		{
			auto locPID = locParticlePair.first;
			if(ParticleCharge(locPID) == 0)
			{
				if(ParticleMass(locPID) > 0.0)
					continue; //massive neutral, can't vote
				auto locNeutralShower = static_cast<const DNeutralShower*>(locParticlePair.second);

				//Do PID cut
				if(!Cut_PhotonPID(locNeutralShower, locVertex, locPropagatedRFTime, false))
				{
					dFullComboTimeCutResults.emplace(locReactionFullCombo, false);
					return false;
				}
			}
			else //charged
			{
				auto locChargedHypo = static_cast<const DChargedTrack*>(locParticlePair.second)->Get_Hypothesis(locParticlePair.first);
				if(!Cut_TrackPID(locChargedHypo, locIsProductionVertex, locVertexPrimaryFullCombo, locVertex, locPropagatedRFTime))
				{
					dFullComboTimeCutResults.emplace(locReactionFullCombo, false);
					return false;
				}
			}
		}
	}

	dFullComboTimeCutResults.emplace(locReactionFullCombo, true);
	return true;
}

void DSourceComboTimeHandler::Fill_Histograms(void)
{
	japp->WriteLock("DSourceComboTimeHandler");
	{
		for(auto& locDeltaTPair : dBeamRFDeltaTs)
			dHist_BeamRFDeltaTVsBeamE->Fill(locDeltaTPair.first, locDeltaTPair.second);

		for(auto& locPIDPair : dSelectedRFDeltaTs)
		{
			auto locPIDIterator = dHistMap_RFDeltaTVsP_BestRF.find(locPIDPair.first);
			if(locPIDIterator == dHistMap_RFDeltaTVsP_BestRF.end())
				continue;
			for(auto& locSystemPair : locPIDPair.second)
			{
				auto locSystemIterator = locPIDIterator->second.find(locSystemPair.first);
				if(locSystemIterator == locPIDIterator->second.end())
					continue;

				auto& locBestHist = locSystemIterator->second;
				for(auto& locVectorPair : locSystemPair.second) //best vector
					locBestHist->Fill(locVectorPair.first, locVectorPair.second);
			}
		}

		for(auto& locPIDPair : dAllRFDeltaTs)
		{
			auto locPIDIterator = dHistMap_RFDeltaTVsP_AllRFs.find(locPIDPair.first);
			if(locPIDIterator == dHistMap_RFDeltaTVsP_AllRFs.end())
				continue;
			for(auto& locSystemPair : locPIDPair.second)
			{
				auto locSystemIterator = locPIDIterator->second.find(locSystemPair.first);
				if(locSystemIterator == locPIDIterator->second.end())
					continue;

				auto& locAllHist = locSystemIterator->second;
				for(auto& locVectorPair : locSystemPair.second) //best vector
					locAllHist->Fill(locVectorPair.first, locVectorPair.second);
			}
		}
	}
	japp->Unlock("DSourceComboTimeHandler");

	//Reset for next event
	dBeamRFDeltaTs.clear();
	for(auto& locPIDPair : dSelectedRFDeltaTs)
	{
		for(auto& locSystemPair : locPIDPair.second)
			locSystemPair.second.clear();
	}
	for(auto& locPIDPair : dAllRFDeltaTs)
	{
		for(auto& locSystemPair : locPIDPair.second)
			locSystemPair.second.clear();
	}
}

bool DSourceComboTimeHandler::Get_RFBunches_ChargedTrack(const DChargedTrackHypothesis* locHypothesis, bool locIsProductionVertex, const DSourceCombo* locVertexPrimaryCombo, DVector3 locVertex, double locTimeOffset, double locPropagatedRFTime, vector<int>& locRFBunches)
{
	locRFBunches.clear();

	auto locPID = locHypothesis->PID();
	auto locSystem = locHypothesis->t1_detector();

	auto locX4 = Get_ChargedPOCAToVertexX4(locHypothesis, locIsProductionVertex, locVertexPrimaryCombo, locVertex);
	auto locVertexTime = locX4.T() - locTimeOffset;

	auto locP = locHypothesis->momentum().Mag();
	auto locCutFunc = Get_TimeCutFunction(locPID, locSystem);
	auto locDeltaTCut = (locCutFunc != nullptr) ? locCutFunc->Eval(locP) : 3.0; //if null, still use for histogramming

//TEMP!
locVertexTime = locHypothesis->time();
locPropagatedRFTime = Calc_PropagatedRFTime(locHypothesis->position().Z(), 0, 0.0);

	locRFBunches = Calc_BeamBunchShifts(locVertexTime, locPropagatedRFTime, locDeltaTCut, false, locPID, locSystem, locP);
	return (locCutFunc != nullptr);
}

//evaluate timing at the POCA to the vertex
DLorentzVector DSourceComboTimeHandler::Get_ChargedPOCAToVertexX4(const DChargedTrackHypothesis* locHypothesis, bool locIsProductionVertex, const DSourceCombo* locVertexPrimaryCombo, DVector3 locVertex)
{
	auto locX4(locHypothesis->x4());
	if(locVertexPrimaryCombo == nullptr)
		return locX4;

	auto locP4(locHypothesis->lorentzMomentum());
	auto locPOCAPair = std::make_pair(locHypothesis, dSourceComboVertexer->Get_ConstrainingParticles(locIsProductionVertex, locVertexPrimaryCombo));
	auto locPOCAIterator = dChargedParticlePOCAToVertexX4.find(locPOCAPair);
	if(locPOCAIterator != dChargedParticlePOCAToVertexX4.end())
		return locPOCAIterator->second;

	//do the propagation
	dAnalysisUtilities->Propagate_Track(locHypothesis->charge(), locVertex, locX4, locP4, nullptr);
	dChargedParticlePOCAToVertexX4.emplace(locPOCAPair, locX4); //save results so we don't have to do it again
	return locX4;
}

}

