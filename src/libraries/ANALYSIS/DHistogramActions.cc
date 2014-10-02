#include "ANALYSIS/DHistogramActions.h"

void DHistogramAction_PID::Initialize(JEventLoop* locEventLoop)
{
	string locHistName, locHistTitle;
	string locParticleName, locParticleName2, locParticleROOTName, locParticleROOTName2;
	Particle_t locPID, locPID2;

	vector<const DParticleID*> locParticleIDs;
	locEventLoop->Get(locParticleIDs);

	deque<Particle_t> locDesiredPIDs;
	Get_Reaction()->Get_DetectedFinalPIDs(locDesiredPIDs);

	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	vector<const DAnalysisUtilities*> locAnalysisUtilitiesVector;
	locEventLoop->Get(locAnalysisUtilitiesVector);

	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		CreateAndChangeTo_ActionDirectory();
		dParticleID = locParticleIDs[0];
		dAnalysisUtilities = locAnalysisUtilitiesVector[0];
		for(size_t loc_i = 0; loc_i < locDesiredPIDs.size(); ++loc_i)
		{
			locPID = locDesiredPIDs[loc_i];
			locParticleName = ParticleType(locPID);
			locParticleROOTName = ParticleName_ROOT(locPID);
			CreateAndChangeTo_Directory(locParticleName, locParticleName);

			// Confidence Level
			locHistName = "PIDConfidenceLevel";
			locHistTitle = locParticleROOTName + string(" PID;PID Confidence Level");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_PIDFOM[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PIDFOM[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumFOMBins, 0.0, 1.0);

			// Confidence Level in BCAL
			locHistName = "TOFConfidenceLevel_BCAL";
			locHistTitle = locParticleROOTName + string(" PID in BCAL;TOF Confidence Level");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TOFFOM_BCAL[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TOFFOM_BCAL[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumFOMBins, 0.0, 1.0);

			// Confidence Level in FCAL
			locHistName = "TOFConfidenceLevel_FCAL";
			locHistTitle = locParticleROOTName + string(" PID in FCAL;TOF Confidence Level");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TOFFOM_FCAL[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TOFFOM_FCAL[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumFOMBins, 0.0, 1.0);

			if(ParticleCharge(locPID) != 0)
			{
				// Confidence Level in TOF
				locHistName = "TOFConfidenceLevel_TOF";
				locHistTitle = locParticleROOTName + string(" PID in TOF;TOF Confidence Level");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TOFFOM_TOF[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TOFFOM_TOF[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumFOMBins, 0.0, 1.0);

				// Confidence Level in CDC
				locHistName = "TOFConfidenceLevel_CDC";
				locHistTitle = locParticleROOTName + string(" PID in CDC;TOF Confidence Level");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TOFFOM_CDC[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TOFFOM_CDC[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumFOMBins, 0.0, 1.0);

				// DC dE/dx Confidence Level
				locHistName = "DCdEdxConfidenceLevel";
				locHistTitle = locParticleROOTName + string(" PID;DC #it{#frac{dE}{dx}} Confidence Level");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_DCdEdxFOM[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_DCdEdxFOM[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumFOMBins, 0.0, 1.0);
			}

			//beta vs p
			locHistName = "BetaVsP";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#beta");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_BetaVsP[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_BetaVsP[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumBetaBins, dMinBeta, dMaxBeta);

			//delta-beta vs p
			locHistName = "DeltaBetaVsP";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta#beta");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaBetaVsP[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaBetaVsP[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaBetaBins, dMinDeltaBeta, dMaxDeltaBeta);

			//TOF Confidence Level vs Delta-Beta
			locHistName = "TOFConfidenceLevelVsDeltaBeta";
			locHistTitle = locParticleROOTName + string(";#Delta#beta;TOF Confidence Level");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TOFFOMVsDeltaBeta[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TOFFOMVsDeltaBeta[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaBetaBins, dMinDeltaBeta, dMaxDeltaBeta, dNumFOMBins, 0.0, 1.0);

			//one per thrown pid:
			if(!locMCThrowns.empty())
			{
				CreateAndChangeTo_Directory("Throwns", "Throwns");
				for(size_t loc_j = 0; loc_j < dThrownPIDs.size(); ++loc_j)
				{
					locPID2 = dThrownPIDs[loc_j];
					if((ParticleCharge(locPID2) != ParticleCharge(locPID)) && (locPID2 != Unknown))
						continue;
					locParticleName2 = ParticleType(locPID2);
					locParticleROOTName2 = ParticleName_ROOT(locPID2);

					//Confidence Level for Thrown PID
					pair<Particle_t, Particle_t> locPIDPair(locPID, locPID2);
					locHistName = string("PIDConfidenceLevel_ForThrownPID_") + locParticleName2;
					locHistTitle = locParticleROOTName + string(" PID, Thrown PID = ") + locParticleROOTName2 + string(";PID Confidence Level");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistMap_PIDFOMForTruePID[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistMap_PIDFOMForTruePID[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumFOMBins, 0.0, 1.0);
				}
				gDirectory->cd("..");
			}

			CreateAndChangeTo_Directory("Pulls", "Pulls");

			//Delta-t Pulls - CDC
			if(ParticleCharge(locPID) != 0) //CDC
			{
				locHistName = "TimePull_CDC";
				locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TimePull_CDC[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePull_CDC[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

				locHistName = "TimePullVsTheta_CDC";
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TimePullVsTheta_CDC[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePullVsTheta_CDC[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);

				locHistName = "TimePullVsP_CDC";
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TimePullVsP_CDC[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePullVsP_CDC[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);
			}

			//Delta-t Pulls - BCAL
			locHistName = "TimePull_BCAL";
			locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TimePull_BCAL[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TimePull_BCAL[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

			locHistName = "TimePullVsTheta_BCAL";
			locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat/#sigma_{#Deltat}");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TimePullVsTheta_BCAL[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TimePullVsTheta_BCAL[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);

			locHistName = "TimePullVsP_BCAL";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TimePullVsP_BCAL[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TimePullVsP_BCAL[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

			//Delta-t Pulls - TOF
			if(ParticleCharge(locPID) != 0) //TOF
			{
				locHistName = "TimePull_TOF";
				locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TimePull_TOF[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePull_TOF[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

				locHistName = "TimePullVsP_TOF";
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TimePullVsP_TOF[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePullVsP_TOF[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);
			}

			//Delta-t Pulls - FCAL
			locHistName = "TimePull_FCAL";
			locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TimePull_FCAL[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TimePull_FCAL[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

			locHistName = "TimePullVsP_FCAL";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TimePullVsP_FCAL[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TimePullVsP_FCAL[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

			gDirectory->cd("..");
			CreateAndChangeTo_Directory("PVsTheta", "PVsTheta");

			// P Vs Theta, PID FOM < 1%
			locHistName = "PVsTheta_LowPIDFOM";
			locHistTitle = locParticleROOTName + string(", PID FOM < 1%;#theta#circ;p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_PVsTheta_LowPIDFOM[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PVsTheta_LowPIDFOM[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			// P Vs Theta, PID FOM = NaN
			locHistName = "PVsTheta_NaNPIDFOM";
			locHistTitle = locParticleROOTName + string(", PID FOM = NaN;#theta#circ;p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_PVsTheta_NaNPIDFOM[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PVsTheta_NaNPIDFOM[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			if(ParticleCharge(locPID) != 0) //no other sources of PID for neutrals
			{
				// P Vs Theta, TOF FOM < 1%
				locHistName = "PVsTheta_LowTOFFOM";
				locHistTitle = locParticleROOTName + string(", TOF FOM < 1%;#theta#circ;p (GeV/c)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_PVsTheta_LowTOFFOM[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_PVsTheta_LowTOFFOM[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

				// P Vs Theta, TOF FOM = NaN
				locHistName = "PVsTheta_NaNTOFFOM";
				locHistTitle = locParticleROOTName + string(", TOF FOM = NaN;#theta#circ;p (GeV/c)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_PVsTheta_NaNTOFFOM[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_PVsTheta_NaNTOFFOM[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

				// P Vs Theta, DC dE/dx FOM < 1%
				locHistName = "PVsTheta_LowDCdEdxFOM";
				locHistTitle = locParticleROOTName + string(", DC #it{#frac{dE}{dx}} FOM < 1%;#theta#circ;p (GeV/c)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_PVsTheta_LowDCdEdxFOM[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_PVsTheta_LowDCdEdxFOM[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

				// P Vs Theta, DC dE/dx FOM = NaN
				locHistName = "PVsTheta_NaNDCdEdxFOM";
				locHistTitle = locParticleROOTName + string(", DC #it{#frac{dE}{dx}} FOM = NaN;#theta#circ;p (GeV/c)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_PVsTheta_NaNDCdEdxFOM[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_PVsTheta_NaNDCdEdxFOM[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

				// P Vs Theta, Beta < 0
				locHistName = "PVsTheta_NegativeBeta";
				locHistTitle = locParticleROOTName + string(", #beta < 0.0;#theta#circ;p (GeV/c)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_PVsTheta_NegativeBeta[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_PVsTheta_NegativeBeta[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);
			}
			gDirectory->cd("..");

			gDirectory->cd("..");
		} //end of PID loop
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_PID::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	if(Get_NumPreviousParticleCombos() == 0)
		dPreviouslyHistogrammedParticles.clear();

	const DMCThrownMatching* locMCThrownMatching = NULL;
	locEventLoop->GetSingle(locMCThrownMatching);
	const DEventRFBunch* locEventRFBunch = locParticleCombo->Get_EventRFBunch();

	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);

		deque<const DKinematicData*> locParticles;
		locParticleComboStep->Get_FinalParticles_Measured(locParticles);
		for(size_t loc_j = 0; loc_j < locParticles.size(); ++loc_j)
		{
			if(!locParticleComboStep->Is_FinalParticleDetected(loc_j))
				continue;

			//check if will be duplicate
			pair<Particle_t, const JObject*> locParticleInfo(locParticles[loc_j]->PID(), locParticleComboStep->Get_FinalParticle_SourceObject(loc_j));
			pair<const DEventRFBunch*, pair<Particle_t, const JObject*> > locHistInfo(locEventRFBunch, locParticleInfo);
			if(dPreviouslyHistogrammedParticles.find(locHistInfo) != dPreviouslyHistogrammedParticles.end())
				continue; //previously histogrammed

			dPreviouslyHistogrammedParticles.insert(locHistInfo);
			if(locParticleComboStep->Is_FinalParticleCharged(loc_j)) //charged
				Fill_ChargedHists(static_cast<const DChargedTrackHypothesis*>(locParticles[loc_j]), locMCThrownMatching, locEventRFBunch);
			else //neutral
				Fill_NeutralHists(static_cast<const DNeutralParticleHypothesis*>(locParticles[loc_j]), locMCThrownMatching, locEventRFBunch);
		}
	}

	return true;
}

void DHistogramAction_PID::Fill_ChargedHists(const DChargedTrackHypothesis* locChargedTrackHypothesis, const DMCThrownMatching* locMCThrownMatching, const DEventRFBunch* locEventRFBunch)
{
	Particle_t locPID = locChargedTrackHypothesis->PID();
	double locBeta_Timing = dAnalysisUtilities->Calc_Beta_Timing(locChargedTrackHypothesis, locEventRFBunch, true);
	double locDeltaBeta = locChargedTrackHypothesis->lorentzMomentum().Beta() - locBeta_Timing;

	double locFOM_Timing = (locChargedTrackHypothesis->dNDF_Timing > 0) ? TMath::Prob(locChargedTrackHypothesis->dChiSq_Timing, locChargedTrackHypothesis->dNDF_Timing) : numeric_limits<double>::quiet_NaN();
	double locFOM_DCdEdx = (locChargedTrackHypothesis->dNDF_DCdEdx > 0) ? TMath::Prob(locChargedTrackHypothesis->dChiSq_DCdEdx, locChargedTrackHypothesis->dNDF_DCdEdx) : numeric_limits<double>::quiet_NaN();

	double locP = locChargedTrackHypothesis->momentum().Mag();
	double locTheta = locChargedTrackHypothesis->momentum().Theta()*180.0/TMath::Pi();
	double locMatchFOM = 0.0;
	const DMCThrown* locMCThrown = (locMCThrownMatching != NULL) ? locMCThrownMatching->Get_MatchingMCThrown(locChargedTrackHypothesis, locMatchFOM) : NULL;

	double locTimePull = 0.0;
	unsigned int locTimeNDF = 0;
	dParticleID->Calc_TimingChiSq(locChargedTrackHypothesis, locEventRFBunch, true, locTimeNDF, locTimePull);

	japp->RootWriteLock();
	{
		dHistMap_PIDFOM[locPID]->Fill(locChargedTrackHypothesis->dFOM);

		if(locChargedTrackHypothesis->t1_detector() == SYS_CDC)
		{
			dHistMap_TOFFOM_CDC[locPID]->Fill(locFOM_Timing);
			dHistMap_TimePull_CDC[locPID]->Fill(locTimePull);
			dHistMap_TimePullVsTheta_CDC[locPID]->Fill(locTheta, locTimePull);
			dHistMap_TimePullVsP_CDC[locPID]->Fill(locP, locTimePull);
		}
		else if(locChargedTrackHypothesis->t1_detector() == SYS_BCAL)
		{
			dHistMap_TOFFOM_BCAL[locPID]->Fill(locFOM_Timing);
			dHistMap_TimePull_BCAL[locPID]->Fill(locTimePull);
			dHistMap_TimePullVsTheta_BCAL[locPID]->Fill(locTheta, locTimePull);
			dHistMap_TimePullVsP_BCAL[locPID]->Fill(locP, locTimePull);
		}
		else if(locChargedTrackHypothesis->t1_detector() == SYS_TOF)
		{
			dHistMap_TOFFOM_TOF[locPID]->Fill(locFOM_Timing);
			dHistMap_TimePull_TOF[locPID]->Fill(locTimePull);
			dHistMap_TimePullVsP_TOF[locPID]->Fill(locP, locTimePull);
		}
		else if(locChargedTrackHypothesis->t1_detector() == SYS_FCAL)
		{
			dHistMap_TOFFOM_FCAL[locPID]->Fill(locFOM_Timing);
			dHistMap_TimePull_FCAL[locPID]->Fill(locTimePull);
			dHistMap_TimePullVsP_FCAL[locPID]->Fill(locP, locTimePull);
		}

		dHistMap_DCdEdxFOM[locPID]->Fill(locFOM_DCdEdx);
		dHistMap_BetaVsP[locPID]->Fill(locP, locBeta_Timing);
		dHistMap_DeltaBetaVsP[locPID]->Fill(locP, locDeltaBeta);
		dHistMap_TOFFOMVsDeltaBeta[locPID]->Fill(locDeltaBeta, locFOM_Timing);
		pair<Particle_t, Particle_t> locPIDPair(locPID, Unknown); //default unless matched
		if((locMCThrown != NULL) && (locMatchFOM >= dMinThrownMatchFOM)) //else bogus track (not matched to any thrown tracks)
			locPIDPair.second = (Particle_t)(locMCThrown->type); //matched
		if(dHistMap_PIDFOMForTruePID.find(locPIDPair) != dHistMap_PIDFOMForTruePID.end()) //else hist not created or PID is weird
			dHistMap_PIDFOMForTruePID[locPIDPair]->Fill(locChargedTrackHypothesis->dFOM);

		if((locChargedTrackHypothesis->dFOM < 0.01) && (locChargedTrackHypothesis->dNDF > 0))
			dHistMap_PVsTheta_LowPIDFOM[locPID]->Fill(locTheta, locP);
		else if(locChargedTrackHypothesis->dNDF == 0) //NaN
			dHistMap_PVsTheta_NaNPIDFOM[locPID]->Fill(locTheta, locP);

		if(locFOM_Timing < 0.01)
			dHistMap_PVsTheta_LowTOFFOM[locPID]->Fill(locTheta, locP);
		else if(locChargedTrackHypothesis->dNDF_Timing == 0) //NaN
			dHistMap_PVsTheta_NaNTOFFOM[locPID]->Fill(locTheta, locP);

		if(locFOM_DCdEdx < 0.01)
			dHistMap_PVsTheta_LowDCdEdxFOM[locPID]->Fill(locTheta, locP);
		else if(locChargedTrackHypothesis->dNDF_DCdEdx == 0) //NaN
			dHistMap_PVsTheta_NaNDCdEdxFOM[locPID]->Fill(locTheta, locP);

		if(locBeta_Timing < 0.0)
			dHistMap_PVsTheta_NegativeBeta[locPID]->Fill(locTheta, locP);
	}
	japp->RootUnLock();
}

void DHistogramAction_PID::Fill_NeutralHists(const DNeutralParticleHypothesis* locNeutralParticleHypothesis, const DMCThrownMatching* locMCThrownMatching, const DEventRFBunch* locEventRFBunch)
{
	Particle_t locPID = locNeutralParticleHypothesis->PID();

	double locBeta_Timing = dAnalysisUtilities->Calc_Beta_Timing(locNeutralParticleHypothesis, locEventRFBunch);
	double locDeltaBeta = locNeutralParticleHypothesis->lorentzMomentum().Beta() - locBeta_Timing;

	double locP = locNeutralParticleHypothesis->momentum().Mag();
	double locTheta = locNeutralParticleHypothesis->momentum().Theta()*180.0/TMath::Pi();
	double locMatchFOM = 0.0;
	const DMCThrown* locMCThrown = (locMCThrownMatching != NULL) ? locMCThrownMatching->Get_MatchingMCThrown(locNeutralParticleHypothesis, locMatchFOM) : NULL;

	double locTimePull = 0.0;
	unsigned int locTimeNDF = 0;
	dParticleID->Calc_TimingChiSq(locNeutralParticleHypothesis, locEventRFBunch, locTimeNDF, locTimePull);

	japp->RootWriteLock();
	{
		dHistMap_PIDFOM[locPID]->Fill(locNeutralParticleHypothesis->dFOM);

		if(locNeutralParticleHypothesis->t1_detector() == SYS_BCAL)
		{
			dHistMap_TOFFOM_BCAL[locPID]->Fill(locNeutralParticleHypothesis->dFOM);
			dHistMap_TimePull_BCAL[locPID]->Fill(locTimePull);
			dHistMap_TimePullVsTheta_BCAL[locPID]->Fill(locTheta, locTimePull);
			dHistMap_TimePullVsP_BCAL[locPID]->Fill(locP, locTimePull);
		}
		else if(locNeutralParticleHypothesis->t1_detector() == SYS_FCAL)
		{
			dHistMap_TOFFOM_FCAL[locPID]->Fill(locNeutralParticleHypothesis->dFOM);
			dHistMap_TimePull_FCAL[locPID]->Fill(locTimePull);
			dHistMap_TimePullVsP_FCAL[locPID]->Fill(locP, locTimePull);
		}

		dHistMap_BetaVsP[locPID]->Fill(locP, locBeta_Timing);
		dHistMap_DeltaBetaVsP[locPID]->Fill(locP, locDeltaBeta);
		dHistMap_TOFFOMVsDeltaBeta[locPID]->Fill(locDeltaBeta, locNeutralParticleHypothesis->dFOM);
		pair<Particle_t, Particle_t> locPIDPair(locPID, Unknown); //default unless matched
		if((locMCThrown != NULL) && (locMatchFOM >= dMinThrownMatchFOM)) //else bogus track (not matched to any thrown tracks)
			locPIDPair.second = (Particle_t)(locMCThrown->type); //matched
		if(dHistMap_PIDFOMForTruePID.find(locPIDPair) != dHistMap_PIDFOMForTruePID.end()) //else hist not created or PID is weird
			dHistMap_PIDFOMForTruePID[locPIDPair]->Fill(locNeutralParticleHypothesis->dFOM);

		if((locNeutralParticleHypothesis->dFOM < 0.01) && (locNeutralParticleHypothesis->dNDF > 0))
			dHistMap_PVsTheta_LowPIDFOM[locPID]->Fill(locTheta, locP);
		else if(locNeutralParticleHypothesis->dNDF == 0) //NaN
			dHistMap_PVsTheta_NaNPIDFOM[locPID]->Fill(locTheta, locP);
	}
	japp->RootUnLock();
}

void DHistogramAction_DetectorStudies::Initialize(JEventLoop* locEventLoop)
{
	//Optional: Create histograms and/or modify member variables.
	//Create any histograms/trees/etc. within a ROOT lock. 
		//This is so that when running multithreaded, only one thread is writing to the ROOT file at a time. 

	//When creating a reaction-independent action, only modify member variables within a ROOT lock. 
		//Objects created within a plugin (such as reaction-independent actions) can be accessed by many threads simultaneously. 

	string locHistName, locHistTitle;

	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry* locGeometry = locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());

	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		// Optional: Useful utility functions.
		// locEventLoop->GetSingle(dAnalysisUtilities);

		locEventLoop->GetSingle(dParticleID);

		double locTargetZCenter = 0.0;
		locGeometry->GetTargetZ(locTargetZCenter);
		dTargetCenter.SetXYZ(0.0, 0.0, locTargetZCenter);

		//Required: Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this action.
			//If another thread has already created the folder, it just changes to it. 
		CreateAndChangeTo_ActionDirectory();

		CreateAndChangeTo_Directory("Reconstruction", "Reconstruction");
		{
			//TAGGER
			locHistName = "TAGHHitTime";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_TAGHHitTime = new TH1I(locHistName.c_str(), ";TAGH Hit Time (ns)", dNumTimeBins, dMinTime, dMaxTime);
			else //already created by another thread
				dHist_TAGHHitTime = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "TAGHHitTimeVsCounter";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_TAGHHitTimeVsCounter = new TH2I(locHistName.c_str(), ";TAGH Hit Counter;TAGH Hit Time (ns)", 274, 0.5, 274.5, dNumTimeBins, dMinTime, dMaxTime);
			else //already created by another thread
				dHist_TAGHHitTimeVsCounter = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "TAGMHitTime";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_TAGMHitTime = new TH1I(locHistName.c_str(), ";TAGM Hit Time (ns)", dNumTimeBins, dMinTime, dMaxTime);
			else //already created by another thread
				dHist_TAGMHitTime = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "TAGMHitTimeVsCounter";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_TAGMHitTimeVsCounter = new TH2I(locHistName.c_str(), ";TAGM Hit Counter;TAGM Hit Time (ns)", 102, 0.5, 102.5, dNumTimeBins, dMinTime, dMaxTime);
			else //already created by another thread
				dHist_TAGMHitTimeVsCounter = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

			//FCAL
			locHistName = "FCALShowerTime";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_FCALShowerTime = new TH1I(locHistName.c_str(), ";FCAL Shower Time (ns)", dNumTimeBins, dMinTime, dMaxTime);
			else //already created by another thread
				dHist_FCALShowerTime = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "FCALShowerEnergy";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_FCALShowerEnergy = new TH1I(locHistName.c_str(), ";FCAL Shower Energy (GeV)", dNumShowerEnergyBins, dMinShowerEnergy, dMaxShowerEnergy);
			else //already created by another thread
				dHist_FCALShowerEnergy = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "FCALShowerEnergyVsTime";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_FCALShowerEnergyVsTime = new TH2I(locHistName.c_str(), ";FCAL Shower Time (ns);FCAL Shower Energy (GeV)", dNum2DTimeBins, dMinTime, dMaxTime, dNum2DShowerEnergyBins, dMinShowerEnergy, dMaxShowerEnergy);
			else //already created by another thread
				dHist_FCALShowerEnergyVsTime = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "FCALShowerYVsX";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_FCALShowerYVsX = new TH2I(locHistName.c_str(), ";FCAL Shower X (cm);FCAL Shower Y (cm)", dNumFCALTOFXYBins, -120.0, 120.0, dNumFCALTOFXYBins, -120.0, 120.0);
			else //already created by another thread
				dHist_FCALShowerYVsX = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

			//BCAL
			locHistName = "BCALShowerPhi";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_BCALShowerPhi = new TH1I(locHistName.c_str(), ";BCAL Shower #phi#circ", dNumPhiBins, dMinPhi, dMaxPhi);
			else //already created by another thread
				dHist_BCALShowerPhi = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "BCALShowerPhiVsZ";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_BCALShowerPhiVsZ = new TH2I(locHistName.c_str(), ";BCAL Shower Z (cm);BCAL Shower #phi#circ", dNum2DBCALZBins, 0.0, 450.0, dNum2DPhiBins, dMinPhi, dMaxPhi);
			else //already created by another thread
				dHist_BCALShowerPhiVsZ = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "BCALShowerTime";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_BCALShowerTime = new TH1I(locHistName.c_str(), ";BCAL Shower Time (ns)", dNumTimeBins, dMinTime, dMaxTime);
			else //already created by another thread
				dHist_BCALShowerTime = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "BCALShowerEnergy";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_BCALShowerEnergy = new TH1I(locHistName.c_str(), ";BCAL Shower Energy (GeV)", dNumShowerEnergyBins, dMinShowerEnergy, dMaxBCALP);
			else //already created by another thread
				dHist_BCALShowerEnergy = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "BCALShowerEnergyVsTime";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_BCALShowerEnergyVsTime = new TH2I(locHistName.c_str(), ";BCAL Shower Time (ns);BCAL Shower Energy (GeV)", dNum2DTimeBins, dMinTime, dMaxTime, dNum2DShowerEnergyBins, dMinShowerEnergy, dMaxBCALP);
			else //already created by another thread
				dHist_BCALShowerEnergyVsTime = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

			//TOF
			locHistName = "TOFPointTime";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_TOFPointTime = new TH1I(locHistName.c_str(), ";TOF Point Time (ns)", dNumTimeBins, dMinTime, dMaxTime);
			else //already created by another thread
				dHist_TOFPointTime = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "TOFPointEnergy";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_TOFPointEnergy = new TH1I(locHistName.c_str(), ";TOF Point Energy (MeV)", dNumHitEnergyBins, dMinHitEnergy, dMaxHitEnergy);
			else //already created by another thread
				dHist_TOFPointEnergy = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "TOFPointYVsX";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_TOFPointYVsX = new TH2I(locHistName.c_str(), ";TOF Point X (cm);TOF Point Y (cm)", dNumFCALTOFXYBins, -120.0, 120.0, dNumFCALTOFXYBins, -120.0, 120.0);
			else //already created by another thread
				dHist_TOFPointYVsX = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

			//SC
			locHistName = "SCHitTime";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_SCHitTime = new TH1I(locHistName.c_str(), ";SC Hit Time (ns)", dNumTimeBins, dMinTime, dMaxTime);
			else //already created by another thread
				dHist_SCHitTime = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "SCHitSector";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_SCHitSector = new TH1I(locHistName.c_str(), ";SC Hit Sector", 24, 0.5, 24.5);
			else //already created by another thread
				dHist_SCHitSector = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "SCHitTimeVsSector";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_SCHitTimeVsSector = new TH2I(locHistName.c_str(), ";SC Hit Sector;SC Hit Time (ns)", 24, 0.5, 24.5, dNum2DTimeBins, dMinTime, dMaxTime);
			else //already created by another thread
				dHist_SCHitTimeVsSector = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "SCHitEnergy";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_SCHitEnergy = new TH1I(locHistName.c_str(), ";SC Hit Energy (MeV)", dNumHitEnergyBins, dMinHitEnergy, dMaxHitEnergy);
			else //already created by another thread
				dHist_SCHitEnergy = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "SCHitEnergyVsSector";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_SCHitEnergyVsSector = new TH2I(locHistName.c_str(), ";SC Hit Sector;SC Hit Energy (MeV)", 24, 0.5, 24.5, dNum2DHitEnergyBins, dMinHitEnergy, dMaxHitEnergy);
			else //already created by another thread
				dHist_SCHitEnergyVsSector = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

			//TRACKING
			locHistName = "NumDCHitsPerTrack";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_NumDCHitsPerTrack = new TH1I(locHistName.c_str(), ";# Track Hits", 50, 0.5, 50.5);
			else //already created by another thread
				dHist_NumDCHitsPerTrack = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "NumDCHitsPerTrackVsTheta";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_NumDCHitsPerTrackVsTheta = new TH2I(locHistName.c_str(), ";#theta#circ;# Track Hits", dNum2DThetaBins, dMinTheta, dMaxTheta, 50, 0.5, 50.5);
			else //already created by another thread
				dHist_NumDCHitsPerTrackVsTheta = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

			//Track Matched to Hit
			for(int locTruePIDFlag = 0; locTruePIDFlag < 2; ++locTruePIDFlag)
			{
				if(locMCThrowns.empty() && (locTruePIDFlag == 1))
					continue; //not a simulated event: don't histogram thrown info!

				string locDirName = (locTruePIDFlag == 1) ? "TruePID" : "ReconstructedPID";
				CreateAndChangeTo_Directory(locDirName.c_str(), locDirName.c_str());
				for(size_t loc_i = 0; loc_i < dTrackingPIDs.size(); ++loc_i)
				{
					Particle_t locPID = dTrackingPIDs[loc_i];
					string locParticleName = ParticleType(locPID);
					string locParticleROOTName = ParticleName_ROOT(locPID);
					CreateAndChangeTo_Directory(locParticleName, locParticleName);
					pair<Particle_t, bool> locPIDPair(locPID, bool(locTruePIDFlag));

					// CDC dE/dx
					locHistName = "dEdx_CDC";
					locHistTitle = locParticleROOTName + string(";CDC dE/dx (keV/cm)");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistMap_CDCdEdx[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistMap_CDCdEdx[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumdEdxBins, dMindEdX, dMaxdEdX);

					// CDC dE/dx Vs P 
					locHistName = "dEdxVsP_CDC";
					locHistTitle = locParticleROOTName + string(";p (GeV/c);CDC dE/dx (keV/cm)");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistMap_CDCdEdxVsP[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistMap_CDCdEdxVsP[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DdEdxBins, dMindEdX, dMaxdEdX);

					// FDC dE/dx
					locHistName = "dEdx_FDC";
					locHistTitle = locParticleROOTName + string(";FDC dE/dx (keV/cm)");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistMap_FDCdEdx[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistMap_FDCdEdx[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumdEdxBins, dMindEdX, dMaxdEdX);

					// FDC dE/dx Vs P 
					locHistName = "dEdxVsP_FDC";
					locHistTitle = locParticleROOTName + string(";p (GeV/c);FDC dE/dx (keV/cm)");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistMap_FDCdEdxVsP[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistMap_FDCdEdxVsP[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DdEdxBins, dMindEdX, dMaxdEdX);

					// Tracking Chisq/NDF
					locHistName = "Tracking_ChiSqPerDF";
					locHistTitle = locParticleROOTName + string(";Tracking #chi^{2}/NDF");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistMap_TrackingChiSqPerDF[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistMap_TrackingChiSqPerDF[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumTrackingChiSqPerDFBins, dMinTrackingChiSqPerDF, dMaxTrackingChiSqPerDF);

					// Tracking Chisq/NDF Vs Theta
					locHistName = "Tracking_ChiSqPerDFVsTheta";
					locHistTitle = locParticleROOTName + string(";#theta#circ;Tracking #chi^{2}/NDF");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistMap_TrackingChiSqPerDFVsTheta[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistMap_TrackingChiSqPerDFVsTheta[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNumTrackingChiSqPerDFBins, dMinTrackingChiSqPerDF, dMaxTrackingChiSqPerDF);

					// Tracking Chisq/NDF Vs P
					locHistName = "Tracking_ChiSqPerDFVsP";
					locHistTitle = locParticleROOTName + string(";p (GeV/c);Tracking #chi^{2}/NDF");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistMap_TrackingChiSqPerDFVsP[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistMap_TrackingChiSqPerDFVsP[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumTrackingChiSqPerDFBins, dMinTrackingChiSqPerDF, dMaxTrackingChiSqPerDF);

					gDirectory->cd("..");
				}
				gDirectory->cd("..");
			}
		}
		gDirectory->cd("..");

		//Showers not matched to tracks, Tracks not matched to hits
		CreateAndChangeTo_Directory("Not-Matched", "Not-Matched");
		{
			//BCAL
			locHistName = "BCALTrackDOCA";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_BCALTrackDOCA = new TH1I(locHistName.c_str(), ";BCAL Shower Distance to Nearest Track (cm)", dNumTrackDOCABins, dMinTrackDOCA, dMaxTrackDOCA);
			else //already created by another thread
				dHist_BCALTrackDOCA = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "BCALNeutralShowerTime";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_BCALNeutralShowerTime = new TH1I(locHistName.c_str(), ";BCAL Neutral Shower Time (ns)", dNumTimeBins, dMinTime, dMaxTime);
			else //already created by another thread
				dHist_BCALNeutralShowerTime = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "BCALNeutralShowerEnergy";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_BCALNeutralShowerEnergy = new TH1I(locHistName.c_str(), ";BCAL Neutral Shower Energy (GeV)", dNumShowerEnergyBins, dMinShowerEnergy, dMaxBCALP);
			else //already created by another thread
				dHist_BCALNeutralShowerEnergy = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "BCALNeutralShowerDeltaT";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_BCALNeutralShowerDeltaT = new TH1I(locHistName.c_str(), ";BCAL Neutral Shower #Deltat (Propagated - RF) (ns)", dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
			else //already created by another thread
				dHist_BCALNeutralShowerDeltaT = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "BCALNeutralShowerDeltaTVsE";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_BCALNeutralShowerDeltaTVsE = new TH2I(locHistName.c_str(), ";BCAL Neutral Shower Energy (GeV);BCAL Neutral Shower #Deltat (ns)", dNum2DShowerEnergyBins, dMinShowerEnergy, dMaxBCALP, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);
			else //already created by another thread
				dHist_BCALNeutralShowerDeltaTVsE = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "BCALNeutralShowerDeltaTVsZ";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_BCALNeutralShowerDeltaTVsZ = new TH2I(locHistName.c_str(), ";BCAL Neutral Shower Z (cm);BCAL Neutral Shower #Deltat (ns)", dNum2DBCALZBins, 0.0, 450.0, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);
			else //already created by another thread
				dHist_BCALNeutralShowerDeltaTVsZ = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

			//FCAL
			locHistName = "FCALTrackDOCA";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_FCALTrackDOCA = new TH1I(locHistName.c_str(), ";FCAL Shower Distance to Nearest Track (cm)", dNumTrackDOCABins, dMinTrackDOCA, dMaxTrackDOCA);
			else //already created by another thread
				dHist_FCALTrackDOCA = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "FCALNeutralShowerTime";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_FCALNeutralShowerTime = new TH1I(locHistName.c_str(), ";FCAL Neutral Shower Time (ns)", dNumTimeBins, dMinTime, dMaxTime);
			else //already created by another thread
				dHist_FCALNeutralShowerTime = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "FCALNeutralShowerEnergy";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_FCALNeutralShowerEnergy = new TH1I(locHistName.c_str(), ";FCAL Neutral Shower Energy (GeV)", dNumShowerEnergyBins, dMinShowerEnergy, dMaxShowerEnergy);
			else //already created by another thread
				dHist_FCALNeutralShowerEnergy = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "FCALNeutralShowerDeltaT";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_FCALNeutralShowerDeltaT = new TH1I(locHistName.c_str(), ";FCAL Neutral Shower #Deltat (Propagated - RF) (ns)", dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
			else //already created by another thread
				dHist_FCALNeutralShowerDeltaT = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

			locHistName = "FCALNeutralShowerDeltaTVsE";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_FCALNeutralShowerDeltaTVsE = new TH2I(locHistName.c_str(), ";FCAL Neutral Shower Energy (GeV);FCAL Neutral Shower #Deltat (ns)", dNum2DShowerEnergyBins, dMinShowerEnergy, dMaxShowerEnergy, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);
			else //already created by another thread
				dHist_FCALNeutralShowerDeltaTVsE = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

			//TRACKING
			locHistName = "TrackPVsTheta_NoHitMatch";
			if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
				dHist_TrackPVsTheta_NoHitMatch = new TH2I(locHistName.c_str(), ";#theta#circ;p (GeV/c)", dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);
			else //already created by another thread
				dHist_TrackPVsTheta_NoHitMatch = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
		}
		gDirectory->cd("..");

		//PID
		CreateAndChangeTo_Directory("dEdxVsPByCharge", "dEdxVsPByCharge");
		{
			//q +/-
			for(int locCharge = -1; locCharge <= 1; locCharge += 2)
			{
				if(locCharge == -1)
					CreateAndChangeTo_Directory("q-", "q-");
				else
					CreateAndChangeTo_Directory("q+", "q+");
				string locParticleROOTName = (locCharge == -1) ? "q^{-}" : "q^{+}";

				locHistName = "TOFdEdXVsP";
				locHistTitle = locParticleROOTName + ";p (GeV/c);TOF dE/dX (MeV/cm)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_QTOFdEdXVsP[locCharge] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DdEdxBins, dMindEdX, dMaxdEdX);
				else //already created by another thread
					dHistMap_QTOFdEdXVsP[locCharge] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "SCdEdXVsP";
				locHistTitle = locParticleROOTName + ";p (GeV/c);SC dE/dX (MeV/cm)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_QSCdEdXVsP[locCharge] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DdEdxBins, dMindEdX, dMaxdEdX);
				else //already created by another thread
					dHistMap_QSCdEdXVsP[locCharge] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "dEdxVsP_CDC";
				locHistTitle = locParticleROOTName + string(";p (GeV/c);CDC dE/dx (keV/cm)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_QCDCdEdXVsP[locCharge] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_QCDCdEdXVsP[locCharge] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DdEdxBins, dMindEdX, dMaxdEdX);

				locHistName = "dEdxVsP_FDC";
				locHistTitle = locParticleROOTName + string(";p (GeV/c);FDC dE/dx (keV/cm)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_QFDCdEdXVsP[locCharge] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_QFDCdEdXVsP[locCharge] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DdEdxBins, dMindEdX, dMaxdEdX);

				gDirectory->cd("..");
			}
		}
		gDirectory->cd("..");

		//Track Matched to Hit
		for(int locTruePIDFlag = 0; locTruePIDFlag < 2; ++locTruePIDFlag)
		{
			if(locMCThrowns.empty() && (locTruePIDFlag == 1))
				continue; //not a simulated event: don't histogram thrown info!

			string locDirName = (locTruePIDFlag == 1) ? "Matched_TruePID" : "Matched_ReconstructedPID";
			CreateAndChangeTo_Directory(locDirName.c_str(), locDirName.c_str());

			//By PID
			for(size_t loc_i = 0; loc_i < dTrackingPIDs.size(); ++loc_i)
			{
				Particle_t locPID = dTrackingPIDs[loc_i];
				string locParticleName = ParticleType(locPID);
				string locParticleROOTName = ParticleName_ROOT(locPID);
				CreateAndChangeTo_Directory(locParticleName, locParticleName);
				pair<Particle_t, bool> locPIDPair(locPID, bool(locTruePIDFlag));

				//BCAL
				locHistName = "BCALTrackDOCA";
				locHistTitle = locParticleROOTName + ";BCAL Shower Distance to Track (cm)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_BCALTrackDOCA[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);
				else //already created by another thread
					dHistMap_BCALTrackDOCA[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "BCALShowerEnergy";
				locHistTitle = locParticleROOTName + ";BCAL Shower Energy (GeV)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_BCALShowerEnergy[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumShowerEnergyBins, dMinShowerEnergy, dMaxBCALP);
				else //already created by another thread
					dHistMap_BCALShowerEnergy[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "BCALShowerTrackDepth";
				locHistTitle = locParticleROOTName + ";BCAL Shower Track Depth (cm)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_BCALShowerTrackDepth[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumShowerDepthBins, dMinShowerDepth, dMaxShowerDepth);
				else //already created by another thread
					dHistMap_BCALShowerTrackDepth[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "BCALShowerTrackDepthVsP";
				locHistTitle = locParticleROOTName + ";p (GeV/c);BCAL Shower Track Depth (cm)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_BCALShowerTrackDepthVsP[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxBCALP, dNumShowerDepthBins, dMinShowerDepth, dMaxShowerDepth);
				else //already created by another thread
					dHistMap_BCALShowerTrackDepthVsP[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "BCALShowerDeltaT";
				locHistTitle = locParticleROOTName + ";BCAL Shower #Deltat (Propagated - RF) (ns)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_BCALShowerDeltaT[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
				else //already created by another thread
					dHistMap_BCALShowerDeltaT[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "BCALShowerDeltaTVsZ";
				locHistTitle = locParticleROOTName + ";BCAL Shower Z (cm);BCAL Shower #Deltat (ns)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_BCALShowerDeltaTVsZ[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DBCALZBins, 0.0, 450.0, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);
				else //already created by another thread
					dHistMap_BCALShowerDeltaTVsZ[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "BCALShowerDeltaTVsP";
				locHistTitle = locParticleROOTName + ";p (GeV/c);BCAL Shower #Deltat (ns)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_BCALShowerDeltaTVsP[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxBCALP, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);
				else //already created by another thread
					dHistMap_BCALShowerDeltaTVsP[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));


				//FCAL
				locHistName = "FCALTrackDOCA";
				locHistTitle = locParticleROOTName + ";FCAL Shower Distance to Track (cm)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_FCALTrackDOCA[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);
				else //already created by another thread
					dHistMap_FCALTrackDOCA[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "FCALShowerEnergy";
				locHistTitle = locParticleROOTName + ";FCAL Shower Energy (GeV)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_FCALShowerEnergy[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumShowerEnergyBins, dMinShowerEnergy, dMaxShowerEnergy);
				else //already created by another thread
					dHistMap_FCALShowerEnergy[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "FCALShowerTrackDepth";
				locHistTitle = locParticleROOTName + ";FCAL Shower Track Depth (cm)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_FCALShowerTrackDepth[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumShowerDepthBins, dMinShowerDepth, dMaxShowerDepth);
				else //already created by another thread
					dHistMap_FCALShowerTrackDepth[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "FCALShowerTrackDepthVsP";
				locHistTitle = locParticleROOTName + ";p (GeV/c);FCAL Shower Track Depth (cm)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_FCALShowerTrackDepthVsP[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumShowerDepthBins, dMinShowerDepth, dMaxShowerDepth);
				else //already created by another thread
					dHistMap_FCALShowerTrackDepthVsP[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "FCALShowerDeltaT";
				locHistTitle = locParticleROOTName + ";FCAL Shower #Deltat (Propagated - RF) (ns)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_FCALShowerDeltaT[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
				else //already created by another thread
					dHistMap_FCALShowerDeltaT[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "FCALShowerDeltaTVsP";
				locHistTitle = locParticleROOTName + ";p (GeV/c);FCAL Shower #Deltat (ns)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_FCALShowerDeltaTVsP[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);
				else //already created by another thread
					dHistMap_FCALShowerDeltaTVsP[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));


				//TOF
				locHistName = "TOFdEdX";
				locHistTitle = locParticleROOTName + ";TOF dE/dX (MeV/cm)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_TOFdEdX[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumdEdxBins, dMindEdX, dMaxdEdX);
				else //already created by another thread
					dHistMap_TOFdEdX[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "TOFdEdXVsP";
				locHistTitle = locParticleROOTName + ";p (GeV/c);TOF dE/dX (MeV/cm)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_TOFdEdXVsP[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DdEdxBins, dMindEdX, dMaxdEdX);
				else //already created by another thread
					dHistMap_TOFdEdXVsP[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "TOFTrackDOCA";
				locHistTitle = locParticleROOTName + ";TOF Point Distance to Track (cm)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_TOFTrackDOCA[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);
				else //already created by another thread
					dHistMap_TOFTrackDOCA[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "TOFDeltaT";
				locHistTitle = locParticleROOTName + ";TOF Point #Deltat (Propagated - RF) (ns)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_TOFDeltaT[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
				else //already created by another thread
					dHistMap_TOFDeltaT[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "TOFDeltaTVsP";
				locHistTitle = locParticleROOTName + ";p (GeV/c);TOF Point #Deltat (ns)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_TOFDeltaTVsP[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);
				else //already created by another thread
					dHistMap_TOFDeltaTVsP[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

				//SC
				locHistName = "SCdEdX";
				locHistTitle = locParticleROOTName + ";SC dE/dX (MeV/cm)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_SCdEdX[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumdEdxBins, dMindEdX, dMaxdEdX);
				else //already created by another thread
					dHistMap_SCdEdX[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "SCdEdXVsP";
				locHistTitle = locParticleROOTName + ";p (GeV/c);SC dE/dX (MeV/cm)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_SCdEdXVsP[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DdEdxBins, dMindEdX, dMaxdEdX);
				else //already created by another thread
					dHistMap_SCdEdXVsP[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "SCTrackDeltaPhi";
				locHistTitle = locParticleROOTName + ";SC Point Track #Delta#phi#circ";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_SCTrackDeltaPhi[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);
				else //already created by another thread
					dHistMap_SCTrackDeltaPhi[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "SCTime";
				locHistTitle = locParticleROOTName + ";#theta#circ;SC Point Time (ns)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_SCTime[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumTimeBins, dMinTime, dMaxTime);
				else //already created by another thread
					dHistMap_SCTime[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "SCTimeVsTheta";
				locHistTitle = locParticleROOTName + ";#theta#circ;SC Point Time (ns)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_SCTimeVsTheta[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DTimeBins, dMinTime, dMaxTime);
				else //already created by another thread
					dHistMap_SCTimeVsTheta[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "SCEnergyVsTheta";
				locHistTitle = locParticleROOTName + ";#theta#circ;SC Point Energy (MeV)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_SCEnergyVsTheta[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DHitEnergyBins, dMinHitEnergy, dMaxHitEnergy);
				else //already created by another thread
					dHistMap_SCEnergyVsTheta[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "SCPhiVsTheta";
				locHistTitle = locParticleROOTName + ";#theta#circ;#phi#circ";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_SCPhiVsTheta[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);
				else //already created by another thread
					dHistMap_SCPhiVsTheta[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "SCDeltaT";
				locHistTitle = locParticleROOTName + ";SC Point #Deltat (Propagated - RF) (ns)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_SCDeltaT[locPIDPair] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
				else //already created by another thread
					dHistMap_SCDeltaT[locPIDPair] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "SCDeltaTVsP";
				locHistTitle = locParticleROOTName + ";p (GeV/c);SC Point #Deltat (Propagated - RF) (ns)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_SCDeltaTVsP[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);
				else //already created by another thread
					dHistMap_SCDeltaTVsP[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "SCDeltaTVsPhi";
				locHistTitle = locParticleROOTName + ";#phi#circ;SC Point #Deltat (ns)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_SCDeltaTVsPhi[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPhiBins, dMinPhi, dMaxPhi, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);
				else //already created by another thread
					dHistMap_SCDeltaTVsPhi[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

				locHistName = "SCDeltaTVsTheta";
				locHistTitle = locParticleROOTName + ";#theta#circ;SC Point #Deltat (ns)";
				if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
					dHistMap_SCDeltaTVsTheta[locPIDPair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);
				else //already created by another thread
					dHistMap_SCDeltaTVsTheta[locPIDPair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));

				gDirectory->cd("..");
			}
			gDirectory->cd("..");
		}
		gDirectory->cd("..");
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_DetectorStudies::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//Expect locParticleCombo to be NULL since this is a reaction-independent action.

	//Optional: Quit the action if it has already been executed this event (else may result in double-counting when filling histograms)
	if(Get_NumPreviousParticleCombos() != 0)
		return true;

	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	Fill_ReconstructionHists(locEventLoop, false);
	Fill_NotMatchedHists(locEventLoop);
	Fill_MatchedHists(locEventLoop, false);
	Fill_PIDHists(locEventLoop);

	if(!locMCThrowns.empty())
	{
		Fill_ReconstructionHists(locEventLoop, true);
		Fill_MatchedHists(locEventLoop, true);
	}

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}

void DHistogramAction_DetectorStudies::Fill_ReconstructionHists(JEventLoop* locEventLoop, bool locUseTruePIDFlag)
{
	if(locUseTruePIDFlag)
	{
		vector<const DChargedTrack*> locChargedTracks;
		locEventLoop->Get(locChargedTracks);

		vector<const DMCThrownMatching*> locMCThrownMatchingVector;
		locEventLoop->Get(locMCThrownMatchingVector);

		japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
		{
			for(size_t loc_i = 0; loc_i < locChargedTracks.size(); ++loc_i)
			{
				double locMatchFOM = 0.0;
				const DMCThrown* locMCThrown = locMCThrownMatchingVector[0]->Get_MatchingMCThrown(locChargedTracks[loc_i], locMatchFOM);
				if((locMCThrown == NULL) || (locMatchFOM < dMinThrownMatchFOM))
					continue;

				//OK, have the thrown. Now, grab the best charged track hypothesis to get the best matching
				const DChargedTrackHypothesis* locChargedTrackHypothesis = locMCThrownMatchingVector[0]->Get_MatchingChargedHypothesis(locMCThrown, locMatchFOM);
				if(locChargedTrackHypothesis->PID() != locMCThrown->PID())
					continue;

				pair<Particle_t, bool> locPIDPair(locChargedTrackHypothesis->PID(), true);
				if(dHistMap_TrackingChiSqPerDF.find(locPIDPair) == dHistMap_TrackingChiSqPerDF.end())
					continue; //disregard PID

				const DTrackTimeBased* locTrackTimeBased = NULL;
				locChargedTrackHypothesis->GetSingle(locTrackTimeBased);
				double locP = locTrackTimeBased->momentum().Mag();
				double locTheta = locTrackTimeBased->momentum().Theta()*180.0/TMath::Pi();

				double locChiSqPerDF = (locTrackTimeBased->Ndof > 0) ? locTrackTimeBased->chisq/locTrackTimeBased->Ndof : numeric_limits<double>::quiet_NaN();
				dHistMap_TrackingChiSqPerDF[locPIDPair]->Fill(locChiSqPerDF);
				dHistMap_TrackingChiSqPerDFVsTheta[locPIDPair]->Fill(locTheta, locChiSqPerDF);
				dHistMap_TrackingChiSqPerDFVsP[locPIDPair]->Fill(locP, locChiSqPerDF);
				if(locTrackTimeBased->dNumHitsUsedFordEdx_CDC > 0)
				{
					dHistMap_CDCdEdx[locPIDPair]->Fill(locTrackTimeBased->ddEdx_CDC*1.0E6);
					dHistMap_CDCdEdxVsP[locPIDPair]->Fill(locP, locTrackTimeBased->ddEdx_CDC*1.0E6);
				}
				if(locTrackTimeBased->dNumHitsUsedFordEdx_FDC > 0)
				{
					dHistMap_FDCdEdx[locPIDPair]->Fill(locTrackTimeBased->ddEdx_FDC*1.0E6);
					dHistMap_FDCdEdxVsP[locPIDPair]->Fill(locP, locTrackTimeBased->ddEdx_FDC*1.0E6);
				}
			}
		}
		japp->RootUnLock(); //RELEASE ROOT LOCK!!
		return;
	}

	vector<const DTAGHHit*> locTAGHHits;
	locEventLoop->Get(locTAGHHits);

	vector<const DTAGMHit*> locTAGMHits;
	locEventLoop->Get(locTAGMHits);

	vector<const DBCALShower*> locBCALShowers;
	locEventLoop->Get(locBCALShowers);

	vector<const DFCALShower*> locFCALShowers;
	locEventLoop->Get(locFCALShowers);

	vector<const DTOFPoint*> locTOFPoints;
	locEventLoop->Get(locTOFPoints);

	vector<const DSCHit*> locSCHits;
	locEventLoop->Get(locSCHits);

	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locEventLoop->Get(locTrackTimeBasedVector);

	//Fill Histograms
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
		{
			dHist_NumDCHitsPerTrack->Fill(locTrackTimeBasedVector[loc_i]->Ndof + 5);
			dHist_NumDCHitsPerTrackVsTheta->Fill(locTrackTimeBasedVector[loc_i]->momentum().Theta()*180.0/TMath::Pi(), locTrackTimeBasedVector[loc_i]->Ndof + 5);

			pair<Particle_t, bool> locPIDPair(locTrackTimeBasedVector[loc_i]->PID(), false);
			if(dHistMap_TrackingChiSqPerDF.find(locPIDPair) == dHistMap_TrackingChiSqPerDF.end())
				continue; //disregard PID

			double locP = locTrackTimeBasedVector[loc_i]->momentum().Mag();
			double locTheta = locTrackTimeBasedVector[loc_i]->momentum().Theta()*180.0/TMath::Pi();

			double locChiSqPerDF = (locTrackTimeBasedVector[loc_i]->Ndof > 0) ? locTrackTimeBasedVector[loc_i]->chisq/locTrackTimeBasedVector[loc_i]->Ndof : numeric_limits<double>::quiet_NaN();
			dHistMap_TrackingChiSqPerDF[locPIDPair]->Fill(locChiSqPerDF);
			dHistMap_TrackingChiSqPerDFVsTheta[locPIDPair]->Fill(locTheta, locChiSqPerDF);
			dHistMap_TrackingChiSqPerDFVsP[locPIDPair]->Fill(locP, locChiSqPerDF);

			if(locTrackTimeBasedVector[loc_i]->dNumHitsUsedFordEdx_CDC > 0)
			{
				dHistMap_CDCdEdx[locPIDPair]->Fill(locTrackTimeBasedVector[loc_i]->ddEdx_CDC*1.0E6);
				dHistMap_CDCdEdxVsP[locPIDPair]->Fill(locP, locTrackTimeBasedVector[loc_i]->ddEdx_CDC*1.0E6);
			}
			if(locTrackTimeBasedVector[loc_i]->dNumHitsUsedFordEdx_FDC > 0)
			{
				dHistMap_FDCdEdx[locPIDPair]->Fill(locTrackTimeBasedVector[loc_i]->ddEdx_FDC*1.0E6);
				dHistMap_FDCdEdxVsP[locPIDPair]->Fill(locP, locTrackTimeBasedVector[loc_i]->ddEdx_FDC*1.0E6);
			}
		}

		for(size_t loc_i = 0; loc_i < locTAGHHits.size(); ++loc_i)
		{
			dHist_TAGHHitTime->Fill(locTAGHHits[loc_i]->t);
			dHist_TAGHHitTimeVsCounter->Fill(locTAGHHits[loc_i]->counter_id, locTAGHHits[loc_i]->t);
		}

		for(size_t loc_i = 0; loc_i < locTAGMHits.size(); ++loc_i)
		{
			dHist_TAGMHitTime->Fill(locTAGMHits[loc_i]->t);
			dHist_TAGMHitTimeVsCounter->Fill(locTAGMHits[loc_i]->column, locTAGMHits[loc_i]->t);
		}

		for(size_t loc_i = 0; loc_i < locFCALShowers.size(); ++loc_i)
		{
			dHist_FCALShowerTime->Fill(locFCALShowers[loc_i]->getTime());
			dHist_FCALShowerEnergy->Fill(locFCALShowers[loc_i]->getEnergy());
			dHist_FCALShowerEnergyVsTime->Fill(locFCALShowers[loc_i]->getTime(), locFCALShowers[loc_i]->getEnergy());
			dHist_FCALShowerYVsX->Fill(locFCALShowers[loc_i]->getPosition().X(), locFCALShowers[loc_i]->getPosition().Y());
		}

		for(size_t loc_i = 0; loc_i < locBCALShowers.size(); ++loc_i)
		{
			dHist_BCALShowerTime->Fill(locBCALShowers[loc_i]->t);
			dHist_BCALShowerEnergy->Fill(locBCALShowers[loc_i]->E);
			dHist_BCALShowerEnergyVsTime->Fill(locBCALShowers[loc_i]->t, locBCALShowers[loc_i]->E);

			DVector3 locBCALPosition(locBCALShowers[loc_i]->x, locBCALShowers[loc_i]->y, locBCALShowers[loc_i]->z);
			double locBCALPhi = locBCALPosition.Phi()*180.0/TMath::Pi();
			dHist_BCALShowerPhi->Fill(locBCALPhi);
			dHist_BCALShowerPhiVsZ->Fill(locBCALPosition.Z(), locBCALPhi);
		}

		for(size_t loc_i = 0; loc_i < locTOFPoints.size(); ++loc_i)
		{
			dHist_TOFPointTime->Fill(locTOFPoints[loc_i]->t);
			dHist_TOFPointEnergy->Fill(locTOFPoints[loc_i]->dE*1.0E3);
			dHist_TOFPointYVsX->Fill(locTOFPoints[loc_i]->pos.X(), locTOFPoints[loc_i]->pos.Y());
		}

		for(size_t loc_i = 0; loc_i < locSCHits.size(); ++loc_i)
		{
			dHist_SCHitSector->Fill(locSCHits[loc_i]->sector);
			dHist_SCHitTime->Fill(locSCHits[loc_i]->t);
			dHist_SCHitTimeVsSector->Fill(locSCHits[loc_i]->sector, locSCHits[loc_i]->t);
			dHist_SCHitEnergy->Fill(locSCHits[loc_i]->dE*1.0E3);
			dHist_SCHitEnergyVsSector->Fill(locSCHits[loc_i]->sector, locSCHits[loc_i]->dE*1.0E3);
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

void DHistogramAction_DetectorStudies::Fill_PIDHists(JEventLoop* locEventLoop)
{
	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locEventLoop->Get(locTrackTimeBasedVector);

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	//Fill Histograms
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
		{
			int locCharge = ParticleCharge(locTrackTimeBasedVector[loc_i]->PID());
			if(dHistMap_QCDCdEdXVsP.find(locCharge) == dHistMap_QCDCdEdXVsP.end())
				continue;

			double locP = locTrackTimeBasedVector[loc_i]->momentum().Mag();

			DSCHitMatchParams locSCHitMatchParams;
			if(dParticleID->Get_BestSCMatchParams(locTrackTimeBasedVector[loc_i], locDetectorMatches, locSCHitMatchParams))
			{
				if(locSCHitMatchParams.dTrackTimeBased != NULL)
					dHistMap_QSCdEdXVsP[locCharge]->Fill(locP, locSCHitMatchParams.dEdx*1.0E3);
			}
			DTOFHitMatchParams locTOFHitMatchParams;
			if(dParticleID->Get_BestTOFMatchParams(locTrackTimeBasedVector[loc_i], locDetectorMatches, locTOFHitMatchParams))
			{
				if(locTOFHitMatchParams.dTrackTimeBased != NULL)
					dHistMap_QTOFdEdXVsP[locCharge]->Fill(locP, locTOFHitMatchParams.dEdx*1.0E3);
			}

			if(locTrackTimeBasedVector[loc_i]->dNumHitsUsedFordEdx_CDC > 0)
				dHistMap_QCDCdEdXVsP[locCharge]->Fill(locP, locTrackTimeBasedVector[loc_i]->ddEdx_CDC*1.0E6);
			if(locTrackTimeBasedVector[loc_i]->dNumHitsUsedFordEdx_FDC > 0)
				dHistMap_QFDCdEdXVsP[locCharge]->Fill(locP, locTrackTimeBasedVector[loc_i]->ddEdx_FDC*1.0E6);
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

void DHistogramAction_DetectorStudies::Fill_NotMatchedHists(JEventLoop* locEventLoop)
{
	vector<const DNeutralShower*> locNeutralShowers;
	locEventLoop->Get(locNeutralShowers);

	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locEventLoop->Get(locTrackTimeBasedVector);

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	vector<const DEventRFBunch*> locEventRFBunches;
	locEventLoop->Get(locEventRFBunches);
	double locStartTime = locEventRFBunches.empty() ? 0.0 : locEventRFBunches[0]->dTime;

	//Fill Histograms
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		for(size_t loc_i = 0; loc_i < locNeutralShowers.size(); ++loc_i)
		{
			//assume is photon
			double locPathLength = (locNeutralShowers[loc_i]->dSpacetimeVertex.Vect() - dTargetCenter).Mag();
			double locDeltaT = locNeutralShowers[loc_i]->dSpacetimeVertex.T() - locPathLength/29.9792458 - locStartTime;

			double locDistance = 9.9E9;
			if(locNeutralShowers[loc_i]->dDetectorSystem == SYS_FCAL)
			{
				const DFCALShower* locFCALShower = NULL;
				locNeutralShowers[loc_i]->GetSingle(locFCALShower);

				if(locDetectorMatches->Get_DistanceToNearestTrack(locFCALShower, locDistance))
					dHist_FCALTrackDOCA->Fill(locDistance);

				dHist_FCALNeutralShowerTime->Fill(locNeutralShowers[loc_i]->dSpacetimeVertex.T());
				dHist_FCALNeutralShowerEnergy->Fill(locNeutralShowers[loc_i]->dEnergy);

				dHist_FCALNeutralShowerDeltaT->Fill(locDeltaT);
				dHist_FCALNeutralShowerDeltaTVsE->Fill(locNeutralShowers[loc_i]->dEnergy, locDeltaT);
			}
			else
			{
				const DBCALShower* locBCALShower = NULL;
				locNeutralShowers[loc_i]->GetSingle(locBCALShower);

				if(locDetectorMatches->Get_DistanceToNearestTrack(locBCALShower, locDistance))
					dHist_BCALTrackDOCA->Fill(locDistance);

				dHist_BCALNeutralShowerTime->Fill(locNeutralShowers[loc_i]->dSpacetimeVertex.T());
				dHist_BCALNeutralShowerEnergy->Fill(locNeutralShowers[loc_i]->dEnergy);

				dHist_BCALNeutralShowerDeltaT->Fill(locDeltaT);
				dHist_BCALNeutralShowerDeltaTVsE->Fill(locNeutralShowers[loc_i]->dEnergy, locDeltaT);
				dHist_BCALNeutralShowerDeltaTVsZ->Fill(locNeutralShowers[loc_i]->dSpacetimeVertex.Z(), locDeltaT);
			}
		}

		for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
		{
			if(!locDetectorMatches->Get_IsMatchedToHit(locTrackTimeBasedVector[loc_i]))
				dHist_TrackPVsTheta_NoHitMatch->Fill(locTrackTimeBasedVector[loc_i]->momentum().Theta()*180.0/TMath::Pi(), locTrackTimeBasedVector[loc_i]->momentum().Mag());
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

void DHistogramAction_DetectorStudies::Fill_MatchedHists(JEventLoop* locEventLoop, bool locUseTruePIDFlag)
{
	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks, "PreSelect");

	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
	locEventLoop->Get(locMCThrownMatchingVector);

	vector<const DEventRFBunch*> locEventRFBunches;
	locEventLoop->Get(locEventRFBunches);
	double locEventStartTime = locEventRFBunches.empty() ? 0.0 : locEventRFBunches[0]->dTime;

	//Fill Histograms
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		for(size_t loc_i = 0; loc_i < locChargedTracks.size(); ++loc_i)
		{
			const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTracks[loc_i]->Get_BestFOM();
			double locStartTime = locEventStartTime + (locChargedTrackHypothesis->position().Z() - dTargetCenter.Z())/29.9792458;

			if(locUseTruePIDFlag && (!locMCThrownMatchingVector.empty()))
			{
				double locMatchFOM = 0.0;
				const DMCThrown* locMCThrown = locMCThrownMatchingVector[0]->Get_MatchingMCThrown(locChargedTrackHypothesis, locMatchFOM);
				if((locMCThrown == NULL) || (locMatchFOM < dMinThrownMatchFOM))
					continue;
				//OK, have the thrown. Now, grab the best charged track hypothesis to get the best matching
				locChargedTrackHypothesis = locMCThrownMatchingVector[0]->Get_MatchingChargedHypothesis(locMCThrown, locMatchFOM);
			}

			Particle_t locPID = locChargedTrackHypothesis->PID();
			pair<Particle_t, bool> locPIDPair(locPID, locUseTruePIDFlag);
			if(dHistMap_BCALTrackDOCA.find(locPIDPair) == dHistMap_BCALTrackDOCA.end())
				continue; //disregard PID

			DVector3 locMomentum = locChargedTrackHypothesis->momentum();
			const DShowerMatchParams& locFCALShowerMatchParams = locChargedTrackHypothesis->dFCALShowerMatchParams;
			const DTOFHitMatchParams& locTOFHitMatchParams = locChargedTrackHypothesis->dTOFHitMatchParams;
			const DSCHitMatchParams& locSCHitMatchParams = locChargedTrackHypothesis->dSCHitMatchParams;
			const DShowerMatchParams& locBCALShowerMatchParams = locChargedTrackHypothesis->dBCALShowerMatchParams;

			//BCAL
			if(locBCALShowerMatchParams.dTrackTimeBased != NULL)
			{
				const DBCALShower* locBCALShower = dynamic_cast<const DBCALShower*>(locBCALShowerMatchParams.dShowerObject);
				dHistMap_BCALTrackDOCA[locPIDPair]->Fill(locBCALShowerMatchParams.dDOCAToShower);
				dHistMap_BCALShowerEnergy[locPIDPair]->Fill(locBCALShower->E);
				dHistMap_BCALShowerTrackDepth[locPIDPair]->Fill(locBCALShowerMatchParams.dx);
				dHistMap_BCALShowerTrackDepthVsP[locPIDPair]->Fill(locMomentum.Mag(), locBCALShowerMatchParams.dx);

				double locDeltaT = locBCALShower->t - locBCALShowerMatchParams.dFlightTime - locStartTime;
				dHistMap_BCALShowerDeltaT[locPIDPair]->Fill(locDeltaT);
				dHistMap_BCALShowerDeltaTVsZ[locPIDPair]->Fill(locBCALShower->z, locDeltaT);
				dHistMap_BCALShowerDeltaTVsP[locPIDPair]->Fill(locMomentum.Mag(), locDeltaT);
			}

			//FCAL
			if(locFCALShowerMatchParams.dTrackTimeBased != NULL)
			{
				const DFCALShower* locFCALShower = dynamic_cast<const DFCALShower*>(locFCALShowerMatchParams.dShowerObject);
				dHistMap_FCALTrackDOCA[locPIDPair]->Fill(locFCALShowerMatchParams.dDOCAToShower);
				dHistMap_FCALShowerEnergy[locPIDPair]->Fill(locFCALShower->getEnergy());
				dHistMap_FCALShowerTrackDepth[locPIDPair]->Fill(locFCALShowerMatchParams.dx);
				dHistMap_FCALShowerTrackDepthVsP[locPIDPair]->Fill(locMomentum.Mag(), locFCALShowerMatchParams.dx);

				double locDeltaT = locFCALShower->getTime() - locFCALShowerMatchParams.dFlightTime - locStartTime;
				dHistMap_FCALShowerDeltaT[locPIDPair]->Fill(locDeltaT);
				dHistMap_FCALShowerDeltaTVsP[locPIDPair]->Fill(locMomentum.Mag(), locDeltaT);
			}

			//TOF
			if(locTOFHitMatchParams.dTrackTimeBased != NULL)
			{
				dHistMap_TOFdEdX[locPIDPair]->Fill(locTOFHitMatchParams.dEdx*1.0E3);
				dHistMap_TOFdEdXVsP[locPIDPair]->Fill(locMomentum.Mag(), locTOFHitMatchParams.dEdx*1.0E3);
				dHistMap_TOFTrackDOCA[locPIDPair]->Fill(locTOFHitMatchParams.dDOCAToHit);

				double locDeltaT = locTOFHitMatchParams.dTOFPoint->t - locTOFHitMatchParams.dFlightTime - locStartTime;
				dHistMap_TOFDeltaT[locPIDPair]->Fill(locDeltaT);
				dHistMap_TOFDeltaTVsP[locPIDPair]->Fill(locMomentum.Mag(), locDeltaT);
			}

			//SC
			if(locSCHitMatchParams.dTrackTimeBased != NULL)
			{
				dHistMap_SCdEdX[locPIDPair]->Fill(locSCHitMatchParams.dEdx*1.0E3);
				dHistMap_SCdEdXVsP[locPIDPair]->Fill(locMomentum.Mag(), locSCHitMatchParams.dEdx*1.0E3);
				dHistMap_SCTrackDeltaPhi[locPIDPair]->Fill(locSCHitMatchParams.dDeltaPhiToHit*180.0/TMath::Pi());

				dHistMap_SCTime[locPIDPair]->Fill(locSCHitMatchParams.dHitTime);
				dHistMap_SCTimeVsTheta[locPIDPair]->Fill(locMomentum.Theta()*180.0/TMath::Pi(), locSCHitMatchParams.dHitTime);
				dHistMap_SCEnergyVsTheta[locPIDPair]->Fill(locMomentum.Theta()*180.0/TMath::Pi(), locSCHitMatchParams.dHitEnergy*1.0E3);

				dHistMap_SCPhiVsTheta[locPIDPair]->Fill(locMomentum.Theta()*180.0/TMath::Pi(), locMomentum.Phi()*180.0/TMath::Pi());

				double locDeltaT = locSCHitMatchParams.dHitTime - locSCHitMatchParams.dFlightTime - locStartTime;
				dHistMap_SCDeltaT[locPIDPair]->Fill(locDeltaT);
				dHistMap_SCDeltaTVsP[locPIDPair]->Fill(locMomentum.Mag(), locDeltaT);
				dHistMap_SCDeltaTVsPhi[locPIDPair]->Fill(locMomentum.Phi()*180.0/TMath::Pi(), locDeltaT);
				dHistMap_SCDeltaTVsTheta[locPIDPair]->Fill(locMomentum.Theta()*180.0/TMath::Pi(), locDeltaT);
			}
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

void DHistogramAction_TrackVertexComparison::Initialize(JEventLoop* locEventLoop)
{
	deque<deque<Particle_t> > locDetectedChargedPIDs;
	Get_Reaction()->Get_DetectedFinalPIDs(locDetectedChargedPIDs, 1);

	deque<deque<Particle_t> > locDetectedChargedPIDs_HasDupes;
	Get_Reaction()->Get_DetectedFinalPIDs(locDetectedChargedPIDs_HasDupes, 1, true);

	string locHistName, locHistTitle, locStepName, locStepROOTName, locParticleName, locParticleROOTName;
	Particle_t locPID, locHigherMassPID, locLowerMassPID;
	string locHigherMassParticleName, locLowerMassParticleName, locHigherMassParticleROOTName, locLowerMassParticleROOTName;

	vector<const DAnalysisUtilities*> locAnalysisUtilitiesVector;
	locEventLoop->Get(locAnalysisUtilitiesVector);

	size_t locNumSteps = Get_Reaction()->Get_NumReactionSteps();
	dHistDeque_TrackZToCommon.resize(locNumSteps);
	dHistDeque_TrackTToCommon.resize(locNumSteps);
	dHistDeque_TrackDOCAToCommon.resize(locNumSteps);
	dHistDeque_MaxTrackDeltaZ.resize(locNumSteps);
	dHistDeque_MaxTrackDeltaT.resize(locNumSteps);
	dHistDeque_MaxTrackDOCA.resize(locNumSteps);
	dHistDeque_TrackDeltaTVsP.resize(locNumSteps);

	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	dAnalysisUtilities = locAnalysisUtilitiesVector[0];
	CreateAndChangeTo_ActionDirectory();
	for(size_t loc_i = 0; loc_i < locNumSteps; ++loc_i)
	{
		if(locDetectedChargedPIDs[loc_i].empty())
			continue;

		const DReactionStep* locReactionStep = Get_Reaction()->Get_ReactionStep(loc_i);
		locStepName = locReactionStep->Get_StepName();
		locStepROOTName = locReactionStep->Get_StepROOTName();
		CreateAndChangeTo_Directory(locStepName, locStepName);

		// Max Track DeltaZ
		locHistName = "MaxTrackDeltaZ";
		locHistTitle = locStepROOTName + string(";Largest Track #DeltaVertex-Z (cm)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistDeque_MaxTrackDeltaZ[loc_i] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistDeque_MaxTrackDeltaZ[loc_i] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaVertexZBins, 0.0, dMaxDeltaVertexZ);

		// Max Track DeltaT
		locHistName = "MaxTrackDeltaT";
		locHistTitle = locStepROOTName + string(";Largest Track #DeltaVertex-T (ns)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistDeque_MaxTrackDeltaT[loc_i] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistDeque_MaxTrackDeltaT[loc_i] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaVertexTBins, 0.0, dMaxDeltaVertexT);

		// Max Track DOCA
		locHistName = "MaxTrackDOCA";
		locHistTitle = locStepROOTName + string(";Largest Track DOCA (cm)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistDeque_MaxTrackDOCA[loc_i] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistDeque_MaxTrackDOCA[loc_i] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDOCABins, 0.0, dMaxDOCA);

		for(size_t loc_j = 0; loc_j < locDetectedChargedPIDs[loc_i].size(); ++loc_j)
		{
			locPID = locDetectedChargedPIDs[loc_i][loc_j];
			if(dHistDeque_TrackZToCommon[loc_i].find(locPID) != dHistDeque_TrackZToCommon[loc_i].end())
				continue; //already created for this pid
			locParticleName = ParticleType(locPID);
			locParticleROOTName = ParticleName_ROOT(locPID);

			// TrackZ To Common
			locHistName = string("TrackZToCommon_") + locParticleName;
			locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";#DeltaVertex-Z (Track, Common) (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistDeque_TrackZToCommon[loc_i][locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistDeque_TrackZToCommon[loc_i][locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaVertexZBins, dMinDeltaVertexZ, dMaxDeltaVertexZ);

			// TrackT To Common
			locHistName = string("TrackTToCommon_") + locParticleName;
			locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";#DeltaVertex-T (Track, Common) (ns)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistDeque_TrackTToCommon[loc_i][locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistDeque_TrackTToCommon[loc_i][locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaVertexTBins, dMinDeltaVertexT, dMaxDeltaVertexT);

			// TrackDOCA To Common
			locHistName = string("TrackDOCAToCommon_") + locParticleName;
			locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";DOCA (Track, Common) (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistDeque_TrackDOCAToCommon[loc_i][locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistDeque_TrackDOCAToCommon[loc_i][locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDOCABins, dMinDOCA, dMaxDOCA);

			// DeltaT Vs P against beam photon
			if((locReactionStep->Get_InitialParticleID() == Gamma) && (dHistMap_BeamTrackDeltaTVsP.find(locPID) == dHistMap_BeamTrackDeltaTVsP.end()))
			{
				locHistName = string("TrackDeltaTVsP_") + ParticleType(locPID) + string("_Beam") + ParticleType(Gamma);
				locHistTitle = locStepROOTName + string(";") + ParticleName_ROOT(locPID) + string(" Momentum (GeV/c);t_{") + ParticleName_ROOT(locPID) + string("} - t_{Beam ") + ParticleName_ROOT(Gamma) + string("} (ns)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_BeamTrackDeltaTVsP[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_BeamTrackDeltaTVsP[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaVertexTBins, dMinDeltaVertexT, dMaxDeltaVertexT);
			}
		}

		//delta-t vs p
		for(int loc_j = 0; loc_j < int(locDetectedChargedPIDs_HasDupes[loc_i].size()) - 1; ++loc_j)
		{
			locPID = locDetectedChargedPIDs_HasDupes[loc_i][loc_j];
			for(size_t loc_k = loc_j + 1; loc_k < locDetectedChargedPIDs_HasDupes[loc_i].size(); ++loc_k)
			{
				if(ParticleMass(locDetectedChargedPIDs_HasDupes[loc_i][loc_k]) > ParticleMass(locPID))
				{
					locHigherMassPID = locDetectedChargedPIDs_HasDupes[loc_i][loc_k];
					locLowerMassPID = locPID;
				}
				else
				{
					locHigherMassPID = locPID;
					locLowerMassPID = locDetectedChargedPIDs_HasDupes[loc_i][loc_k];
				}
				pair<Particle_t, Particle_t> locParticlePair(locHigherMassPID, locLowerMassPID);
				if(dHistDeque_TrackDeltaTVsP[loc_i].find(locParticlePair) != dHistDeque_TrackDeltaTVsP[loc_i].end())
					continue; //already created for this pair

				locHigherMassParticleName = ParticleType(locHigherMassPID);
				locHigherMassParticleROOTName = ParticleName_ROOT(locHigherMassPID);
				locLowerMassParticleName = ParticleType(locLowerMassPID);
				locLowerMassParticleROOTName = ParticleName_ROOT(locLowerMassPID);

				// DeltaT Vs P
				locHistName = string("TrackDeltaTVsP_") + locHigherMassParticleName + string("_") + locLowerMassParticleName;
				locHistTitle = locStepROOTName + string(";") + locHigherMassParticleROOTName + string(" Momentum (GeV/c);t_{") + locHigherMassParticleROOTName + string("} - t_{") + locLowerMassParticleROOTName + string("} (ns)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_TrackDeltaTVsP[loc_i][locParticlePair] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_TrackDeltaTVsP[loc_i][locParticlePair] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaVertexTBins, dMinDeltaVertexT, dMaxDeltaVertexT);
			}
		}
		gDirectory->cd("..");
	}

	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_TrackVertexComparison::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	Particle_t locPID;
	double locDOCA, locDeltaZ, locDeltaT, locMaxDOCA, locMaxDeltaZ, locMaxDeltaT;

	//note: very difficult to tell when results will be duplicate: just histogram all combos
	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);
		deque<const DKinematicData*> locParticles;
		locParticleComboStep->Get_DetectedFinalChargedParticles_Measured(locParticles);
		if(locParticles.empty())
			continue;

		//Grab/Find common vertex & time
		DVector3 locVertex = locParticleComboStep->Get_Position();
		double locVertexTime = locParticleComboStep->Get_Time();

		for(size_t loc_j = 0; loc_j < locParticles.size(); ++loc_j)
		{
			locPID = locParticles[loc_j]->PID();

			//find max's
			locMaxDOCA = -1.0;
			locMaxDeltaZ = -1.0;
			locMaxDeltaT = -1.0;
			for(size_t loc_k = loc_j + 1; loc_k < locParticles.size(); ++loc_k)
			{
				locDeltaZ = fabs(locParticles[loc_j]->position().Z() - locParticles[loc_k]->position().Z());
				if(locDeltaZ > locMaxDeltaZ)
					locMaxDeltaZ = locDeltaZ;

				locDeltaT = fabs(locParticles[loc_j]->time() - locParticles[loc_k]->time());
				if(locDeltaT > locMaxDeltaT)
					locMaxDeltaT = locDeltaT;

				locDOCA = dAnalysisUtilities->Calc_DOCA(locParticles[loc_j], locParticles[loc_k]);
				if(locDOCA > locMaxDOCA)
					locMaxDOCA = locDOCA;
			}

			//delta-t vs p
			deque<pair<const DKinematicData*, size_t> > locParticlePairs;
			size_t locHigherMassParticleIndex, locLowerMassParticleIndex;
			//minimize number of locks by keeping track of the results and saving them at the end
			deque<pair<Particle_t, Particle_t> > locPIDPairs;
			deque<double> locPs;
			deque<double> locDeltaTs;
			for(size_t loc_k = loc_j + 1; loc_k < locParticles.size(); ++loc_k)
			{
				locParticlePairs.clear();
				locParticlePairs.push_back(pair<const DKinematicData*, size_t>(locParticles[loc_j], loc_i));
				locParticlePairs.push_back(pair<const DKinematicData*, size_t>(locParticles[loc_k], loc_i));

				if(locParticles[loc_k]->mass() > locParticles[loc_j]->mass())
				{
					locHigherMassParticleIndex = loc_k;
					locLowerMassParticleIndex = loc_j;
				}
				else
				{
					locHigherMassParticleIndex = loc_j;
					locLowerMassParticleIndex = loc_k;
				}

				locDeltaTs.push_back(locParticles[locHigherMassParticleIndex]->time() - locParticles[locLowerMassParticleIndex]->time());
				locPs.push_back(locParticles[locHigherMassParticleIndex]->momentum().Mag());
				locPIDPairs.push_back(pair<Particle_t, Particle_t>(locParticles[locHigherMassParticleIndex]->PID(), locParticles[locLowerMassParticleIndex]->PID()));
			}

			const DKinematicData* locBeamParticle = (locParticleComboStep->Get_InitialParticleID() == Gamma) ? locParticleComboStep->Get_InitialParticle() : NULL;
			double locBeamDeltaT = (locBeamParticle != NULL) ? locParticles[loc_j]->time() - locBeamParticle->time() : numeric_limits<double>::quiet_NaN();
			locDOCA = dAnalysisUtilities->Calc_DOCAToVertex(locParticles[loc_j], locVertex);

			//HISTOGRAM //do all at once to reduce #locks & amount of time within the lock
			japp->RootWriteLock();
			{
				//comparison to common vertex/time
				dHistDeque_TrackZToCommon[loc_i][locPID]->Fill(locParticles[loc_j]->position().Z() - locVertex.Z());
				dHistDeque_TrackTToCommon[loc_i][locPID]->Fill(locParticles[loc_j]->time() - locVertexTime);
				dHistDeque_TrackDOCAToCommon[loc_i][locPID]->Fill(locDOCA);
				//hist max's
				if(locMaxDeltaZ > 0.0) //else none found (e.g. only 1 detected charged track)
				{
					dHistDeque_MaxTrackDeltaZ[loc_i]->Fill(locMaxDeltaZ);
					dHistDeque_MaxTrackDeltaT[loc_i]->Fill(locMaxDeltaT);
					dHistDeque_MaxTrackDOCA[loc_i]->Fill(locMaxDOCA);
				}
				//delta-t's
				if(locBeamParticle != NULL)
					dHistMap_BeamTrackDeltaTVsP[locPID]->Fill(locParticles[loc_j]->momentum().Mag(), locBeamDeltaT);
				for(size_t loc_k = 0; loc_k < locPIDPairs.size(); ++loc_k)
				{
					if(dHistDeque_TrackDeltaTVsP[loc_i].find(locPIDPairs[loc_k]) == dHistDeque_TrackDeltaTVsP[loc_i].end())
					{
						//pair not found: equal masses and order switched somehow //e.g. mass set differently between REST and reconstruction
						pair<Particle_t, Particle_t> locTempPIDPair(locPIDPairs[loc_k]);
						locPIDPairs[loc_k].first = locTempPIDPair.second;
						locPIDPairs[loc_k].second = locTempPIDPair.first;
					}
					dHistDeque_TrackDeltaTVsP[loc_i][locPIDPairs[loc_k]]->Fill(locPs[loc_k], locDeltaTs[loc_k]);
				}
			}
			japp->RootUnLock();
		} //end of particle loop
	} //end of step loop
	return true;
}

void DHistogramAction_ParticleComboKinematics::Initialize(JEventLoop* locEventLoop)
{
	if(Get_UseKinFitResultsFlag() && (Get_Reaction()->Get_KinFitType() == d_NoFit))
	{
		cout << "WARNING: REQUESTED HISTOGRAM OF KINEMAITIC FIT RESULTS WHEN NO KINEMATIC FIT!!!" << endl;
		return; //no fit performed, but kinfit data requested!!
	}

	vector<const DParticleID*> locParticleIDs;
	locEventLoop->Get(locParticleIDs);

	string locHistName, locHistTitle, locStepName, locStepROOTName, locParticleName, locParticleROOTName;
	Particle_t locPID;

	size_t locNumSteps = Get_Reaction()->Get_NumReactionSteps();
	dHistDeque_PVsTheta.resize(locNumSteps);
	dHistDeque_PhiVsTheta.resize(locNumSteps);
	dHistDeque_BetaVsP.resize(locNumSteps);
	dHistDeque_DeltaBetaVsP.resize(locNumSteps);
	dHistDeque_P.resize(locNumSteps);
	dHistDeque_Theta.resize(locNumSteps);
	dHistDeque_Phi.resize(locNumSteps);
	dHistDeque_VertexZ.resize(locNumSteps);
	dHistDeque_VertexT.resize(locNumSteps);
	dHistDeque_VertexYVsX.resize(locNumSteps);

	deque<deque<Particle_t> > locDetectedFinalPIDs;
	Get_Reaction()->Get_DetectedFinalPIDs(locDetectedFinalPIDs);

	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry *locGeometry = locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());
	vector<const DAnalysisUtilities*> locAnalysisUtilitiesVector;
	locEventLoop->Get(locAnalysisUtilitiesVector);

	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		CreateAndChangeTo_ActionDirectory();

		dParticleID = locParticleIDs[0];
		dAnalysisUtilities = locAnalysisUtilitiesVector[0];
		locGeometry->GetTargetZ(dTargetZCenter);

		//beam particle
		locPID = Get_Reaction()->Get_ReactionStep(0)->Get_InitialParticleID();
		bool locBeamParticleUsed = (locPID == Gamma);
		if(locBeamParticleUsed)
		{
			locParticleName = string("Beam_") + ParticleType(locPID);
			CreateAndChangeTo_Directory(locParticleName, locParticleName);
			locParticleROOTName = ParticleName_ROOT(locPID);

			// Momentum
			locHistName = "Momentum";
			locHistTitle = string("Beam ") + locParticleROOTName + string(";p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dBeamParticleHist_P = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dBeamParticleHist_P = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPBins, dMinP, dMaxP);

			// Theta
			locHistName = "Theta";
			locHistTitle = string("Beam ") + locParticleROOTName + string(";#theta#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dBeamParticleHist_Theta = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dBeamParticleHist_Theta = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumThetaBins, dMinTheta, dMaxTheta);

			// Phi
			locHistName = "Phi";
			locHistTitle = string("Beam ") + locParticleROOTName + string(";#phi#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dBeamParticleHist_Phi = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dBeamParticleHist_Phi = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPhiBins, dMinPhi, dMaxPhi);

			// P Vs Theta
			locHistName = "PVsTheta";
			locHistTitle = string("Beam ") + locParticleROOTName + string(";#theta#circ;p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dBeamParticleHist_PVsTheta = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dBeamParticleHist_PVsTheta = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			// Phi Vs Theta
			locHistName = "PhiVsTheta";
			locHistTitle = string("Beam ") + locParticleROOTName + string(";#theta#circ;#phi#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dBeamParticleHist_PhiVsTheta = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dBeamParticleHist_PhiVsTheta = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

			// Vertex-Z
			locHistName = "VertexZ";
			locHistTitle = string("Beam ") + locParticleROOTName + string(";Vertex-Z (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dBeamParticleHist_VertexZ = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dBeamParticleHist_VertexZ = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumVertexZBins, dMinVertexZ, dMaxVertexZ);

			// Vertex-Y Vs Vertex-X
			locHistName = "VertexYVsX";
			locHistTitle = string("Beam ") + locParticleROOTName + string(";Vertex-X (cm);Vertex-Y (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dBeamParticleHist_VertexYVsX = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dBeamParticleHist_VertexYVsX = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNumVertexXYBins, dMinVertexXY, dMaxVertexXY, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY);

			// Vertex-T
			locHistName = "VertexT";
			locHistTitle = string("Beam ") + locParticleROOTName + string(";Vertex-T (ns)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dBeamParticleHist_VertexT = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dBeamParticleHist_VertexT = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumTBins, dMinT, dMaxT);

			// Delta-T (Beam, RF)
			locHistName = "DeltaTRF";
			locHistTitle = string("Beam ") + locParticleROOTName + string(";#Deltat_{Beam - RF} (ns)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dBeamParticleHist_DeltaTRF = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dBeamParticleHist_DeltaTRF = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTRFBins, dMinDeltaTRF, dMaxDeltaTRF);

			// Delta-T (Beam, RF) Vs Beam E
			locHistName = "DeltaTRFVsBeamE";
			locHistTitle = string("Beam ") + locParticleROOTName + string(";E (GeV);#Deltat_{Beam - RF} (ns)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dBeamParticleHist_DeltaTRFVsBeamE = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dBeamParticleHist_DeltaTRFVsBeamE = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaTRFBins, dMinDeltaTRF, dMaxDeltaTRF);

			gDirectory->cd("..");
		}

		//Steps
		deque<Particle_t> locPIDs;
		for(size_t loc_i = 0; loc_i < locNumSteps; ++loc_i)
		{
			const DReactionStep* locReactionStep = Get_Reaction()->Get_ReactionStep(loc_i);
			locStepName = locReactionStep->Get_StepName();
			locStepROOTName = locReactionStep->Get_StepROOTName();

			Particle_t locInitialPID = locReactionStep->Get_InitialParticleID();

			//get PIDs
			if(!Get_UseKinFitResultsFlag()) //measured, ignore missing & decaying particles (ignore target anyway)
				locPIDs = locDetectedFinalPIDs[loc_i];
			else //kinematic fit: decaying & missing particles are reconstructed
			{
				locReactionStep->Get_FinalParticleIDs(locPIDs);
				if((!locBeamParticleUsed) || (loc_i != 0)) //add decaying parent particle //skip if on beam particle!
					locPIDs.insert(locPIDs.begin(), locInitialPID);
			}

			if(locPIDs.empty())
				continue;

			CreateAndChangeTo_Directory(locStepName, locStepName);
			for(size_t loc_j = 0; loc_j < locPIDs.size(); ++loc_j)
			{
				locPID = locPIDs[loc_j];
				locParticleName = ParticleType(locPID);
				locParticleROOTName = ParticleName_ROOT(locPID);
				if(dHistDeque_P[loc_i].find(locPID) != dHistDeque_P[loc_i].end())
					continue; //pid already done

				CreateAndChangeTo_Directory(locParticleName, locParticleName);

				// Momentum
				locHistName = "Momentum";
				locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";p (GeV/c)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_P[loc_i][locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_P[loc_i][locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPBins, dMinP, dMaxP);

				// Theta
				locHistName = "Theta";
				locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";#theta#circ");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_Theta[loc_i][locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_Theta[loc_i][locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumThetaBins, dMinTheta, dMaxTheta);

				// Phi
				locHistName = "Phi";
				locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";#phi#circ");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_Phi[loc_i][locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_Phi[loc_i][locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPhiBins, dMinPhi, dMaxPhi);

				// P Vs Theta
				locHistName = "PVsTheta";
				locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";#theta#circ;p (GeV/c)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_PVsTheta[loc_i][locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_PVsTheta[loc_i][locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

				// Phi Vs Theta
				locHistName = "PhiVsTheta";
				locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";#theta#circ;#phi#circ");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_PhiVsTheta[loc_i][locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_PhiVsTheta[loc_i][locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

				//beta vs p
				locHistName = "BetaVsP";
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#beta");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_BetaVsP[loc_i][locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_BetaVsP[loc_i][locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumBetaBins, dMinBeta, dMaxBeta);

				//delta-beta vs p
				locHistName = "DeltaBetaVsP";
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta#beta");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaBetaVsP[loc_i][locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaBetaVsP[loc_i][locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaBetaBins, dMinDeltaBeta, dMaxDeltaBeta);

				// Vertex-Z
				locHistName = "VertexZ";
				locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";Vertex-Z (cm)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_VertexZ[loc_i][locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_VertexZ[loc_i][locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumVertexZBins, dMinVertexZ, dMaxVertexZ);

				// Vertex-Y Vs Vertex-X
				locHistName = "VertexYVsX";
				locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";Vertex-X (cm);Vertex-Y (cm)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_VertexYVsX[loc_i][locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_VertexYVsX[loc_i][locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNumVertexXYBins, dMinVertexXY, dMaxVertexXY, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY);

				// Vertex-T
				locHistName = "VertexT";
				locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";Vertex-T (ns)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_VertexT[loc_i][locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_VertexT[loc_i][locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumTBins, dMinT, dMaxT);

				gDirectory->cd("..");
			} //end of particle loop
			gDirectory->cd("..");
		} //end of step loop
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_ParticleComboKinematics::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	if(Get_UseKinFitResultsFlag() && (Get_Reaction()->Get_KinFitType() == d_NoFit))
	{
		cout << "WARNING: REQUESTED HISTOGRAM OF KINEMAITIC FIT RESULTS WHEN NO KINEMATIC FIT!!! Skipping histogram." << endl;
		return true; //no fit performed, but kinfit data requested!!
	}

	if(Get_NumPreviousParticleCombos() == 0)
	{
		dPreviouslyHistogrammedParticles.clear();
		dPreviouslyHistogrammedBeamParticles.clear();
	}

	const DEventRFBunch* locEventRFBunch = locParticleCombo->Get_EventRFBunch();

	const DKinematicData* locKinematicData;
	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);

		//initial particle
		if(Get_UseKinFitResultsFlag())
			locKinematicData = locParticleComboStep->Get_InitialParticle();
		else
			locKinematicData = locParticleComboStep->Get_InitialParticle_Measured();
		if(locKinematicData != NULL)
		{
			if(locKinematicData->PID() == Gamma)
			{
				//check if will be duplicate
				const JObject* locSourceObject = locParticleComboStep->Get_InitialParticle_Measured();
				if(dPreviouslyHistogrammedBeamParticles.find(locSourceObject) == dPreviouslyHistogrammedBeamParticles.end())
				{
					dPreviouslyHistogrammedBeamParticles.insert(locSourceObject);
					Fill_BeamHists(locKinematicData, locEventRFBunch);
				}
			}
			else if(Get_UseKinFitResultsFlag()) //decaying particle, but kinfit so can hist
				Fill_Hists(locEventLoop, locKinematicData, locEventRFBunch, loc_i); //has many source object, and is unique to this combo: no dupes to check against: let it ride
		}

		//final particles
		for(size_t loc_j = 0; loc_j < locParticleComboStep->Get_NumFinalParticles(); ++loc_j)
		{
			if(Get_UseKinFitResultsFlag())
				locKinematicData = locParticleComboStep->Get_FinalParticle(loc_j);
			else
				locKinematicData = locParticleComboStep->Get_FinalParticle_Measured(loc_j);
			if(locKinematicData == NULL)
				continue; //e.g. a decaying or missing particle whose params aren't set yet

			//check if duplicate
			const JObject* locSourceObject = locParticleComboStep->Get_FinalParticle_SourceObject(loc_j);
			if(locSourceObject != NULL) //else is reconstructed missing/decaying particle: has many source object, and is unique to this combo: no dupes to check against: let it ride
			{
				pair<Particle_t, const JObject*> locParticleInfo(locKinematicData->PID(), locSourceObject);
				pair<size_t, pair<Particle_t, const JObject*> > locHistInfo(loc_i, locParticleInfo);
				if(dPreviouslyHistogrammedParticles.find(locHistInfo) != dPreviouslyHistogrammedParticles.end())
					continue; //previously histogrammed
				dPreviouslyHistogrammedParticles.insert(locHistInfo);
			}

			Fill_Hists(locEventLoop, locKinematicData, locEventRFBunch, loc_i);
		} //end of particle loop
	} //end of step loop
	return true;
}

void DHistogramAction_ParticleComboKinematics::Fill_Hists(JEventLoop* locEventLoop, const DKinematicData* locKinematicData, const DEventRFBunch* locEventRFBunch, size_t locStepIndex)
{
	Particle_t locPID = locKinematicData->PID();
	DVector3 locMomentum = locKinematicData->momentum();
	DVector3 locPosition = locKinematicData->position();
	double locPhi = locMomentum.Phi()*180.0/TMath::Pi();
	double locTheta = locMomentum.Theta()*180.0/TMath::Pi();
	double locP = locMomentum.Mag();

	double locBeta_Timing, locDeltaBeta;
	if(ParticleCharge(locPID) == 0)
	{
		const DNeutralParticleHypothesis* locNeutralParticleHypothesis = static_cast<const DNeutralParticleHypothesis*>(locKinematicData);
		locBeta_Timing = dAnalysisUtilities->Calc_Beta_Timing(locNeutralParticleHypothesis, locEventRFBunch);
		locDeltaBeta = locNeutralParticleHypothesis->lorentzMomentum().Beta() - locBeta_Timing;
	}
	else
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locKinematicData);
		locBeta_Timing = dAnalysisUtilities->Calc_Beta_Timing(locChargedTrackHypothesis, locEventRFBunch, true);
		locDeltaBeta = locChargedTrackHypothesis->lorentzMomentum().Beta() - locBeta_Timing;
	}

	japp->RootWriteLock();
	{
		dHistDeque_P[locStepIndex][locPID]->Fill(locP);
		dHistDeque_Phi[locStepIndex][locPID]->Fill(locPhi);
		dHistDeque_Theta[locStepIndex][locPID]->Fill(locTheta);
		dHistDeque_PVsTheta[locStepIndex][locPID]->Fill(locTheta, locP);
		dHistDeque_PhiVsTheta[locStepIndex][locPID]->Fill(locTheta, locPhi);
		dHistDeque_BetaVsP[locStepIndex][locPID]->Fill(locP, locBeta_Timing);
		dHistDeque_DeltaBetaVsP[locStepIndex][locPID]->Fill(locP, locDeltaBeta);
		dHistDeque_VertexZ[locStepIndex][locPID]->Fill(locKinematicData->position().Z());
		dHistDeque_VertexYVsX[locStepIndex][locPID]->Fill(locKinematicData->position().X(), locKinematicData->position().Y());
		dHistDeque_VertexT[locStepIndex][locPID]->Fill(locKinematicData->time());
	}
	japp->RootUnLock();
}

void DHistogramAction_ParticleComboKinematics::Fill_BeamHists(const DKinematicData* locKinematicData, const DEventRFBunch* locEventRFBunch)
{
	DVector3 locMomentum = locKinematicData->momentum();
	double locPhi = locMomentum.Phi()*180.0/TMath::Pi();
	double locTheta = locMomentum.Theta()*180.0/TMath::Pi();
	double locP = locMomentum.Mag();
	double locDeltaTRF = locKinematicData->time() - (locEventRFBunch->dTime + (locKinematicData->z() - dTargetZCenter)/29.9792458);

	japp->RootWriteLock();
	{
		dBeamParticleHist_P->Fill(locP);
		dBeamParticleHist_Phi->Fill(locPhi);
		dBeamParticleHist_Theta->Fill(locTheta);
		dBeamParticleHist_PVsTheta->Fill(locTheta, locP);
		dBeamParticleHist_PhiVsTheta->Fill(locTheta, locPhi);
		dBeamParticleHist_VertexZ->Fill(locKinematicData->position().Z());
		dBeamParticleHist_VertexYVsX->Fill(locKinematicData->position().X(), locKinematicData->position().Y());
		dBeamParticleHist_VertexT->Fill(locKinematicData->time());
		dBeamParticleHist_DeltaTRF->Fill(locDeltaTRF);
		dBeamParticleHist_DeltaTRFVsBeamE->Fill(locKinematicData->energy(), locDeltaTRF);
	}
	japp->RootUnLock();
}

void DHistogramAction_ParticleComboGenReconComparison::Initialize(JEventLoop* locEventLoop)
{
	if(Get_UseKinFitResultsFlag() && (Get_Reaction()->Get_KinFitType() == d_NoFit))
	{
		cout << "WARNING: REQUESTED HISTOGRAM OF KINEMAITIC FIT RESULTS WHEN NO KINEMATIC FIT!!!" << endl;
		return; //no fit performed, but kinfit data requested!!
	}

	vector<const DParticleID*> locParticleIDs;
	locEventLoop->Get(locParticleIDs);

	string locHistName, locHistTitle, locStepName, locStepROOTName, locParticleName, locParticleROOTName;
	Particle_t locPID;

	size_t locNumSteps = Get_Reaction()->Get_NumReactionSteps();

	dHistDeque_DeltaPOverP.resize(locNumSteps);
	dHistDeque_DeltaTheta.resize(locNumSteps);
	dHistDeque_DeltaPhi.resize(locNumSteps);
	dHistDeque_DeltaT.resize(locNumSteps);
	dHistDeque_DeltaT_TOF.resize(locNumSteps);
	dHistDeque_DeltaT_BCAL.resize(locNumSteps);
	dHistDeque_DeltaT_FCAL.resize(locNumSteps);
	dHistDeque_DeltaVertexZ.resize(locNumSteps);
	dHistDeque_DeltaPOverPVsP.resize(locNumSteps);
	dHistDeque_DeltaPOverPVsTheta.resize(locNumSteps);
	dHistDeque_DeltaThetaVsP.resize(locNumSteps);
	dHistDeque_DeltaThetaVsTheta.resize(locNumSteps);
	dHistDeque_DeltaPhiVsP.resize(locNumSteps);
	dHistDeque_DeltaPhiVsTheta.resize(locNumSteps);
	dHistDeque_DeltaTVsTheta.resize(locNumSteps);
	dHistDeque_DeltaTVsP.resize(locNumSteps);
	dHistDeque_DeltaVertexZVsTheta.resize(locNumSteps);

	dHistDeque_Pulls.resize(locNumSteps);
	dHistDeque_PullsVsP.resize(locNumSteps);
	dHistDeque_PullsVsTheta.resize(locNumSteps);

	dHistDeque_TimePull_CDC.resize(locNumSteps);
	dHistDeque_TimePull_ST.resize(locNumSteps);
	dHistDeque_TimePull_BCAL.resize(locNumSteps);
	dHistDeque_TimePull_TOF.resize(locNumSteps);
	dHistDeque_TimePull_FCAL.resize(locNumSteps);

	dHistDeque_TimePullVsTheta_CDC.resize(locNumSteps);
	dHistDeque_TimePullVsTheta_BCAL.resize(locNumSteps);
	dHistDeque_TimePullVsTheta_ST.resize(locNumSteps);

	dHistDeque_TimePullVsP_CDC.resize(locNumSteps);
	dHistDeque_TimePullVsP_BCAL.resize(locNumSteps);
	dHistDeque_TimePullVsP_ST.resize(locNumSteps);
	dHistDeque_TimePullVsP_TOF.resize(locNumSteps);
	dHistDeque_TimePullVsP_FCAL.resize(locNumSteps);

	deque<deque<Particle_t> > locDetectedFinalPIDs;
	Get_Reaction()->Get_DetectedFinalPIDs(locDetectedFinalPIDs);

	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry *locGeometry = locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());

	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		CreateAndChangeTo_ActionDirectory();

		locGeometry->GetTargetZ(dTargetZCenter);

		//RF
		locHistName = "DeltaT_RFBeamBunch";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dRFBeamBunchDeltaT_Hist = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dRFBeamBunchDeltaT_Hist = new TH1I(locHistName.c_str(), ";RF #Deltat (Reconstructed - Thrown)", dNumRFDeltaTBins, dMinRFDeltaT, dMaxRFDeltaT);

		//beam particle
		locPID = Get_Reaction()->Get_ReactionStep(0)->Get_InitialParticleID();
		bool locBeamParticleUsed = (locPID == Gamma);
		if(locBeamParticleUsed)
		{
			locParticleName = string("Beam_") + ParticleType(locPID);
			CreateAndChangeTo_Directory(locParticleName, locParticleName);
			locParticleROOTName = ParticleName_ROOT(locPID);

			// DeltaP/P
			locHistName = string("DeltaPOverP");
			locHistTitle = locParticleROOTName + string(";#Deltap/p (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dBeamParticleHist_DeltaPOverP = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dBeamParticleHist_DeltaPOverP = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

			// DeltaP/P Vs P
			locHistName = string("DeltaPOverPVsP");
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltap/p (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dBeamParticleHist_DeltaPOverPVsP = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dBeamParticleHist_DeltaPOverPVsP = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

			// DeltaT
			locHistName = string("DeltaT");
			locHistTitle = locParticleROOTName + string(";#Deltat (ns) (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dBeamParticleHist_DeltaT = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dBeamParticleHist_DeltaT = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

			gDirectory->cd("..");
		}

		deque<string> locPullNames(8, "");
		locPullNames[0] = "E";  locPullNames[1] = "Px";  locPullNames[2] = "Py";  locPullNames[3] = "Pz";
		locPullNames[4] = "Xx";  locPullNames[5] = "Xy";  locPullNames[6] = "Xz";  locPullNames[7] = "T";

		deque<string> locPullTitles(8, "");
		locPullTitles[0] = "E";  locPullTitles[1] = "p_{x}";  locPullTitles[2] = "p_{y}";  locPullTitles[3] = "p_{z}";
		locPullTitles[4] = "x_{x}";  locPullTitles[5] = "x_{y}";  locPullTitles[6] = "x_{z}";  locPullTitles[7] = "t";

		//CREATE THE HISTOGRAMS
		deque<Particle_t> locPIDs;
		for(size_t loc_i = 0; loc_i < locNumSteps; ++loc_i)
		{
			const DReactionStep* locReactionStep = Get_Reaction()->Get_ReactionStep(loc_i);
			locStepName = locReactionStep->Get_StepName();
			locStepROOTName = locReactionStep->Get_StepROOTName();

			Particle_t locInitialPID = locReactionStep->Get_InitialParticleID();

			//get PIDs
			if(!Get_UseKinFitResultsFlag()) //measured, ignore missing & decaying particles (ignore target anyway)
				locPIDs = locDetectedFinalPIDs[loc_i];
			else //kinematic fit: decaying & missing particles are reconstructed
			{
				locReactionStep->Get_FinalParticleIDs(locPIDs);
				if((!locBeamParticleUsed) || (loc_i != 0)) //add decaying parent particle //skip if on beam particle!
					locPIDs.insert(locPIDs.begin(), locInitialPID);
			}

			if(locPIDs.empty())
				continue;

			CreateAndChangeTo_Directory(locStepName, locStepName);
			for(size_t loc_j = 0; loc_j < locPIDs.size(); ++loc_j)
			{
				locPID = locPIDs[loc_j];
				locParticleName = ParticleType(locPID);
				locParticleROOTName = ParticleName_ROOT(locPID);
				if(dHistDeque_DeltaPOverP[loc_i].find(locPID) != dHistDeque_DeltaPOverP[loc_i].end())
					continue; //pid already done

				CreateAndChangeTo_Directory(locParticleName, locParticleName);

				// DeltaP/P
				locHistName = string("DeltaPOverP");
				locHistTitle = locParticleROOTName + string(";#Deltap/p (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaPOverP[loc_i][locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaPOverP[loc_i][locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

				// DeltaTheta
				locHistName = string("DeltaTheta");
				locHistTitle = locParticleROOTName + string(";#Delta#theta#circ (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaTheta[loc_i][locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaTheta[loc_i][locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaThetaBins, dMinDeltaTheta, dMaxDeltaTheta);

				// DeltaPhi
				locHistName = string("DeltaPhi");
				locHistTitle = locParticleROOTName + string(";#Delta#phi#circ (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaPhi[loc_i][locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaPhi[loc_i][locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

				// DeltaT
				locHistName = string("DeltaT");
				locHistTitle = locParticleROOTName + string(";#Deltat (ns) (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaT[loc_i][locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaT[loc_i][locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

				// DeltaT - BCAL
				locHistName = string("DeltaT_BCAL");
				locHistTitle = locParticleROOTName + string(" in BCAL;#Deltat (ns) (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaT_BCAL[loc_i][locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaT_BCAL[loc_i][locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

				// DeltaT - TOF (charged only)
				if(ParticleCharge(locPID) != 0)
				{
					locHistName = string("DeltaT_TOF");
					locHistTitle = locParticleROOTName + string(" in TOF;#Deltat (ns) (Reconstructed - Thrown)");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistDeque_DeltaT_TOF[loc_i][locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistDeque_DeltaT_TOF[loc_i][locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
				}

				// DeltaT - FCAL (neutral only)
				if(ParticleCharge(locPID) == 0)
				{
					locHistName = string("DeltaT_FCAL");
					locHistTitle = locParticleROOTName + string(" in FCAL;#Deltat (ns) (Reconstructed - Thrown)");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistDeque_DeltaT_FCAL[loc_i][locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistDeque_DeltaT_FCAL[loc_i][locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
				}

				// DeltaVertexZ
				locHistName = string("DeltaVertexZ");
				locHistTitle = locParticleROOTName + string(";#DeltaVertex-Z (cm) (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaVertexZ[loc_i][locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaVertexZ[loc_i][locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaVertexZBins, dMinDeltaVertexZ, dMaxDeltaVertexZ);

				// DeltaP/P Vs P
				locHistName = string("DeltaPOverPVsP");
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltap/p (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaPOverPVsP[loc_i][locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaPOverPVsP[loc_i][locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

				// DeltaP/P Vs Theta
				locHistName = string("DeltaPOverPVsTheta");
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltap/p (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaPOverPVsTheta[loc_i][locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaPOverPVsTheta[loc_i][locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

				// DeltaTheta Vs P
				locHistName = string("DeltaThetaVsP");
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta#theta#circ (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaThetaVsP[loc_i][locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaThetaVsP[loc_i][locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaThetaBins, dMinDeltaTheta, dMaxDeltaTheta);

				// DeltaTheta Vs Theta
				locHistName = string("DeltaThetaVsTheta");
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Delta#theta#circ (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaThetaVsTheta[loc_i][locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaThetaVsTheta[loc_i][locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaThetaBins, dMinDeltaTheta, dMaxDeltaTheta);

				// DeltaPhi Vs P
				locHistName = string("DeltaPhiVsP");
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta#phi#circ (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaPhiVsP[loc_i][locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaPhiVsP[loc_i][locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

				// DeltaPhi Vs Theta
				locHistName = string("DeltaPhiVsTheta");
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Delta#phi#circ (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaPhiVsTheta[loc_i][locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaPhiVsTheta[loc_i][locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

				// DeltaT Vs Theta
				locHistName = string("DeltaTVsTheta");
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat (ns) (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaTVsTheta[loc_i][locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaTVsTheta[loc_i][locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

				// DeltaT Vs P
				locHistName = string("DeltaTVsP");
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat (ns) (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaTVsP[loc_i][locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaTVsP[loc_i][locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

				// DeltaVertexZ Vs Theta
				locHistName = string("DeltaVertexZVsTheta");
				locHistTitle = locParticleROOTName + string(";#theta#circ;#DeltaVertex-Z (cm) (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaVertexZVsTheta[loc_i][locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaVertexZVsTheta[loc_i][locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaVertexZBins, dMinDeltaVertexZ, dMaxDeltaVertexZ);

				/************************************************************************ Pulls ************************************************************************/

				CreateAndChangeTo_Directory("Pulls", "Pulls");
				for(size_t loc_j = 0; loc_j < dPullTypes.size(); ++loc_j)
				{
					if((ParticleCharge(locPID) != 0) && (dPullTypes[loc_j] == d_EPull))
						continue;
					if((ParticleCharge(locPID) == 0) && ((dPullTypes[loc_j] >= d_PxPull) && (dPullTypes[loc_j] <= d_PzPull)))
						continue;

					//Pull 1D
					locHistName = locPullNames[loc_j] + string("Pull");
					locHistTitle = locParticleROOTName + string(";#Delta") + locPullTitles[loc_j] + string("/#sigma_{") + locPullTitles[loc_j] + string("} (Reconstructed - Thrown)");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						(dHistDeque_Pulls[loc_i][locPID])[dPullTypes[loc_j]] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
					else
						(dHistDeque_Pulls[loc_i][locPID])[dPullTypes[loc_j]] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

					//Pull vs P
					locHistName = locPullNames[loc_j] + string("PullVsP");
					locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta") + locPullTitles[loc_j] + string("/#sigma_{") + locPullTitles[loc_j] + string("} (Reconstructed - Thrown)");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						(dHistDeque_PullsVsP[loc_i][locPID])[dPullTypes[loc_j]] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
					else
						(dHistDeque_PullsVsP[loc_i][locPID])[dPullTypes[loc_j]] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

					//Pull vs Theta
					locHistName = locPullNames[loc_j] + string("PullVsTheta");
					locHistTitle = locParticleROOTName + string(";#theta#circ;#Delta") + locPullTitles[loc_j] + string("/#sigma_{") + locPullTitles[loc_j] + string("} (Reconstructed - Thrown)");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						(dHistDeque_PullsVsTheta[loc_i][locPID])[dPullTypes[loc_j]] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
					else
						(dHistDeque_PullsVsTheta[loc_i][locPID])[dPullTypes[loc_j]] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);
				}

				//Delta-t Pulls - CDC & ST
				if(ParticleCharge(locPID) != 0)
				{
					//CDC
					locHistName = "TimePull_CDC";
					locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistDeque_TimePull_CDC[loc_i][locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistDeque_TimePull_CDC[loc_i][locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

					locHistName = "TimePullVsTheta_CDC";
					locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat/#sigma_{#Deltat}");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistDeque_TimePullVsTheta_CDC[loc_i][locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistDeque_TimePullVsTheta_CDC[loc_i][locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);

					locHistName = "TimePullVsP_CDC";
					locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistDeque_TimePullVsP_CDC[loc_i][locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistDeque_TimePullVsP_CDC[loc_i][locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

					//ST
					locHistName = "TimePull_ST";
					locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistDeque_TimePull_ST[loc_i][locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistDeque_TimePull_ST[loc_i][locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

					locHistName = "TimePullVsTheta_ST";
					locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat/#sigma_{#Deltat}");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistDeque_TimePullVsTheta_ST[loc_i][locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistDeque_TimePullVsTheta_ST[loc_i][locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);

					locHistName = "TimePullVsP_ST";
					locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistDeque_TimePullVsP_ST[loc_i][locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistDeque_TimePullVsP_ST[loc_i][locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);
				}

				//Delta-t Pulls - BCAL
				locHistName = "TimePull_BCAL";
				locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_TimePull_BCAL[loc_i][locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_TimePull_BCAL[loc_i][locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

				locHistName = "TimePullVsTheta_BCAL";
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_TimePullVsTheta_BCAL[loc_i][locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_TimePullVsTheta_BCAL[loc_i][locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);

				locHistName = "TimePullVsP_BCAL";
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_TimePullVsP_BCAL[loc_i][locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_TimePullVsP_BCAL[loc_i][locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

				//Delta-t Pulls - TOF
				if(ParticleCharge(locPID) != 0) //TOF
				{
					locHistName = "TimePull_TOF";
					locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistDeque_TimePull_TOF[loc_i][locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistDeque_TimePull_TOF[loc_i][locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

					locHistName = "TimePullVsP_TOF";
					locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistDeque_TimePullVsP_TOF[loc_i][locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistDeque_TimePullVsP_TOF[loc_i][locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);
				}

				//Delta-t Pulls - FCAL
				locHistName = "TimePull_FCAL";
				locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_TimePull_FCAL[loc_i][locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_TimePull_FCAL[loc_i][locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

				locHistName = "TimePullVsP_FCAL";
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_TimePullVsP_FCAL[loc_i][locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_TimePullVsP_FCAL[loc_i][locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

				gDirectory->cd("..");

				gDirectory->cd("..");
			} //end of particle loop
			gDirectory->cd("..");
		} //end of step loop
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_ParticleComboGenReconComparison::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	if(Get_UseKinFitResultsFlag() && (Get_Reaction()->Get_KinFitType() == d_NoFit))
	{
		cout << "WARNING: REQUESTED HISTOGRAM OF KINEMAITIC FIT RESULTS WHEN NO KINEMATIC FIT!!! Skipping histogram." << endl;
		return true; //no fit performed, but kinfit data requested!!
	}

	if(Get_NumPreviousParticleCombos() == 0)
	{
		dPreviouslyHistogrammedParticles.clear();
		dPreviouslyHistogrammedBeamParticles.clear();
	}

	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);
	if(locMCThrowns.empty())
		return true; //e.g. non-simulated event

	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons, "MCGEN");
	const DBeamPhoton* locThrownBeamPhoton = locBeamPhotons.empty() ? NULL : locBeamPhotons[0];

	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
	locEventLoop->Get(locMCThrownMatchingVector);
	if(locMCThrownMatchingVector.empty())
		return true;
	const DMCThrownMatching* locMCThrownMatching = locMCThrownMatchingVector[0];

	const DEventRFBunch* locEventRFBunch = locParticleCombo->Get_EventRFBunch();
	const DEventRFBunch* locThrownEventRFBunch = NULL;
	locEventLoop->GetSingle(locThrownEventRFBunch, "Thrown");

	//RF time difference
	double locRFTime = locEventRFBunch->dTime;
	double locRFDeltaT = locRFTime - locThrownEventRFBunch->dTime;
	japp->RootWriteLock();
	{
		dRFBeamBunchDeltaT_Hist->Fill(locRFDeltaT);
	}
	japp->RootUnLock();

	const DKinematicData* locKinematicData;
	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);

		//initial particle
		if(Get_UseKinFitResultsFlag())
			locKinematicData = locParticleComboStep->Get_InitialParticle();
		else
			locKinematicData = locParticleComboStep->Get_InitialParticle_Measured();
		if(locKinematicData != NULL)
		{
			if(locKinematicData->PID() == Gamma)
			{
				//check if will be duplicate
				const JObject* locSourceObject = locParticleComboStep->Get_InitialParticle_Measured();
				if(dPreviouslyHistogrammedBeamParticles.find(locSourceObject) == dPreviouslyHistogrammedBeamParticles.end())
				{
					dPreviouslyHistogrammedBeamParticles.insert(locSourceObject);
					Fill_BeamHists(locKinematicData, locThrownBeamPhoton);
				}
			}
		}

		//final particles
		for(size_t loc_j = 0; loc_j < locParticleComboStep->Get_NumFinalParticles(); ++loc_j)
		{
			if(Get_UseKinFitResultsFlag())
				locKinematicData = locParticleComboStep->Get_FinalParticle(loc_j);
			else
				locKinematicData = locParticleComboStep->Get_FinalParticle_Measured(loc_j);
			if(locKinematicData == NULL)
				continue; //e.g. a decaying or missing particle whose params aren't set yet

			//check if duplicate
			const JObject* locSourceObject = locParticleComboStep->Get_FinalParticle_SourceObject(loc_j);
			if(locSourceObject != NULL) //else is reconstructed missing/decaying particle: has many source object, and is unique to this combo: no dupes to check against: let it ride
			{
				pair<Particle_t, const JObject*> locParticleInfo(locKinematicData->PID(), locSourceObject);
				pair<size_t, pair<Particle_t, const JObject*> > locHistInfo(loc_i, locParticleInfo);
				if(dPreviouslyHistogrammedParticles.find(locHistInfo) != dPreviouslyHistogrammedParticles.end())
					continue; //previously histogrammed
				dPreviouslyHistogrammedParticles.insert(locHistInfo);
			}

			if(ParticleCharge(locKinematicData->PID()) != 0)
			{
				const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locKinematicData);
				double locMatchFOM = 0.0;
				const DMCThrown* locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locChargedTrackHypothesis, locMatchFOM);
				if(locMCThrown != NULL)
					Fill_ChargedHists(locChargedTrackHypothesis, locMCThrown, locThrownEventRFBunch, loc_i);
			}
			else
			{
				const DNeutralParticleHypothesis* locNeutralParticleHypothesis = static_cast<const DNeutralParticleHypothesis*>(locKinematicData);
				double locMatchFOM = 0.0;
				const DMCThrown* locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locNeutralParticleHypothesis, locMatchFOM);
				if(locMCThrown != NULL)
					Fill_NeutralHists(locNeutralParticleHypothesis, locMCThrown, locThrownEventRFBunch, loc_i);
			}
		} //end of particle loop
	} //end of step loop
	return true;
}

void DHistogramAction_ParticleComboGenReconComparison::Fill_BeamHists(const DKinematicData* locKinematicData, const DKinematicData* locThrownKinematicData)
{
	if(locThrownKinematicData == NULL)
		return;

	DVector3 locMomentum = locKinematicData->momentum();
	DVector3 locThrownMomentum = locThrownKinematicData->momentum();

	double locThrownP = locThrownMomentum.Mag();
	double locDeltaPOverP = (locMomentum.Mag() - locThrownP)/locThrownP;
	double locDeltaT = locKinematicData->time() - locThrownKinematicData->time();

	japp->RootWriteLock();
	{
		dBeamParticleHist_DeltaPOverP->Fill(locDeltaPOverP);
		dBeamParticleHist_DeltaPOverPVsP->Fill(locThrownP, locDeltaPOverP);
		dBeamParticleHist_DeltaT->Fill(locDeltaT);
	}
	japp->RootUnLock();
}

void DHistogramAction_ParticleComboGenReconComparison::Fill_ChargedHists(const DChargedTrackHypothesis* locChargedTrackHypothesis, const DMCThrown* locMCThrown, const DEventRFBunch* locThrownEventRFBunch, size_t locStepIndex)
{
	Particle_t locPID = locChargedTrackHypothesis->PID();
	double locThrownP = locMCThrown->momentum().Mag();
	double locThrownTheta = locMCThrown->momentum().Theta()*180.0/TMath::Pi();
	double locDeltaPOverP = (locChargedTrackHypothesis->momentum().Mag() - locThrownP)/locThrownP;
	double locDeltaTheta = locChargedTrackHypothesis->momentum().Theta()*180.0/TMath::Pi() - locThrownTheta;
	double locDeltaPhi = locChargedTrackHypothesis->momentum().Phi()*180.0/TMath::Pi() - locMCThrown->momentum().Phi()*180.0/TMath::Pi();
	double locDeltaT = locChargedTrackHypothesis->time() - locMCThrown->time(); //time comparison isn't fair if track comes from a detached vertex!!!
	double locDeltaVertexZ = locChargedTrackHypothesis->position().Z() - locMCThrown->position().Z();
	const TMatrixDSym& locCovarianceMatrix = locChargedTrackHypothesis->errorMatrix();

	double locStartTime = locThrownEventRFBunch->dTime + (locMCThrown->z() - dTargetZCenter)/29.9792458;
	double locTimePull = (locStartTime - locChargedTrackHypothesis->time())/sqrt(locCovarianceMatrix(6, 6));
	double locT0Pull = (locStartTime - locChargedTrackHypothesis->t0())/locChargedTrackHypothesis->t0_err();
	japp->RootWriteLock();
	{
		dHistDeque_DeltaPOverP[locStepIndex][locPID]->Fill(locDeltaPOverP);
		dHistDeque_DeltaTheta[locStepIndex][locPID]->Fill(locDeltaTheta);
		dHistDeque_DeltaPhi[locStepIndex][locPID]->Fill(locDeltaPhi);
		dHistDeque_DeltaT[locStepIndex][locPID]->Fill(locDeltaT);
		if(locChargedTrackHypothesis->t0_detector() == SYS_START)
		{
			dHistDeque_TimePull_ST[locStepIndex][locPID]->Fill(locT0Pull);
			dHistDeque_TimePullVsTheta_ST[locStepIndex][locPID]->Fill(locThrownTheta, locT0Pull);
			dHistDeque_TimePullVsP_ST[locStepIndex][locPID]->Fill(locThrownP, locT0Pull);
		}
		if(locChargedTrackHypothesis->t0_detector() == SYS_CDC)
		{
			dHistDeque_TimePull_CDC[locStepIndex][locPID]->Fill(locT0Pull);
			dHistDeque_TimePullVsTheta_CDC[locStepIndex][locPID]->Fill(locThrownTheta, locT0Pull);
			dHistDeque_TimePullVsP_CDC[locStepIndex][locPID]->Fill(locThrownP, locT0Pull);
		}
		else if(locChargedTrackHypothesis->t1_detector() == SYS_CDC)
		{
			dHistDeque_TimePull_CDC[locStepIndex][locPID]->Fill(locTimePull);
			dHistDeque_TimePullVsTheta_CDC[locStepIndex][locPID]->Fill(locThrownTheta, locTimePull);
			dHistDeque_TimePullVsP_CDC[locStepIndex][locPID]->Fill(locThrownP, locTimePull);
		}
		if(locChargedTrackHypothesis->t1_detector() == SYS_BCAL)
		{
			dHistDeque_DeltaT_BCAL[locStepIndex][locPID]->Fill(locDeltaT);
			dHistDeque_TimePull_BCAL[locStepIndex][locPID]->Fill(locTimePull);
			dHistDeque_TimePullVsTheta_BCAL[locStepIndex][locPID]->Fill(locThrownTheta, locTimePull);
			dHistDeque_TimePullVsP_BCAL[locStepIndex][locPID]->Fill(locThrownP, locTimePull);
		}
		else if(locChargedTrackHypothesis->t1_detector() == SYS_TOF)
		{
			dHistDeque_DeltaT_TOF[locStepIndex][locPID]->Fill(locDeltaT);
			dHistDeque_TimePull_TOF[locStepIndex][locPID]->Fill(locTimePull);
			dHistDeque_TimePullVsP_TOF[locStepIndex][locPID]->Fill(locThrownP, locTimePull);
		}
		else if(locChargedTrackHypothesis->t1_detector() == SYS_FCAL)
		{
			dHistDeque_TimePull_FCAL[locStepIndex][locPID]->Fill(locTimePull);
			dHistDeque_TimePullVsP_FCAL[locStepIndex][locPID]->Fill(locThrownP, locTimePull);
		}
		dHistDeque_DeltaVertexZ[locStepIndex][locPID]->Fill(locDeltaVertexZ);
		dHistDeque_DeltaPOverPVsP[locStepIndex][locPID]->Fill(locThrownP, locDeltaPOverP);
		dHistDeque_DeltaPOverPVsTheta[locStepIndex][locPID]->Fill(locThrownTheta, locDeltaPOverP);
		dHistDeque_DeltaThetaVsP[locStepIndex][locPID]->Fill(locThrownP, locDeltaTheta);
		dHistDeque_DeltaThetaVsTheta[locStepIndex][locPID]->Fill(locThrownTheta, locDeltaTheta);
		dHistDeque_DeltaPhiVsP[locStepIndex][locPID]->Fill(locThrownP, locDeltaPhi);
		dHistDeque_DeltaPhiVsTheta[locStepIndex][locPID]->Fill(locThrownTheta, locDeltaPhi);
		dHistDeque_DeltaTVsTheta[locStepIndex][locPID]->Fill(locThrownTheta, locDeltaT);
		dHistDeque_DeltaTVsP[locStepIndex][locPID]->Fill(locThrownP, locDeltaT);
		dHistDeque_DeltaVertexZVsTheta[locStepIndex][locPID]->Fill(locThrownTheta, locDeltaVertexZ);

		for(size_t loc_j = 0; loc_j < dPullTypes.size(); ++loc_j)
		{
			if(dPullTypes[loc_j] == d_EPull)
				continue;
			double locPull = 0.0;
			if((dPullTypes[loc_j] >= d_PxPull) && (dPullTypes[loc_j] <= d_PzPull))
			{
				int locIndex = int(dPullTypes[loc_j] - d_PxPull);
				locPull = (locChargedTrackHypothesis->momentum()(locIndex) - locMCThrown->momentum()(locIndex))/sqrt(locCovarianceMatrix(locIndex, locIndex));
			}
			else if((dPullTypes[loc_j] >= d_XxPull) && (dPullTypes[loc_j] <= d_XzPull))
			{
				int locIndex = int(dPullTypes[loc_j] - d_XxPull);
				locPull = (locChargedTrackHypothesis->position()(locIndex) - locMCThrown->position()(locIndex))/sqrt(locCovarianceMatrix(locIndex + 3, locIndex + 3));
			}
			else if(dPullTypes[loc_j] == d_TPull)
				locPull = (locChargedTrackHypothesis->time() - locMCThrown->time())/sqrt(locCovarianceMatrix(6, 6));
			(dHistDeque_Pulls[locStepIndex][locPID])[dPullTypes[loc_j]]->Fill(locPull);
			(dHistDeque_PullsVsP[locStepIndex][locPID])[dPullTypes[loc_j]]->Fill(locThrownP, locPull);
			(dHistDeque_PullsVsTheta[locStepIndex][locPID])[dPullTypes[loc_j]]->Fill(locThrownTheta, locPull);
		}
	}
	japp->RootUnLock();
}

void DHistogramAction_ParticleComboGenReconComparison::Fill_NeutralHists(const DNeutralParticleHypothesis* locNeutralParticleHypothesis, const DMCThrown* locMCThrown, const DEventRFBunch* locThrownEventRFBunch, size_t locStepIndex)
{
	Particle_t locPID = locNeutralParticleHypothesis->PID();

	const DNeutralShower* locNeutralShower = NULL;
	locNeutralParticleHypothesis->GetSingle(locNeutralShower);
	if(locNeutralShower == NULL)
		return; //shouldn't be possible ...

	double locThrownP = locMCThrown->momentum().Mag();
	double locThrownTheta = locMCThrown->momentum().Theta()*180.0/TMath::Pi();
	double locDeltaPOverP = (locNeutralParticleHypothesis->momentum().Mag() - locThrownP)/locThrownP;
	double locDeltaTheta = locNeutralParticleHypothesis->momentum().Theta()*180.0/TMath::Pi() - locThrownTheta;
	double locDeltaPhi = locNeutralParticleHypothesis->momentum().Phi()*180.0/TMath::Pi() - locMCThrown->momentum().Phi()*180.0/TMath::Pi();
	double locDeltaT = locNeutralParticleHypothesis->time() - locMCThrown->time(); //time comparison isn't fair if track comes from a detached vertex!!!
	double locDeltaVertexZ = locNeutralParticleHypothesis->position().Z() - locMCThrown->position().Z();
	const TMatrixDSym& locCovarianceMatrix = locNeutralParticleHypothesis->errorMatrix();

	double locStartTime = locThrownEventRFBunch->dTime + (locMCThrown->z() - dTargetZCenter)/29.9792458;
	double locTimePull = (locStartTime - locNeutralParticleHypothesis->time())/sqrt(locCovarianceMatrix(6, 6));

	japp->RootWriteLock();
	{
		dHistDeque_DeltaPOverP[locStepIndex][locPID]->Fill(locDeltaPOverP);
		dHistDeque_DeltaTheta[locStepIndex][locPID]->Fill(locDeltaTheta);
		dHistDeque_DeltaPhi[locStepIndex][locPID]->Fill(locDeltaPhi);
		dHistDeque_DeltaT[locStepIndex][locPID]->Fill(locDeltaT);
		if(locNeutralParticleHypothesis->t1_detector() == SYS_BCAL)
		{
			dHistDeque_DeltaT_BCAL[locStepIndex][locPID]->Fill(locDeltaT);
			dHistDeque_TimePull_BCAL[locStepIndex][locPID]->Fill(locTimePull);
			dHistDeque_TimePullVsTheta_BCAL[locStepIndex][locPID]->Fill(locThrownTheta, locTimePull);
			dHistDeque_TimePullVsP_BCAL[locStepIndex][locPID]->Fill(locThrownP, locTimePull);
		}
		else if(locNeutralParticleHypothesis->t1_detector() == SYS_FCAL)
		{
			dHistDeque_DeltaT_FCAL[locStepIndex][locPID]->Fill(locDeltaT);
			dHistDeque_TimePull_FCAL[locStepIndex][locPID]->Fill(locTimePull);
			dHistDeque_TimePullVsP_FCAL[locStepIndex][locPID]->Fill(locThrownP, locTimePull);
		}

		dHistDeque_DeltaVertexZ[locStepIndex][locPID]->Fill(locDeltaVertexZ);
		dHistDeque_DeltaPOverPVsP[locStepIndex][locPID]->Fill(locThrownP, locDeltaPOverP);
		dHistDeque_DeltaPOverPVsTheta[locStepIndex][locPID]->Fill(locThrownTheta, locDeltaPOverP);
		dHistDeque_DeltaThetaVsP[locStepIndex][locPID]->Fill(locThrownP, locDeltaTheta);
		dHistDeque_DeltaThetaVsTheta[locStepIndex][locPID]->Fill(locThrownTheta, locDeltaTheta);
		dHistDeque_DeltaPhiVsP[locStepIndex][locPID]->Fill(locThrownP, locDeltaPhi);
		dHistDeque_DeltaPhiVsTheta[locStepIndex][locPID]->Fill(locThrownTheta, locDeltaPhi);
		dHistDeque_DeltaTVsTheta[locStepIndex][locPID]->Fill(locThrownTheta, locDeltaT);
		dHistDeque_DeltaTVsP[locStepIndex][locPID]->Fill(locThrownP, locDeltaT);
		dHistDeque_DeltaVertexZVsTheta[locStepIndex][locPID]->Fill(locThrownTheta, locDeltaVertexZ);

		for(size_t loc_j = 0; loc_j < dPullTypes.size(); ++loc_j)
		{
			if((dPullTypes[loc_j] >= d_PxPull) && (dPullTypes[loc_j] <= d_PzPull))
				continue;
			double locPull = 0.0;
			if(dPullTypes[loc_j] == d_EPull)
				locPull = (locNeutralShower->dEnergy - locMCThrown->energy())/sqrt(locNeutralShower->dCovarianceMatrix(0, 0));
			else if((dPullTypes[loc_j] >= d_XxPull) && (dPullTypes[loc_j] <= d_XzPull))
			{
				int locIndex = int(dPullTypes[loc_j] - d_XxPull);
				locPull = (locNeutralParticleHypothesis->position()(locIndex) - locMCThrown->position()(locIndex))/sqrt(locCovarianceMatrix(locIndex + 3, locIndex + 3));
			}
			else if(dPullTypes[loc_j] == d_TPull)
				locPull = (locNeutralParticleHypothesis->time() - locMCThrown->time())/sqrt(locCovarianceMatrix(6, 6));
			(dHistDeque_Pulls[locStepIndex][locPID])[dPullTypes[loc_j]]->Fill(locPull);
			(dHistDeque_PullsVsP[locStepIndex][locPID])[dPullTypes[loc_j]]->Fill(locThrownP, locPull);
			(dHistDeque_PullsVsTheta[locStepIndex][locPID])[dPullTypes[loc_j]]->Fill(locThrownTheta, locPull);
		}

	}
	japp->RootUnLock();
}

void DHistogramAction_ThrownParticleKinematics::Initialize(JEventLoop* locEventLoop)
{
	string locHistName, locHistTitle, locParticleName, locParticleROOTName;
	Particle_t locPID;

	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		CreateAndChangeTo_ActionDirectory();

		// Beam Particle: MCGEN
		{
			locPID = Gamma;
			locParticleName = string("MCGENBeamParticle_") + ParticleType(locPID);
			locParticleROOTName = ParticleName_ROOT(locPID);
			CreateAndChangeTo_Directory(locParticleName, locParticleName);

			locHistName = "Momentum";
			locHistTitle = string("MCGEN Thrown Beam ") + locParticleROOTName + string(";p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dMCGENBeamParticle_P = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dMCGENBeamParticle_P = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPBins, dMinP, dMaxP);

			locHistName = "Time";
			locHistTitle = string("MCGEN Thrown Beam ") + locParticleROOTName + string(";t (ns)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dMCGENBeamParticle_Time = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dMCGENBeamParticle_Time = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumTBins, dMinT, dMaxT);

			gDirectory->cd("..");
		}

		// Beam Particle: All
		{
			locPID = Gamma;
			locParticleName = string("TRUTHBeamParticles_") + ParticleType(locPID);
			locParticleROOTName = ParticleName_ROOT(locPID);
			CreateAndChangeTo_Directory(locParticleName, locParticleName);

			locHistName = "Momentum";
			locHistTitle = string("TRUTH Thrown Beam ") + locParticleROOTName + string(";p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dAllBeamParticle_P = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dAllBeamParticle_P = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPBins, dMinP, dMaxP);

			locHistName = "Time";
			locHistTitle = string("TRUTH Thrown Beam ") + locParticleROOTName + string(";t (ns)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dAllBeamParticle_Time = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dAllBeamParticle_Time = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumTBins, dMinT, dMaxT);

			gDirectory->cd("..");
		}

		for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
		{
			locPID = dFinalStatePIDs[loc_i];
			locParticleName = ParticleType(locPID);
			locParticleROOTName = ParticleName_ROOT(locPID);
			CreateAndChangeTo_Directory(locParticleName, locParticleName);

			// Momentum
			locHistName = "Momentum";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_P[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_P[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPBins, dMinP, dMaxP);

			// Theta
			locHistName = "Theta";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";#theta#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_Theta[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_Theta[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumThetaBins, dMinTheta, dMaxTheta);

			// Phi
			locHistName = "Phi";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";#phi#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_Phi[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_Phi[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPhiBins, dMinPhi, dMaxPhi);

			// P Vs Theta
			locHistName = "PVsTheta";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";#theta#circ;p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_PVsTheta[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PVsTheta[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			// Phi Vs Theta
			locHistName = "PhiVsTheta";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";#theta#circ;#phi#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_PhiVsTheta[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PhiVsTheta[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

			// Vertex-Z
			locHistName = "VertexZ";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";Vertex-Z (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_VertexZ[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_VertexZ[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumVertexZBins, dMinVertexZ, dMaxVertexZ);

			// Vertex-Y Vs Vertex-X
			locHistName = "VertexYVsX";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";Vertex-X (cm);Vertex-Y (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_VertexYVsX[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_VertexYVsX[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNumVertexXYBins, dMinVertexXY, dMaxVertexXY, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY);

			// Vertex-T
			locHistName = "VertexT";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";Vertex-T (ns)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_VertexT[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_VertexT[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumTBins, dMinT, dMaxT);

			gDirectory->cd("..");
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_ThrownParticleKinematics::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	vector<const DMCThrown*> locMCThrowns, locMCThrowns_Decaying;
	locEventLoop->Get(locMCThrowns, "FinalState");
	locEventLoop->Get(locMCThrowns_Decaying, "Decaying");
	locMCThrowns.insert(locMCThrowns.begin(), locMCThrowns_Decaying.begin(), locMCThrowns_Decaying.end());
	if(locMCThrowns.empty())
		return true; //e.g. non-simulated event

	if(Get_NumPreviousParticleCombos() != 0)
		return true; //else double-counting!

	Particle_t locPID;
	const DMCThrown* locMCThrown;

	vector<const DBeamPhoton*> locMCGENBeamPhotons;
	locEventLoop->Get(locMCGENBeamPhotons, "MCGEN");

	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons, "TRUTH");
	japp->RootWriteLock();
	{
		for(size_t loc_i = 0; loc_i < locMCGENBeamPhotons.size(); ++loc_i)
		{
			dMCGENBeamParticle_P->Fill(locMCGENBeamPhotons[loc_i]->energy());
			dMCGENBeamParticle_Time->Fill(locMCGENBeamPhotons[loc_i]->time());
		}
		for(size_t loc_i = 0; loc_i < locBeamPhotons.size(); ++loc_i)
		{
			dAllBeamParticle_P->Fill(locBeamPhotons[loc_i]->energy());
			dAllBeamParticle_Time->Fill(locBeamPhotons[loc_i]->time());
		}
	}
	japp->RootUnLock();

	for(size_t loc_i = 0; loc_i < locMCThrowns.size(); ++loc_i)
	{
		locMCThrown = locMCThrowns[loc_i];
		locPID = (Particle_t)locMCThrown->type;
		if(dHistMap_P.find(locPID) == dHistMap_P.end())
			continue; //not interested in histogramming

		DVector3 locMomentum = locMCThrown->momentum();
		double locPhi = locMomentum.Phi()*180.0/TMath::Pi();
		double locTheta = locMomentum.Theta()*180.0/TMath::Pi();
		double locP = locMomentum.Mag();
		japp->RootWriteLock();
		{
			dHistMap_P[locPID]->Fill(locP);
			dHistMap_Phi[locPID]->Fill(locPhi);
			dHistMap_Theta[locPID]->Fill(locTheta);
			dHistMap_PVsTheta[locPID]->Fill(locTheta, locP);
			dHistMap_PhiVsTheta[locPID]->Fill(locTheta, locPhi);
			dHistMap_VertexZ[locPID]->Fill(locMCThrown->position().Z());
			dHistMap_VertexYVsX[locPID]->Fill(locMCThrown->position().X(), locMCThrown->position().Y());
			dHistMap_VertexT[locPID]->Fill(locMCThrown->time());
		}
		japp->RootUnLock();
	}
	return true;
}

void DHistogramAction_ReconnedThrownKinematics::Initialize(JEventLoop* locEventLoop)
{
	string locHistName, locHistTitle, locParticleName, locParticleROOTName;
	Particle_t locPID;

	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		locEventLoop->GetSingle(dAnalysisUtilities);

		CreateAndChangeTo_ActionDirectory();

		// Beam Particle
		locPID = Gamma;
		locParticleName = string("Beam_") + ParticleType(locPID);
		locParticleROOTName = ParticleName_ROOT(locPID);
		CreateAndChangeTo_Directory(locParticleName, locParticleName);
		locHistName = "Momentum";
		locHistTitle = string("Thrown Beam ") + locParticleROOTName + string(";p (GeV/c)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dBeamParticle_P = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dBeamParticle_P = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPBins, dMinP, dMaxP);
		gDirectory->cd("..");

		//PID
		CreateAndChangeTo_Directory("PID", "PID");
		{
			//beta vs p
			locHistName = "BetaVsP_Q+";
			locHistTitle = "q^{+};p (GeV/c);#beta";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_QBetaVsP[1] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_QBetaVsP[1] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumBetaBins, dMinBeta, dMaxBeta);

			locHistName = "BetaVsP_Q-";
			locHistTitle = "q^{-};p (GeV/c);#beta";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_QBetaVsP[-1] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_QBetaVsP[-1] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumBetaBins, dMinBeta, dMaxBeta);
		}
		gDirectory->cd("..");

		for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
		{
			locPID = dFinalStatePIDs[loc_i];
			locParticleName = ParticleType(locPID);
			locParticleROOTName = ParticleName_ROOT(locPID);
			CreateAndChangeTo_Directory(locParticleName, locParticleName);

			// Momentum
			locHistName = "Momentum";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_P[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_P[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPBins, dMinP, dMaxP);

			// Theta
			locHistName = "Theta";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";#theta#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_Theta[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_Theta[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumThetaBins, dMinTheta, dMaxTheta);

			// Phi
			locHistName = "Phi";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";#phi#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_Phi[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_Phi[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPhiBins, dMinPhi, dMaxPhi);

			// P Vs Theta
			locHistName = "PVsTheta";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";#theta#circ;p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_PVsTheta[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PVsTheta[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			// Phi Vs Theta
			locHistName = "PhiVsTheta";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";#theta#circ;#phi#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_PhiVsTheta[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PhiVsTheta[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

			// Vertex-Z
			locHistName = "VertexZ";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";Vertex-Z (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_VertexZ[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_VertexZ[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumVertexZBins, dMinVertexZ, dMaxVertexZ);

			// Vertex-Y Vs Vertex-X
			locHistName = "VertexYVsX";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";Vertex-X (cm);Vertex-Y (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_VertexYVsX[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_VertexYVsX[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNumVertexXYBins, dMinVertexXY, dMaxVertexXY, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY);

			// Vertex-T
			locHistName = "VertexT";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";Vertex-T (ns)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_VertexT[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_VertexT[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumTBins, dMinT, dMaxT);

			gDirectory->cd("..");
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_ReconnedThrownKinematics::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	vector<const DMCThrown*> locMCThrowns, locMCThrowns_Decaying;
	locEventLoop->Get(locMCThrowns, "FinalState");
	locEventLoop->Get(locMCThrowns_Decaying, "Decaying");
	locMCThrowns.insert(locMCThrowns.begin(), locMCThrowns_Decaying.begin(), locMCThrowns_Decaying.end());
	if(locMCThrowns.empty())
		return true; //e.g. non-simulated event

	if(Get_NumPreviousParticleCombos() != 0)
		return true; //else double-counting!

	const DEventRFBunch* locEventRFBunch = NULL;
	locEventLoop->GetSingle(locEventRFBunch);
	if(locParticleCombo != NULL)
		locEventRFBunch = locParticleCombo->Get_EventRFBunch();

	const DMCThrownMatching* locMCThrownMatching = NULL;
	locEventLoop->GetSingle(locMCThrownMatching);

	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons);
	japp->RootWriteLock();
	{
		for(size_t loc_i = 0; loc_i < locBeamPhotons.size(); ++loc_i)
			dBeamParticle_P->Fill(locBeamPhotons[loc_i]->energy());
	}
	japp->RootUnLock();

	for(size_t loc_i = 0; loc_i < locMCThrowns.size(); ++loc_i)
	{
		const DMCThrown* locMCThrown = locMCThrowns[loc_i];
		Particle_t locPID = (Particle_t)locMCThrown->type;

		double locMatchFOM = 0.0;
		double locBeta_Timing = 0.0;
		if(ParticleCharge(locPID) != 0)
		{
			const DChargedTrackHypothesis* locChargedTrackHypothesis = locMCThrownMatching->Get_MatchingChargedHypothesis(locMCThrown, locMatchFOM);
			if(locChargedTrackHypothesis == NULL)
				continue; //not reconstructed
			locBeta_Timing = dAnalysisUtilities->Calc_Beta_Timing(locChargedTrackHypothesis, locEventRFBunch, true);
		}
		else
		{
			const DNeutralParticleHypothesis* locNeutralParticleHypothesis = locMCThrownMatching->Get_MatchingNeutralHypothesis(locMCThrown, locMatchFOM);
			if(locNeutralParticleHypothesis == NULL)
				continue; //not reconstructed
			locBeta_Timing = dAnalysisUtilities->Calc_Beta_Timing(locNeutralParticleHypothesis, locEventRFBunch);
		}

		if(locMatchFOM < dMinThrownMatchFOM)
			continue; //not well reconstructed

		DVector3 locMomentum = locMCThrown->momentum();
		double locPhi = locMomentum.Phi()*180.0/TMath::Pi();
		double locTheta = locMomentum.Theta()*180.0/TMath::Pi();
		double locP = locMomentum.Mag();
		int locCharge = ParticleCharge(locPID);

		japp->RootWriteLock();
		{
			if(dHistMap_QBetaVsP.find(locCharge) != dHistMap_QBetaVsP.end())
				dHistMap_QBetaVsP[locCharge]->Fill(locP, locBeta_Timing);
			if(dHistMap_P.find(locPID) != dHistMap_P.end())
			{
				dHistMap_P[locPID]->Fill(locP);
				dHistMap_Phi[locPID]->Fill(locPhi);
				dHistMap_Theta[locPID]->Fill(locTheta);
				dHistMap_PVsTheta[locPID]->Fill(locTheta, locP);
				dHistMap_PhiVsTheta[locPID]->Fill(locTheta, locPhi);
				dHistMap_VertexZ[locPID]->Fill(locMCThrown->position().Z());
				dHistMap_VertexYVsX[locPID]->Fill(locMCThrown->position().X(), locMCThrown->position().Y());
				dHistMap_VertexT[locPID]->Fill(locMCThrown->time());
			}
		}
		japp->RootUnLock();
	}
	return true;
}

void DHistogramAction_DetectedParticleKinematics::Initialize(JEventLoop* locEventLoop)
{
	string locHistName, locHistTitle, locParticleName, locParticleROOTName;
	Particle_t locPID;

	vector<const DParticleID*> locParticleIDs;
	locEventLoop->Get(locParticleIDs);

	const DAnalysisUtilities* locAnalysisUtilities = NULL;
	locEventLoop->GetSingle(locAnalysisUtilities);

	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		CreateAndChangeTo_ActionDirectory();

		dParticleID = locParticleIDs[0];
		dAnalysisUtilities = locAnalysisUtilities;

		// Event Vertex-Z
		locHistName = "EventVertexZ";
		locHistTitle = ";Event Vertex-Z (cm)";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dEventVertexZ = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dEventVertexZ = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumVertexZBins, dMinVertexZ, dMaxVertexZ);

		// Event Vertex-Y Vs Vertex-X
		locHistName = "EventVertexYVsX";
		locHistTitle = ";Event Vertex-X (cm);Event Vertex-Y (cm)";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dEventVertexYVsX = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
		else
			dEventVertexYVsX = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNumVertexXYBins, dMinVertexXY, dMaxVertexXY, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY);

		// Event Vertex-T
		locHistName = "EventVertexT";
		locHistTitle = ";Event Vertex Time (ns)";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dEventVertexT = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dEventVertexT = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumTBins, dMinT, dMaxT);

		// Beam Particle
		locPID = Gamma;
		locParticleName = string("Beam_") + ParticleType(locPID);
		locParticleROOTName = ParticleName_ROOT(locPID);
		CreateAndChangeTo_Directory(locParticleName, locParticleName);
		locHistName = "Momentum";
		locHistTitle = string("Beam ") + locParticleROOTName + string(";p (GeV/c)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dBeamParticle_P = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dBeamParticle_P = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPBins, dMinP, dMaxP);
		gDirectory->cd("..");

		//PID
		CreateAndChangeTo_Directory("PID", "PID");
		{
			//beta vs p
			locHistName = "BetaVsP_Q+";
			locHistTitle = "q^{+};p (GeV/c);#beta";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_QBetaVsP[1] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_QBetaVsP[1] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumBetaBins, dMinBeta, dMaxBeta);

			locHistName = "BetaVsP_Q-";
			locHistTitle = "q^{-};p (GeV/c);#beta";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_QBetaVsP[-1] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_QBetaVsP[-1] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumBetaBins, dMinBeta, dMaxBeta);
		}
		gDirectory->cd("..");

		for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
		{
			locPID = dFinalStatePIDs[loc_i];
			locParticleName = ParticleType(locPID);
			locParticleROOTName = ParticleName_ROOT(locPID);
			CreateAndChangeTo_Directory(locParticleName, locParticleName);

			// Momentum
			locHistName = "Momentum";
			locHistTitle = locParticleROOTName + string(";p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_P[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_P[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPBins, dMinP, dMaxP);

			// Theta
			locHistName = "Theta";
			locHistTitle = locParticleROOTName + string(";#theta#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_Theta[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_Theta[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumThetaBins, dMinTheta, dMaxTheta);

			// Phi
			locHistName = "Phi";
			locHistTitle = locParticleROOTName + string(";#phi#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_Phi[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_Phi[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPhiBins, dMinPhi, dMaxPhi);

			// P Vs Theta
			locHistName = "PVsTheta";
			locHistTitle = locParticleROOTName + string(";#theta#circ;p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_PVsTheta[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PVsTheta[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			// Phi Vs Theta
			locHistName = "PhiVsTheta";
			locHistTitle = locParticleROOTName + string(";#theta#circ;#phi#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_PhiVsTheta[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PhiVsTheta[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

			//beta vs p
			locHistName = "BetaVsP";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#beta");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_BetaVsP[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_BetaVsP[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumBetaBins, dMinBeta, dMaxBeta);

			//delta-beta vs p
			locHistName = "DeltaBetaVsP";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta#beta");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaBetaVsP[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaBetaVsP[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaBetaBins, dMinDeltaBeta, dMaxDeltaBeta);

			// Vertex-Z
			locHistName = "VertexZ";
			locHistTitle = locParticleROOTName + string(";Vertex-Z (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_VertexZ[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_VertexZ[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumVertexZBins, dMinVertexZ, dMaxVertexZ);

			// Vertex-Y Vs Vertex-X
			locHistName = "VertexYVsX";
			locHistTitle = locParticleROOTName + string(";Vertex-X (cm);Vertex-Y (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_VertexYVsX[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_VertexYVsX[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNumVertexXYBins, dMinVertexXY, dMaxVertexXY, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY);

			// Vertex-T
			locHistName = "VertexT";
			locHistTitle = locParticleROOTName + string(";Vertex-T (ns)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_VertexT[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_VertexT[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumTBins, dMinT, dMaxT);

			gDirectory->cd("..");
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_DetectedParticleKinematics::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	if(Get_NumPreviousParticleCombos() != 0)
		return true; //else double-counting!

	const DEventRFBunch* locEventRFBunch = NULL;
	locEventLoop->GetSingle(locEventRFBunch);
	if(locParticleCombo != NULL)
		locEventRFBunch = locParticleCombo->Get_EventRFBunch();

	const DVertex* locVertex = NULL;
	locEventLoop->GetSingle(locVertex);
	DLorentzVector locEventVertex = locVertex->dSpacetimeVertex;
	if(locParticleCombo != NULL)
		locEventVertex = locParticleCombo->Get_ParticleComboStep(0)->Get_SpacetimeVertex();

	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons);
	japp->RootWriteLock();
	{
		dEventVertexZ->Fill(locEventVertex.Z());
		dEventVertexYVsX->Fill(locEventVertex.X(), locEventVertex.Y());
		dEventVertexT->Fill(locEventVertex.T());
		for(size_t loc_i = 0; loc_i < locBeamPhotons.size(); ++loc_i)
			dBeamParticle_P->Fill(locBeamPhotons[loc_i]->energy());
	}
	japp->RootUnLock();

	vector<const DChargedTrack*> locPreSelectChargedTracks;
	locEventLoop->Get(locPreSelectChargedTracks, "PreSelect");

	for(size_t loc_i = 0; loc_i < locPreSelectChargedTracks.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locPreSelectChargedTracks[loc_i]->Get_BestFOM();
		int locCharge = ParticleCharge(locChargedTrackHypothesis->PID());

		DVector3 locMomentum = locChargedTrackHypothesis->momentum();
		double locPhi = locMomentum.Phi()*180.0/TMath::Pi();
		double locTheta = locMomentum.Theta()*180.0/TMath::Pi();
		double locP = locMomentum.Mag();
		double locBeta_Timing = dAnalysisUtilities->Calc_Beta_Timing(locChargedTrackHypothesis, locEventRFBunch, true);
		double locDeltaBeta = locChargedTrackHypothesis->lorentzMomentum().Beta() - locBeta_Timing;

		if(dHistMap_QBetaVsP.find(locCharge) != dHistMap_QBetaVsP.end())
			dHistMap_QBetaVsP[locCharge]->Fill(locP, locBeta_Timing);

		Particle_t locPID = (locChargedTrackHypothesis->dFOM < dMinPIDFOM) ? Unknown : locChargedTrackHypothesis->PID();
		if(dHistMap_P.find(locPID) == dHistMap_P.end())
			continue; //not interested in histogramming

		japp->RootWriteLock();
		{
			dHistMap_P[locPID]->Fill(locP);
			dHistMap_Phi[locPID]->Fill(locPhi);
			dHistMap_Theta[locPID]->Fill(locTheta);
			dHistMap_PVsTheta[locPID]->Fill(locTheta, locP);
			dHistMap_PhiVsTheta[locPID]->Fill(locTheta, locPhi);
			dHistMap_BetaVsP[locPID]->Fill(locP, locBeta_Timing);
			dHistMap_DeltaBetaVsP[locPID]->Fill(locP, locDeltaBeta);
			dHistMap_VertexZ[locPID]->Fill(locChargedTrackHypothesis->position().Z());
			dHistMap_VertexYVsX[locPID]->Fill(locChargedTrackHypothesis->position().X(), locChargedTrackHypothesis->position().Y());
			dHistMap_VertexT[locPID]->Fill(locChargedTrackHypothesis->time());
		}
		japp->RootUnLock();
	}

	vector<const DNeutralParticle*> locNeutralParticles;
	locEventLoop->Get(locNeutralParticles);

	vector<const DNeutralShower*> locNeutralShowers;
	locEventLoop->Get(locNeutralShowers, "PreSelect");

	for(size_t loc_i = 0; loc_i < locNeutralShowers.size(); ++loc_i)
	{
		const DNeutralParticle* locNeutralParticle = NULL;
		for(size_t loc_j = 0; loc_j < locNeutralParticles.size(); ++loc_j)
		{
			if(locNeutralParticles[loc_j]->dNeutralShower != locNeutralShowers[loc_i])
				continue;
			locNeutralParticle = locNeutralParticles[loc_j];
			break;
		}
		if(locNeutralParticle == NULL)
			continue;
		const DNeutralParticleHypothesis* locNeutralParticleHypothesis = locNeutralParticles[loc_i]->Get_Hypothesis(Gamma);
		if(locNeutralParticleHypothesis->dFOM < dMinPIDFOM)
			continue;

		Particle_t locPID = locNeutralParticleHypothesis->PID();

		if(dHistMap_P.find(locPID) == dHistMap_P.end())
			continue; //e.g. a decaying particle, or not interested in histogramming

		DVector3 locMomentum = locNeutralParticleHypothesis->momentum();
		double locPhi = locMomentum.Phi()*180.0/TMath::Pi();
		double locTheta = locMomentum.Theta()*180.0/TMath::Pi();
		double locP = locMomentum.Mag();

		double locBeta_Timing = dAnalysisUtilities->Calc_Beta_Timing(locNeutralParticleHypothesis, locEventRFBunch);
		double locDeltaBeta = locNeutralParticleHypothesis->lorentzMomentum().Beta() - locBeta_Timing;

		japp->RootWriteLock();
		{
			dHistMap_P[locPID]->Fill(locP);
			dHistMap_Phi[locPID]->Fill(locPhi);
			dHistMap_Theta[locPID]->Fill(locTheta);
			dHistMap_PVsTheta[locPID]->Fill(locTheta, locP);
			dHistMap_PhiVsTheta[locPID]->Fill(locTheta, locPhi);
			dHistMap_BetaVsP[locPID]->Fill(locP, locBeta_Timing);
			dHistMap_DeltaBetaVsP[locPID]->Fill(locP, locDeltaBeta);
			dHistMap_VertexZ[locPID]->Fill(locNeutralParticleHypothesis->position().Z());
			dHistMap_VertexYVsX[locPID]->Fill(locNeutralParticleHypothesis->position().X(), locNeutralParticleHypothesis->position().Y());
			dHistMap_VertexT[locPID]->Fill(locNeutralParticleHypothesis->time());
		}
		japp->RootUnLock();
	}
	return true;
}

void DHistogramAction_GenReconTrackComparison::Initialize(JEventLoop* locEventLoop)
{
	string locHistName, locHistTitle, locParticleName, locParticleROOTName;
	Particle_t locPID;

	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry *locGeometry = locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());

	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	CreateAndChangeTo_ActionDirectory();
	{
		locGeometry->GetTargetZ(dTargetZCenter);

		locHistName = "DeltaT_RFBeamBunch";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dRFBeamBunchDeltaT_Hist = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dRFBeamBunchDeltaT_Hist = new TH1I(locHistName.c_str(), ";RF #Deltat (Reconstructed - Thrown)", dNumRFDeltaTBins, dMinRFDeltaT, dMaxRFDeltaT);

		deque<string> locPullNames(8, "");
		locPullNames[0] = "E";  locPullNames[1] = "Px";  locPullNames[2] = "Py";  locPullNames[3] = "Pz";
		locPullNames[4] = "Xx";  locPullNames[5] = "Xy";  locPullNames[6] = "Xz";  locPullNames[7] = "T";

		deque<string> locPullTitles(8, "");
		locPullTitles[0] = "E";  locPullTitles[1] = "p_{x}";  locPullTitles[2] = "p_{y}";  locPullTitles[3] = "p_{z}";
		locPullTitles[4] = "x_{x}";  locPullTitles[5] = "x_{y}";  locPullTitles[6] = "x_{z}";  locPullTitles[7] = "t";

		for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
		{
			locPID = dFinalStatePIDs[loc_i];
			locParticleName = ParticleType(locPID);
			locParticleROOTName = ParticleName_ROOT(locPID);
			CreateAndChangeTo_Directory(locParticleName, locParticleName);

			// MatchChiSqPerDF
			locHistName = string("MatchFOM");
			locHistTitle = locParticleROOTName + string(";Thrown/Reconstructed Matching FOM");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_MatchFOM[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_MatchFOM[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumMCMatchingFOMBins, 0.0, 1.0);

			// DeltaP/P
			locHistName = string("DeltaPOverP");
			locHistTitle = locParticleROOTName + string(";#Deltap/p (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaPOverP[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaPOverP[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

			// DeltaTheta
			locHistName = string("DeltaTheta");
			locHistTitle = locParticleROOTName + string(";#Delta#theta#circ (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaTheta[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaTheta[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaThetaBins, dMinDeltaTheta, dMaxDeltaTheta);

			// DeltaPhi
			locHistName = string("DeltaPhi");
			locHistTitle = locParticleROOTName + string(";#Delta#phi#circ (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaPhi[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaPhi[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

			// DeltaT
			locHistName = string("DeltaT");
			locHistTitle = locParticleROOTName + string(";#Deltat (ns) (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaT[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaT[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

			// DeltaT - BCAL
			locHistName = string("DeltaT_BCAL");
			locHistTitle = locParticleROOTName + string(" in BCAL;#Deltat (ns) (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaT_BCAL[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaT_BCAL[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

			// DeltaT - TOF (charged only)
			if(ParticleCharge(locPID) != 0)
			{
				locHistName = string("DeltaT_TOF");
				locHistTitle = locParticleROOTName + string(" in TOF;#Deltat (ns) (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_DeltaT_TOF[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_DeltaT_TOF[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
			}

			// DeltaT - FCAL (neutral only)
			if(ParticleCharge(locPID) == 0)
			{
				locHistName = string("DeltaT_FCAL");
				locHistTitle = locParticleROOTName + string(" in FCAL;#Deltat (ns) (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_DeltaT_FCAL[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_DeltaT_FCAL[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
			}

			// DeltaVertexZ
			locHistName = string("DeltaVertexZ");
			locHistTitle = locParticleROOTName + string(";#DeltaVertex-Z (cm) (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaVertexZ[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaVertexZ[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaVertexZBins, dMinDeltaVertexZ, dMaxDeltaVertexZ);

			// DeltaP/P Vs P
			locHistName = string("DeltaPOverPVsP");
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltap/p (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaPOverPVsP[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaPOverPVsP[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

			// DeltaP/P Vs Theta
			locHistName = string("DeltaPOverPVsTheta");
			locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltap/p (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaPOverPVsTheta[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaPOverPVsTheta[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

			// DeltaTheta Vs P
			locHistName = string("DeltaThetaVsP");
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta#theta#circ (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaThetaVsP[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaThetaVsP[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaThetaBins, dMinDeltaTheta, dMaxDeltaTheta);

			// DeltaTheta Vs Theta
			locHistName = string("DeltaThetaVsTheta");
			locHistTitle = locParticleROOTName + string(";#theta#circ;#Delta#theta#circ (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaThetaVsTheta[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaThetaVsTheta[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaThetaBins, dMinDeltaTheta, dMaxDeltaTheta);

			// DeltaPhi Vs P
			locHistName = string("DeltaPhiVsP");
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta#phi#circ (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaPhiVsP[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaPhiVsP[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

			// DeltaPhi Vs Theta
			locHistName = string("DeltaPhiVsTheta");
			locHistTitle = locParticleROOTName + string(";#theta#circ;#Delta#phi#circ (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaPhiVsTheta[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaPhiVsTheta[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

			// DeltaT Vs Theta
			locHistName = string("DeltaTVsTheta");
			locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat (ns) (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaTVsTheta[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaTVsTheta[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

			// DeltaT Vs P
			locHistName = string("DeltaTVsP");
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat (ns) (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaTVsP[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaTVsP[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

			// DeltaVertexZ Vs Theta
			locHistName = string("DeltaVertexZVsTheta");
			locHistTitle = locParticleROOTName + string(";#theta#circ;#DeltaVertex-Z (cm) (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaVertexZVsTheta[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaVertexZVsTheta[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaVertexZBins, dMinDeltaVertexZ, dMaxDeltaVertexZ);

			// P Vs Theta
			locHistName = "PVsTheta_LargeDeltaT";
			locHistTitle = locParticleROOTName + string(";#theta#circ;p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_PVsTheta_LargeDeltaT[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PVsTheta_LargeDeltaT[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			/************************************************************************ Pulls ************************************************************************/

			CreateAndChangeTo_Directory("Pulls", "Pulls");
			for(size_t loc_j = 0; loc_j < dPullTypes.size(); ++loc_j)
			{
				if((ParticleCharge(locPID) != 0) && (dPullTypes[loc_j] == d_EPull))
					continue;
				if((ParticleCharge(locPID) == 0) && ((dPullTypes[loc_j] >= d_PxPull) && (dPullTypes[loc_j] <= d_PzPull)))
					continue;

				//Pull 1D
				locHistName = locPullNames[loc_j] + string("Pull");
				locHistTitle = locParticleROOTName + string(";#Delta") + locPullTitles[loc_j] + string("/#sigma_{") + locPullTitles[loc_j] + string("} (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					(dHistMap_Pulls[locPID])[dPullTypes[loc_j]] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					(dHistMap_Pulls[locPID])[dPullTypes[loc_j]] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

				//Pull vs P
				locHistName = locPullNames[loc_j] + string("PullVsP");
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta") + locPullTitles[loc_j] + string("/#sigma_{") + locPullTitles[loc_j] + string("} (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					(dHistMap_PullsVsP[locPID])[dPullTypes[loc_j]] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					(dHistMap_PullsVsP[locPID])[dPullTypes[loc_j]] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

				//Pull vs Theta
				locHistName = locPullNames[loc_j] + string("PullVsTheta");
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Delta") + locPullTitles[loc_j] + string("/#sigma_{") + locPullTitles[loc_j] + string("} (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					(dHistMap_PullsVsTheta[locPID])[dPullTypes[loc_j]] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					(dHistMap_PullsVsTheta[locPID])[dPullTypes[loc_j]] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);
			}

			//Delta-t Pulls - CDC & ST
			if(ParticleCharge(locPID) != 0)
			{
				//CDC
				locHistName = "TimePull_CDC";
				locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TimePull_CDC[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePull_CDC[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

				locHistName = "TimePullVsTheta_CDC";
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TimePullVsTheta_CDC[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePullVsTheta_CDC[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);

				locHistName = "TimePullVsP_CDC";
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TimePullVsP_CDC[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePullVsP_CDC[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

				//ST
				locHistName = "TimePull_ST";
				locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TimePull_ST[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePull_ST[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

				locHistName = "TimePullVsTheta_ST";
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TimePullVsTheta_ST[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePullVsTheta_ST[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);

				locHistName = "TimePullVsP_ST";
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TimePullVsP_ST[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePullVsP_ST[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);
			}

			//Delta-t Pulls - BCAL
			locHistName = "TimePull_BCAL";
			locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TimePull_BCAL[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TimePull_BCAL[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

			locHistName = "TimePullVsTheta_BCAL";
			locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat/#sigma_{#Deltat}");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TimePullVsTheta_BCAL[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TimePullVsTheta_BCAL[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);

			locHistName = "TimePullVsP_BCAL";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TimePullVsP_BCAL[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TimePullVsP_BCAL[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

			//Delta-t Pulls - TOF
			if(ParticleCharge(locPID) != 0) //TOF
			{
				locHistName = "TimePull_TOF";
				locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TimePull_TOF[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePull_TOF[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

				locHistName = "TimePullVsP_TOF";
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TimePullVsP_TOF[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePullVsP_TOF[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);
			}

			//Delta-t Pulls - FCAL
			locHistName = "TimePull_FCAL";
			locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TimePull_FCAL[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TimePull_FCAL[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

			locHistName = "TimePullVsP_FCAL";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TimePullVsP_FCAL[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TimePullVsP_FCAL[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

			gDirectory->cd("..");

			gDirectory->cd("..");
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_GenReconTrackComparison::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);
	if(locMCThrowns.empty())
		return true; //e.g. non-simulated event

	if(Get_NumPreviousParticleCombos() != 0)
		return true; //else double-counting!

	Particle_t locPID;
	double locDeltaPOverP, locDeltaTheta, locDeltaPhi, locDeltaVertexZ;
	double locThrownP, locThrownTheta, locDeltaT;

	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
	locEventLoop->Get(locMCThrownMatchingVector);
	if(locMCThrownMatchingVector.empty())
		return true;
	const DMCThrownMatching* locMCThrownMatching = locMCThrownMatchingVector[0];

	const DEventRFBunch* locThrownEventRFBunch = NULL;
	locEventLoop->GetSingle(locThrownEventRFBunch, "Thrown");

	//RF time difference
	vector<const DEventRFBunch*> locEventRFBunches;
	locEventLoop->Get(locEventRFBunches);
	const DEventRFBunch* locEventRFBunch = locEventRFBunches[0];
	double locRFTime = locEventRFBunch->dTime;
	double locRFDeltaT = locRFTime - locThrownEventRFBunch->dTime;
	japp->RootWriteLock();
	{
		dRFBeamBunchDeltaT_Hist->Fill(locRFDeltaT);
	}
	japp->RootUnLock();

	//charged particles
	map<const DMCThrown*, pair<const DChargedTrack*, double> > locThrownToChargedMap;
	locMCThrownMatching->Get_ThrownToChargedMap(locThrownToChargedMap);
	map<const DMCThrown*, pair<const DChargedTrack*, double> >::iterator locChargedIterator = locThrownToChargedMap.begin();
	for(; locChargedIterator != locThrownToChargedMap.end(); ++locChargedIterator)
	{
		const DMCThrown* locMCThrown = locChargedIterator->first;
		locPID = (Particle_t)locMCThrown->type;
		if(dHistMap_DeltaPOverP.find(locPID) == dHistMap_DeltaPOverP.end())
			continue; //e.g. not interested in histogramming

		double locMatchFOM = locChargedIterator->second.second;
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedIterator->second.first->Get_Hypothesis(locPID);
		if(locChargedTrackHypothesis == NULL)
			locChargedTrackHypothesis = locChargedIterator->second.first->Get_BestFOM();

		locThrownP = locMCThrown->momentum().Mag();
		locThrownTheta = locMCThrown->momentum().Theta()*180.0/TMath::Pi();
		locDeltaPOverP = (locChargedTrackHypothesis->momentum().Mag() - locThrownP)/locThrownP;
		locDeltaTheta = locChargedTrackHypothesis->momentum().Theta()*180.0/TMath::Pi() - locThrownTheta;
		locDeltaPhi = locChargedTrackHypothesis->momentum().Phi()*180.0/TMath::Pi() - locMCThrown->momentum().Phi()*180.0/TMath::Pi();
		locDeltaT = locChargedTrackHypothesis->time() - locMCThrown->time(); //time comparison isn't fair if track comes from a detached vertex!!!
		locDeltaVertexZ = locChargedTrackHypothesis->position().Z() - locMCThrown->position().Z();
		const TMatrixDSym& locCovarianceMatrix = locChargedTrackHypothesis->errorMatrix();

		vector<const DTrackTimeBased*> locTrackTimeBasedVector;
		locChargedTrackHypothesis->Get(locTrackTimeBasedVector);
		const DTrackTimeBased* locTrackTimeBased = locTrackTimeBasedVector[0];

		double locStartTime = locThrownEventRFBunch->dTime + (locMCThrown->z() - dTargetZCenter)/29.9792458;
		double locTimePull = (locStartTime - locChargedTrackHypothesis->time())/sqrt(locCovarianceMatrix(6, 6));
		double locT0Pull = (locStartTime - locChargedTrackHypothesis->t0())/locChargedTrackHypothesis->t0_err();
		japp->RootWriteLock();
		{
			dHistMap_MatchFOM[locPID]->Fill(locMatchFOM);
			dHistMap_DeltaPOverP[locPID]->Fill(locDeltaPOverP);
			dHistMap_DeltaTheta[locPID]->Fill(locDeltaTheta);
			dHistMap_DeltaPhi[locPID]->Fill(locDeltaPhi);
			dHistMap_DeltaT[locPID]->Fill(locDeltaT);
			if(locChargedTrackHypothesis->t0_detector() == SYS_START)
			{
				dHistMap_TimePull_ST[locPID]->Fill(locT0Pull);
				dHistMap_TimePullVsTheta_ST[locPID]->Fill(locThrownTheta, locT0Pull);
				dHistMap_TimePullVsP_ST[locPID]->Fill(locThrownP, locT0Pull);
			}
			if(locChargedTrackHypothesis->t0_detector() == SYS_CDC)
			{
				dHistMap_TimePull_CDC[locPID]->Fill(locT0Pull);
				dHistMap_TimePullVsTheta_CDC[locPID]->Fill(locThrownTheta, locT0Pull);
				dHistMap_TimePullVsP_CDC[locPID]->Fill(locThrownP, locT0Pull);
			}
			else if(locChargedTrackHypothesis->t1_detector() == SYS_CDC)
			{
				dHistMap_TimePull_CDC[locPID]->Fill(locTimePull);
				dHistMap_TimePullVsTheta_CDC[locPID]->Fill(locThrownTheta, locTimePull);
				dHistMap_TimePullVsP_CDC[locPID]->Fill(locThrownP, locTimePull);
			}
			if(locChargedTrackHypothesis->t1_detector() == SYS_BCAL)
			{
				dHistMap_DeltaT_BCAL[locPID]->Fill(locDeltaT);
				dHistMap_TimePull_BCAL[locPID]->Fill(locTimePull);
				dHistMap_TimePullVsTheta_BCAL[locPID]->Fill(locThrownTheta, locTimePull);
				dHistMap_TimePullVsP_BCAL[locPID]->Fill(locThrownP, locTimePull);
			}
			else if(locChargedTrackHypothesis->t1_detector() == SYS_TOF)
			{
				dHistMap_DeltaT_TOF[locPID]->Fill(locDeltaT);
				dHistMap_TimePull_TOF[locPID]->Fill(locTimePull);
				dHistMap_TimePullVsP_TOF[locPID]->Fill(locThrownP, locTimePull);
			}
			else if(locChargedTrackHypothesis->t1_detector() == SYS_FCAL)
			{
				dHistMap_TimePull_FCAL[locPID]->Fill(locTimePull);
				dHistMap_TimePullVsP_FCAL[locPID]->Fill(locThrownP, locTimePull);
			}
			dHistMap_DeltaVertexZ[locPID]->Fill(locDeltaVertexZ);
			dHistMap_DeltaPOverPVsP[locPID]->Fill(locThrownP, locDeltaPOverP);
			dHistMap_DeltaPOverPVsTheta[locPID]->Fill(locThrownTheta, locDeltaPOverP);
			dHistMap_DeltaThetaVsP[locPID]->Fill(locThrownP, locDeltaTheta);
			dHistMap_DeltaThetaVsTheta[locPID]->Fill(locThrownTheta, locDeltaTheta);
			dHistMap_DeltaPhiVsP[locPID]->Fill(locThrownP, locDeltaPhi);
			dHistMap_DeltaPhiVsTheta[locPID]->Fill(locThrownTheta, locDeltaPhi);
			dHistMap_DeltaTVsTheta[locPID]->Fill(locThrownTheta, locDeltaT);
			dHistMap_DeltaTVsP[locPID]->Fill(locThrownP, locDeltaT);
			dHistMap_DeltaVertexZVsTheta[locPID]->Fill(locThrownTheta, locDeltaVertexZ);
			if((locTrackTimeBased->FOM > 0.01) && (locDeltaT >= 1.002))
				dHistMap_PVsTheta_LargeDeltaT[locPID]->Fill(locThrownTheta, locThrownP);

			for(size_t loc_j = 0; loc_j < dPullTypes.size(); ++loc_j)
			{
				if(dPullTypes[loc_j] == d_EPull)
					continue;
				double locPull = 0.0;
				if((dPullTypes[loc_j] >= d_PxPull) && (dPullTypes[loc_j] <= d_PzPull))
				{
					int locIndex = int(dPullTypes[loc_j] - d_PxPull);
					locPull = (locChargedTrackHypothesis->momentum()(locIndex) - locMCThrown->momentum()(locIndex))/sqrt(locCovarianceMatrix(locIndex, locIndex));
				}
				else if((dPullTypes[loc_j] >= d_XxPull) && (dPullTypes[loc_j] <= d_XzPull))
				{
					int locIndex = int(dPullTypes[loc_j] - d_XxPull);
					locPull = (locChargedTrackHypothesis->position()(locIndex) - locMCThrown->position()(locIndex))/sqrt(locCovarianceMatrix(locIndex + 3, locIndex + 3));
				}
				else if(dPullTypes[loc_j] == d_TPull)
					locPull = (locChargedTrackHypothesis->time() - locMCThrown->time())/sqrt(locCovarianceMatrix(6, 6));
				(dHistMap_Pulls[locPID])[dPullTypes[loc_j]]->Fill(locPull);
				(dHistMap_PullsVsP[locPID])[dPullTypes[loc_j]]->Fill(locThrownP, locPull);
				(dHistMap_PullsVsTheta[locPID])[dPullTypes[loc_j]]->Fill(locThrownTheta, locPull);
			}
		}
		japp->RootUnLock();
	}

	//neutral particles
	map<const DMCThrown*, pair<const DNeutralParticle*, double> > locThrownToNeutralMap;
	locMCThrownMatching->Get_ThrownToNeutralMap(locThrownToNeutralMap);
	map<const DMCThrown*, pair<const DNeutralParticle*, double> >::iterator locNeutralIterator = locThrownToNeutralMap.begin();
	for(; locNeutralIterator != locThrownToNeutralMap.end(); ++locNeutralIterator)
	{
		const DMCThrown* locMCThrown = locNeutralIterator->first;
		locPID = (Particle_t)locMCThrown->type;
		if(dHistMap_DeltaPOverP.find(locPID) == dHistMap_DeltaPOverP.end())
			continue; //e.g. not interested in histogramming

		double locMatchFOM = locNeutralIterator->second.second;
		const DNeutralParticleHypothesis* locNeutralParticleHypothesis = locNeutralIterator->second.first->Get_Hypothesis(locPID);
		if(locNeutralParticleHypothesis == NULL)
			locNeutralParticleHypothesis = locNeutralIterator->second.first->Get_BestFOM();

		const DNeutralShower* locNeutralShower = NULL;
		locNeutralParticleHypothesis->GetSingle(locNeutralShower);
		if(locNeutralShower == NULL)
			continue; //shouldn't be possible ...

		locThrownP = locMCThrown->momentum().Mag();
		locThrownTheta = locMCThrown->momentum().Theta()*180.0/TMath::Pi();
		locDeltaPOverP = (locNeutralParticleHypothesis->momentum().Mag() - locThrownP)/locThrownP;
		locDeltaTheta = locNeutralParticleHypothesis->momentum().Theta()*180.0/TMath::Pi() - locThrownTheta;
		locDeltaPhi = locNeutralParticleHypothesis->momentum().Phi()*180.0/TMath::Pi() - locMCThrown->momentum().Phi()*180.0/TMath::Pi();
		locDeltaT = locNeutralParticleHypothesis->time() - locMCThrown->time(); //time comparison isn't fair if track comes from a detached vertex!!!
		locDeltaVertexZ = locNeutralParticleHypothesis->position().Z() - locMCThrown->position().Z();
		const TMatrixDSym& locCovarianceMatrix = locNeutralParticleHypothesis->errorMatrix();

		double locStartTime = locThrownEventRFBunch->dTime + (locMCThrown->z() - dTargetZCenter)/29.9792458;
		double locTimePull = (locStartTime - locNeutralParticleHypothesis->time())/sqrt(locCovarianceMatrix(6, 6));

		japp->RootWriteLock();
		{
			dHistMap_MatchFOM[locPID]->Fill(locMatchFOM);
			dHistMap_DeltaPOverP[locPID]->Fill(locDeltaPOverP);
			dHistMap_DeltaTheta[locPID]->Fill(locDeltaTheta);
			dHistMap_DeltaPhi[locPID]->Fill(locDeltaPhi);
			dHistMap_DeltaT[locPID]->Fill(locDeltaT);
			if(locNeutralParticleHypothesis->t1_detector() == SYS_BCAL)
			{
				dHistMap_DeltaT_BCAL[locPID]->Fill(locDeltaT);
				dHistMap_TimePull_BCAL[locPID]->Fill(locTimePull);
				dHistMap_TimePullVsTheta_BCAL[locPID]->Fill(locThrownTheta, locTimePull);
				dHistMap_TimePullVsP_BCAL[locPID]->Fill(locThrownP, locTimePull);
			}
			else if(locNeutralParticleHypothesis->t1_detector() == SYS_FCAL)
			{
				dHistMap_DeltaT_FCAL[locPID]->Fill(locDeltaT);
				dHistMap_TimePull_FCAL[locPID]->Fill(locTimePull);
				dHistMap_TimePullVsP_FCAL[locPID]->Fill(locThrownP, locTimePull);
			}

			dHistMap_DeltaVertexZ[locPID]->Fill(locDeltaVertexZ);
			dHistMap_DeltaPOverPVsP[locPID]->Fill(locThrownP, locDeltaPOverP);
			dHistMap_DeltaPOverPVsTheta[locPID]->Fill(locThrownTheta, locDeltaPOverP);
			dHistMap_DeltaThetaVsP[locPID]->Fill(locThrownP, locDeltaTheta);
			dHistMap_DeltaThetaVsTheta[locPID]->Fill(locThrownTheta, locDeltaTheta);
			dHistMap_DeltaPhiVsP[locPID]->Fill(locThrownP, locDeltaPhi);
			dHistMap_DeltaPhiVsTheta[locPID]->Fill(locThrownTheta, locDeltaPhi);
			dHistMap_DeltaTVsTheta[locPID]->Fill(locThrownTheta, locDeltaT);
			dHistMap_DeltaTVsP[locPID]->Fill(locThrownP, locDeltaT);
			dHistMap_DeltaVertexZVsTheta[locPID]->Fill(locThrownTheta, locDeltaVertexZ);
			if(locDeltaT >= 1.002)
				dHistMap_PVsTheta_LargeDeltaT[locPID]->Fill(locThrownTheta, locThrownP);

			for(size_t loc_j = 0; loc_j < dPullTypes.size(); ++loc_j)
			{
				if((dPullTypes[loc_j] >= d_PxPull) && (dPullTypes[loc_j] <= d_PzPull))
					continue;
				double locPull = 0.0;
				if(dPullTypes[loc_j] == d_EPull)
					locPull = (locNeutralShower->dEnergy - locMCThrown->energy())/sqrt(locNeutralShower->dCovarianceMatrix(0, 0));
				else if((dPullTypes[loc_j] >= d_XxPull) && (dPullTypes[loc_j] <= d_XzPull))
				{
					int locIndex = int(dPullTypes[loc_j] - d_XxPull);
					locPull = (locNeutralParticleHypothesis->position()(locIndex) - locMCThrown->position()(locIndex))/sqrt(locCovarianceMatrix(locIndex + 3, locIndex + 3));
				}
				else if(dPullTypes[loc_j] == d_TPull)
					locPull = (locNeutralParticleHypothesis->time() - locMCThrown->time())/sqrt(locCovarianceMatrix(6, 6));
				(dHistMap_Pulls[locPID])[dPullTypes[loc_j]]->Fill(locPull);
				(dHistMap_PullsVsP[locPID])[dPullTypes[loc_j]]->Fill(locThrownP, locPull);
				(dHistMap_PullsVsTheta[locPID])[dPullTypes[loc_j]]->Fill(locThrownTheta, locPull);
			}

		}
		japp->RootUnLock();
	}
	return true;
}

void DHistogramAction_TOFHitStudy::Initialize(JEventLoop* locEventLoop)
{
	string locHistName, locHistTitle, locParticleName, locParticleROOTName;
	Particle_t locPID;

	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	CreateAndChangeTo_ActionDirectory();

	for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
	{
		locPID = dFinalStatePIDs[loc_i];
		locParticleName = ParticleType(locPID);
		locParticleROOTName = ParticleName_ROOT(locPID);
		CreateAndChangeTo_Directory(locParticleName, locParticleName);

		// DeltaT
		locHistName = string("DeltaT_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";#Deltat (ns) (Reconstructed - Thrown)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_DeltaT[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_DeltaT[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

		// DeltaX
		locHistName = string("DeltaX_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";#Deltax (cm) (Reconstructed - Thrown)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_DeltaX[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_DeltaX[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaXBins, dMinDeltaX, dMaxDeltaX);

		// DeltaY
		locHistName = string("DeltaY_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";#Deltay (cm) (Reconstructed - Thrown)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_DeltaY[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_DeltaY[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaXBins, dMinDeltaX, dMaxDeltaX);

		// dE
		locHistName = string("dE_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";dE (MeV)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_dE[locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_dE[locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumdEBins, dMindE, dMaxdE);

		// DeltaT Vs P
		locHistName = string("DeltaTVsP_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat (ns) (Reconstructed - Thrown)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_DeltaTVsP[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_DeltaTVsP[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

		// DeltaX Vs P
		locHistName = string("DeltaXVsP_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltax (cm) (Reconstructed - Thrown)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_DeltaXVsP[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_DeltaXVsP[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaXBins, dMinDeltaX, dMaxDeltaX);

		// DeltaY Vs P
		locHistName = string("DeltaYVsP_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltay (cm) (Reconstructed - Thrown)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_DeltaYVsP[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_DeltaYVsP[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaXBins, dMinDeltaX, dMaxDeltaX);

		// dE Vs P
		locHistName = string("dEVsP_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";p (GeV/c);dE (GeV)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_dEVsP[locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_dEVsP[locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumdEBins, dMindE, dMaxdE);

		gDirectory->cd("..");
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_TOFHitStudy::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	if(Get_NumPreviousParticleCombos() != 0)
		return true; //else double-counting!

	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
	locEventLoop->Get(locMCThrownMatchingVector);
	if(locMCThrownMatchingVector.empty())
		return true;
	const DMCThrownMatching* locMCThrownMatching = locMCThrownMatchingVector[0];

	vector<const DMCThrown*> locMCThrownVector;
	locEventLoop->Get(locMCThrownVector);

	map<const DTOFTruth*, pair<const DTOFPoint*, double> > locTOFTruthToPointMap;
	locMCThrownMatching->Get_TOFTruthToPointMap(locTOFTruthToPointMap);

	map<const DTOFTruth*, pair< const DTOFPoint*, double> >::iterator locTOFIterator;
	for(locTOFIterator = locTOFTruthToPointMap.begin(); locTOFIterator != locTOFTruthToPointMap.end(); ++locTOFIterator)
	{
		const DTOFTruth* locTOFTruth = locTOFIterator->first;
		const DTOFPoint* locTOFPoint = locTOFIterator->second.first;
		const DMCThrown* locMCThrown = NULL;
		for(size_t loc_i = 0; loc_i < locMCThrownVector.size(); ++loc_i)
		{
			if(locMCThrownVector[loc_i]->myid != locTOFTruth->track)
				continue;
			locMCThrown = locMCThrownVector[loc_i];
			break;
		}

		Particle_t locPID = (locMCThrown == NULL) ? Unknown : locMCThrown->PID();
		if(dHistMap_DeltaT.find(locPID) == dHistMap_DeltaT.end())
			continue;

		DVector3 locMomentumAtTOF(locTOFTruth->px, locTOFTruth->py, locTOFTruth->pz);
		DVector3 locThrownMomentum = (locMCThrown == NULL) ? locMomentumAtTOF : locMCThrown->momentum();
		double locThrownPMag = locThrownMomentum.Mag();

		//DTOFPoint and DTOFTruth reported at different z's (I think center vs. detector face): propagate truth information to the reconstructed z
		double locDeltaZ = locTOFPoint->pos.Z() - locTOFTruth->z;
		double locDeltaPathLength = locDeltaZ/cos(locMomentumAtTOF.Theta());
		double locPropagatedTrueX = locTOFTruth->x + locDeltaPathLength*sin(locMomentumAtTOF.Theta())*cos(locMomentumAtTOF.Phi());
		double locPropagatedTrueY = locTOFTruth->y + locDeltaPathLength*sin(locMomentumAtTOF.Theta())*sin(locMomentumAtTOF.Phi());
		double locVelocity = 29.9792458*locMomentumAtTOF.Mag()/locTOFTruth->E;
		double locPropagatedTrueT = locTOFTruth->t + locDeltaPathLength/locVelocity;

		double locDeltaT = locTOFPoint->t - locPropagatedTrueT;
		double locDeltaX = locTOFPoint->pos.X() - locPropagatedTrueX;
		double locDeltaY = locTOFPoint->pos.Y() - locPropagatedTrueY;

		double locdE_MeV = locTOFPoint->dE*1000.0;

		japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
		{
			dHistMap_DeltaT[locPID]->Fill(locDeltaT);
			dHistMap_DeltaX[locPID]->Fill(locDeltaX);
			dHistMap_DeltaY[locPID]->Fill(locDeltaY);
			dHistMap_dE[locPID]->Fill(locdE_MeV);
			dHistMap_DeltaTVsP[locPID]->Fill(locThrownPMag, locDeltaT);
			dHistMap_DeltaXVsP[locPID]->Fill(locThrownPMag, locDeltaX);
			dHistMap_DeltaYVsP[locPID]->Fill(locThrownPMag, locDeltaY);
			dHistMap_dEVsP[locPID]->Fill(locThrownPMag, locdE_MeV);
		}
		japp->RootUnLock(); //RELEASE ROOT LOCK!!
	}

	return true;
}

void DHistogramAction_NumReconstructedObjects::Initialize(JEventLoop* locEventLoop)
{
	string locHistName;

	bool locIsRESTEvent = (string(locEventLoop->GetJEvent().GetJEventSource()->className()) == string("DEventSourceREST"));

	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		CreateAndChangeTo_ActionDirectory();

		//2D Summary
		locHistName = "NumHighLevelObjects";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumHighLevelObjects = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
		else
		{
			dHist_NumHighLevelObjects = new TH2D(locHistName.c_str(), ";;# Objects / Event", 14, 0.5, 14.5, dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(1, "DTAGMHit");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(2, "DTAGHHit");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(3, "DSCHit");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(4, "DTOFPoint");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(5, "DBCALShower");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(6, "DFCALShower");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(7, "DTimeBasedTrack");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(8, "TrackSCMatches");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(9, "TrackTOFMatches");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(10, "TrackBCALMatches");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(11, "TrackFCALMatches");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(12, "DBeamPhoton");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(13, "DChargedTrack");
			dHist_NumHighLevelObjects->GetXaxis()->SetBinLabel(14, "DNeutralShower");
		}

		//Charged
		locHistName = "NumChargedTracks";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumChargedTracks = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumChargedTracks = new TH1D(locHistName.c_str(), ";# DChargedTrack", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumPosChargedTracks";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumPosChargedTracks = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumPosChargedTracks = new TH1D(locHistName.c_str(), ";# #it{q}^{+} DChargedTrack", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumNegChargedTracks";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumNegChargedTracks = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumNegChargedTracks = new TH1D(locHistName.c_str(), ";# #it{q}^{-} DChargedTrack", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		//TBT
		locHistName = "NumTimeBasedTracks";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumTimeBasedTracks = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumTimeBasedTracks = new TH1D(locHistName.c_str(), ";# Time-Based Tracks", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumPosTimeBasedTracks";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumPosTimeBasedTracks = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumPosTimeBasedTracks = new TH1D(locHistName.c_str(), ";# #it{q}^{+} Time-Based Tracks", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumNegTimeBasedTracks";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumNegTimeBasedTracks = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumNegTimeBasedTracks = new TH1D(locHistName.c_str(), ";# #it{q}^{-} Time-Based Tracks", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		if(!locIsRESTEvent)
		{
			//WBT
			locHistName = "NumWireBasedTracks";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumWireBasedTracks = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumWireBasedTracks = new TH1D(locHistName.c_str(), ";# Wire-Based Tracks", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			locHistName = "NumPosWireBasedTracks";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumPosWireBasedTracks = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumPosWireBasedTracks = new TH1D(locHistName.c_str(), ";# #it{q}^{-} Wire-Based Tracks", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			locHistName = "NumNegWireBasedTracks";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumNegWireBasedTracks = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumNegWireBasedTracks = new TH1D(locHistName.c_str(), ";# #it{q}^{-} Wire-Based Tracks", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			//Track Candidates
			locHistName = "NumTrackCandidates";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumTrackCandidates = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumTrackCandidates = new TH1D(locHistName.c_str(), ";# Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			locHistName = "NumPosTrackCandidates";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumPosTrackCandidates = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumPosTrackCandidates = new TH1D(locHistName.c_str(), ";# #it{q}^{+} Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			locHistName = "NumNegTrackCandidates";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumNegTrackCandidates = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumNegTrackCandidates = new TH1D(locHistName.c_str(), ";# #it{q}^{-} Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			//CDC Track Candidates
			locHistName = "NumPosTrackCandidates_CDC";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumPosTrackCandidates_CDC = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumPosTrackCandidates_CDC = new TH1D(locHistName.c_str(), ";# #it{q}^{+} CDC Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			locHistName = "NumNegTrackCandidates_CDC";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumNegTrackCandidates_CDC = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumNegTrackCandidates_CDC = new TH1D(locHistName.c_str(), ";# #it{q}^{-} CDC Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			//FDC Track Candidates
			locHistName = "NumPosTrackCandidates_FDC";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumPosTrackCandidates_FDC = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumPosTrackCandidates_FDC = new TH1D(locHistName.c_str(), ";# #it{q}^{+} FDC Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			locHistName = "NumNegTrackCandidates_FDC";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumNegTrackCandidates_FDC = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumNegTrackCandidates_FDC = new TH1D(locHistName.c_str(), ";# #it{q}^{-} FDC Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);
		}

		//Beam Photons
		locHistName = "NumBeamPhotons";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumBeamPhotons = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumBeamPhotons = new TH1D(locHistName.c_str(), ";# DBeamPhoton", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		//Showers / Neutrals / TOF / SC
		locHistName = "NumFCALShowers";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumFCALShowers = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumFCALShowers = new TH1D(locHistName.c_str(), ";# DFCALShower", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumBCALShowers";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumBCALShowers = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumBCALShowers = new TH1D(locHistName.c_str(), ";# DBCALShower", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumNeutralShowers";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumNeutralShowers = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumNeutralShowers = new TH1D(locHistName.c_str(), ";# DNeutralShower", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumTOFPoints";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumTOFPoints = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumTOFPoints = new TH1D(locHistName.c_str(), ";# DTOFPoint", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumSCHits";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumSCHits = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumSCHits = new TH1D(locHistName.c_str(), ";# DSCHit", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumTAGMHits";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumTAGMHits = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumTAGMHits = new TH1D(locHistName.c_str(), ";# DTAGMHit", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumTAGHHits";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumTAGHHits = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumTAGHHits = new TH1D(locHistName.c_str(), ";# DTAGHHit", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		//Matches
		locHistName = "NumTrackBCALMatches";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumTrackBCALMatches = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumTrackBCALMatches = new TH1D(locHistName.c_str(), ";# Track-BCAL Matches", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumTrackFCALMatches";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumTrackFCALMatches = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumTrackFCALMatches = new TH1D(locHistName.c_str(), ";# Track-FCAL Matches", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumTrackTOFMatches";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumTrackTOFMatches = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumTrackTOFMatches = new TH1D(locHistName.c_str(), ";# Track-TOF Matches", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumTrackSCMatches";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumTrackSCMatches = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumTrackSCMatches = new TH1D(locHistName.c_str(), ";# Track-SC Matches", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		if(!locIsRESTEvent)
		{
			//Hits
			locHistName = "NumCDCHits";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumCDCHits = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumCDCHits = new TH1I(locHistName.c_str(), ";# DCDCHit", dMaxNumCDCHits + 1, -0.5, (float)dMaxNumCDCHits + 0.5);

			locHistName = "NumFDCWireHits";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumFDCWireHits = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumFDCWireHits = new TH1I(locHistName.c_str(), ";# Wire DFDCHit", dMaxNumFDCHits/2 + 1, -0.5, (float)dMaxNumFDCHits - 0.5 + 2.0);

			locHistName = "NumFDCCathodeHits";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumFDCCathodeHits = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumFDCCathodeHits = new TH1I(locHistName.c_str(), ";# Cathode DFDCHit", dMaxNumFDCHits/2 + 1, -0.5, (float)dMaxNumFDCHits - 0.5 + 2.0);

			locHistName = "NumTOFHits";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumTOFHits = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumTOFHits = new TH1I(locHistName.c_str(), ";# DTOFHit", dMaxNumTOFCalorimeterHits + 1, -0.5, (float)dMaxNumTOFCalorimeterHits + 0.5);

			locHistName = "NumBCALHits";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumBCALHits = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumBCALHits = new TH1I(locHistName.c_str(), ";# DBCALHit", dMaxNumTOFCalorimeterHits + 1, -0.5, (float)dMaxNumTOFCalorimeterHits + 0.5);

			locHistName = "NumFCALHits";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumFCALHits = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumFCALHits = new TH1I(locHistName.c_str(), ";# DFCALHit", dMaxNumTOFCalorimeterHits + 1, -0.5, (float)dMaxNumTOFCalorimeterHits + 0.5);
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_NumReconstructedObjects::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	if(Get_NumPreviousParticleCombos() != 0)
		return true; //else double-counting!

	bool locIsRESTEvent = (string(locEventLoop->GetJEvent().GetJEventSource()->className()) == string("DEventSourceREST"));

	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locEventLoop->Get(locTrackTimeBasedVector);

	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons);

	vector<const DFCALShower*> locFCALShowers;
	locEventLoop->Get(locFCALShowers);

	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks);

	vector<const DBCALShower*> locBCALShowers;
	locEventLoop->Get(locBCALShowers);

	vector<const DNeutralShower*> locNeutralShowers;
	locEventLoop->Get(locNeutralShowers);

	vector<const DTOFPoint*> locTOFPoints;
	locEventLoop->Get(locTOFPoints);

	vector<const DSCHit*> locSCHits;
	locEventLoop->Get(locSCHits);

	vector<const DTAGHHit*> locTAGHHits;
	locEventLoop->Get(locTAGHHits);

	vector<const DTAGMHit*> locTAGMHits;
	locEventLoop->Get(locTAGMHits);

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	//if not REST
	vector<const DTrackWireBased*> locTrackWireBasedVector;
	vector<const DTrackCandidate*> locTrackCandidates;
	vector<const DTrackCandidate*> locTrackCandidates_CDC;
	vector<const DTrackCandidate*> locTrackCandidates_FDC;
	vector<const DCDCHit*> locCDCHits;
	vector<const DFDCHit*> locFDCHits;
	vector<const DTOFHit*> locTOFHits;
	vector<const DBCALHit*> locBCALHits;
	vector<const DFCALHit*> locFCALHits;

	size_t locNumFDCWireHits = 0, locNumFDCCathodeHits = 0;
	if(!locIsRESTEvent)
	{
		locEventLoop->Get(locTrackWireBasedVector);
		locEventLoop->Get(locTrackCandidates);
		locEventLoop->Get(locTrackCandidates_CDC, "CDC");
		locEventLoop->Get(locTrackCandidates_FDC, "FDCCathodes");
		locEventLoop->Get(locCDCHits);
		locEventLoop->Get(locFDCHits);
		locEventLoop->Get(locTOFHits);
		locEventLoop->Get(locBCALHits);
		locEventLoop->Get(locFCALHits);

		for(size_t loc_i = 0; loc_i < locFDCHits.size(); ++loc_i)
		{
			if(locFDCHits[loc_i]->type == DFDCHit::AnodeWire)
				++locNumFDCWireHits;
			else
				++locNumFDCCathodeHits;
		}
	}

	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//# High-Level Objects
		dHist_NumHighLevelObjects->Fill(1, (Double_t)locTAGMHits.size());
		dHist_NumHighLevelObjects->Fill(2, (Double_t)locTAGHHits.size());
		dHist_NumHighLevelObjects->Fill(3, (Double_t)locSCHits.size());
		dHist_NumHighLevelObjects->Fill(4, (Double_t)locTOFPoints.size());
		dHist_NumHighLevelObjects->Fill(5, (Double_t)locBCALShowers.size());
		dHist_NumHighLevelObjects->Fill(6, (Double_t)locFCALShowers.size());
		dHist_NumHighLevelObjects->Fill(7, (Double_t)locTrackTimeBasedVector.size());
		dHist_NumHighLevelObjects->Fill(8, (Double_t)locDetectorMatches->Get_NumTrackSCMatches());
		dHist_NumHighLevelObjects->Fill(9, (Double_t)locDetectorMatches->Get_NumTrackTOFMatches());
		dHist_NumHighLevelObjects->Fill(10, (Double_t)locDetectorMatches->Get_NumTrackBCALMatches());
		dHist_NumHighLevelObjects->Fill(11, (Double_t)locDetectorMatches->Get_NumTrackFCALMatches());
		dHist_NumHighLevelObjects->Fill(12, (Double_t)locBeamPhotons.size());
		dHist_NumHighLevelObjects->Fill(13, (Double_t)locChargedTracks.size());
		dHist_NumHighLevelObjects->Fill(14, (Double_t)locNeutralShowers.size());

		//Charged
		unsigned int locNumPos = 0, locNumNeg = 0;
		for(size_t loc_i = 0; loc_i < locChargedTracks.size(); ++loc_i)
		{
			if(ParticleCharge(locChargedTracks[loc_i]->Get_BestFOM()->PID()) > 0)
				++locNumPos;
			else
				++locNumNeg;
		}
		dHist_NumChargedTracks->Fill(locChargedTracks.size());
		dHist_NumPosChargedTracks->Fill(locNumPos);
		dHist_NumNegChargedTracks->Fill(locNumNeg);

		//TBT
		locNumPos = 0;  locNumNeg = 0;
		for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
		{
			if(ParticleCharge(locTrackTimeBasedVector[loc_i]->PID()) > 0)
				++locNumPos;
			else
				++locNumNeg;
		}
		dHist_NumTimeBasedTracks->Fill(locTrackTimeBasedVector.size());
		dHist_NumPosTimeBasedTracks->Fill(locNumPos);
		dHist_NumNegTimeBasedTracks->Fill(locNumNeg);

		if(!locIsRESTEvent)
		{
			//WBT
			locNumPos = 0;  locNumNeg = 0;
			for(size_t loc_i = 0; loc_i < locTrackWireBasedVector.size(); ++loc_i)
			{
				if(ParticleCharge(locTrackWireBasedVector[loc_i]->PID()) > 0)
					++locNumPos;
				else
					++locNumNeg;
			}
			dHist_NumWireBasedTracks->Fill(locTrackWireBasedVector.size());
			dHist_NumPosWireBasedTracks->Fill(locNumPos);
			dHist_NumNegWireBasedTracks->Fill(locNumNeg);

			//Candidates
			locNumPos = 0;  locNumNeg = 0;
			for(size_t loc_i = 0; loc_i < locTrackCandidates.size(); ++loc_i)
			{
				if(locTrackCandidates[loc_i]->charge() > 0.0)
					++locNumPos;
				else
					++locNumNeg;
			}
			dHist_NumTrackCandidates->Fill(locTrackCandidates.size());
			dHist_NumPosTrackCandidates->Fill(locNumPos);
			dHist_NumNegTrackCandidates->Fill(locNumNeg);

			//CDC Candidates
			locNumPos = 0;  locNumNeg = 0;
			for(size_t loc_i = 0; loc_i < locTrackCandidates_CDC.size(); ++loc_i)
			{
				if(locTrackCandidates_CDC[loc_i]->charge() > 0.0)
					++locNumPos;
				else
					++locNumNeg;
			}
			dHist_NumPosTrackCandidates_CDC->Fill(locNumPos);
			dHist_NumNegTrackCandidates_CDC->Fill(locNumNeg);

			//FDC Candidates
			locNumPos = 0;  locNumNeg = 0;
			for(size_t loc_i = 0; loc_i < locTrackCandidates_FDC.size(); ++loc_i)
			{
				if(locTrackCandidates_FDC[loc_i]->charge() > 0.0)
					++locNumPos;
				else
					++locNumNeg;
			}
			dHist_NumPosTrackCandidates_FDC->Fill(locNumPos);
			dHist_NumNegTrackCandidates_FDC->Fill(locNumNeg);
		}

		//Beam Photons
		dHist_NumBeamPhotons->Fill((Double_t)locBeamPhotons.size());

		//Showers
		dHist_NumFCALShowers->Fill((Double_t)locFCALShowers.size());
		dHist_NumBCALShowers->Fill((Double_t)locBCALShowers.size());
		dHist_NumNeutralShowers->Fill((Double_t)locNeutralShowers.size());

		//TOF & SC
		dHist_NumTOFPoints->Fill((Double_t)locTOFPoints.size());
		dHist_NumSCHits->Fill((Double_t)locSCHits.size());

		//TAGGER
		dHist_NumTAGMHits->Fill((Double_t)locTAGMHits.size());
		dHist_NumTAGHHits->Fill((Double_t)locTAGHHits.size());

		//Matches
		dHist_NumTrackBCALMatches->Fill((Double_t)locDetectorMatches->Get_NumTrackBCALMatches());
		dHist_NumTrackFCALMatches->Fill((Double_t)locDetectorMatches->Get_NumTrackFCALMatches());
		dHist_NumTrackTOFMatches->Fill((Double_t)locDetectorMatches->Get_NumTrackTOFMatches());
		dHist_NumTrackSCMatches->Fill((Double_t)locDetectorMatches->Get_NumTrackSCMatches());

		//Hits
		if(!locIsRESTEvent)
		{
			dHist_NumCDCHits->Fill((Double_t)locCDCHits.size());
			dHist_NumFDCWireHits->Fill((Double_t)locNumFDCWireHits);
			dHist_NumFDCCathodeHits->Fill((Double_t)locNumFDCCathodeHits);
			dHist_NumTOFHits->Fill((Double_t)locTOFHits.size());
			dHist_NumBCALHits->Fill((Double_t)locBCALHits.size());
			dHist_NumFCALHits->Fill((Double_t)locFCALHits.size());
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!

	return true;
}

void DHistogramAction_TrackMultiplicity::Initialize(JEventLoop* locEventLoop)
{
	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		CreateAndChangeTo_ActionDirectory();

		string locHistName("NumReconstructedParticles");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumReconstructedParticles = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
		else
		{
			dHist_NumReconstructedParticles = new TH2D("NumReconstructedParticles", ";Particle Type;Num Particles / Event", 5 + dFinalStatePIDs.size(), -0.5, 4.5 + dFinalStatePIDs.size(), dMaxNumTracks + 1, -0.5, (float)dMaxNumTracks + 0.5);
			dHist_NumReconstructedParticles->GetXaxis()->SetBinLabel(1, "# Total");
			dHist_NumReconstructedParticles->GetXaxis()->SetBinLabel(2, "# q != 0");
			dHist_NumReconstructedParticles->GetXaxis()->SetBinLabel(3, "# q = 0");
			dHist_NumReconstructedParticles->GetXaxis()->SetBinLabel(4, "# q = +");
			dHist_NumReconstructedParticles->GetXaxis()->SetBinLabel(5, "# q = -");
			for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
			{
				string locLabelName = string("# ") + string(ParticleName_ROOT(dFinalStatePIDs[loc_i]));
				dHist_NumReconstructedParticles->GetXaxis()->SetBinLabel(6 + loc_i, locLabelName.c_str());
			}
		}

		locHistName = "NumGoodReconstructedParticles";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumGoodReconstructedParticles = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
		else
		{
			dHist_NumGoodReconstructedParticles = new TH2D("NumGoodReconstructedParticles", ";Particle Type;Num Particles / Event", 5 + dFinalStatePIDs.size(), -0.5, 4.5 + dFinalStatePIDs.size(), dMaxNumTracks + 1, -0.5, (float)dMaxNumTracks + 0.5);
			dHist_NumGoodReconstructedParticles->GetXaxis()->SetBinLabel(1, "# Total");
			dHist_NumGoodReconstructedParticles->GetXaxis()->SetBinLabel(2, "# q != 0");
			dHist_NumGoodReconstructedParticles->GetXaxis()->SetBinLabel(3, "# q = 0");
			dHist_NumGoodReconstructedParticles->GetXaxis()->SetBinLabel(4, "# q = +");
			dHist_NumGoodReconstructedParticles->GetXaxis()->SetBinLabel(5, "# q = -");
			for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
			{
				string locLabelName = string("# ") + string(ParticleName_ROOT(dFinalStatePIDs[loc_i]));
				dHist_NumGoodReconstructedParticles->GetXaxis()->SetBinLabel(6 + loc_i, locLabelName.c_str());
			}
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_TrackMultiplicity::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	if(Get_NumPreviousParticleCombos() != 0)
		return true; //else double-counting!

	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks);

	vector<const DChargedTrack*> locGoodChargedTracks;
	locEventLoop->Get(locGoodChargedTracks, "PreSelect");

	// get #tracks by PID/q type 
	size_t locNumPositiveTracks = 0, locNumNegativeTracks = 0; 
	map<Particle_t, size_t> locNumTracksByPID;
	for(size_t loc_i = 0; loc_i < locChargedTracks.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTracks[loc_i]->Get_BestFOM();
		Particle_t locPID = locChargedTrackHypothesis->PID();

		if(locChargedTrackHypothesis->charge() > 0.0)
			++locNumPositiveTracks;
		else
			++locNumNegativeTracks;

		if(locNumTracksByPID.find(locPID) != locNumTracksByPID.end())
			++locNumTracksByPID[locPID];
		else
			locNumTracksByPID[locPID] = 1;
	}

	// get # good tracks by PID/q type 
	size_t locNumGoodPositiveTracks = 0, locNumGoodNegativeTracks = 0;
	map<Particle_t, size_t> locNumGoodTracksByPID;
	for(size_t loc_i = 0; loc_i < locGoodChargedTracks.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locGoodChargedTracks[loc_i]->Get_BestFOM();
		Particle_t locPID = locChargedTrackHypothesis->PID();

		double locPIDFOM = locChargedTrackHypothesis->dFOM;

		if(locChargedTrackHypothesis->charge() > 0.0)
			++locNumGoodPositiveTracks;
		else
			++locNumGoodNegativeTracks;

		if(locPIDFOM < dMinPIDFOM)
			continue;

		if(locNumGoodTracksByPID.find(locPID) != locNumGoodTracksByPID.end())
			++locNumGoodTracksByPID[locPID];
		else
			locNumGoodTracksByPID[locPID] = 1;
	}

	vector<const DNeutralParticle*> locNeutralParticles;
	locEventLoop->Get(locNeutralParticles);

	vector<const DNeutralShower*> locGoodNeutralShowers;
	locEventLoop->Get(locGoodNeutralShowers, "PreSelect");

	// neutrals by pid
	for(size_t loc_i = 0; loc_i < locNeutralParticles.size(); ++loc_i)
	{
		const DNeutralParticleHypothesis* locNeutralParticleHypothesis = locNeutralParticles[loc_i]->Get_BestFOM();
		if(locNeutralParticleHypothesis->dFOM < dMinPIDFOM)
			continue;

		Particle_t locPID = locNeutralParticleHypothesis->PID();
		if(locNumTracksByPID.find(locPID) != locNumTracksByPID.end())
			++locNumTracksByPID[locPID];
		else
			locNumTracksByPID[locPID] = 1;
	}

	// good neutrals
	size_t locNumGoodNeutrals = 0;
	for(size_t loc_i = 0; loc_i < locGoodNeutralShowers.size(); ++loc_i)
	{
		const DNeutralParticle* locNeutralParticle = NULL;
		for(size_t loc_j = 0; loc_j < locNeutralParticles.size(); ++loc_j)
		{
			if(locNeutralParticles[loc_j]->dNeutralShower != locGoodNeutralShowers[loc_i])
				continue;
			locNeutralParticle = locNeutralParticles[loc_j];
			break;
		}
		if(locNeutralParticle == NULL)
			continue;
		++locNumGoodNeutrals;

		const DNeutralParticleHypothesis* locNeutralParticleHypothesis = locNeutralParticles[loc_i]->Get_BestFOM();
		if(locNeutralParticleHypothesis->dFOM < dMinPIDFOM)
			continue;

		Particle_t locPID = locNeutralParticleHypothesis->PID();
		if(locNumGoodTracksByPID.find(locPID) != locNumGoodTracksByPID.end())
			++locNumGoodTracksByPID[locPID];
		else
			locNumGoodTracksByPID[locPID] = 1;
	}

	size_t locNumGoodTracks = locNumGoodPositiveTracks + locNumGoodNegativeTracks;
	japp->RootWriteLock();
	{
		dHist_NumReconstructedParticles->Fill(0.0, (Double_t)(locChargedTracks.size() + locNeutralParticles.size()));
		dHist_NumReconstructedParticles->Fill(1.0, (Double_t)locChargedTracks.size());
		dHist_NumReconstructedParticles->Fill(2.0, (Double_t)locNeutralParticles.size());
		dHist_NumReconstructedParticles->Fill(3.0, (Double_t)locNumPositiveTracks);
		dHist_NumReconstructedParticles->Fill(4.0, (Double_t)locNumNegativeTracks);
		for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
			dHist_NumReconstructedParticles->Fill(5.0 + (Double_t)loc_i, (Double_t)locNumTracksByPID[dFinalStatePIDs[loc_i]]);

		dHist_NumGoodReconstructedParticles->Fill(0.0, (Double_t)(locNumGoodTracks + locNumGoodNeutrals));
		dHist_NumGoodReconstructedParticles->Fill(1.0, (Double_t)locNumGoodTracks);
		dHist_NumGoodReconstructedParticles->Fill(2.0, (Double_t)locNumGoodNeutrals);
		dHist_NumGoodReconstructedParticles->Fill(3.0, (Double_t)locNumGoodPositiveTracks);
		dHist_NumGoodReconstructedParticles->Fill(4.0, (Double_t)locNumGoodNegativeTracks);
		for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
			dHist_NumGoodReconstructedParticles->Fill(5.0 + (Double_t)loc_i, (Double_t)locNumGoodTracksByPID[dFinalStatePIDs[loc_i]]);
	}
	japp->RootUnLock();

	return true;
}

void DHistogramAction_TruePID::Initialize(JEventLoop* locEventLoop)
{
	string locStepName, locStepROOTName, locHistTitle, locHistName, locParticleName, locParticleROOTName;
	Particle_t locPID;

	size_t locNumSteps = Get_Reaction()->Get_NumReactionSteps();
	dHistDeque_P_CorrectID.resize(locNumSteps);
	dHistDeque_P_IncorrectID.resize(locNumSteps);
	dHistDeque_PVsTheta_CorrectID.resize(locNumSteps);
	dHistDeque_PVsTheta_IncorrectID.resize(locNumSteps);

	deque<deque<Particle_t> > locDetectedPIDs;
	Get_Reaction()->Get_DetectedFinalPIDs(locDetectedPIDs);

	vector<const DAnalysisUtilities*> locAnalysisUtilitiesVector;
	locEventLoop->Get(locAnalysisUtilitiesVector);

	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		CreateAndChangeTo_ActionDirectory();
		dAnalysisUtilities = locAnalysisUtilitiesVector[0];
		for(size_t loc_i = 0; loc_i < locNumSteps; ++loc_i)
		{
			if(locDetectedPIDs[loc_i].empty())
				continue;

			const DReactionStep* locReactionStep = Get_Reaction()->Get_ReactionStep(loc_i);
			locStepName = locReactionStep->Get_StepName();
			locStepROOTName = locReactionStep->Get_StepROOTName();
			CreateAndChangeTo_Directory(locStepName, locStepName);

			for(size_t loc_j = 0; loc_j < locDetectedPIDs[loc_i].size(); ++loc_j)
			{
				locPID = locDetectedPIDs[loc_i][loc_j];
				locParticleName = ParticleType(locPID);
				locParticleROOTName = ParticleName_ROOT(locPID);

				if(dHistDeque_P_CorrectID[loc_i].find(locPID) != dHistDeque_P_CorrectID[loc_i].end())
					continue; //hists already created for this pid

				//P of Correct ID
				locHistName = string("Momentum_CorrectID_") + locParticleName;
				locHistTitle = string("Correct ") + locParticleROOTName + string(" ID, ") + locStepROOTName + string(";p (GeV/c)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_P_CorrectID[loc_i][locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_P_CorrectID[loc_i][locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPBins, dMinP, dMaxP);

				//P of Incorrect ID
				locHistName = string("Momentum_IncorrectID_") + locParticleName;
				locHistTitle = string("Incorrect ") + locParticleROOTName + string(" ID, ") + locStepROOTName + string(";p (GeV/c)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_P_IncorrectID[loc_i][locPID] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_P_IncorrectID[loc_i][locPID] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPBins, dMinP, dMaxP);

				//P Vs Theta of Correct ID
				locHistName = string("PVsTheta_CorrectID_") + locParticleName;
				locHistTitle = string("Correct ") + locParticleROOTName + string(" ID, ") + locStepROOTName + string(";#theta#circ;p (GeV/c)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_PVsTheta_CorrectID[loc_i][locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_PVsTheta_CorrectID[loc_i][locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNumThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

				//P Vs Theta of Incorrect ID
				locHistName = string("PVsTheta_IncorrectID_") + locParticleName;
				locHistTitle = string("Incorrect ") + locParticleROOTName + string(" ID, ") + locStepROOTName + string(";#theta#circ;p (GeV/c)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_PVsTheta_IncorrectID[loc_i][locPID] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_PVsTheta_IncorrectID[loc_i][locPID] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNumThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);
			}
			gDirectory->cd("..");
		} //end of step loop

		//# Combos Pass/Fail All True PID
		locHistName = "Combo_TruePIDStatus";
		locHistTitle = Get_Reaction()->Get_ReactionName() + string(";# Combos;All Combo Particles True PID Status");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_TruePIDStatus = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_TruePIDStatus = new TH1I(locHistName.c_str(), locHistTitle.c_str(), 2, -0.5, 1.5);
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_TruePID::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
	locEventLoop->Get(locMCThrownMatchingVector);
	if(locMCThrownMatchingVector.empty())
		return true;
	const DMCThrownMatching* locMCThrownMatching = locMCThrownMatchingVector[0];
	double locP, locTheta;
	const DMCThrown* locMCThrown;
	Particle_t locPID;

	if(Get_NumPreviousParticleCombos() == 0)
		dPreviouslyHistogrammedParticles.clear();

	deque<const DKinematicData*> locParticles;
	int locComboTruePIDStatus = 1;
	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);
		locParticleComboStep->Get_FinalParticles_Measured(locParticles);

		for(size_t loc_j = 0; loc_j < locParticles.size(); ++loc_j)
		{
			if(!locParticleComboStep->Is_FinalParticleDetected(loc_j))
				continue;
			locPID = locParticles[loc_j]->PID();

			double locMatchFOM = 0.0;
			if(ParticleCharge(locPID) == 0)
				locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(static_cast<const DNeutralParticleHypothesis*>(locParticles[loc_j]), locMatchFOM);
			else
				locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(static_cast<const DChargedTrackHypothesis*>(locParticles[loc_j]), locMatchFOM);

			bool locCutResult = ((locMCThrown == NULL) || (locMatchFOM < dMinThrownMatchFOM)) ? false : (((Particle_t)locMCThrown->type) == locPID);
			if(!locCutResult)
				locComboTruePIDStatus = 0;

			//check if duplicate
			const JObject* locSourceObject = locParticleComboStep->Get_FinalParticle_SourceObject(loc_j);
			pair<Particle_t, const JObject*> locParticleInfo(locParticles[loc_j]->PID(), locSourceObject);
			pair<size_t, pair<Particle_t, const JObject*> > locHistInfo(loc_i, locParticleInfo);
			if(dPreviouslyHistogrammedParticles.find(locHistInfo) != dPreviouslyHistogrammedParticles.end())
				continue; //previously histogrammed
			dPreviouslyHistogrammedParticles.insert(locHistInfo);

			locP = locParticles[loc_j]->momentum().Mag();
			locTheta = locParticles[loc_j]->momentum().Theta()*180.0/TMath::Pi();

			japp->RootWriteLock();
			{
				if(locCutResult)
				{
					dHistDeque_P_CorrectID[loc_i][locPID]->Fill(locP);
					dHistDeque_PVsTheta_CorrectID[loc_i][locPID]->Fill(locTheta, locP);
				}
				else
				{
					dHistDeque_P_IncorrectID[loc_i][locPID]->Fill(locP);
					dHistDeque_PVsTheta_IncorrectID[loc_i][locPID]->Fill(locTheta, locP);
				}
			}
			japp->RootUnLock();
		}
	}
	japp->RootWriteLock();
	{
		dHist_TruePIDStatus->Fill(locComboTruePIDStatus);
	}
	japp->RootUnLock();

	return true;
}

void DHistogramAction_InvariantMass::Initialize(JEventLoop* locEventLoop)
{
	string locHistName, locHistTitle;
	double locMassPerBin = 1000.0*(dMaxMass - dMinMass)/((double)dNumMassBins);

	string locParticleNamesForHist = Get_Reaction()->Get_DecayChainFinalParticlesROOTNames(dInitialPID, Get_UseKinFitResultsFlag());

	vector<const DAnalysisUtilities*> locAnalysisUtilitiesVector;
	locEventLoop->Get(locAnalysisUtilitiesVector);

	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		dAnalysisUtilities = locAnalysisUtilitiesVector[0];
		CreateAndChangeTo_ActionDirectory();

		locHistName = "InvariantMass";
		ostringstream locStream;
		locStream << locMassPerBin;
		locHistTitle = string(";") + locParticleNamesForHist + string(" Invariant Mass (GeV/c^{2});# Combos / ") + locStream.str() + string(" MeV/c^{2}");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_InvaraintMass = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_InvaraintMass = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumMassBins, dMinMass, dMaxMass);
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_InvariantMass::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	if(Get_NumPreviousParticleCombos() == 0)
		dPreviousSourceObjects.clear();

	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);
		if(locParticleComboStep->Get_InitialParticleID() != dInitialPID)
			continue;

		set<pair<const JObject*, Particle_t> > locSourceObjects;
		DLorentzVector locFinalStateP4 = dAnalysisUtilities->Calc_FinalStateP4(locParticleCombo, loc_i, locSourceObjects, Get_UseKinFitResultsFlag());

		if(!dEnableDoubleCounting)
		{
			if(dPreviousSourceObjects.find(locSourceObjects) != dPreviousSourceObjects.end())
				return true; //dupe: already histed!
			dPreviousSourceObjects.insert(locSourceObjects);
		}

		//get particle names so can select the correct histogram
		double locInvariantMass = locFinalStateP4.M();
		japp->RootWriteLock();
		{
			dHist_InvaraintMass->Fill(locInvariantMass);
		}
		japp->RootUnLock();
		//don't break: e.g. if multiple pi0's, histogram invariant mass of each one
	}
	return true;
}

void DHistogramAction_MissingMass::Initialize(JEventLoop* locEventLoop)
{
	double locMassPerBin = 1000.0*(dMaxMass - dMinMass)/((double)dNumMassBins);
	string locInitialParticlesROOTName = Get_Reaction()->Get_InitialParticlesROOTName();
	string locFinalParticlesROOTName = Get_Reaction()->Get_DecayChainFinalParticlesROOTNames(0, dMissingMassOffOfStepIndex, dMissingMassOffOfPIDs, Get_UseKinFitResultsFlag(), true);

	vector<const DAnalysisUtilities*> locAnalysisUtilitiesVector;
	locEventLoop->Get(locAnalysisUtilitiesVector);

	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		dAnalysisUtilities = locAnalysisUtilitiesVector[0];
		CreateAndChangeTo_ActionDirectory();
		string locHistName = "MissingMass";
		ostringstream locStream;
		locStream << locMassPerBin;
		string locHistTitle = string(";") + locInitialParticlesROOTName + string("#rightarrow") + locFinalParticlesROOTName + string(" Missing Mass (GeV/c^{2});# Combos / ") + locStream.str() + string(" MeV/c^{2}");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_MissingMass = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_MissingMass = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumMassBins, dMinMass, dMaxMass);
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_MissingMass::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	if(Get_NumPreviousParticleCombos() == 0)
		dPreviousSourceObjects.clear();

	set<pair<const JObject*, Particle_t> > locSourceObjects;
	DLorentzVector locMissingP4;

	if(dMissingMassOffOfStepIndex == -1)
		locMissingP4 = dAnalysisUtilities->Calc_MissingP4(locParticleCombo, locSourceObjects, Get_UseKinFitResultsFlag());
	else
		locMissingP4 = dAnalysisUtilities->Calc_MissingP4(locParticleCombo, 0, dMissingMassOffOfStepIndex, dMissingMassOffOfPIDs, locSourceObjects, Get_UseKinFitResultsFlag());

	if(!dEnableDoubleCounting)
	{
		if(dPreviousSourceObjects.find(locSourceObjects) != dPreviousSourceObjects.end())
			return true; //dupe: already histed!
		dPreviousSourceObjects.insert(locSourceObjects);
	}

	double locMissingMass = locMissingP4.M();
	japp->RootWriteLock();
	{
		dHist_MissingMass->Fill(locMissingMass);
	}
	japp->RootUnLock();

	return true;
}

void DHistogramAction_MissingMassSquared::Initialize(JEventLoop* locEventLoop)
{
	double locMassSqPerBin = 1000.0*1000.0*(dMaxMassSq - dMinMassSq)/((double)dNumMassBins);
	string locInitialParticlesROOTName = Get_Reaction()->Get_InitialParticlesROOTName();
	string locFinalParticlesROOTName = Get_Reaction()->Get_DecayChainFinalParticlesROOTNames(0, dMissingMassOffOfStepIndex, dMissingMassOffOfPIDs, Get_UseKinFitResultsFlag(), true);

	vector<const DAnalysisUtilities*> locAnalysisUtilitiesVector;
	locEventLoop->Get(locAnalysisUtilitiesVector);

	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		dAnalysisUtilities = locAnalysisUtilitiesVector[0];
		CreateAndChangeTo_ActionDirectory();
		string locHistName = "MissingMassSquared";
		ostringstream locStream;
		locStream << locMassSqPerBin;
		string locHistTitle = string(";") + locInitialParticlesROOTName + string("#rightarrow") + locFinalParticlesROOTName + string(" Missing Mass Squared (GeV/c^{2})^{2};# Combos / ") + locStream.str() + string(" (MeV/c^{2})^{2}");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_MissingMassSquared = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_MissingMassSquared = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumMassBins, dMinMassSq, dMaxMassSq);
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_MissingMassSquared::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	if(Get_NumPreviousParticleCombos() == 0)
		dPreviousSourceObjects.clear();

	set<pair<const JObject*, Particle_t> > locSourceObjects;
	DLorentzVector locMissingP4;

	if(dMissingMassOffOfStepIndex == -1)
		locMissingP4 = dAnalysisUtilities->Calc_MissingP4(locParticleCombo, locSourceObjects, Get_UseKinFitResultsFlag());
	else
		locMissingP4 = dAnalysisUtilities->Calc_MissingP4(locParticleCombo, 0, dMissingMassOffOfStepIndex, dMissingMassOffOfPIDs, locSourceObjects, Get_UseKinFitResultsFlag());

	if(!dEnableDoubleCounting)
	{
		if(dPreviousSourceObjects.find(locSourceObjects) != dPreviousSourceObjects.end())
			return true; //dupe: already histed!
		dPreviousSourceObjects.insert(locSourceObjects);
	}

	double locMissingMassSquared = locMissingP4.M2();
	japp->RootWriteLock();
	{
		dHist_MissingMassSquared->Fill(locMissingMassSquared);
	}
	japp->RootUnLock();

	return true;
}

void DHistogramAction_KinFitResults::Initialize(JEventLoop* locEventLoop)
{
	string locHistName, locHistTitle, locParticleName, locParticleROOTName, locStepName, locStepROOTName;
	Particle_t locPID;

	DKinFitType locKinFitType = Get_Reaction()->Get_KinFitType();
	if(locKinFitType == d_NoFit)
		return;

	deque<deque<Particle_t> > locDetectedPIDs;
	Get_Reaction()->Get_DetectedFinalPIDs(locDetectedPIDs);

	vector<const DAnalysisUtilities*> locAnalysisUtilitiesVector;
	locEventLoop->Get(locAnalysisUtilitiesVector);

	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	dAnalysisUtilities = locAnalysisUtilitiesVector[0];
	CreateAndChangeTo_ActionDirectory();

	string locKinFitTypeString;
	if(locKinFitType == d_P4Fit)
		locKinFitTypeString = "P4";
	else if(locKinFitType == d_VertexFit)
		locKinFitTypeString = "Vertex";
	else if(locKinFitType == d_SpacetimeFit)
		locKinFitTypeString = "Spacetime";
	else if(locKinFitType == d_P4AndVertexFit)
		locKinFitTypeString = "P4 & Vertex";
	else if(locKinFitType == d_P4AndSpacetimeFit)
		locKinFitTypeString = "P4 & Spacetime";

	// Confidence Level
	locHistName = "ConfidenceLevel";
	locHistTitle = locKinFitTypeString + string(" Kinematic Fit;Confidence Level;# Combos");
	if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
		dHist_ConfidenceLevel = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
	else
		dHist_ConfidenceLevel = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumConfidenceLevelBins, 0.0, 1.0);

	// Pulls
	map<DKinFitPullType, TH1I*> locParticlePulls;

	//beam pulls
	bool locBeamFlag = (Get_Reaction()->Get_ReactionStep(0)->Get_InitialParticleID() == Gamma);
	if(locBeamFlag)
	{
		CreateAndChangeTo_Directory("Beam", "Beam");
		Create_ParticlePulls(true, "Beam", Gamma, dHistMap_BeamPulls, locKinFitTypeString);
		gDirectory->cd("..");
	}

	//final particle pulls
	for(size_t loc_i = 0; loc_i < Get_Reaction()->Get_NumReactionSteps(); ++loc_i)
	{
		const DReactionStep* locReactionStep = Get_Reaction()->Get_ReactionStep(loc_i);
		locStepName = locReactionStep->Get_StepName();
		locStepROOTName = locReactionStep->Get_StepROOTName();

		if(locDetectedPIDs[loc_i].empty())
			continue;

		CreateAndChangeTo_Directory(locStepName, locStepName);

		for(size_t loc_j = 0; loc_j < locDetectedPIDs[loc_i].size(); ++loc_j)
		{
			locPID = locDetectedPIDs[loc_i][loc_j];
			locParticleName = ParticleType(locPID);
			CreateAndChangeTo_Directory(locParticleName, locParticleName);

			Create_ParticlePulls(false, locStepROOTName, locPID, locParticlePulls, locKinFitTypeString);
			dHistMap_Pulls[pair<size_t, Particle_t>(loc_i, locPID)] = locParticlePulls;

			gDirectory->cd("..");
		} //end of particle loop
		gDirectory->cd("..");
	} //end of step loop

	//RF Time Pull
	if(locBeamFlag && ((locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndSpacetimeFit)))
	{
		CreateAndChangeTo_Directory("RF", "RF");

		//T Pull
		locHistName = "Pull_RF_T";
		locHistTitle = string("RF Bunch, ") + locKinFitTypeString + string(" Fit;t Pull;# Combos");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_RFTimePull = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_RFTimePull = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, dMinPull, dMaxPull);

		gDirectory->cd("..");
	}

	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

void DHistogramAction_KinFitResults::Create_ParticlePulls(bool locIsBeamFlag, string locStepROOTName, Particle_t locPID, map<DKinFitPullType, TH1I*>& locParticlePulls, const string& locKinFitTypeString)
{
	string locHistName, locHistTitle, locParticleName, locParticleROOTName;
	locParticleName = ParticleType(locPID);
	locParticleROOTName = ParticleName_ROOT(locPID);

	DKinFitType locKinFitType = Get_Reaction()->Get_KinFitType();

	locParticlePulls.clear();

	bool locNeutralShowerFlag = ((ParticleCharge(locPID) == 0) && (!locIsBeamFlag) && (locKinFitType != d_P4Fit));
	//p4 pulls:
	if(locNeutralShowerFlag)
	{
		//neutral shower not in a p4-only fit
		//E Pull
		locHistName = "Pull_E";
		locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(", ") + locKinFitTypeString + string(" Fit;E Pull;# Combos");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			locParticlePulls[d_EPull] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			locParticlePulls[d_EPull] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, dMinPull, dMaxPull);
	}
	else
	{
		//Px Pull
		locHistName = "Pull_Px";
		locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(", ") + locKinFitTypeString + string(" Fit;p_{x} Pull;# Combos");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			locParticlePulls[d_PxPull] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			locParticlePulls[d_PxPull] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, dMinPull, dMaxPull);

		//Py Pull
		locHistName = "Pull_Py";
		locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(", ") + locKinFitTypeString + string(" Fit;p_{y} Pull;# Combos");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			locParticlePulls[d_PyPull] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			locParticlePulls[d_PyPull] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, dMinPull, dMaxPull);

		//Pz Pull
		locHistName = "Pull_Pz";
		locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(", ") + locKinFitTypeString + string(" Fit;p_{z} Pull;# Combos");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			locParticlePulls[d_PzPull] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			locParticlePulls[d_PzPull] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, dMinPull, dMaxPull);
	}

	//vertex pulls:
	if(locNeutralShowerFlag || (locKinFitType == d_VertexFit) || (locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit))
	{
		//Xx Pull
		locHistName = "Pull_Xx";
		locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(", ") + locKinFitTypeString + string(" Fit;x_{x} Pull;# Combos");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			locParticlePulls[d_XxPull] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			locParticlePulls[d_XxPull] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, dMinPull, dMaxPull);

		//Xy Pull
		locHistName = "Pull_Xy";
		locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(", ") + locKinFitTypeString + string(" Fit;x_{y} Pull;# Combos");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			locParticlePulls[d_XyPull] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			locParticlePulls[d_XyPull] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, dMinPull, dMaxPull);

		//Xz Pull
		locHistName = "Pull_Xz";
		locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(", ") + locKinFitTypeString + string(" Fit;x_{z} Pull;# Combos");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			locParticlePulls[d_XzPull] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			locParticlePulls[d_XzPull] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, dMinPull, dMaxPull);
	}

	//time pulls:
	if((locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndSpacetimeFit))
	{
		//T Pull
		locHistName = "Pull_T";
		locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(", ") + locKinFitTypeString + string(" Fit;t Pull;# Combos");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			locParticlePulls[d_TPull] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
		else
			locParticlePulls[d_TPull] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, dMinPull, dMaxPull);
	}
}

bool DHistogramAction_KinFitResults::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//kinfit results are unique for each DParticleCombo: no need to check for duplicates
	const DKinFitResults* locKinFitResults = locParticleCombo->Get_KinFitResults();
	if(locKinFitResults == NULL)
		return true;

	// Confidence Level
	double locConfidenceLevel = locKinFitResults->Get_ConfidenceLevel();
	japp->RootWriteLock();
	{
		if(string(dHist_ConfidenceLevel->GetXaxis()->GetTitle()) == string("Confidence Level"))
		{

			deque<const DKinFitConstraint*> locKinFitConstraints;
			locKinFitResults->Get_KinFitConstraints(locKinFitConstraints);

			string locHistTitle = "Kinematic Fit Constraints: ";
			bool locFirstConstraintFlag = true;
			for(size_t loc_i = 0; loc_i < locKinFitConstraints.size(); ++loc_i)
			{
				string locConstraintString = locKinFitConstraints[loc_i]->Get_ConstraintString();
				if(locConstraintString == "")
					continue;
				if(!locFirstConstraintFlag)
					locHistTitle += ", ";
				locFirstConstraintFlag = false;
				locHistTitle += locConstraintString;
			}
			dHist_ConfidenceLevel->SetTitle(locHistTitle.c_str());
			ostringstream locAxisTitle;
			locAxisTitle << "Confidence Level (" << locKinFitResults->Get_NumConstraints() << " Constraints, " << locKinFitResults->Get_NumUnknowns();
			locAxisTitle << " Unknowns: " << locKinFitResults->Get_NDF() << "-C Fit)";
			dHist_ConfidenceLevel->GetXaxis()->SetTitle(locAxisTitle.str().c_str());
		}
		dHist_ConfidenceLevel->Fill(locConfidenceLevel);
	}
	japp->RootUnLock();
	if(locConfidenceLevel < dPullHistConfidenceLevelCut)
		return true; //don't histogram pulls

	// Pulls
	map<const DKinematicData*, map<DKinFitPullType, double> > locPulls; //DKinematicData is the MEASURED particle
	locKinFitResults->Get_Pulls(locPulls);
	deque<const DKinematicData*> locParticles;
	map<DKinFitPullType, double> locParticlePulls;
	map<DKinFitPullType, double>::iterator locIterator;
	const DKinematicData* locKinematicData;

	// beam pulls
	bool locBeamFlag = (Get_Reaction()->Get_ReactionStep(0)->Get_InitialParticleID() == Gamma);
	if(locBeamFlag)
	{
		locKinematicData = locParticleCombo->Get_ParticleComboStep(0)->Get_InitialParticle_Measured();
		locParticlePulls = locPulls[locKinematicData];
		japp->RootWriteLock();
		{
			for(locIterator = locParticlePulls.begin(); locIterator != locParticlePulls.end(); ++locIterator)
				dHistMap_BeamPulls[locIterator->first]->Fill(locIterator->second);
		}
		japp->RootUnLock();
	}

	// final particle pulls
	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);
		locParticleComboStep->Get_DetectedFinalParticles_Measured(locParticles);
		for(size_t loc_j = 0; loc_j < locParticles.size(); ++loc_j)
		{
			locParticlePulls = locPulls[locParticles[loc_j]];
			pair<size_t, Particle_t> locParticlePair(loc_i, locParticles[loc_j]->PID());
			japp->RootWriteLock();
			{
				for(locIterator = locParticlePulls.begin(); locIterator != locParticlePulls.end(); ++locIterator)
					(dHistMap_Pulls[locParticlePair])[locIterator->first]->Fill(locIterator->second);
			}
			japp->RootUnLock();
		}
	}

	//rf time pull
	DKinFitType locKinFitType = Get_Reaction()->Get_KinFitType();
	if(locBeamFlag && ((locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndSpacetimeFit)))
	{
		locParticlePulls = locPulls[NULL];
		japp->RootWriteLock();
		{
			dHist_RFTimePull->Fill(locParticlePulls[d_TPull]);
		}
		japp->RootUnLock();
	}

	return true;
}

