#include "ANALYSIS/DHistogramActions.h"

void DHistogramAction_NumParticleCombos::Initialize(JEventLoop* locEventLoop)
{
	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		CreateAndChangeTo_ActionDirectory();
		if(gDirectory->Get("NumParticleCombos") != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumParticleCombos = static_cast<TH1D*>(gDirectory->Get("NumParticleCombos"));
		else
			dHist_NumParticleCombos = new TH1D("NumParticleCombos", ";# Particle Combos;Events", dNumComboBins, 0.5, 0.5 + float(dNumComboBins));
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_NumParticleCombos::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	size_t locNumPreviousParticleCombos = Get_NumPreviousParticleCombos();

	//currently not an ideal/clean way to do this, but this works
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		Double_t locNumEntries = dHist_NumParticleCombos->GetEntries();
		if(locNumPreviousParticleCombos != 0) //undo previous fill
		{
			Int_t locPreviousBin = dHist_NumParticleCombos->FindBin(locNumPreviousParticleCombos);
			Double_t locPreviousContent = dHist_NumParticleCombos->GetBinContent(locPreviousBin);
			dHist_NumParticleCombos->SetBinContent(locPreviousBin, locPreviousContent - 1.0);
		}
		dHist_NumParticleCombos->Fill(locNumPreviousParticleCombos + 1);
		dHist_NumParticleCombos->SetEntries(locNumEntries + 1);
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!

	return true;
}

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
				dHistMap_PIDFOM[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PIDFOM[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumFOMBins, 0.0, 1.0);

			// Confidence Level in BCAL
			locHistName = "TOFConfidenceLevel_BCAL";
			locHistTitle = locParticleROOTName + string(" PID in BCAL;TOF Confidence Level");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TOFFOM_BCAL[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TOFFOM_BCAL[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumFOMBins, 0.0, 1.0);

			// Confidence Level in FCAL
			locHistName = "TOFConfidenceLevel_FCAL";
			locHistTitle = locParticleROOTName + string(" PID in FCAL;TOF Confidence Level");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TOFFOM_FCAL[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TOFFOM_FCAL[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumFOMBins, 0.0, 1.0);

			if(ParticleCharge(locPID) != 0)
			{
				// Confidence Level in TOF
				locHistName = "TOFConfidenceLevel_TOF";
				locHistTitle = locParticleROOTName + string(" PID in TOF;TOF Confidence Level");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TOFFOM_TOF[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TOFFOM_TOF[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumFOMBins, 0.0, 1.0);

				// Confidence Level in CDC
				locHistName = "TOFConfidenceLevel_CDC";
				locHistTitle = locParticleROOTName + string(" PID in CDC;TOF Confidence Level");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TOFFOM_CDC[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TOFFOM_CDC[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumFOMBins, 0.0, 1.0);

				// DC dE/dx Confidence Level
				locHistName = "DCdEdxConfidenceLevel";
				locHistTitle = locParticleROOTName + string(" PID;DC #it{#frac{dE}{dx}} Confidence Level");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_DCdEdxFOM[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_DCdEdxFOM[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumFOMBins, 0.0, 1.0);
			}

			//beta vs p
			locHistName = "BetaVsP";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#beta");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_BetaVsP[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_BetaVsP[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumBetaBins, dMinBeta, dMaxBeta);

			//delta-beta vs p
			locHistName = "DeltaBetaVsP";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta#beta");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaBetaVsP[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaBetaVsP[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaBetaBins, dMinDeltaBeta, dMaxDeltaBeta);

			//TOF Confidence Level vs Delta-Beta
			locHistName = "TOFConfidenceLevelVsDeltaBeta";
			locHistTitle = locParticleROOTName + string(";#Delta#beta;TOF Confidence Level");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TOFFOMVsDeltaBeta[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TOFFOMVsDeltaBeta[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaBetaBins, dMinDeltaBeta, dMaxDeltaBeta, dNumFOMBins, 0.0, 1.0);

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
						dHistMap_PIDFOMForTruePID[locPIDPair] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistMap_PIDFOMForTruePID[locPIDPair] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumFOMBins, 0.0, 1.0);
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
					dHistMap_TimePull_CDC[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePull_CDC[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

				locHistName = "TimePullVsTheta_CDC";
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TimePullVsTheta_CDC[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePullVsTheta_CDC[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);

				locHistName = "TimePullVsP_CDC";
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TimePullVsP_CDC[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePullVsP_CDC[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);
			}

			//Delta-t Pulls - BCAL
			locHistName = "TimePull_BCAL";
			locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TimePull_BCAL[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TimePull_BCAL[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

			locHistName = "TimePullVsTheta_BCAL";
			locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat/#sigma_{#Deltat}");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TimePullVsTheta_BCAL[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TimePullVsTheta_BCAL[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);

			locHistName = "TimePullVsP_BCAL";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TimePullVsP_BCAL[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TimePullVsP_BCAL[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

			//Delta-t Pulls - TOF
			if(ParticleCharge(locPID) != 0) //TOF
			{
				locHistName = "TimePull_TOF";
				locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TimePull_TOF[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePull_TOF[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

				locHistName = "TimePullVsP_TOF";
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TimePullVsP_TOF[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePullVsP_TOF[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);
			}

			//Delta-t Pulls - FCAL
			locHistName = "TimePull_FCAL";
			locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TimePull_FCAL[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TimePull_FCAL[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

			locHistName = "TimePullVsP_FCAL";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TimePullVsP_FCAL[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TimePullVsP_FCAL[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

			gDirectory->cd("..");
			CreateAndChangeTo_Directory("PVsTheta", "PVsTheta");

			// P Vs Theta, PID FOM < 1%
			locHistName = "PVsTheta_LowPIDFOM";
			locHistTitle = locParticleROOTName + string(", PID FOM < 1%;#theta#circ;p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_PVsTheta_LowPIDFOM[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PVsTheta_LowPIDFOM[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			// P Vs Theta, PID FOM = NaN
			locHistName = "PVsTheta_NaNPIDFOM";
			locHistTitle = locParticleROOTName + string(", PID FOM = NaN;#theta#circ;p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_PVsTheta_NaNPIDFOM[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PVsTheta_NaNPIDFOM[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			if(ParticleCharge(locPID) != 0) //no other sources of PID for neutrals
			{
				// P Vs Theta, TOF FOM < 1%
				locHistName = "PVsTheta_LowTOFFOM";
				locHistTitle = locParticleROOTName + string(", TOF FOM < 1%;#theta#circ;p (GeV/c)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_PVsTheta_LowTOFFOM[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_PVsTheta_LowTOFFOM[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

				// P Vs Theta, TOF FOM = NaN
				locHistName = "PVsTheta_NaNTOFFOM";
				locHistTitle = locParticleROOTName + string(", TOF FOM = NaN;#theta#circ;p (GeV/c)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_PVsTheta_NaNTOFFOM[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_PVsTheta_NaNTOFFOM[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

				// P Vs Theta, DC dE/dx FOM < 1%
				locHistName = "PVsTheta_LowDCdEdxFOM";
				locHistTitle = locParticleROOTName + string(", DC #it{#frac{dE}{dx}} FOM < 1%;#theta#circ;p (GeV/c)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_PVsTheta_LowDCdEdxFOM[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_PVsTheta_LowDCdEdxFOM[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

				// P Vs Theta, DC dE/dx FOM = NaN
				locHistName = "PVsTheta_NaNDCdEdxFOM";
				locHistTitle = locParticleROOTName + string(", DC #it{#frac{dE}{dx}} FOM = NaN;#theta#circ;p (GeV/c)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_PVsTheta_NaNDCdEdxFOM[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_PVsTheta_NaNDCdEdxFOM[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

				// P Vs Theta, Beta < 0
				locHistName = "PVsTheta_NegativeBeta";
				locHistTitle = locParticleROOTName + string(", #beta < 0.0;#theta#circ;p (GeV/c)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_PVsTheta_NegativeBeta[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_PVsTheta_NegativeBeta[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);
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

void DHistogramAction_TrackVertexComparison::Initialize(JEventLoop* locEventLoop)
{
	deque<deque<Particle_t> > locDetectedChargedPIDs;
	Get_Reaction()->Get_DetectedFinalChargedPIDs(locDetectedChargedPIDs);

	deque<deque<Particle_t> > locDetectedChargedPIDs_HasDupes;
	Get_Reaction()->Get_DetectedFinalChargedPIDs(locDetectedChargedPIDs_HasDupes, true);

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
			dHistDeque_MaxTrackDeltaZ[loc_i] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistDeque_MaxTrackDeltaZ[loc_i] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaVertexZBins, 0.0, dMaxDeltaVertexZ);

		// Max Track DeltaT
		locHistName = "MaxTrackDeltaT";
		locHistTitle = locStepROOTName + string(";Largest Track #DeltaVertex-T (ns)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistDeque_MaxTrackDeltaT[loc_i] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistDeque_MaxTrackDeltaT[loc_i] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaVertexTBins, 0.0, dMaxDeltaVertexT);

		// Max Track DOCA
		locHistName = "MaxTrackDOCA";
		locHistTitle = locStepROOTName + string(";Largest Track DOCA (cm)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistDeque_MaxTrackDOCA[loc_i] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistDeque_MaxTrackDOCA[loc_i] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDOCABins, 0.0, dMaxDOCA);

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
				dHistDeque_TrackZToCommon[loc_i][locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistDeque_TrackZToCommon[loc_i][locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaVertexZBins, dMinDeltaVertexZ, dMaxDeltaVertexZ);

			// TrackT To Common
			locHistName = string("TrackTToCommon_") + locParticleName;
			locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";#DeltaVertex-T (Track, Common) (ns)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistDeque_TrackTToCommon[loc_i][locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistDeque_TrackTToCommon[loc_i][locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaVertexTBins, dMinDeltaVertexT, dMaxDeltaVertexT);

			// TrackDOCA To Common
			locHistName = string("TrackDOCAToCommon_") + locParticleName;
			locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";DOCA (Track, Common) (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistDeque_TrackDOCAToCommon[loc_i][locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistDeque_TrackDOCAToCommon[loc_i][locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDOCABins, dMinDOCA, dMaxDOCA);

			// DeltaT Vs P against beam photon
			if((locReactionStep->Get_InitialParticleID() == Gamma) && (dHistMap_BeamTrackDeltaTVsP.find(locPID) == dHistMap_BeamTrackDeltaTVsP.end()))
			{
				locHistName = string("TrackDeltaTVsP_") + ParticleType(locPID) + string("_Beam") + ParticleType(Gamma);
				locHistTitle = locStepROOTName + string(";") + ParticleName_ROOT(locPID) + string(" Momentum (GeV/c);t_{") + ParticleName_ROOT(locPID) + string("} - t_{Beam ") + ParticleName_ROOT(Gamma) + string("} (ns)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_BeamTrackDeltaTVsP[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_BeamTrackDeltaTVsP[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaVertexTBins, dMinDeltaVertexT, dMaxDeltaVertexT);
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
					dHistDeque_TrackDeltaTVsP[loc_i][locParticlePair] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_TrackDeltaTVsP[loc_i][locParticlePair] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaVertexTBins, dMinDeltaVertexT, dMaxDeltaVertexT);
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
			dBeamParticleHist_P = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dBeamParticleHist_P = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPBins, dMinP, dMaxP);

		// Theta
		locHistName = "Theta";
		locHistTitle = string("Beam ") + locParticleROOTName + string(";#theta#circ");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dBeamParticleHist_Theta = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dBeamParticleHist_Theta = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumThetaBins, dMinTheta, dMaxTheta);

		// Phi
		locHistName = "Phi";
		locHistTitle = string("Beam ") + locParticleROOTName + string(";#phi#circ");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dBeamParticleHist_Phi = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dBeamParticleHist_Phi = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPhiBins, dMinPhi, dMaxPhi);

		// P Vs Theta
		locHistName = "PVsTheta";
		locHistTitle = string("Beam ") + locParticleROOTName + string(";#theta#circ;p (GeV/c)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dBeamParticleHist_PVsTheta = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
		else
			dBeamParticleHist_PVsTheta = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

		// Phi Vs Theta
		locHistName = "PhiVsTheta";
		locHistTitle = string("Beam ") + locParticleROOTName + string(";#theta#circ;#phi#circ");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dBeamParticleHist_PhiVsTheta = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
		else
			dBeamParticleHist_PhiVsTheta = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

		// Vertex-Z
		locHistName = "VertexZ";
		locHistTitle = string("Beam ") + locParticleROOTName + string(";Vertex-Z (cm)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dBeamParticleHist_VertexZ = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dBeamParticleHist_VertexZ = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumVertexZBins, dMinVertexZ, dMaxVertexZ);

		// Vertex-Y Vs Vertex-X
		locHistName = "VertexYVsX";
		locHistTitle = string("Beam ") + locParticleROOTName + string(";Vertex-X (cm);Vertex-Y (cm)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dBeamParticleHist_VertexYVsX = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
		else
			dBeamParticleHist_VertexYVsX = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNumVertexXYBins, dMinVertexXY, dMaxVertexXY, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY);

		// Vertex-T
		locHistName = "VertexT";
		locHistTitle = string("Beam ") + locParticleROOTName + string(";Vertex-T (ns)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dBeamParticleHist_VertexT = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dBeamParticleHist_VertexT = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumTBins, dMinT, dMaxT);

		// Delta-T (Beam, RF)
		locHistName = "DeltaTRF";
		locHistTitle = string("Beam ") + locParticleROOTName + string(";#Deltat_{Beam - RF} (ns)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dBeamParticleHist_DeltaTRF = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dBeamParticleHist_DeltaTRF = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTRFBins, dMinDeltaTRF, dMaxDeltaTRF);

		// Delta-T (Beam, RF) Vs Beam E
		locHistName = "DeltaTRFVsBeamE";
		locHistTitle = string("Beam ") + locParticleROOTName + string(";E (GeV);#Deltat_{Beam - RF} (ns)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dBeamParticleHist_DeltaTRFVsBeamE = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
		else
			dBeamParticleHist_DeltaTRFVsBeamE = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaTRFBins, dMinDeltaTRF, dMaxDeltaTRF);

		gDirectory->cd("..");
	}

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
			if(dHistDeque_P[loc_i].find(locPID) != dHistDeque_P[loc_i].end())
				continue; //pid already done

			CreateAndChangeTo_Directory(locParticleName, locParticleName);

			// Momentum
			locHistName = "Momentum";
			locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistDeque_P[loc_i][locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistDeque_P[loc_i][locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPBins, dMinP, dMaxP);

			// Theta
			locHistName = "Theta";
			locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";#theta#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistDeque_Theta[loc_i][locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistDeque_Theta[loc_i][locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumThetaBins, dMinTheta, dMaxTheta);

			// Phi
			locHistName = "Phi";
			locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";#phi#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistDeque_Phi[loc_i][locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistDeque_Phi[loc_i][locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPhiBins, dMinPhi, dMaxPhi);

			// P Vs Theta
			locHistName = "PVsTheta";
			locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";#theta#circ;p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistDeque_PVsTheta[loc_i][locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistDeque_PVsTheta[loc_i][locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			// Phi Vs Theta
			locHistName = "PhiVsTheta";
			locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";#theta#circ;#phi#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistDeque_PhiVsTheta[loc_i][locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistDeque_PhiVsTheta[loc_i][locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

			//beta vs p
			locHistName = "BetaVsP";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#beta");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistDeque_BetaVsP[loc_i][locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistDeque_BetaVsP[loc_i][locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumBetaBins, dMinBeta, dMaxBeta);

			//delta-beta vs p
			locHistName = "DeltaBetaVsP";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta#beta");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistDeque_DeltaBetaVsP[loc_i][locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistDeque_DeltaBetaVsP[loc_i][locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaBetaBins, dMinDeltaBeta, dMaxDeltaBeta);

			// Vertex-Z
			locHistName = "VertexZ";
			locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";Vertex-Z (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistDeque_VertexZ[loc_i][locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistDeque_VertexZ[loc_i][locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumVertexZBins, dMinVertexZ, dMaxVertexZ);

			// Vertex-Y Vs Vertex-X
			locHistName = "VertexYVsX";
			locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";Vertex-X (cm);Vertex-Y (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistDeque_VertexYVsX[loc_i][locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistDeque_VertexYVsX[loc_i][locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNumVertexXYBins, dMinVertexXY, dMaxVertexXY, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY);

			// Vertex-T
			locHistName = "VertexT";
			locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";Vertex-T (ns)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistDeque_VertexT[loc_i][locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistDeque_VertexT[loc_i][locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumTBins, dMinT, dMaxT);

			gDirectory->cd("..");
		} //end of particle loop
		gDirectory->cd("..");
	} //end of step loop
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
			dRFBeamBunchDeltaT_Hist = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dRFBeamBunchDeltaT_Hist = new TH1D(locHistName.c_str(), ";RF #Deltat (Reconstructed - Thrown)", dNumRFDeltaTBins, dMinRFDeltaT, dMaxRFDeltaT);

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
				dBeamParticleHist_DeltaPOverP = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dBeamParticleHist_DeltaPOverP = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

			// DeltaP/P Vs P
			locHistName = string("DeltaPOverPVsP");
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltap/p (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dBeamParticleHist_DeltaPOverPVsP = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dBeamParticleHist_DeltaPOverPVsP = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

			// DeltaT
			locHistName = string("DeltaT");
			locHistTitle = locParticleROOTName + string(";#Deltat (ns) (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dBeamParticleHist_DeltaT = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dBeamParticleHist_DeltaT = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

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
					dHistDeque_DeltaPOverP[loc_i][locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaPOverP[loc_i][locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

				// DeltaTheta
				locHistName = string("DeltaTheta");
				locHistTitle = locParticleROOTName + string(";#Delta#theta#circ (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaTheta[loc_i][locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaTheta[loc_i][locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaThetaBins, dMinDeltaTheta, dMaxDeltaTheta);

				// DeltaPhi
				locHistName = string("DeltaPhi");
				locHistTitle = locParticleROOTName + string(";#Delta#phi#circ (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaPhi[loc_i][locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaPhi[loc_i][locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

				// DeltaT
				locHistName = string("DeltaT");
				locHistTitle = locParticleROOTName + string(";#Deltat (ns) (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaT[loc_i][locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaT[loc_i][locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

				// DeltaT - BCAL
				locHistName = string("DeltaT_BCAL");
				locHistTitle = locParticleROOTName + string(" in BCAL;#Deltat (ns) (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaT_BCAL[loc_i][locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaT_BCAL[loc_i][locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

				// DeltaT - TOF (charged only)
				if(ParticleCharge(locPID) != 0)
				{
					locHistName = string("DeltaT_TOF");
					locHistTitle = locParticleROOTName + string(" in TOF;#Deltat (ns) (Reconstructed - Thrown)");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistDeque_DeltaT_TOF[loc_i][locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistDeque_DeltaT_TOF[loc_i][locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
				}

				// DeltaT - FCAL (neutral only)
				if(ParticleCharge(locPID) == 0)
				{
					locHistName = string("DeltaT_FCAL");
					locHistTitle = locParticleROOTName + string(" in FCAL;#Deltat (ns) (Reconstructed - Thrown)");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistDeque_DeltaT_FCAL[loc_i][locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistDeque_DeltaT_FCAL[loc_i][locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
				}

				// DeltaVertexZ
				locHistName = string("DeltaVertexZ");
				locHistTitle = locParticleROOTName + string(";#DeltaVertex-Z (cm) (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaVertexZ[loc_i][locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaVertexZ[loc_i][locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaVertexZBins, dMinDeltaVertexZ, dMaxDeltaVertexZ);

				// DeltaP/P Vs P
				locHistName = string("DeltaPOverPVsP");
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltap/p (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaPOverPVsP[loc_i][locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaPOverPVsP[loc_i][locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

				// DeltaP/P Vs Theta
				locHistName = string("DeltaPOverPVsTheta");
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltap/p (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaPOverPVsTheta[loc_i][locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaPOverPVsTheta[loc_i][locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

				// DeltaTheta Vs P
				locHistName = string("DeltaThetaVsP");
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta#theta#circ (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaThetaVsP[loc_i][locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaThetaVsP[loc_i][locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaThetaBins, dMinDeltaTheta, dMaxDeltaTheta);

				// DeltaTheta Vs Theta
				locHistName = string("DeltaThetaVsTheta");
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Delta#theta#circ (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaThetaVsTheta[loc_i][locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaThetaVsTheta[loc_i][locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaThetaBins, dMinDeltaTheta, dMaxDeltaTheta);

				// DeltaPhi Vs P
				locHistName = string("DeltaPhiVsP");
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta#phi#circ (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaPhiVsP[loc_i][locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaPhiVsP[loc_i][locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

				// DeltaPhi Vs Theta
				locHistName = string("DeltaPhiVsTheta");
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Delta#phi#circ (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaPhiVsTheta[loc_i][locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaPhiVsTheta[loc_i][locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

				// DeltaT Vs Theta
				locHistName = string("DeltaTVsTheta");
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat (ns) (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaTVsTheta[loc_i][locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaTVsTheta[loc_i][locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

				// DeltaT Vs P
				locHistName = string("DeltaTVsP");
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat (ns) (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaTVsP[loc_i][locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaTVsP[loc_i][locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

				// DeltaVertexZ Vs Theta
				locHistName = string("DeltaVertexZVsTheta");
				locHistTitle = locParticleROOTName + string(";#theta#circ;#DeltaVertex-Z (cm) (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_DeltaVertexZVsTheta[loc_i][locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_DeltaVertexZVsTheta[loc_i][locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaVertexZBins, dMinDeltaVertexZ, dMaxDeltaVertexZ);

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
						(dHistDeque_Pulls[loc_i][locPID])[dPullTypes[loc_j]] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
					else
						(dHistDeque_Pulls[loc_i][locPID])[dPullTypes[loc_j]] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

					//Pull vs P
					locHistName = locPullNames[loc_j] + string("PullVsP");
					locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta") + locPullTitles[loc_j] + string("/#sigma_{") + locPullTitles[loc_j] + string("} (Reconstructed - Thrown)");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						(dHistDeque_PullsVsP[loc_i][locPID])[dPullTypes[loc_j]] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
					else
						(dHistDeque_PullsVsP[loc_i][locPID])[dPullTypes[loc_j]] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

					//Pull vs Theta
					locHistName = locPullNames[loc_j] + string("PullVsTheta");
					locHistTitle = locParticleROOTName + string(";#theta#circ;#Delta") + locPullTitles[loc_j] + string("/#sigma_{") + locPullTitles[loc_j] + string("} (Reconstructed - Thrown)");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						(dHistDeque_PullsVsTheta[loc_i][locPID])[dPullTypes[loc_j]] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
					else
						(dHistDeque_PullsVsTheta[loc_i][locPID])[dPullTypes[loc_j]] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);
				}

				//Delta-t Pulls - CDC & ST
				if(ParticleCharge(locPID) != 0)
				{
					//CDC
					locHistName = "TimePull_CDC";
					locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistDeque_TimePull_CDC[loc_i][locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistDeque_TimePull_CDC[loc_i][locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

					locHistName = "TimePullVsTheta_CDC";
					locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat/#sigma_{#Deltat}");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistDeque_TimePullVsTheta_CDC[loc_i][locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistDeque_TimePullVsTheta_CDC[loc_i][locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);

					locHistName = "TimePullVsP_CDC";
					locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistDeque_TimePullVsP_CDC[loc_i][locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistDeque_TimePullVsP_CDC[loc_i][locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

					//ST
					locHistName = "TimePull_ST";
					locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistDeque_TimePull_ST[loc_i][locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistDeque_TimePull_ST[loc_i][locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

					locHistName = "TimePullVsTheta_ST";
					locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat/#sigma_{#Deltat}");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistDeque_TimePullVsTheta_ST[loc_i][locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistDeque_TimePullVsTheta_ST[loc_i][locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);

					locHistName = "TimePullVsP_ST";
					locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistDeque_TimePullVsP_ST[loc_i][locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistDeque_TimePullVsP_ST[loc_i][locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);
				}

				//Delta-t Pulls - BCAL
				locHistName = "TimePull_BCAL";
				locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_TimePull_BCAL[loc_i][locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_TimePull_BCAL[loc_i][locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

				locHistName = "TimePullVsTheta_BCAL";
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_TimePullVsTheta_BCAL[loc_i][locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_TimePullVsTheta_BCAL[loc_i][locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);

				locHistName = "TimePullVsP_BCAL";
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_TimePullVsP_BCAL[loc_i][locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_TimePullVsP_BCAL[loc_i][locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

				//Delta-t Pulls - TOF
				if(ParticleCharge(locPID) != 0) //TOF
				{
					locHistName = "TimePull_TOF";
					locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistDeque_TimePull_TOF[loc_i][locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistDeque_TimePull_TOF[loc_i][locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

					locHistName = "TimePullVsP_TOF";
					locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
					if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
						dHistDeque_TimePullVsP_TOF[loc_i][locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
					else
						dHistDeque_TimePullVsP_TOF[loc_i][locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);
				}

				//Delta-t Pulls - FCAL
				locHistName = "TimePull_FCAL";
				locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_TimePull_FCAL[loc_i][locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_TimePull_FCAL[loc_i][locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

				locHistName = "TimePullVsP_FCAL";
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistDeque_TimePullVsP_FCAL[loc_i][locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistDeque_TimePullVsP_FCAL[loc_i][locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

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
	double locRFTime = locEventRFBunch->dMatchedToTracksFlag ? locEventRFBunch->dTime : numeric_limits<double>::quiet_NaN();
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

		// Beam Particle
		locPID = Gamma;
		locParticleName = string("Beam_") + ParticleType(locPID);
		locParticleROOTName = ParticleName_ROOT(locPID);
		CreateAndChangeTo_Directory(locParticleName, locParticleName);
		locHistName = "Momentum";
		locHistTitle = string("Thrown Beam ") + locParticleROOTName + string(";p (GeV/c)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dBeamParticle_P = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dBeamParticle_P = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPBins, dMinP, dMaxP);
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
				dHistMap_P[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_P[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPBins, dMinP, dMaxP);

			// Theta
			locHistName = "Theta";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";#theta#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_Theta[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_Theta[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumThetaBins, dMinTheta, dMaxTheta);

			// Phi
			locHistName = "Phi";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";#phi#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_Phi[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_Phi[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPhiBins, dMinPhi, dMaxPhi);

			// P Vs Theta
			locHistName = "PVsTheta";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";#theta#circ;p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_PVsTheta[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PVsTheta[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			// Phi Vs Theta
			locHistName = "PhiVsTheta";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";#theta#circ;#phi#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_PhiVsTheta[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PhiVsTheta[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

			// Vertex-Z
			locHistName = "VertexZ";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";Vertex-Z (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_VertexZ[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_VertexZ[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumVertexZBins, dMinVertexZ, dMaxVertexZ);

			// Vertex-Y Vs Vertex-X
			locHistName = "VertexYVsX";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";Vertex-X (cm);Vertex-Y (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_VertexYVsX[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_VertexYVsX[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNumVertexXYBins, dMinVertexXY, dMaxVertexXY, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY);

			// Vertex-T
			locHistName = "VertexT";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";Vertex-T (ns)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_VertexT[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_VertexT[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumTBins, dMinT, dMaxT);

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

	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons, "MCGEN");
	japp->RootWriteLock();
	{
		for(size_t loc_i = 0; loc_i < locBeamPhotons.size(); ++loc_i)
			dBeamParticle_P->Fill(locBeamPhotons[loc_i]->energy());
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
		CreateAndChangeTo_ActionDirectory();

		// Beam Particle
		locPID = Gamma;
		locParticleName = string("Beam_") + ParticleType(locPID);
		locParticleROOTName = ParticleName_ROOT(locPID);
		CreateAndChangeTo_Directory(locParticleName, locParticleName);
		locHistName = "Momentum";
		locHistTitle = string("Thrown Beam ") + locParticleROOTName + string(";p (GeV/c)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dBeamParticle_P = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dBeamParticle_P = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPBins, dMinP, dMaxP);
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
				dHistMap_P[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_P[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPBins, dMinP, dMaxP);

			// Theta
			locHistName = "Theta";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";#theta#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_Theta[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_Theta[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumThetaBins, dMinTheta, dMaxTheta);

			// Phi
			locHistName = "Phi";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";#phi#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_Phi[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_Phi[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPhiBins, dMinPhi, dMaxPhi);

			// P Vs Theta
			locHistName = "PVsTheta";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";#theta#circ;p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_PVsTheta[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PVsTheta[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			// Phi Vs Theta
			locHistName = "PhiVsTheta";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";#theta#circ;#phi#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_PhiVsTheta[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PhiVsTheta[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

			// Vertex-Z
			locHistName = "VertexZ";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";Vertex-Z (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_VertexZ[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_VertexZ[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumVertexZBins, dMinVertexZ, dMaxVertexZ);

			// Vertex-Y Vs Vertex-X
			locHistName = "VertexYVsX";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";Vertex-X (cm);Vertex-Y (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_VertexYVsX[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_VertexYVsX[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNumVertexXYBins, dMinVertexXY, dMaxVertexXY, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY);

			// Vertex-T
			locHistName = "VertexT";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";Vertex-T (ns)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_VertexT[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_VertexT[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumTBins, dMinT, dMaxT);

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
		if(dHistMap_P.find(locPID) == dHistMap_P.end())
			continue; //not interested in histogramming

		double locMatchFOM = 0.0;
		if(ParticleCharge(locPID) != 0)
		{
			if(locMCThrownMatching->Get_MatchingChargedHypothesis(locMCThrown, locMatchFOM) == NULL)
				continue; //not reconstructed
		}
		else if(locMCThrownMatching->Get_MatchingNeutralHypothesis(locMCThrown, locMatchFOM) == NULL)
			continue; //not reconstructed

		if(locMatchFOM < dMinThrownMatchFOM)
			continue; //not well reconstructed

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

		// Beam Particle
		locPID = Gamma;
		locParticleName = string("Beam_") + ParticleType(locPID);
		locParticleROOTName = ParticleName_ROOT(locPID);
		CreateAndChangeTo_Directory(locParticleName, locParticleName);
		locHistName = "Momentum";
		locHistTitle = string("Beam ") + locParticleROOTName + string(";p (GeV/c)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dBeamParticle_P = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dBeamParticle_P = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPBins, dMinP, dMaxP);
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
				dHistMap_P[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_P[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPBins, dMinP, dMaxP);

			// Theta
			locHistName = "Theta";
			locHistTitle = locParticleROOTName + string(";#theta#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_Theta[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_Theta[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumThetaBins, dMinTheta, dMaxTheta);

			// Phi
			locHistName = "Phi";
			locHistTitle = locParticleROOTName + string(";#phi#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_Phi[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_Phi[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPhiBins, dMinPhi, dMaxPhi);

			// P Vs Theta
			locHistName = "PVsTheta";
			locHistTitle = locParticleROOTName + string(";#theta#circ;p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_PVsTheta[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PVsTheta[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			// Phi Vs Theta
			locHistName = "PhiVsTheta";
			locHistTitle = locParticleROOTName + string(";#theta#circ;#phi#circ");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_PhiVsTheta[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PhiVsTheta[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

			//beta vs p
			locHistName = "BetaVsP";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#beta");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_BetaVsP[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_BetaVsP[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumBetaBins, dMinBeta, dMaxBeta);

			//delta-beta vs p
			locHistName = "DeltaBetaVsP";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta#beta");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaBetaVsP[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaBetaVsP[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaBetaBins, dMinDeltaBeta, dMaxDeltaBeta);

			// Vertex-Z
			locHistName = "VertexZ";
			locHistTitle = locParticleROOTName + string(";Vertex-Z (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_VertexZ[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_VertexZ[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumVertexZBins, dMinVertexZ, dMaxVertexZ);

			if(ParticleCharge(locPID) != 0)
			{
				// Tracking FOM Vs Vertex-Z
				locHistName = "TrackingFOMVsVertexZ";
				locHistTitle = locParticleROOTName + string(";Vertex-Z (cm); Time-Based Track Reconstruction FOM");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TrackingFOMVsVertexZ[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TrackingFOMVsVertexZ[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DVertexZBins, dMinVertexZ, dMaxVertexZ, dNumTrackFOMBins, 0.0, 1.0);
			}

			// Vertex-Z Vs Theta
			locHistName = "VertexZVsTheta";
			locHistTitle = locParticleROOTName + string(";#theta#circ;Vertex-Z (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_VertexZVsTheta[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_VertexZVsTheta[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DVertexZBins, dMinVertexZ, dMaxVertexZ);

			// Vertex-Y Vs Vertex-X
			locHistName = "VertexYVsX";
			locHistTitle = locParticleROOTName + string(";Vertex-X (cm);Vertex-Y (cm)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_VertexYVsX[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_VertexYVsX[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNumVertexXYBins, dMinVertexXY, dMaxVertexXY, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY);

			// Vertex-T
			locHistName = "VertexT";
			locHistTitle = locParticleROOTName + string(";Vertex-T (ns)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_VertexT[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_VertexT[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumTBins, dMinT, dMaxT);

			gDirectory->cd("..");
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_DetectedParticleKinematics::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	if(Get_NumPreviousParticleCombos() != 0)
		return true; //else double-counting!

	Particle_t locPID;

	vector<const DEventRFBunch*> locEventRFBunches;
	locEventLoop->Get(locEventRFBunches);
	const DEventRFBunch* locEventRFBunch = locEventRFBunches[0];
	if(locParticleCombo != NULL)
		locEventRFBunch = locParticleCombo->Get_EventRFBunch();

	double locBeta_Timing, locDeltaBeta;

	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons);
	for(size_t loc_i = 0; loc_i < locBeamPhotons.size(); ++loc_i)
		dBeamParticle_P->Fill(locBeamPhotons[loc_i]->energy());

	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks);
	const DChargedTrackHypothesis* locChargedTrackHypothesis = NULL;

	for(size_t loc_i = 0; loc_i < locChargedTracks.size(); ++loc_i)
	{
		locChargedTrackHypothesis = locChargedTracks[loc_i]->Get_BestFOM();
		locPID = (locChargedTrackHypothesis->dFOM < dMinimumPIDFOM) ? Unknown : locChargedTrackHypothesis->PID();

		double locTrackingFOM = TMath::Prob(locChargedTrackHypothesis->dChiSq_Track, locChargedTrackHypothesis->dNDF_Track);
		if(locTrackingFOM < dMinimumTrackingFOM)
			continue;

		if(dHistMap_P.find(locPID) == dHistMap_P.end())
			continue; //e.g. a decaying particle, or not interested in histogramming

		DVector3 locMomentum = locChargedTrackHypothesis->momentum();
		double locPhi = locMomentum.Phi()*180.0/TMath::Pi();
		double locTheta = locMomentum.Theta()*180.0/TMath::Pi();
		double locP = locMomentum.Mag();
		locBeta_Timing = dAnalysisUtilities->Calc_Beta_Timing(locChargedTrackHypothesis, locEventRFBunch, true);
		locDeltaBeta = locChargedTrackHypothesis->lorentzMomentum().Beta() - locBeta_Timing;

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
			dHistMap_TrackingFOMVsVertexZ[locPID]->Fill(locChargedTrackHypothesis->position().Z(), locTrackingFOM);
			dHistMap_VertexZVsTheta[locPID]->Fill(locTheta, locChargedTrackHypothesis->position().Z());
			dHistMap_VertexYVsX[locPID]->Fill(locChargedTrackHypothesis->position().X(), locChargedTrackHypothesis->position().Y());
			dHistMap_VertexT[locPID]->Fill(locChargedTrackHypothesis->time());
		}
		japp->RootUnLock();
	}

	vector<const DNeutralParticle*> locNeutralParticles;
	locEventLoop->Get(locNeutralParticles);
	const DNeutralParticleHypothesis* locNeutralParticleHypothesis;

	for(size_t loc_i = 0; loc_i < locNeutralParticles.size(); ++loc_i)
	{
		locNeutralParticleHypothesis = locNeutralParticles[loc_i]->Get_Hypothesis(Gamma);
		if(locNeutralParticleHypothesis->dFOM < dMinimumPIDFOM)
			locNeutralParticleHypothesis = locNeutralParticles[loc_i]->Get_Hypothesis(Neutron);

		locPID = locNeutralParticleHypothesis->PID();

		if(dHistMap_P.find(locPID) == dHistMap_P.end())
			continue; //e.g. a decaying particle, or not interested in histogramming

		DVector3 locMomentum = locNeutralParticleHypothesis->momentum();
		double locPhi = locMomentum.Phi()*180.0/TMath::Pi();
		double locTheta = locMomentum.Theta()*180.0/TMath::Pi();
		double locP = locMomentum.Mag();

		locBeta_Timing = dAnalysisUtilities->Calc_Beta_Timing(locNeutralParticleHypothesis, locEventRFBunch);
		locDeltaBeta = locNeutralParticleHypothesis->lorentzMomentum().Beta() - locBeta_Timing;

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
			dHistMap_VertexZVsTheta[locPID]->Fill(locTheta, locNeutralParticleHypothesis->position().Z());
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
			dRFBeamBunchDeltaT_Hist = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dRFBeamBunchDeltaT_Hist = new TH1D(locHistName.c_str(), ";RF #Deltat (Reconstructed - Thrown)", dNumRFDeltaTBins, dMinRFDeltaT, dMaxRFDeltaT);

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

			// DeltaP/P
			locHistName = string("DeltaPOverP");
			locHistTitle = locParticleROOTName + string(";#Deltap/p (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaPOverP[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaPOverP[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

			// DeltaTheta
			locHistName = string("DeltaTheta");
			locHistTitle = locParticleROOTName + string(";#Delta#theta#circ (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaTheta[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaTheta[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaThetaBins, dMinDeltaTheta, dMaxDeltaTheta);

			// DeltaPhi
			locHistName = string("DeltaPhi");
			locHistTitle = locParticleROOTName + string(";#Delta#phi#circ (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaPhi[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaPhi[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

			// DeltaT
			locHistName = string("DeltaT");
			locHistTitle = locParticleROOTName + string(";#Deltat (ns) (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaT[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaT[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

			// DeltaT - BCAL
			locHistName = string("DeltaT_BCAL");
			locHistTitle = locParticleROOTName + string(" in BCAL;#Deltat (ns) (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaT_BCAL[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaT_BCAL[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

			// DeltaT - TOF (charged only)
			if(ParticleCharge(locPID) != 0)
			{
				locHistName = string("DeltaT_TOF");
				locHistTitle = locParticleROOTName + string(" in TOF;#Deltat (ns) (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_DeltaT_TOF[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_DeltaT_TOF[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
			}

			// DeltaT - FCAL (neutral only)
			if(ParticleCharge(locPID) == 0)
			{
				locHistName = string("DeltaT_FCAL");
				locHistTitle = locParticleROOTName + string(" in FCAL;#Deltat (ns) (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_DeltaT_FCAL[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_DeltaT_FCAL[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
			}

			// DeltaVertexZ
			locHistName = string("DeltaVertexZ");
			locHistTitle = locParticleROOTName + string(";#DeltaVertex-Z (cm) (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaVertexZ[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaVertexZ[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaVertexZBins, dMinDeltaVertexZ, dMaxDeltaVertexZ);

			// DeltaP/P Vs P
			locHistName = string("DeltaPOverPVsP");
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltap/p (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaPOverPVsP[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaPOverPVsP[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

			// DeltaP/P Vs Theta
			locHistName = string("DeltaPOverPVsTheta");
			locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltap/p (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaPOverPVsTheta[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaPOverPVsTheta[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

			// DeltaTheta Vs P
			locHistName = string("DeltaThetaVsP");
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta#theta#circ (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaThetaVsP[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaThetaVsP[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaThetaBins, dMinDeltaTheta, dMaxDeltaTheta);

			// DeltaTheta Vs Theta
			locHistName = string("DeltaThetaVsTheta");
			locHistTitle = locParticleROOTName + string(";#theta#circ;#Delta#theta#circ (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaThetaVsTheta[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaThetaVsTheta[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaThetaBins, dMinDeltaTheta, dMaxDeltaTheta);

			// DeltaPhi Vs P
			locHistName = string("DeltaPhiVsP");
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta#phi#circ (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaPhiVsP[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaPhiVsP[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

			// DeltaPhi Vs Theta
			locHistName = string("DeltaPhiVsTheta");
			locHistTitle = locParticleROOTName + string(";#theta#circ;#Delta#phi#circ (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaPhiVsTheta[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaPhiVsTheta[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

			// DeltaT Vs Theta
			locHistName = string("DeltaTVsTheta");
			locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat (ns) (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaTVsTheta[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaTVsTheta[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

			// DeltaT Vs P
			locHistName = string("DeltaTVsP");
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat (ns) (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaTVsP[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaTVsP[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

			// DeltaVertexZ Vs Theta
			locHistName = string("DeltaVertexZVsTheta");
			locHistTitle = locParticleROOTName + string(";#theta#circ;#DeltaVertex-Z (cm) (Reconstructed - Thrown)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DeltaVertexZVsTheta[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DeltaVertexZVsTheta[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaVertexZBins, dMinDeltaVertexZ, dMaxDeltaVertexZ);

			// P Vs Theta
			locHistName = "PVsTheta_LargeDeltaT";
			locHistTitle = locParticleROOTName + string(";#theta#circ;p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_PVsTheta_LargeDeltaT[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_PVsTheta_LargeDeltaT[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

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
					(dHistMap_Pulls[locPID])[dPullTypes[loc_j]] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
				else
					(dHistMap_Pulls[locPID])[dPullTypes[loc_j]] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

				//Pull vs P
				locHistName = locPullNames[loc_j] + string("PullVsP");
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta") + locPullTitles[loc_j] + string("/#sigma_{") + locPullTitles[loc_j] + string("} (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					(dHistMap_PullsVsP[locPID])[dPullTypes[loc_j]] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					(dHistMap_PullsVsP[locPID])[dPullTypes[loc_j]] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

				//Pull vs Theta
				locHistName = locPullNames[loc_j] + string("PullVsTheta");
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Delta") + locPullTitles[loc_j] + string("/#sigma_{") + locPullTitles[loc_j] + string("} (Reconstructed - Thrown)");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					(dHistMap_PullsVsTheta[locPID])[dPullTypes[loc_j]] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					(dHistMap_PullsVsTheta[locPID])[dPullTypes[loc_j]] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);
			}

			//Delta-t Pulls - CDC & ST
			if(ParticleCharge(locPID) != 0)
			{
				//CDC
				locHistName = "TimePull_CDC";
				locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TimePull_CDC[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePull_CDC[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

				locHistName = "TimePullVsTheta_CDC";
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TimePullVsTheta_CDC[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePullVsTheta_CDC[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);

				locHistName = "TimePullVsP_CDC";
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TimePullVsP_CDC[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePullVsP_CDC[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

				//ST
				locHistName = "TimePull_ST";
				locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TimePull_ST[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePull_ST[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

				locHistName = "TimePullVsTheta_ST";
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TimePullVsTheta_ST[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePullVsTheta_ST[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);

				locHistName = "TimePullVsP_ST";
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TimePullVsP_ST[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePullVsP_ST[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);
			}

			//Delta-t Pulls - BCAL
			locHistName = "TimePull_BCAL";
			locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TimePull_BCAL[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TimePull_BCAL[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

			locHistName = "TimePullVsTheta_BCAL";
			locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat/#sigma_{#Deltat}");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TimePullVsTheta_BCAL[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TimePullVsTheta_BCAL[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);

			locHistName = "TimePullVsP_BCAL";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TimePullVsP_BCAL[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TimePullVsP_BCAL[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

			//Delta-t Pulls - TOF
			if(ParticleCharge(locPID) != 0) //TOF
			{
				locHistName = "TimePull_TOF";
				locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TimePull_TOF[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePull_TOF[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

				locHistName = "TimePullVsP_TOF";
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
				if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
					dHistMap_TimePullVsP_TOF[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
				else
					dHistMap_TimePullVsP_TOF[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);
			}

			//Delta-t Pulls - FCAL
			locHistName = "TimePull_FCAL";
			locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TimePull_FCAL[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TimePull_FCAL[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, -10.0, 10.0);

			locHistName = "TimePullVsP_FCAL";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TimePullVsP_FCAL[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TimePullVsP_FCAL[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

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
	double locRFTime = locEventRFBunch->dMatchedToTracksFlag ? locEventRFBunch->dTime : numeric_limits<double>::quiet_NaN();
	double locRFDeltaT = locRFTime - locThrownEventRFBunch->dTime;
	japp->RootWriteLock();
	{
		dRFBeamBunchDeltaT_Hist->Fill(locRFDeltaT);
	}
	japp->RootUnLock();

	//charged particles
	map<const DMCThrown*, pair<const DChargedTrack*, double> > locThrownToChargedMap;
	locMCThrownMatching->Get_ThrownToChargedMap(locThrownToChargedMap);
	for(map<const DMCThrown*, pair<const DChargedTrack*, double> >::iterator locIterator = locThrownToChargedMap.begin(); locIterator != locThrownToChargedMap.end(); ++locIterator)
	{
		const DMCThrown* locMCThrown = locIterator->first;
		locPID = (Particle_t)locMCThrown->type;
		if(dHistMap_DeltaPOverP.find(locPID) == dHistMap_DeltaPOverP.end())
			continue; //e.g. not interested in histogramming

		const DChargedTrackHypothesis* locChargedTrackHypothesis = locIterator->second.first->Get_Hypothesis(locPID);
		if(locChargedTrackHypothesis == NULL)
			locChargedTrackHypothesis = locIterator->second.first->Get_BestFOM();

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
	for(map<const DMCThrown*, pair<const DNeutralParticle*, double> >::iterator locIterator = locThrownToNeutralMap.begin(); locIterator != locThrownToNeutralMap.end(); ++locIterator)
	{
		const DMCThrown* locMCThrown = locIterator->first;
		locPID = (Particle_t)locMCThrown->type;
		if(dHistMap_DeltaPOverP.find(locPID) == dHistMap_DeltaPOverP.end())
			continue; //e.g. not interested in histogramming

		const DNeutralParticleHypothesis* locNeutralParticleHypothesis = locIterator->second.first->Get_Hypothesis(locPID);
		if(locNeutralParticleHypothesis == NULL)
			locNeutralParticleHypothesis = locIterator->second.first->Get_BestFOM();

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
			dHistMap_DeltaT[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_DeltaT[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

		// DeltaX
		locHistName = string("DeltaX_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";#Deltax (cm) (Reconstructed - Thrown)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_DeltaX[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_DeltaX[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaXBins, dMinDeltaX, dMaxDeltaX);

		// DeltaY
		locHistName = string("DeltaY_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";#Deltay (cm) (Reconstructed - Thrown)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_DeltaY[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_DeltaY[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaXBins, dMinDeltaX, dMaxDeltaX);

		// dE
		locHistName = string("dE_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";dE (MeV)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_dE[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_dE[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumdEBins, dMindE, dMaxdE);

		// DeltaT Vs P
		locHistName = string("DeltaTVsP_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat (ns) (Reconstructed - Thrown)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_DeltaTVsP[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_DeltaTVsP[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

		// DeltaX Vs P
		locHistName = string("DeltaXVsP_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltax (cm) (Reconstructed - Thrown)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_DeltaXVsP[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_DeltaXVsP[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaXBins, dMinDeltaX, dMaxDeltaX);

		// DeltaY Vs P
		locHistName = string("DeltaYVsP_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltay (cm) (Reconstructed - Thrown)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_DeltaYVsP[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_DeltaYVsP[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaXBins, dMinDeltaX, dMaxDeltaX);

		// dE Vs P
		locHistName = string("dEVsP_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";p (GeV/c);dE (GeV)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_dEVsP[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_dEVsP[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumdEBins, dMindE, dMaxdE);

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

		//TBT
		locHistName = "NumPosTimeBasedTracks";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumPosTimeBasedTracks = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumPosTimeBasedTracks = new TH1D("NumPosTimeBasedTracks", ";# #it{q}^{+} Time-Based Tracks", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumNegTimeBasedTracks";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumNegTimeBasedTracks = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumNegTimeBasedTracks = new TH1D("NumNegTimeBasedTracks", ";# #it{q}^{-} Time-Based Tracks", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		if(!locIsRESTEvent)
		{
			//WBT
			locHistName = "NumPosWireBasedTracks";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumPosWireBasedTracks = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumPosWireBasedTracks = new TH1D("NumPosWireBasedTracks", ";# #it{q}^{-} Wire-Based Tracks", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			locHistName = "NumNegWireBasedTracks";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumNegWireBasedTracks = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumNegWireBasedTracks = new TH1D("NumNegWireBasedTracks", ";# #it{q}^{-} Wire-Based Tracks", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			//Track Candidates
			locHistName = "NumPosTrackCandidates";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumPosTrackCandidates = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumPosTrackCandidates = new TH1D("NumPosTrackCandidates", ";# #it{q}^{+} Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			locHistName = "NumNegTrackCandidates";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumNegTrackCandidates = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumNegTrackCandidates = new TH1D("NumNegTrackCandidates", ";# #it{q}^{-} Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			//CDC Track Candidates
			locHistName = "NumPosTrackCandidates_CDC";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumPosTrackCandidates_CDC = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumPosTrackCandidates_CDC = new TH1D("NumPosTrackCandidates_CDC", ";# #it{q}^{+} CDC Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			locHistName = "NumNegTrackCandidates_CDC";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumNegTrackCandidates_CDC = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumNegTrackCandidates_CDC = new TH1D("NumNegTrackCandidates_CDC", ";# #it{q}^{-} CDC Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			//FDC Track Candidates
			locHistName = "NumPosTrackCandidates_FDC";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumPosTrackCandidates_FDC = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumPosTrackCandidates_FDC = new TH1D("NumPosTrackCandidates_FDC", ";# #it{q}^{+} FDC Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

			locHistName = "NumNegTrackCandidates_FDC";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumNegTrackCandidates_FDC = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumNegTrackCandidates_FDC = new TH1D("NumNegTrackCandidates_FDC", ";# #it{q}^{-} FDC Track Candidates", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);
		}

		//Beam Photons
		locHistName = "NumBeamPhotons";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumBeamPhotons = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumBeamPhotons = new TH1D("NumBeamPhotons", ";# Beam Photons", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		//Showers / Neutrals / TOF / SC
		locHistName = "NumFCALShowers";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumFCALShowers = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumFCALShowers = new TH1D("NumFCALShowers", ";# FCAL Showers", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumBCALShowers";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumBCALShowers = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumBCALShowers = new TH1D("NumBCALShowers", ";# BCAL Showers", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumNeutralShowers";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumNeutralShowers = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumNeutralShowers = new TH1D("NumNeutralShowers", ";# Neutral Showers", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumTOFPoints";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumTOFPoints = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumTOFPoints = new TH1D("NumTOFPoints", ";# TOF Points", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumSCHits";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumSCHits = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumSCHits = new TH1D("NumSCHits", ";# SC Hits", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		//Matches
		locHistName = "NumTrackBCALMatches";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumTrackBCALMatches = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumTrackBCALMatches = new TH1D("NumTrackBCALMatches", ";# Track-BCAL Matches", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumTrackFCALMatches";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumTrackFCALMatches = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumTrackFCALMatches = new TH1D("NumTrackFCALMatches", ";# Track-FCAL Matches", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumTrackTOFMatches";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumTrackTOFMatches = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumTrackTOFMatches = new TH1D("NumTrackTOFMatches", ";# Track-TOF Matches", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		locHistName = "NumTrackSCMatches";
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumTrackSCMatches = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_NumTrackSCMatches = new TH1D("NumTrackSCMatches", ";# Track-SC Matches", dMaxNumObjects + 1, -0.5, (float)dMaxNumObjects + 0.5);

		if(!locIsRESTEvent)
		{
			//Hits
			locHistName = "NumCDCHits";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumCDCHits = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumCDCHits = new TH1D("NumCDCHits", ";# CDC Hits", dMaxNumCDCHits + 1, -0.5, (float)dMaxNumCDCHits + 0.5);

			locHistName = "NumFDCHits";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumFDCHits = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumFDCHits = new TH1D("NumFDCHits", ";# FDC Hits", dMaxNumFDCHits + 1, -0.5, (float)dMaxNumFDCHits + 0.5);

			locHistName = "NumTOFHits";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumTOFHits = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumTOFHits = new TH1D("NumTOFHits", ";# TOF Hits", dMaxNumTOFCalorimeterHits + 1, -0.5, (float)dMaxNumTOFCalorimeterHits + 0.5);

			locHistName = "NumBCALHits";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumBCALHits = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumBCALHits = new TH1D("NumBCALHits", ";# BCAL Hits", dMaxNumTOFCalorimeterHits + 1, -0.5, (float)dMaxNumTOFCalorimeterHits + 0.5);

			locHistName = "NumFCALHits";
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHist_NumFCALHits = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHist_NumFCALHits = new TH1D("NumFCALHits", ";# FCAL Hits", dMaxNumTOFCalorimeterHits + 1, -0.5, (float)dMaxNumTOFCalorimeterHits + 0.5);
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_NumReconstructedObjects::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	if(Get_NumPreviousParticleCombos() != 0)
		return true; //else double-counting!

	bool locIsRESTEvent = (string(locEventLoop->GetJEvent().GetJEventSource()->className()) == string("DEventSourceREST"));

	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//TBT
		vector<const DTrackTimeBased*> locTrackTimeBasedVector;
		locEventLoop->Get(locTrackTimeBasedVector);
		unsigned int locNumPos = 0, locNumNeg = 0;
		for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
		{
			if(ParticleCharge(locTrackTimeBasedVector[loc_i]->PID()) > 0)
				++locNumPos;
			else
				++locNumNeg;
		}
		dHist_NumPosTimeBasedTracks->Fill(locNumPos);
		dHist_NumPosTimeBasedTracks->SetEntries(dHist_NumPosTimeBasedTracks->GetEntries() - 1 + (Double_t)locNumPos);
		dHist_NumNegTimeBasedTracks->Fill(locNumNeg);
		dHist_NumNegTimeBasedTracks->SetEntries(dHist_NumNegTimeBasedTracks->GetEntries() - 1 + (Double_t)locNumNeg);

		if(!locIsRESTEvent)
		{
			//WBT
			vector<const DTrackWireBased*> locTrackWireBasedVector;
			locEventLoop->Get(locTrackWireBasedVector);
			locNumPos = 0;  locNumNeg = 0;
			for(size_t loc_i = 0; loc_i < locTrackWireBasedVector.size(); ++loc_i)
			{
				if(ParticleCharge(locTrackWireBasedVector[loc_i]->PID()) > 0)
					++locNumPos;
				else
					++locNumNeg;
			}
			dHist_NumPosWireBasedTracks->Fill(locNumPos);
			dHist_NumPosWireBasedTracks->SetEntries(dHist_NumPosWireBasedTracks->GetEntries() - 1 + (Double_t)locNumPos);
			dHist_NumNegWireBasedTracks->Fill(locNumNeg);
			dHist_NumNegWireBasedTracks->SetEntries(dHist_NumNegWireBasedTracks->GetEntries() - 1 + (Double_t)locNumNeg);

			//Candidates
			vector<const DTrackCandidate*> locTrackCandidates;
			locEventLoop->Get(locTrackCandidates);
			locNumPos = 0;  locNumNeg = 0;
			for(size_t loc_i = 0; loc_i < locTrackCandidates.size(); ++loc_i)
			{
				if(locTrackCandidates[loc_i]->charge() > 0.0)
					++locNumPos;
				else
					++locNumNeg;
			}
			dHist_NumPosTrackCandidates->Fill(locNumPos);
			dHist_NumPosTrackCandidates->SetEntries(dHist_NumPosTrackCandidates->GetEntries() - 1 + (Double_t)locNumPos);
			dHist_NumNegTrackCandidates->Fill(locNumNeg);
			dHist_NumNegTrackCandidates->SetEntries(dHist_NumNegTrackCandidates->GetEntries() - 1 + (Double_t)locNumNeg);

			//CDC Candidates
			locEventLoop->Get(locTrackCandidates, "CDC");
			locNumPos = 0;  locNumNeg = 0;
			for(size_t loc_i = 0; loc_i < locTrackCandidates.size(); ++loc_i)
			{
				if(locTrackCandidates[loc_i]->charge() > 0.0)
					++locNumPos;
				else
					++locNumNeg;
			}
			dHist_NumPosTrackCandidates_CDC->Fill(locNumPos);
			dHist_NumPosTrackCandidates_CDC->SetEntries(dHist_NumPosTrackCandidates_CDC->GetEntries() - 1 + (Double_t)locNumPos);
			dHist_NumNegTrackCandidates_CDC->Fill(locNumNeg);
			dHist_NumNegTrackCandidates_CDC->SetEntries(dHist_NumNegTrackCandidates_CDC->GetEntries() - 1 + (Double_t)locNumNeg);

			//FDC Candidates
			locEventLoop->Get(locTrackCandidates, "FDCCathodes");
			locNumPos = 0;  locNumNeg = 0;
			for(size_t loc_i = 0; loc_i < locTrackCandidates.size(); ++loc_i)
			{
				if(locTrackCandidates[loc_i]->charge() > 0.0)
					++locNumPos;
				else
					++locNumNeg;
			}
			dHist_NumPosTrackCandidates_FDC->Fill(locNumPos);
			dHist_NumPosTrackCandidates_FDC->SetEntries(dHist_NumPosTrackCandidates_FDC->GetEntries() - 1 + (Double_t)locNumPos);
			dHist_NumNegTrackCandidates_FDC->Fill(locNumNeg);
			dHist_NumNegTrackCandidates_FDC->SetEntries(dHist_NumNegTrackCandidates_FDC->GetEntries() - 1 + (Double_t)locNumNeg);
		}

		//Beam Photons
		vector<const DBeamPhoton*> locBeamPhotons;
		locEventLoop->Get(locBeamPhotons);
		dHist_NumBeamPhotons->Fill((Double_t)locBeamPhotons.size());
		dHist_NumBeamPhotons->SetEntries(dHist_NumBeamPhotons->GetEntries() - 1 + (Double_t)locBeamPhotons.size());

		//Showers
		vector<const DFCALShower*> locFCALShowers;
		locEventLoop->Get(locFCALShowers);
		dHist_NumFCALShowers->Fill((Double_t)locFCALShowers.size());
		dHist_NumFCALShowers->SetEntries(dHist_NumFCALShowers->GetEntries() - 1 + (Double_t)locFCALShowers.size());

		vector<const DBCALShower*> locBCALShowers;
		locEventLoop->Get(locBCALShowers);
		dHist_NumBCALShowers->Fill((Double_t)locBCALShowers.size());
		dHist_NumBCALShowers->SetEntries(dHist_NumBCALShowers->GetEntries() - 1 + (Double_t)locBCALShowers.size());

		vector<const DNeutralShower*> locNeutralShowers;
		locEventLoop->Get(locNeutralShowers);
		dHist_NumNeutralShowers->Fill((Double_t)locNeutralShowers.size());
		dHist_NumNeutralShowers->SetEntries(dHist_NumNeutralShowers->GetEntries() - 1 + (Double_t)locNeutralShowers.size());

		//TOF & SC
		vector<const DTOFPoint*> locTOFPoints;
		locEventLoop->Get(locTOFPoints);
		dHist_NumTOFPoints->Fill((Double_t)locTOFPoints.size());
		dHist_NumTOFPoints->SetEntries(dHist_NumTOFPoints->GetEntries() - 1 + (Double_t)locTOFPoints.size());

		vector<const DSCHit*> locSCHits;
		locEventLoop->Get(locSCHits);
		dHist_NumSCHits->Fill((Double_t)locSCHits.size());
		dHist_NumSCHits->SetEntries(dHist_NumSCHits->GetEntries() - 1 + (Double_t)locSCHits.size());

		//Matches
		const DDetectorMatches* locDetectorMatches = NULL;
		locEventLoop->GetSingle(locDetectorMatches);
		if(locDetectorMatches != NULL)
		{
			dHist_NumTrackBCALMatches->Fill((Double_t)locDetectorMatches->Get_NumTrackBCALMatches());
			dHist_NumTrackBCALMatches->SetEntries(dHist_NumTrackBCALMatches->GetEntries() - 1 + (Double_t)locDetectorMatches->Get_NumTrackBCALMatches());
			dHist_NumTrackFCALMatches->Fill((Double_t)locDetectorMatches->Get_NumTrackFCALMatches());
			dHist_NumTrackFCALMatches->SetEntries(dHist_NumTrackFCALMatches->GetEntries() - 1 + (Double_t)locDetectorMatches->Get_NumTrackFCALMatches());
			dHist_NumTrackTOFMatches->Fill((Double_t)locDetectorMatches->Get_NumTrackTOFMatches());
			dHist_NumTrackTOFMatches->SetEntries(dHist_NumTrackTOFMatches->GetEntries() - 1 + (Double_t)locDetectorMatches->Get_NumTrackTOFMatches());
			dHist_NumTrackSCMatches->Fill((Double_t)locDetectorMatches->Get_NumTrackSCMatches());
			dHist_NumTrackSCMatches->SetEntries(dHist_NumTrackSCMatches->GetEntries() - 1 + (Double_t)locDetectorMatches->Get_NumTrackSCMatches());
		}
		else
		{
			dHist_NumTrackBCALMatches->Fill(0.0);
			dHist_NumTrackBCALMatches->SetEntries(dHist_NumTrackBCALMatches->GetEntries() - 1);
			dHist_NumTrackFCALMatches->Fill(0.0);
			dHist_NumTrackFCALMatches->SetEntries(dHist_NumTrackFCALMatches->GetEntries() - 1);
			dHist_NumTrackTOFMatches->Fill(0.0);
			dHist_NumTrackTOFMatches->SetEntries(dHist_NumTrackTOFMatches->GetEntries() - 1);
			dHist_NumTrackSCMatches->Fill(0.0);
			dHist_NumTrackSCMatches->SetEntries(dHist_NumTrackSCMatches->GetEntries() - 1);
		}

		//Hits
		if(!locIsRESTEvent)
		{
			vector<const DCDCHit*> locCDCHits;
			locEventLoop->Get(locCDCHits);
			dHist_NumCDCHits->Fill((Double_t)locCDCHits.size());
			dHist_NumCDCHits->SetEntries(dHist_NumCDCHits->GetEntries() - 1 + (Double_t)locCDCHits.size());

			vector<const DFDCHit*> locFDCHits;
			locEventLoop->Get(locFDCHits);
			dHist_NumFDCHits->Fill((Double_t)locFDCHits.size());
			dHist_NumFDCHits->SetEntries(dHist_NumFDCHits->GetEntries() - 1 + (Double_t)locFDCHits.size());

			vector<const DTOFHit*> locTOFHits;
			locEventLoop->Get(locTOFHits);
			dHist_NumTOFHits->Fill((Double_t)locTOFHits.size());
			dHist_NumTOFHits->SetEntries(dHist_NumTOFHits->GetEntries() - 1 + (Double_t)locTOFHits.size());

			vector<const DBCALHit*> locBCALHits;
			locEventLoop->Get(locBCALHits);
			dHist_NumBCALHits->Fill((Double_t)locBCALHits.size());
			dHist_NumBCALHits->SetEntries(dHist_NumBCALHits->GetEntries() - 1 + (Double_t)locBCALHits.size());

			vector<const DFCALHit*> locFCALHits;
			locEventLoop->Get(locFCALHits);
			dHist_NumFCALHits->Fill((Double_t)locFCALHits.size());
			dHist_NumFCALHits->SetEntries(dHist_NumFCALHits->GetEntries() - 1 + (Double_t)locFCALHits.size());
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

		string locHistName("dHist_NumReconstructedTracks");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_NumReconstructedTracks = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
		else
		{
			dHist_NumReconstructedTracks = new TH2D("NumReconstructedTracks", ";Track Type;Num Tracks / Event", 5 + dFinalStatePIDs.size(), -0.5, 4.5 + dFinalStatePIDs.size(), dMaxNumTracks + 1, -0.5, (float)dMaxNumTracks + 0.5);
			dHist_NumReconstructedTracks->GetXaxis()->SetBinLabel(1, "# Total");
			dHist_NumReconstructedTracks->GetXaxis()->SetBinLabel(2, "# q = +");
			dHist_NumReconstructedTracks->GetXaxis()->SetBinLabel(3, "# q = -");
			dHist_NumReconstructedTracks->GetXaxis()->SetBinLabel(4, "# q = 0");
			dHist_NumReconstructedTracks->GetXaxis()->SetBinLabel(5, "# q != 0");
			for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
			{
				string locLabelName = string("# ") + string(ParticleName_ROOT(dFinalStatePIDs[loc_i]));
				dHist_NumReconstructedTracks->GetXaxis()->SetBinLabel(6 + loc_i, locLabelName.c_str());
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

	vector<const DNeutralParticle*> locNeutralParticles;
	locEventLoop->Get(locNeutralParticles);

	const DMCThrown* locMCThrown;
	Particle_t locPID;
	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
	locEventLoop->Get(locMCThrownMatchingVector);
	const DMCThrownMatching* locMCThrownMatching = locMCThrownMatchingVector.empty() ? NULL : locMCThrownMatchingVector[0];

	size_t locNumPositiveTracks = 0;
	size_t locNumNegativeTracks = 0;
	for(size_t loc_i = 0; loc_i < locChargedTracks.size(); ++loc_i)
	{
		if(locChargedTracks[loc_i]->Get_BestFOM()->charge() > 0.0)
			++locNumPositiveTracks;
		else
			++locNumNegativeTracks;
	}

	// get #tracks by pid type //USES MC PID IF EXISTS!!!
	map<Particle_t, size_t> locNumTracksByPID;

	// charged by pid
	for(size_t loc_i = 0; loc_i < locChargedTracks.size(); ++loc_i)
	{
		if(locMCThrownMatching != NULL)
		{
			double locMatchFOM = 0.0;
			locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locChargedTracks[loc_i], locMatchFOM);
			if((locMCThrown == NULL) || (locMatchFOM < dMinThrownMatchFOM))
				locPID = Unknown;
			else
				locPID = (Particle_t)locMCThrown->type;
		}
		else
			locPID = locChargedTracks[loc_i]->Get_BestFOM()->PID();

		if(locNumTracksByPID.find(locPID) != locNumTracksByPID.end())
			++locNumTracksByPID[locPID];
		else
			locNumTracksByPID[locPID] = 1;
	}

	// neutrals by pid
	for(size_t loc_i = 0; loc_i < locNeutralParticles.size(); ++loc_i)
	{
		if(locMCThrownMatching != NULL)
		{
			double locMatchFOM = 0.0;
			locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locNeutralParticles[loc_i], locMatchFOM);
			if((locMCThrown == NULL) || (locMatchFOM < dMinThrownMatchFOM))
				locPID = Unknown;
			else
				locPID = (Particle_t)locMCThrown->type;
		}
		else
			locPID = locNeutralParticles[loc_i]->Get_BestFOM()->PID();

		if(locNumTracksByPID.find(locPID) != locNumTracksByPID.end())
			++locNumTracksByPID[locPID];
		else
			locNumTracksByPID[locPID] = 1;
	}

	japp->RootWriteLock();
	{
		dHist_NumReconstructedTracks->Fill(0.0, (Double_t)(locChargedTracks.size() + locNeutralParticles.size()));
		dHist_NumReconstructedTracks->Fill(1.0, (Double_t)locNumPositiveTracks);
		dHist_NumReconstructedTracks->Fill(2.0, (Double_t)locNumNegativeTracks);
		dHist_NumReconstructedTracks->Fill(3.0, (Double_t)locNeutralParticles.size());
		dHist_NumReconstructedTracks->Fill(4.0, (Double_t)locChargedTracks.size());
		for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
			dHist_NumReconstructedTracks->Fill(5.0 + (Double_t)loc_i, (Double_t)locNumTracksByPID[dFinalStatePIDs[loc_i]]);
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
				dHistDeque_P_CorrectID[loc_i][locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistDeque_P_CorrectID[loc_i][locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPBins, dMinP, dMaxP);

			//P of Incorrect ID
			locHistName = string("Momentum_IncorrectID_") + locParticleName;
			locHistTitle = string("Incorrect ") + locParticleROOTName + string(" ID, ") + locStepROOTName + string(";p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistDeque_P_IncorrectID[loc_i][locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistDeque_P_IncorrectID[loc_i][locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPBins, dMinP, dMaxP);

			//P Vs Theta of Correct ID
			locHistName = string("PVsTheta_CorrectID_") + locParticleName;
			locHistTitle = string("Correct ") + locParticleROOTName + string(" ID, ") + locStepROOTName + string(";#theta#circ;p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistDeque_PVsTheta_CorrectID[loc_i][locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistDeque_PVsTheta_CorrectID[loc_i][locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNumThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			//P Vs Theta of Incorrect ID
			locHistName = string("PVsTheta_IncorrectID_") + locParticleName;
			locHistTitle = string("Incorrect ") + locParticleROOTName + string(" ID, ") + locStepROOTName + string(";#theta#circ;p (GeV/c)");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistDeque_PVsTheta_IncorrectID[loc_i][locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistDeque_PVsTheta_IncorrectID[loc_i][locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNumThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);
		}
		gDirectory->cd("..");
	} //end of step loop

	//# Combos Pass/Fail All True PID
	locHistName = "Combo_TruePIDStatus";
	locHistTitle = Get_Reaction()->Get_ReactionName() + string(";# Combos;All Combo Particles True PID Status");
	if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
		dHist_TruePIDStatus = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
	else
		dHist_TruePIDStatus = new TH1D(locHistName.c_str(), locHistTitle.c_str(), 2, -0.5, 1.5);

	//# Combos in Signal Region Pass/Fail All True PID
	if(dMinMassSq < dMaxMassSq)
	{
		locHistName = "Combo_TruePIDStatus_SignalRegion";
		locHistTitle = Get_Reaction()->Get_ReactionName() + string(";# Combos in ") + ParticleName_ROOT(dInitialPID) + string(" Signal Region;All Combo Particles True PID Status");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHist_TruePIDStatus_SignalRegion = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_TruePIDStatus_SignalRegion = new TH1D(locHistName.c_str(), locHistTitle.c_str(), 2, -0.5, 1.5);
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
	dHist_TruePIDStatus->Fill(locComboTruePIDStatus);

	if(dMinMassSq < dMaxMassSq)
	{
		bool locInSignalRegionFlag = true; //all possibilities have to be in the signal region
		for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
		{
			const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);
			if(locParticleComboStep->Get_InitialParticleID() != dInitialPID)
				continue;

			DLorentzVector locFinalStateP4 = dAnalysisUtilities->Calc_FinalStateP4(locParticleCombo, loc_i, Get_UseKinFitResultsFlag());
			if((locFinalStateP4.M2() < dMinMassSq) || (locFinalStateP4.M2() > dMaxMassSq))
			{
				locInSignalRegionFlag = false;
				break;
			}
		}
		if(locInSignalRegionFlag)
			dHist_TruePIDStatus_SignalRegion->Fill(locComboTruePIDStatus);
	}
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
			dHist_InvaraintMass = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_InvaraintMass = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumMassBins, dMinMass, dMaxMass);
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
			dHist_MissingMass = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_MissingMass = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumMassBins, dMinMass, dMaxMass);
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
			dHist_MissingMassSquared = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_MissingMassSquared = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumMassBins, dMinMassSq, dMaxMassSq);
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
		dHist_ConfidenceLevel = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
	else
		dHist_ConfidenceLevel = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumConfidenceLevelBins, 0.0, 1.0);

	// Pulls
	map<DKinFitPullType, TH1D*> locParticlePulls;

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
			dHist_RFTimePull = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHist_RFTimePull = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, dMinPull, dMaxPull);

		gDirectory->cd("..");
	}

	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

void DHistogramAction_KinFitResults::Create_ParticlePulls(bool locIsBeamFlag, string locStepROOTName, Particle_t locPID, map<DKinFitPullType, TH1D*>& locParticlePulls, const string& locKinFitTypeString)
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
			locParticlePulls[d_EPull] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			locParticlePulls[d_EPull] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, dMinPull, dMaxPull);
	}
	else
	{
		//Px Pull
		locHistName = "Pull_Px";
		locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(", ") + locKinFitTypeString + string(" Fit;p_{x} Pull;# Combos");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			locParticlePulls[d_PxPull] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			locParticlePulls[d_PxPull] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, dMinPull, dMaxPull);

		//Py Pull
		locHistName = "Pull_Py";
		locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(", ") + locKinFitTypeString + string(" Fit;p_{y} Pull;# Combos");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			locParticlePulls[d_PyPull] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			locParticlePulls[d_PyPull] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, dMinPull, dMaxPull);

		//Pz Pull
		locHistName = "Pull_Pz";
		locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(", ") + locKinFitTypeString + string(" Fit;p_{z} Pull;# Combos");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			locParticlePulls[d_PzPull] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			locParticlePulls[d_PzPull] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, dMinPull, dMaxPull);
	}

	//vertex pulls:
	if(locNeutralShowerFlag || (locKinFitType == d_VertexFit) || (locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit))
	{
		//Xx Pull
		locHistName = "Pull_Xx";
		locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(", ") + locKinFitTypeString + string(" Fit;x_{x} Pull;# Combos");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			locParticlePulls[d_XxPull] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			locParticlePulls[d_XxPull] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, dMinPull, dMaxPull);

		//Xy Pull
		locHistName = "Pull_Xy";
		locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(", ") + locKinFitTypeString + string(" Fit;x_{y} Pull;# Combos");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			locParticlePulls[d_XyPull] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			locParticlePulls[d_XyPull] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, dMinPull, dMaxPull);

		//Xz Pull
		locHistName = "Pull_Xz";
		locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(", ") + locKinFitTypeString + string(" Fit;x_{z} Pull;# Combos");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			locParticlePulls[d_XzPull] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			locParticlePulls[d_XzPull] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, dMinPull, dMaxPull);
	}

	//time pulls:
	if((locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndSpacetimeFit))
	{
		//T Pull
		locHistName = "Pull_T";
		locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(", ") + locKinFitTypeString + string(" Fit;t Pull;# Combos");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			locParticlePulls[d_TPull] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			locParticlePulls[d_TPull] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumPullBins, dMinPull, dMaxPull);
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
		for(locIterator = locParticlePulls.begin(); locIterator != locParticlePulls.end(); ++locIterator)
			dHistMap_BeamPulls[locIterator->first]->Fill(locIterator->second);
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
			for(locIterator = locParticlePulls.begin(); locIterator != locParticlePulls.end(); ++locIterator)
				(dHistMap_Pulls[locParticlePair])[locIterator->first]->Fill(locIterator->second);
			japp->RootUnLock();
		}
	}

	//rf time pull
	DKinFitType locKinFitType = Get_Reaction()->Get_KinFitType();
	if(locBeamFlag && ((locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndSpacetimeFit)))
	{
		locParticlePulls = locPulls[NULL];
		japp->RootWriteLock();
		dHist_RFTimePull->Fill(locParticlePulls[d_TPull]);
		japp->RootUnLock();
	}

	return true;
}

