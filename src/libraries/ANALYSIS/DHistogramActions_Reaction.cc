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

			if(ParticleCharge(locPID) == 0)
			{
				//BCAL
				CreateAndChangeTo_Directory("BCAL", "BCAL");
				locHistName = "Beta";
				locHistTitle =  string("BCAL ") + locParticleROOTName + string(" Candidates;#beta");
				dHistMap_Beta[locPID][SYS_BCAL] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumBetaBins, dMinBeta, dMaxBeta);
				gDirectory->cd("..");

				//FCAL
				CreateAndChangeTo_Directory("FCAL", "FCAL");
				locHistName = "Beta";
				locHistTitle =  string("FCAL ") + locParticleROOTName + string(" Candidates;#beta");
				dHistMap_Beta[locPID][SYS_FCAL] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumBetaBins, dMinBeta, dMaxBeta);
				gDirectory->cd("..");
			}

			//q = 0
			if(locPID == Gamma)
			{
				//BCAL
				CreateAndChangeTo_Directory("BCAL", "BCAL");

				locHistName = "BetaVsP";
				locHistTitle =  string("BCAL ") + locParticleROOTName + string(" Candidates;Shower Energy (GeV);#beta");
				dHistMap_BetaVsP[locPID][SYS_BCAL] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxBCALP, dNum2DBetaBins, dMinBeta, dMaxBeta);

				locHistName = "DeltaTVsShowerE_Photon";
				locHistTitle = string("BCAL ") + locParticleROOTName + string(" Candidates;Shower Energy (GeV);#Deltat_{BCAL - RF}");
				dHistMap_DeltaTVsP[locPID][SYS_BCAL] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);

				locHistName = "TimePullVsShowerE_Photon";
				locHistTitle = string("BCAL ") + locParticleROOTName + string(";Shower Energy (GeV);#Deltat/#sigma_{#Deltat}");
				dHistMap_TimePullVsP[Gamma][SYS_BCAL] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, dMinPull, dMaxPull);

				locHistName = "TimeFOMVsShowerE_Photon";
				locHistTitle = string("BCAL ") + locParticleROOTName + string(";Shower Energy (GeV);Timing PID Confidence Level");
				dHistMap_TimeFOMVsP[Gamma][SYS_BCAL] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DFOMBins, 0.0, 1.0);

				gDirectory->cd("..");

				//FCAL
				CreateAndChangeTo_Directory("FCAL", "FCAL");

				locHistName = "BetaVsP";
				locHistTitle =  string("FCAL ") + locParticleROOTName + string(" Candidates;Shower Energy (GeV);#beta");
				dHistMap_BetaVsP[locPID][SYS_FCAL] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DBetaBins, dMinBeta, dMaxBeta);

				locHistName = "DeltaTVsShowerE_Photon";
				locHistTitle = string("FCAL ") + locParticleROOTName + string(";p (GeV/c);#Deltat_{FCAL - RF}");
				dHistMap_DeltaTVsP[locPID][SYS_FCAL] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);

				locHistName = "TimePullVsShowerE_Photon";
				locHistTitle = string("FCAL ") + locParticleROOTName + string(";Shower Energy (GeV);#Deltat/#sigma_{#Deltat}");
				dHistMap_TimePullVsP[Gamma][SYS_FCAL] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, dMinPull, dMaxPull);

				locHistName = "TimeFOMVsShowerE_Photon";
				locHistTitle = string("FCAL ") + locParticleROOTName + string(";Shower Energy (GeV);Timing PID Confidence Level");
				dHistMap_TimeFOMVsP[Gamma][SYS_FCAL] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DFOMBins, 0.0, 1.0);

				gDirectory->cd("..");
			}
			else if(ParticleCharge(locPID) != 0)
			{
				//SC
				CreateAndChangeTo_Directory("SC", "SC");

				locHistName = "dEdXVsP";
				locHistTitle = locParticleROOTName + string(" Candidates;p (GeV/c);SC dE/dX (MeV/cm)");
				dHistMap_dEdXVsP[locPID][SYS_START] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DdEdxBins, dMindEdX, dMaxdEdX);

				locHistName = "DeltadEdXVsP";
				locHistTitle = locParticleROOTName + string(" Candidates;p (GeV/c);SC #Delta(dE/dX) (MeV/cm)");
				dHistMap_DeltadEdXVsP[locPID][SYS_START] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltadEdxBins, dMinDeltadEdx, dMaxDeltadEdx);

				gDirectory->cd("..");

				//TOF
				CreateAndChangeTo_Directory("TOF", "TOF");

				locHistName = "dEdXVsP";
				locHistTitle = locParticleROOTName + string(" Candidates;p (GeV/c);TOF dE/dX (MeV/cm)");
				dHistMap_dEdXVsP[locPID][SYS_TOF] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DdEdxBins, dMindEdX, dMaxdEdX);

				locHistName = "DeltadEdXVsP";
				locHistTitle = locParticleROOTName + string(" Candidates;p (GeV/c);TOF #Delta(dE/dX) (MeV/cm)");
				dHistMap_DeltadEdXVsP[locPID][SYS_TOF] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltadEdxBins, dMinDeltadEdx, dMaxDeltadEdx);

				locHistName = "DeltaBetaVsP";
				locHistTitle = locParticleROOTName + string(" Candidates;p (GeV/c);TOF #Delta#beta");
				dHistMap_DeltaBetaVsP[locPID][SYS_TOF] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltaBetaBins, dMinDeltaBeta, dMaxDeltaBeta);

				locHistName = "BetaVsP";
				locHistTitle = locParticleROOTName + string(" Candidates;p (GeV/c);TOF #beta");
				dHistMap_BetaVsP[locPID][SYS_TOF] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DBetaBins, dMinBeta, dMaxBeta);

				locHistName = "DeltaTVsP";
				locHistTitle = locParticleROOTName + string(" Candidates;p (GeV/c);#Deltat_{TOF - RF}");
				dHistMap_DeltaTVsP[locPID][SYS_TOF] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);

				locHistName = string("TimePullVsP_") + locParticleName;
				locHistTitle = string("TOF ") + locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
				dHistMap_TimePullVsP[locPID][SYS_TOF] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, dMinPull, dMaxPull);

				locHistName = string("TimeFOMVsP_") + locParticleName;
				locHistTitle = string("TOF ") + locParticleROOTName + string(";p (GeV/c);Timing PID Confidence Level");
				dHistMap_TimeFOMVsP[locPID][SYS_TOF] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DFOMBins, 0.0, 1.0);

				gDirectory->cd("..");

				//BCAL
				CreateAndChangeTo_Directory("BCAL", "BCAL");

				locHistName = "BetaVsP";
				locHistTitle = locParticleROOTName + string(" Candidates;p (GeV/c);BCAL #beta");
				dHistMap_BetaVsP[locPID][SYS_BCAL] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxBCALP, dNum2DBetaBins, dMinBeta, dMaxBeta);

				locHistName = "DeltaBetaVsP";
				locHistTitle = locParticleROOTName + string(" Candidates;p (GeV/c);BCAL #Delta#beta");
				dHistMap_DeltaBetaVsP[locPID][SYS_BCAL] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxBCALP, dNum2DDeltaBetaBins, dMinDeltaBeta, dMaxDeltaBeta);

				locHistName = "DeltaTVsP";
				locHistTitle = locParticleROOTName + string(" Candidates;p (GeV/c);#Deltat_{BCAL - RF}");
				dHistMap_DeltaTVsP[locPID][SYS_BCAL] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);

				locHistName = string("TimePullVsP_") + locParticleName;
				locHistTitle = string("BCAL ") + locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
				dHistMap_TimePullVsP[locPID][SYS_BCAL] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, dMinPull, dMaxPull);

				locHistName = string("TimeFOMVsP_") + locParticleName;
				locHistTitle = string("BCAL ") + locParticleROOTName + string(";p (GeV/c);Timing PID Confidence Level");
				dHistMap_TimeFOMVsP[locPID][SYS_BCAL] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DFOMBins, 0.0, 1.0);

				locHistName = "EOverPVsP";
				locHistTitle = locParticleROOTName + string(" Candidates;p (GeV/c);BCAL E_{Shower}/p_{Track} (c);");
				dHistMap_EOverPVsP[locPID][SYS_BCAL] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxBCALP, dNum2DEOverPBins, dMinEOverP, dMaxEOverP);

				locHistName = "EOverPVsTheta";
				locHistTitle = locParticleROOTName + string(" Candidates;#theta#circ;BCAL E_{Shower}/p_{Track} (c);");
				dHistMap_EOverPVsTheta[locPID][SYS_BCAL] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DBCALThetaBins, dMinBCALTheta, dMaxBCALTheta, dNum2DEOverPBins, dMinEOverP, dMaxEOverP);

				gDirectory->cd("..");

				//FCAL
				CreateAndChangeTo_Directory("FCAL", "FCAL");

				locHistName = "BetaVsP";
				locHistTitle = locParticleROOTName + string(" Candidates;p (GeV/c);FCAL #beta");
				dHistMap_BetaVsP[locPID][SYS_FCAL] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DBetaBins, dMinBeta, dMaxBeta);

				locHistName = "DeltaBetaVsP";
				locHistTitle = locParticleROOTName + string(" Candidates;p (GeV/c);FCAL #Delta#beta");
				dHistMap_DeltaBetaVsP[locPID][SYS_FCAL] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltaBetaBins, dMinDeltaBeta, dMaxDeltaBeta);

				locHistName = "DeltaTVsP";
				locHistTitle = locParticleROOTName + string(" Candidates;p (GeV/c);#Deltat_{FCAL - RF}");
				dHistMap_DeltaTVsP[locPID][SYS_FCAL] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltaTBins, dMinDeltaT, dMaxDeltaT);

				locHistName = string("TimePullVsP_") + locParticleName;
				locHistTitle = string("FCAL ") + locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
				dHistMap_TimePullVsP[locPID][SYS_FCAL] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, dMinPull, dMaxPull);

				locHistName = string("TimeFOMVsP_") + locParticleName;
				locHistTitle = string("FCAL ") + locParticleROOTName + string(";p (GeV/c);Timing PID Confidence Level");
				dHistMap_TimeFOMVsP[locPID][SYS_FCAL] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DFOMBins, 0.0, 1.0);

				locHistName = "EOverPVsP";
				locHistTitle = locParticleROOTName + string(" Candidates;p (GeV/c);FCAL E_{Shower}/p_{Track} (c);");
				dHistMap_EOverPVsP[locPID][SYS_FCAL] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DEOverPBins, dMinEOverP, dMaxEOverP);

				locHistName = "EOverPVsTheta";
				locHistTitle = locParticleROOTName + string(" Candidates;#theta#circ;FCAL E_{Shower}/p_{Track} (c);");
				dHistMap_EOverPVsTheta[locPID][SYS_FCAL] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DFCALThetaBins, dMinFCALTheta, dMaxFCALTheta, dNum2DEOverPBins, dMinEOverP, dMaxEOverP);

				gDirectory->cd("..");

				//CDC
				CreateAndChangeTo_Directory("CDC", "CDC");

				locHistName = "dEdXVsP";
				locHistTitle = locParticleROOTName + string(" Candidates;p (GeV/c);CDC dE/dx (keV/cm)");
				dHistMap_dEdXVsP[locPID][SYS_CDC] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DdEdxBins, dMindEdX, dMaxdEdX);

				locHistName = "DeltadEdXVsP";
				locHistTitle = locParticleROOTName + string(" Candidates;p (GeV/c);CDC #Delta(dE/dX) (keV/cm)");
				dHistMap_DeltadEdXVsP[locPID][SYS_CDC] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltadEdxBins, dMinDeltadEdx, dMaxDeltadEdx);

				locHistName = string("dEdXPullVsP_") + locParticleName;
				locHistTitle = locParticleROOTName + string(";p (GeV/c);CDC #Delta(dE/dX)/#sigma_{#Delta(dE/dX)}");
				dHistMap_dEdXPullVsP[locPID][SYS_CDC] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, dMinPull, dMaxPull);

				locHistName = string("dEdXFOMVsP_") + locParticleName;
				locHistTitle = locParticleROOTName + string(";p (GeV/c);CDC dE/dx PID Confidence Level");
				dHistMap_dEdXFOMVsP[locPID][SYS_CDC] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DFOMBins, 0.0, 1.0);

				gDirectory->cd("..");

				//FDC
				CreateAndChangeTo_Directory("FDC", "FDC");

				locHistName = "dEdXVsP";
				locHistTitle = locParticleROOTName + string(" Candidates;p (GeV/c);FDC dE/dx (keV/cm)");
				dHistMap_dEdXVsP[locPID][SYS_FDC] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DdEdxBins, dMindEdX, dMaxdEdX);

				locHistName = "DeltadEdXVsP";
				locHistTitle = locParticleROOTName + string(" Candidates;p (GeV/c);FDC #Delta(dE/dX) (keV/cm)");
				dHistMap_DeltadEdXVsP[locPID][SYS_FDC] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltadEdxBins, dMinDeltadEdx, dMaxDeltadEdx);

				locHistName = string("dEdXPullVsP_") + locParticleName;
				locHistTitle = locParticleROOTName + string(";p (GeV/c);FDC #Delta(dE/dX)/#sigma_{#Delta(dE/dX)}");
				dHistMap_dEdXPullVsP[locPID][SYS_FDC] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, dMinPull, dMaxPull);

				locHistName = string("dEdXFOMVsP_") + locParticleName;
				locHistTitle = locParticleROOTName + string(";p (GeV/c);FDC dE/dx PID Confidence Level");
				dHistMap_dEdXFOMVsP[locPID][SYS_FDC] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DFOMBins, 0.0, 1.0);

				gDirectory->cd("..");
			} //end of charged

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
					dHistMap_PIDFOMForTruePID[locPIDPair] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumFOMBins, 0.0, 1.0);
				}
				gDirectory->cd("..");
			}

			//no FOM for massive neutrals
			if((ParticleCharge(locPID) != 0) || (locPID == Gamma))
			{
				// Overall Confidence Level
				locHistName = "PIDConfidenceLevel";
				locHistTitle = locParticleROOTName + string(" PID;PID Confidence Level");
				dHistMap_PIDFOM[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumFOMBins, 0.0, 1.0);

				// P Vs Theta, PID FOM = NaN
				locHistName = "PVsTheta_NaNPIDFOM";
				locHistTitle = locParticleROOTName + string(", PID FOM = NaN;#theta#circ;p (GeV/c)");
				dHistMap_PVsTheta_NaNPIDFOM[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);
			}

			gDirectory->cd("..");
		} //end of PID loop

		//Return to the base directory
		ChangeTo_BaseDirectory();
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_PID::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	if(Get_NumPreviousParticleCombos() == 0)
		dPreviouslyHistogrammedParticles.clear();

	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
	locEventLoop->Get(locMCThrownMatchingVector);
	const DMCThrownMatching* locMCThrownMatching = locMCThrownMatchingVector.empty() ? NULL : locMCThrownMatchingVector[0];

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
	double locBeta_Timing = locChargedTrackHypothesis->measuredBeta();
	double locDeltaBeta = locChargedTrackHypothesis->deltaBeta();
	double locDeltaT = (locChargedTrackHypothesis->time() - locChargedTrackHypothesis->t0());
	double locFOM_Timing = (locChargedTrackHypothesis->dNDF_Timing > 0) ? TMath::Prob(locChargedTrackHypothesis->dChiSq_Timing, locChargedTrackHypothesis->dNDF_Timing) : numeric_limits<double>::quiet_NaN();

	double locP = locChargedTrackHypothesis->momentum().Mag();
	double locTheta = locChargedTrackHypothesis->momentum().Theta()*180.0/TMath::Pi();
	double locMatchFOM = 0.0;
	const DMCThrown* locMCThrown = (locMCThrownMatching != NULL) ? locMCThrownMatching->Get_MatchingMCThrown(locChargedTrackHypothesis, locMatchFOM) : NULL;

	const DBCALShowerMatchParams* locBCALShowerMatchParams = locChargedTrackHypothesis->Get_BCALShowerMatchParams();
	const DFCALShowerMatchParams* locFCALShowerMatchParams = locChargedTrackHypothesis->Get_FCALShowerMatchParams();
	const DTOFHitMatchParams* locTOFHitMatchParams = locChargedTrackHypothesis->Get_TOFHitMatchParams();
	const DSCHitMatchParams* locSCHitMatchParams = locChargedTrackHypothesis->Get_SCHitMatchParams();

	const DTrackTimeBased* locTrackTimeBased = NULL;
	locChargedTrackHypothesis->GetSingle(locTrackTimeBased);

	double locTimePull = 0.0;
	unsigned int locTimeNDF = 0;
	dParticleID->Calc_TimingChiSq(locChargedTrackHypothesis, locTimeNDF, locTimePull);
	DetectorSystem_t locTimeDetector = locChargedTrackHypothesis->t1_detector();

	japp->RootWriteLock();
	{
		dHistMap_PIDFOM[locPID]->Fill(locChargedTrackHypothesis->dFOM);

		//SC dE/dx
		if(locSCHitMatchParams != NULL)
		{
			dHistMap_dEdXVsP[locPID][SYS_START]->Fill(locP, locSCHitMatchParams->dEdx*1.0E3);
			double locdx = locSCHitMatchParams->dHitEnergy/locSCHitMatchParams->dEdx;
			double locProbabledEdx = 0.0, locSigmadEdx = 0.0;
			dParticleID->GetScintMPdEandSigma(locP, locChargedTrackHypothesis->mass(), locdx, locProbabledEdx, locSigmadEdx);
			dHistMap_DeltadEdXVsP[locPID][SYS_START]->Fill(locP, (locSCHitMatchParams->dEdx - locProbabledEdx)*1.0E3);
		}

		//TOF dE/dx
		if(locTOFHitMatchParams != NULL)
		{
			dHistMap_dEdXVsP[locPID][SYS_TOF]->Fill(locP, locTOFHitMatchParams->dEdx*1.0E3);
			double locdx = locTOFHitMatchParams->dHitEnergy/locTOFHitMatchParams->dEdx;
			double locProbabledEdx = 0.0, locSigmadEdx = 0.0;
			dParticleID->GetScintMPdEandSigma(locP, locChargedTrackHypothesis->mass(), locdx, locProbabledEdx, locSigmadEdx);
			dHistMap_DeltadEdXVsP[locPID][SYS_TOF]->Fill(locP, (locTOFHitMatchParams->dEdx - locProbabledEdx)*1.0E3);
		}

		//BCAL E/p
		if(locBCALShowerMatchParams != NULL)
		{
			const DBCALShower* locBCALShower = locBCALShowerMatchParams->dBCALShower;
			double locEOverP = locBCALShower->E/locP;
			dHistMap_EOverPVsP[locPID][SYS_BCAL]->Fill(locP, locEOverP);
			dHistMap_EOverPVsTheta[locPID][SYS_BCAL]->Fill(locTheta, locEOverP);
		}

		//FCAL E/p
		if(locFCALShowerMatchParams != NULL)
		{
			const DFCALShower* locFCALShower = locFCALShowerMatchParams->dFCALShower;
			double locEOverP = locFCALShower->getEnergy()/locP;
			dHistMap_EOverPVsP[locPID][SYS_FCAL]->Fill(locP, locEOverP);
			dHistMap_EOverPVsTheta[locPID][SYS_FCAL]->Fill(locTheta, locEOverP);
		}

		//Timing
		if((locTimeDetector == SYS_TOF) || (locTimeDetector == SYS_BCAL) || (locTimeDetector == SYS_FCAL))
		{
			dHistMap_BetaVsP[locPID][locTimeDetector]->Fill(locP, locBeta_Timing);
			dHistMap_DeltaBetaVsP[locPID][locTimeDetector]->Fill(locP, locDeltaBeta);
			dHistMap_DeltaTVsP[locPID][locTimeDetector]->Fill(locP, locDeltaT);
			dHistMap_TimePullVsP[locPID][locTimeDetector]->Fill(locP, locTimePull);
			dHistMap_TimeFOMVsP[locPID][locTimeDetector]->Fill(locP, locFOM_Timing);
		}

		//CDC dE/dx
		if(locTrackTimeBased->dNumHitsUsedFordEdx_CDC > 0)
		{
			dHistMap_dEdXVsP[locPID][SYS_CDC]->Fill(locP, locTrackTimeBased->ddEdx_CDC*1.0E6);
			double locProbabledEdx = dParticleID->GetMostProbabledEdx_DC(locP, locChargedTrackHypothesis->mass(), locTrackTimeBased->ddx_CDC, true);
			double locDeltadEdx = locTrackTimeBased->ddEdx_CDC - locProbabledEdx;
			dHistMap_DeltadEdXVsP[locPID][SYS_CDC]->Fill(locP, 1.0E6*locDeltadEdx);
			double locMeandx = locTrackTimeBased->ddx_CDC/locTrackTimeBased->dNumHitsUsedFordEdx_CDC;
			double locSigmadEdx = dParticleID->GetdEdxSigma_DC(locTrackTimeBased->dNumHitsUsedFordEdx_CDC, locP, locChargedTrackHypothesis->mass(), locMeandx, true);
			double locdEdXPull = locDeltadEdx/locSigmadEdx;
			dHistMap_dEdXPullVsP[locPID][SYS_CDC]->Fill(locP, locDeltadEdx/locSigmadEdx);
			double locdEdXChiSq = locdEdXPull*locdEdXPull;
			double locdEdXFOM = TMath::Prob(locdEdXChiSq, locTrackTimeBased->dNumHitsUsedFordEdx_CDC);
			dHistMap_dEdXFOMVsP[locPID][SYS_CDC]->Fill(locP, locdEdXFOM);
		}

		//FDC dE/dx
		if(locTrackTimeBased->dNumHitsUsedFordEdx_FDC > 0)
		{
			dHistMap_dEdXVsP[locPID][SYS_FDC]->Fill(locP, locTrackTimeBased->ddEdx_FDC*1.0E6);
			double locProbabledEdx = dParticleID->GetMostProbabledEdx_DC(locP, locChargedTrackHypothesis->mass(), locTrackTimeBased->ddx_FDC, false);
			double locDeltadEdx = locTrackTimeBased->ddEdx_FDC - locProbabledEdx;
			dHistMap_DeltadEdXVsP[locPID][SYS_FDC]->Fill(locP, 1.0E6*locDeltadEdx);
			double locMeandx = locTrackTimeBased->ddx_FDC/locTrackTimeBased->dNumHitsUsedFordEdx_FDC;
			double locSigmadEdx = dParticleID->GetdEdxSigma_DC(locTrackTimeBased->dNumHitsUsedFordEdx_FDC, locP, locChargedTrackHypothesis->mass(), locMeandx, false);
			double locdEdXPull = locDeltadEdx/locSigmadEdx;
			dHistMap_dEdXPullVsP[locPID][SYS_FDC]->Fill(locP, locDeltadEdx/locSigmadEdx);
			double locdEdXChiSq = locdEdXPull*locdEdXPull;
			double locdEdXFOM = TMath::Prob(locdEdXChiSq, locTrackTimeBased->dNumHitsUsedFordEdx_FDC);
			dHistMap_dEdXFOMVsP[locPID][SYS_FDC]->Fill(locP, locdEdXFOM);
		}

		pair<Particle_t, Particle_t> locPIDPair(locPID, Unknown); //default unless matched
		if(locMCThrown != NULL) //else bogus track (not matched to any thrown tracks)
			locPIDPair.second = (Particle_t)(locMCThrown->type); //matched
		if(dHistMap_PIDFOMForTruePID.find(locPIDPair) != dHistMap_PIDFOMForTruePID.end()) //else hist not created or PID is weird
			dHistMap_PIDFOMForTruePID[locPIDPair]->Fill(locChargedTrackHypothesis->dFOM);

		if(locChargedTrackHypothesis->dNDF == 0) //NaN
			dHistMap_PVsTheta_NaNPIDFOM[locPID]->Fill(locTheta, locP);
	}
	japp->RootUnLock();
}

void DHistogramAction_PID::Fill_NeutralHists(const DNeutralParticleHypothesis* locNeutralParticleHypothesis, const DMCThrownMatching* locMCThrownMatching, const DEventRFBunch* locEventRFBunch)
{
	Particle_t locPID = locNeutralParticleHypothesis->PID();

	double locBeta_Timing = locNeutralParticleHypothesis->measuredBeta();
	double locDeltaT = (locNeutralParticleHypothesis->t0() - locNeutralParticleHypothesis->time());

	double locP = locNeutralParticleHypothesis->momentum().Mag();
	double locMatchFOM = 0.0;
	const DMCThrown* locMCThrown = (locMCThrownMatching != NULL) ? locMCThrownMatching->Get_MatchingMCThrown(locNeutralParticleHypothesis, locMatchFOM) : NULL;

	double locTimePull = 0.0;
	unsigned int locTimeNDF = 0;
	dParticleID->Calc_TimingChiSq(locNeutralParticleHypothesis, locTimeNDF, locTimePull);

	DetectorSystem_t locSystem = locNeutralParticleHypothesis->t1_detector();
	japp->RootWriteLock();
	{
		//Beta (good for all PIDs)
		dHistMap_Beta[locPID][locSystem]->Fill(locBeta_Timing);
		if(locPID != Gamma)
		{
			japp->RootUnLock();
			return;
		}

		dHistMap_PIDFOM[locPID]->Fill(locNeutralParticleHypothesis->dFOM);
		dHistMap_BetaVsP[locPID][locSystem]->Fill(locP, locBeta_Timing);
		dHistMap_DeltaTVsP[locPID][locSystem]->Fill(locP, locDeltaT);
		dHistMap_TimePullVsP[locPID][locSystem]->Fill(locP, locTimePull);
		dHistMap_TimeFOMVsP[locPID][locSystem]->Fill(locP, locNeutralParticleHypothesis->dFOM);

		pair<Particle_t, Particle_t> locPIDPair(locPID, Unknown); //default unless matched
		if(locMCThrown != NULL) //else bogus track (not matched to any thrown tracks)
			locPIDPair.second = (Particle_t)(locMCThrown->type); //matched
		if(dHistMap_PIDFOMForTruePID.find(locPIDPair) != dHistMap_PIDFOMForTruePID.end()) //else hist not created or PID is weird
			dHistMap_PIDFOMForTruePID[locPIDPair]->Fill(locNeutralParticleHypothesis->dFOM);
	}
	japp->RootUnLock();
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
	{
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
			dHistDeque_MaxTrackDeltaZ[loc_i] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaVertexZBins, 0.0, dMaxDeltaVertexZ);

			// Max Track DeltaT
			locHistName = "MaxTrackDeltaT";
			locHistTitle = locStepROOTName + string(";Largest Track #DeltaVertex-T (ns)");
			dHistDeque_MaxTrackDeltaT[loc_i] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaVertexTBins, 0.0, dMaxDeltaVertexT);

			// Max Track DOCA
			locHistName = "MaxTrackDOCA";
			locHistTitle = locStepROOTName + string(";Largest Track DOCA (cm)");
			dHistDeque_MaxTrackDOCA[loc_i] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDOCABins, 0.0, dMaxDOCA);

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
				dHistDeque_TrackZToCommon[loc_i][locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaVertexZBins, dMinDeltaVertexZ, dMaxDeltaVertexZ);

				// TrackT To Common
				locHistName = string("TrackTToCommon_") + locParticleName;
				locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";#DeltaVertex-T (Track, Common) (ns)");
				dHistDeque_TrackTToCommon[loc_i][locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaVertexTBins, dMinDeltaVertexT, dMaxDeltaVertexT);

				// TrackDOCA To Common
				locHistName = string("TrackDOCAToCommon_") + locParticleName;
				locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";DOCA (Track, Common) (cm)");
				dHistDeque_TrackDOCAToCommon[loc_i][locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDOCABins, dMinDOCA, dMaxDOCA);

				// DeltaT Vs P against beam photon
				if((locReactionStep->Get_InitialParticleID() == Gamma) && (dHistMap_BeamTrackDeltaTVsP.find(locPID) == dHistMap_BeamTrackDeltaTVsP.end()))
				{
					locHistName = string("TrackDeltaTVsP_") + ParticleType(locPID) + string("_Beam") + ParticleType(Gamma);
					locHistTitle = locStepROOTName + string(";") + ParticleName_ROOT(locPID) + string(" Momentum (GeV/c);t_{") + ParticleName_ROOT(locPID) + string("} - t_{Beam ") + ParticleName_ROOT(Gamma) + string("} (ns)");
					dHistMap_BeamTrackDeltaTVsP[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNumDeltaVertexTBins, dMinDeltaVertexT, dMaxDeltaVertexT);
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
					dHistDeque_TrackDeltaTVsP[loc_i][locParticlePair] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNumDeltaVertexTBins, dMinDeltaVertexT, dMaxDeltaVertexT);
				}
			}
			gDirectory->cd("..");
		}
		//Return to the base directory
		ChangeTo_BaseDirectory();
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
			dBeamParticleHist_P = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPBins, dMinP, dMaxP);

			// Theta
			locHistName = "Theta";
			locHistTitle = string("Beam ") + locParticleROOTName + string(";#theta#circ");
			dBeamParticleHist_Theta = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumThetaBins, dMinTheta, dMaxTheta);

			// Phi
			locHistName = "Phi";
			locHistTitle = string("Beam ") + locParticleROOTName + string(";#phi#circ");
			dBeamParticleHist_Phi = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPhiBins, dMinPhi, dMaxPhi);

			// P Vs Theta
			locHistName = "PVsTheta";
			locHistTitle = string("Beam ") + locParticleROOTName + string(";#theta#circ;p (GeV/c)");
			dBeamParticleHist_PVsTheta = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			// Phi Vs Theta
			locHistName = "PhiVsTheta";
			locHistTitle = string("Beam ") + locParticleROOTName + string(";#theta#circ;#phi#circ");
			dBeamParticleHist_PhiVsTheta = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

			// Vertex-Z
			locHistName = "VertexZ";
			locHistTitle = string("Beam ") + locParticleROOTName + string(";Vertex-Z (cm)");
			dBeamParticleHist_VertexZ = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumVertexZBins, dMinVertexZ, dMaxVertexZ);

			// Vertex-Y Vs Vertex-X
			locHistName = "VertexYVsX";
			locHistTitle = string("Beam ") + locParticleROOTName + string(";Vertex-X (cm);Vertex-Y (cm)");
			dBeamParticleHist_VertexYVsX = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY);

			// Vertex-T
			locHistName = "VertexT";
			locHistTitle = string("Beam ") + locParticleROOTName + string(";Vertex-T (ns)");
			dBeamParticleHist_VertexT = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumTBins, dMinT, dMaxT);

			// Delta-T (Beam, RF)
			locHistName = "DeltaTRF";
			locHistTitle = string("Beam ") + locParticleROOTName + string(";#Deltat_{Beam - RF} (ns)");
			dBeamParticleHist_DeltaTRF = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaTRFBins, dMinDeltaTRF, dMaxDeltaTRF);

			// Delta-T (Beam, RF) Vs Beam E
			locHistName = "DeltaTRFVsBeamE";
			locHistTitle = string("Beam ") + locParticleROOTName + string(";E (GeV);#Deltat_{Beam - RF} (ns)");
			dBeamParticleHist_DeltaTRFVsBeamE = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNumDeltaTRFBins, dMinDeltaTRF, dMaxDeltaTRF);

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
				dHistDeque_P[loc_i][locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPBins, dMinP, dMaxP);

				// Theta
				locHistName = "Theta";
				locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";#theta#circ");
				dHistDeque_Theta[loc_i][locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumThetaBins, dMinTheta, dMaxTheta);

				// Phi
				locHistName = "Phi";
				locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";#phi#circ");
				dHistDeque_Phi[loc_i][locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPhiBins, dMinPhi, dMaxPhi);

				// P Vs Theta
				locHistName = "PVsTheta";
				locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";#theta#circ;p (GeV/c)");
				dHistDeque_PVsTheta[loc_i][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

				// Phi Vs Theta
				locHistName = "PhiVsTheta";
				locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";#theta#circ;#phi#circ");
				dHistDeque_PhiVsTheta[loc_i][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

				//beta vs p
				locHistName = "BetaVsP";
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#beta");
				dHistDeque_BetaVsP[loc_i][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNumBetaBins, dMinBeta, dMaxBeta);

				//delta-beta vs p
				locHistName = "DeltaBetaVsP";
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta#beta");
				dHistDeque_DeltaBetaVsP[loc_i][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNumDeltaBetaBins, dMinDeltaBeta, dMaxDeltaBeta);

				// Vertex-Z
				locHistName = "VertexZ";
				locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";Vertex-Z (cm)");
				dHistDeque_VertexZ[loc_i][locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumVertexZBins, dMinVertexZ, dMaxVertexZ);

				// Vertex-Y Vs Vertex-X
				locHistName = "VertexYVsX";
				locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";Vertex-X (cm);Vertex-Y (cm)");
				dHistDeque_VertexYVsX[loc_i][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY);

				// Vertex-T
				locHistName = "VertexT";
				locHistTitle = locParticleROOTName + string(", ") + locStepROOTName + string(";Vertex-T (ns)");
				dHistDeque_VertexT[loc_i][locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumTBins, dMinT, dMaxT);

				gDirectory->cd("..");
			} //end of particle loop
			gDirectory->cd("..");
		} //end of step loop

		//Return to the base directory
		ChangeTo_BaseDirectory();
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

	double locBeta_Timing = locKinematicData->measuredBeta();
	double locDeltaBeta = locKinematicData->deltaBeta();

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

void DHistogramAction_InvariantMass::Initialize(JEventLoop* locEventLoop)
{
	string locHistName, locHistTitle;
	double locMassPerBin = 1000.0*(dMaxMass - dMinMass)/((double)dNumMassBins);

	string locParticleNamesForHist = Get_Reaction()->Get_DecayChainFinalParticlesROOTNames(dInitialPID, Get_UseKinFitResultsFlag());

	locEventLoop->GetSingle(dAnalysisUtilities);

	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		CreateAndChangeTo_ActionDirectory();

		locHistName = "InvariantMass";
		ostringstream locStream;
		locStream << locMassPerBin;
		locHistTitle = string(";") + locParticleNamesForHist + string(" Invariant Mass (GeV/c^{2});# Combos / ") + locStream.str() + string(" MeV/c^{2}");
		dHist_InvaraintMass = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumMassBins, dMinMass, dMaxMass);

		//Return to the base directory
		ChangeTo_BaseDirectory();
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

	locEventLoop->GetSingle(dAnalysisUtilities);

	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		CreateAndChangeTo_ActionDirectory();

		string locHistName = "MissingMass";
		ostringstream locStream;
		locStream << locMassPerBin;
		string locHistTitle = string(";") + locInitialParticlesROOTName + string("#rightarrow") + locFinalParticlesROOTName + string(" Missing Mass (GeV/c^{2});# Combos / ") + locStream.str() + string(" MeV/c^{2}");
		dHist_MissingMass = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumMassBins, dMinMass, dMaxMass);

		locHistName = "MissingMassVsBeamE";
		locMassPerBin *= ((double)dNumMassBins)/((double)dNum2DMassBins);
		locStream.str("");
		locStream << locMassPerBin;
		locHistTitle = string(";Beam Energy (GeV);") + locInitialParticlesROOTName + string("#rightarrow") + locFinalParticlesROOTName + string(" Missing Mass (GeV/c^{2})");
		dHist_MissingMassVsBeamE = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DBeamEBins, dMinBeamE, dMaxBeamE, dNum2DMassBins, dMinMass, dMaxMass);

		//Return to the base directory
		ChangeTo_BaseDirectory();
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
	double locBeamEnergy = locParticleCombo->Get_ParticleComboStep(0)->Get_InitialParticle()->energy();
	japp->RootWriteLock();
	{
		dHist_MissingMass->Fill(locMissingMass);
		dHist_MissingMassVsBeamE->Fill(locBeamEnergy, locMissingMass);
	}
	japp->RootUnLock();

	return true;
}

void DHistogramAction_MissingMassSquared::Initialize(JEventLoop* locEventLoop)
{
	double locMassSqPerBin = 1000.0*1000.0*(dMaxMassSq - dMinMassSq)/((double)dNumMassBins);
	string locInitialParticlesROOTName = Get_Reaction()->Get_InitialParticlesROOTName();
	string locFinalParticlesROOTName = Get_Reaction()->Get_DecayChainFinalParticlesROOTNames(0, dMissingMassOffOfStepIndex, dMissingMassOffOfPIDs, Get_UseKinFitResultsFlag(), true);

	locEventLoop->GetSingle(dAnalysisUtilities);

	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		CreateAndChangeTo_ActionDirectory();

		string locHistName = "MissingMassSquared";
		ostringstream locStream;
		locStream << locMassSqPerBin;
		string locHistTitle = string(";") + locInitialParticlesROOTName + string("#rightarrow") + locFinalParticlesROOTName + string(" Missing Mass Squared (GeV/c^{2})^{2};# Combos / ") + locStream.str() + string(" (MeV/c^{2})^{2}");
		dHist_MissingMassSquared = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumMassBins, dMinMassSq, dMaxMassSq);

		locHistName = "MissingMassSquaredVsBeamE";
		locMassSqPerBin *= ((double)dNumMassBins)/((double)dNum2DMassBins);
		locStream.str();
		locStream << locMassSqPerBin;
		locHistTitle = string(";Beam Energy (GeV);") + locInitialParticlesROOTName + string("#rightarrow") + locFinalParticlesROOTName + string(" Missing Mass Squared (GeV/c^{2})^{2};");
		dHist_MissingMassSquaredVsBeamE = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DBeamEBins, dMinBeamE, dMaxBeamE, dNum2DMassBins, dMinMassSq, dMaxMassSq);

		//Return to the base directory
		ChangeTo_BaseDirectory();
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
	double locBeamEnergy = locParticleCombo->Get_ParticleComboStep(0)->Get_InitialParticle()->energy();
	japp->RootWriteLock();
	{
		dHist_MissingMassSquared->Fill(locMissingMassSquared);
		dHist_MissingMassSquaredVsBeamE->Fill(locBeamEnergy, locMissingMassSquared);
	}
	japp->RootUnLock();

	return true;
}

void DHistogramAction_KinFitResults::Initialize(JEventLoop* locEventLoop)
{
	DKinFitType locKinFitType = Get_Reaction()->Get_KinFitType();
	if(locKinFitType == d_NoFit)
		return;

	locEventLoop->GetSingle(dAnalysisUtilities);
	dKinFitUtils = new DKinFitUtils_GlueX(locEventLoop);

	set<pair<int, int> > locVertexParticles = dKinFitUtils->Get_KinFitVertexParticles(Get_Reaction());

	size_t locNumConstraints = 0, locNumUnknowns = 0;
	string locConstraintString = dKinFitUtils->Get_ConstraintInfo(Get_Reaction(), locKinFitType, locNumConstraints, locNumUnknowns);

	size_t locNDF = locNumConstraints - locNumUnknowns;
	bool locIncludeBeamlineInVertexFitFlag = dKinFitUtils->Get_IncludeBeamlineInVertexFitFlag();

	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		CreateAndChangeTo_ActionDirectory();

		// Confidence Level
		string locHistName = "ConfidenceLevel";
		ostringstream locHistTitle;
		locHistTitle << "Kinematic Fit Constraints: " << locConstraintString << ";Confidence Level (" << locNumConstraints;
		locHistTitle << " Constraints, " << locNumUnknowns << " Unknowns: " << locNDF << "-C Fit);# Combos";
		dHist_ConfidenceLevel = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle.str(), dNumConfidenceLevelBins, 0.0, 1.0);

		// Pulls
		map<DKinFitPullType, TH1I*> locParticlePulls;

		//beam pulls
		Particle_t locInitialPID = Get_Reaction()->Get_ReactionStep(0)->Get_InitialParticleID();
		bool locBeamFlag = (locInitialPID == Gamma);
		if(locBeamFlag)
		{
			string locFullROOTName = string("Beam ") + ParticleName_ROOT(locInitialPID);
			CreateAndChangeTo_Directory("Beam", "Beam");

			pair<int, int> locParticlePair(0, -1);
			bool locIsInVertexFitFlag = (locVertexParticles.find(locParticlePair) != locVertexParticles.end());
			if(!locIncludeBeamlineInVertexFitFlag)
				locIsInVertexFitFlag = false;

			Create_ParticlePulls(locFullROOTName, locIsInVertexFitFlag, false, dHistMap_BeamPulls);
			gDirectory->cd("..");
		}

		//final particle pulls
		for(size_t loc_i = 0; loc_i < Get_Reaction()->Get_NumReactionSteps(); ++loc_i)
		{
			const DReactionStep* locReactionStep = Get_Reaction()->Get_ReactionStep(loc_i);
			string locStepName = locReactionStep->Get_StepName();
			string locStepROOTName = locReactionStep->Get_StepROOTName();

			deque<Particle_t> locPIDs;
			locReactionStep->Get_FinalParticleIDs(locPIDs);
			set<Particle_t> locPIDSet;

			for(size_t loc_j = 0; loc_j < locPIDs.size(); ++loc_j)
			{
				if(Get_Reaction()->Get_DecayStepIndex(loc_i, loc_j) >= 0)
					continue; //decaying particle
				if(locReactionStep->Get_MissingParticleIndex() == int(loc_j))
					continue; //missing particle

				Particle_t locPID = locPIDs[loc_j];
				if(locPIDSet.find(locPID) != locPIDSet.end())
					continue; //histograms already created for this pid

				if(locPIDSet.empty()) //first call
					CreateAndChangeTo_Directory(locStepName, locStepName);

				string locParticleName = ParticleType(locPID);
				CreateAndChangeTo_Directory(locParticleName, locParticleName);

				string locParticleROOTName = ParticleName_ROOT(locPID);
				string locFullROOTName = locParticleROOTName + string(", ") + locStepROOTName;

				pair<int, int> locParticlePair(loc_i, loc_j);
				bool locIsInVertexFitFlag = (locVertexParticles.find(locParticlePair) != locVertexParticles.end());

				bool locIsNeutralShowerFlag = (locIsInVertexFitFlag && (ParticleCharge(locPID) == 0));

				Create_ParticlePulls(locFullROOTName, locIsInVertexFitFlag, locIsNeutralShowerFlag, locParticlePulls);
				dHistMap_Pulls[pair<size_t, Particle_t>(loc_i, locPID)] = locParticlePulls;
				locPIDSet.insert(locPID);

				gDirectory->cd("..");
			} //end of particle loop

			if(!locPIDSet.empty())
				gDirectory->cd("..");
		} //end of step loop

		//Return to the base directory
		ChangeTo_BaseDirectory();
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

void DHistogramAction_KinFitResults::Create_ParticlePulls(string locFullROOTName, bool locIsInVertexFitFlag, bool locIsNeutralShowerFlag, map<DKinFitPullType, TH1I*>& locParticlePulls)
{
	locParticlePulls.clear();

	DKinFitType locKinFitType = Get_Reaction()->Get_KinFitType();
	bool locP4IsFit = ((locKinFitType != d_VertexFit) && (locKinFitType != d_SpacetimeFit));

	//p4 pulls:
	string locHistName, locHistTitle;
	if(locP4IsFit && locIsNeutralShowerFlag)
	{
		//neutral shower not in a p4-only fit
		//E Pull
		locHistName = "Pull_E";
		locHistTitle = locFullROOTName + string(";E Pull;# Combos");
		locParticlePulls[d_EPull] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, dMinPull, dMaxPull);
	}
	else if(locP4IsFit || locIsInVertexFitFlag)
	{
		//all detected particles (except neutral showers) have p3 used in p4 fits and vertex fits
		//however, don't include if the particles aren't actually in the fit (vertex-only, and too few particles at that vertex to constrain)

		//Px Pull
		locHistName = "Pull_Px";
		locHistTitle = locFullROOTName + string(";p_{x} Pull;# Combos");
		locParticlePulls[d_PxPull] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, dMinPull, dMaxPull);

		//Py Pull
		locHistName = "Pull_Py";
		locHistTitle = locFullROOTName + string(";p_{y} Pull;# Combos");
		locParticlePulls[d_PyPull] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, dMinPull, dMaxPull);

		//Pz Pull
		locHistName = "Pull_Pz";
		locHistTitle = locFullROOTName + string(";p_{z} Pull;# Combos");
		locParticlePulls[d_PzPull] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, dMinPull, dMaxPull);
	}

	//vertex pulls:
	if((locIsNeutralShowerFlag && locP4IsFit) || (!locIsNeutralShowerFlag && locIsInVertexFitFlag))
	{
		//Xx Pull
		locHistName = "Pull_Xx";
		locHistTitle = locFullROOTName + string(";x_{x} Pull;# Combos");
		locParticlePulls[d_XxPull] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, dMinPull, dMaxPull);

		//Xy Pull
		locHistName = "Pull_Xy";
		locHistTitle = locFullROOTName + string(";x_{y} Pull;# Combos");
		locParticlePulls[d_XyPull] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, dMinPull, dMaxPull);

		//Xz Pull
		locHistName = "Pull_Xz";
		locHistTitle = locFullROOTName + string(";x_{z} Pull;# Combos");
		locParticlePulls[d_XzPull] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, dMinPull, dMaxPull);
	}

	//time pulls:
	if(locIsInVertexFitFlag && ((locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndSpacetimeFit)))
	{
		//T Pull
		locHistName = "Pull_T";
		locHistTitle = locFullROOTName + string(";t Pull;# Combos");
		locParticlePulls[d_TPull] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, dMinPull, dMaxPull);
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
		dHist_ConfidenceLevel->Fill(locConfidenceLevel);
	}
	japp->RootUnLock();
	if(locConfidenceLevel < dPullHistConfidenceLevelCut)
		return true; //don't histogram pulls

	// Pulls
	map<const JObject*, map<DKinFitPullType, double> > locPulls; //DKinematicData is the MEASURED particle
	locKinFitResults->Get_Pulls(locPulls);

	bool locBeamFlag = (Get_Reaction()->Get_ReactionStep(0)->Get_InitialParticleID() == Gamma);
	japp->RootWriteLock();
	{
		// beam pulls
		if(locBeamFlag)
		{
			const DKinematicData* locKinematicData = locParticleCombo->Get_ParticleComboStep(0)->Get_InitialParticle_Measured();
			map<DKinFitPullType, double> locParticlePulls = locPulls[locKinematicData];
			map<DKinFitPullType, double>::iterator locIterator = locParticlePulls.begin();
			for(; locIterator != locParticlePulls.end(); ++locIterator)
				dHistMap_BeamPulls[locIterator->first]->Fill(locIterator->second);
		}

		// final particle pulls
		for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
		{
			const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);

			deque<const DKinematicData*> locParticles;
			locParticleComboStep->Get_FinalParticles_Measured(locParticles);
			for(size_t loc_j = 0; loc_j < locParticles.size(); ++loc_j)
			{
				if(!locParticleComboStep->Is_FinalParticleDetected(loc_j))
					continue;

				//get pulls for this particle
				map<DKinFitPullType, double> locParticlePulls;
				map<const JObject*, map<DKinFitPullType, double> >::iterator locParticleIterator = locPulls.find(locParticles[loc_j]);
				if(locParticleIterator != locPulls.end())
					locParticlePulls = locParticleIterator->second;
				else //is neutral shower
					locParticlePulls = locPulls[locParticleComboStep->Get_FinalParticle_SourceObject(loc_j)];

				//fill histograms
				pair<size_t, Particle_t> locParticlePair(loc_i, locParticles[loc_j]->PID());
				map<DKinFitPullType, double>::iterator locIterator = locParticlePulls.begin();
				for(; locIterator != locParticlePulls.end(); ++locIterator)
					(dHistMap_Pulls[locParticlePair])[locIterator->first]->Fill(locIterator->second);
			}
		}
	}
	japp->RootUnLock();

	return true;
}

void DHistogramAction_MissingTransverseMomentum::Initialize(JEventLoop* locEventLoop)
{
	string locHistName, locHistTitle;
	double locPtPerBin = 1000.0*(dMaxPt - dMinPt)/((double)dNumPtBins);

	vector<const DAnalysisUtilities*> locAnalysisUtilitiesVector;
	locEventLoop->Get(locAnalysisUtilitiesVector);

	//CREATE THE HISTOGRAMS
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		dAnalysisUtilities = locAnalysisUtilitiesVector[0];
		CreateAndChangeTo_ActionDirectory();

		locHistName = "MissingTransverseMomentum";
		ostringstream locStream;
		locStream << locPtPerBin;
		locHistTitle = string(";") + string(" Missing Transverse Momentum (GeV/c);# Combos / ") + locStream.str() + string(" MeV/c");
		dHist_MissingTransverseMomentum = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPtBins, dMinPt, dMaxPt);

		//Return to the base directory
		ChangeTo_BaseDirectory();
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_MissingTransverseMomentum::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	if(Get_NumPreviousParticleCombos() == 0)
		dPreviousSourceObjects.clear();

	set<pair<const JObject*, Particle_t> > locSourceObjects;
	DLorentzVector locFinalStateP4 = dAnalysisUtilities->Calc_FinalStateP4(locParticleCombo, 0, locSourceObjects, Get_UseKinFitResultsFlag()); // Use step '0'

	if(!dEnableDoubleCounting)
	{
		if(dPreviousSourceObjects.find(locSourceObjects) != dPreviousSourceObjects.end())
			return true; //dupe: already histed!
		dPreviousSourceObjects.insert(locSourceObjects);
	}

	double locMissingTransverseMomentum = locFinalStateP4.Pt();

	japp->RootWriteLock();
	{
		dHist_MissingTransverseMomentum->Fill(locMissingTransverseMomentum);
	}
	japp->RootUnLock();

	return true;
}
