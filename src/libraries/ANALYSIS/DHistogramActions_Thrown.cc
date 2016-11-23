#include "ANALYSIS/DHistogramActions.h"

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
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		CreateAndChangeTo_ActionDirectory();

		locGeometry->GetTargetZ(dTargetZCenter);

		//RF
		locHistName = "DeltaT_RFBeamBunch";
		dRFBeamBunchDeltaT_Hist = GetOrCreate_Histogram<TH1I>(locHistName, ";RF #Deltat (Reconstructed - Thrown)", dNumRFDeltaTBins, dMinRFDeltaT, dMaxRFDeltaT);

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
			dBeamParticleHist_DeltaPOverP = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

			// DeltaP/P Vs P
			locHistName = string("DeltaPOverPVsP");
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltap/p (Reconstructed - Thrown)");
			dBeamParticleHist_DeltaPOverPVsP = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

			// DeltaT
			locHistName = string("DeltaT");
			locHistTitle = locParticleROOTName + string(";#Deltat (ns) (Reconstructed - Thrown)");
			dBeamParticleHist_DeltaT = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

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
				dHistDeque_DeltaPOverP[loc_i][locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

				// DeltaTheta
				locHistName = string("DeltaTheta");
				locHistTitle = locParticleROOTName + string(";#Delta#theta#circ (Reconstructed - Thrown)");
				dHistDeque_DeltaTheta[loc_i][locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaThetaBins, dMinDeltaTheta, dMaxDeltaTheta);

				// DeltaPhi
				locHistName = string("DeltaPhi");
				locHistTitle = locParticleROOTName + string(";#Delta#phi#circ (Reconstructed - Thrown)");
				dHistDeque_DeltaPhi[loc_i][locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

				// DeltaT
				locHistName = string("DeltaT");
				locHistTitle = locParticleROOTName + string(";#Deltat (ns) (Reconstructed - Thrown)");
				dHistDeque_DeltaT[loc_i][locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

				// DeltaT - BCAL
				locHistName = string("DeltaT_BCAL");
				locHistTitle = locParticleROOTName + string(" in BCAL;#Deltat (ns) (Reconstructed - Thrown)");
				dHistDeque_DeltaT_BCAL[loc_i][locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

				// DeltaT - TOF (charged only)
				if(ParticleCharge(locPID) != 0)
				{
					locHistName = string("DeltaT_TOF");
					locHistTitle = locParticleROOTName + string(" in TOF;#Deltat (ns) (Reconstructed - Thrown)");
					dHistDeque_DeltaT_TOF[loc_i][locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
				}

				// DeltaT - FCAL (neutral only)
				if(ParticleCharge(locPID) == 0)
				{
					locHistName = string("DeltaT_FCAL");
					locHistTitle = locParticleROOTName + string(" in FCAL;#Deltat (ns) (Reconstructed - Thrown)");
					dHistDeque_DeltaT_FCAL[loc_i][locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
				}

				// DeltaVertexZ
				locHistName = string("DeltaVertexZ");
				locHistTitle = locParticleROOTName + string(";#DeltaVertex-Z (cm) (Reconstructed - Thrown)");
				dHistDeque_DeltaVertexZ[loc_i][locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaVertexZBins, dMinDeltaVertexZ, dMaxDeltaVertexZ);

				// DeltaP/P Vs P
				locHistName = string("DeltaPOverPVsP");
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltap/p (Reconstructed - Thrown)");
				dHistDeque_DeltaPOverPVsP[loc_i][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

				// DeltaP/P Vs Theta
				locHistName = string("DeltaPOverPVsTheta");
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltap/p (Reconstructed - Thrown)");
				dHistDeque_DeltaPOverPVsTheta[loc_i][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

				// DeltaTheta Vs P
				locHistName = string("DeltaThetaVsP");
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta#theta#circ (Reconstructed - Thrown)");
				dHistDeque_DeltaThetaVsP[loc_i][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNumDeltaThetaBins, dMinDeltaTheta, dMaxDeltaTheta);

				// DeltaTheta Vs Theta
				locHistName = string("DeltaThetaVsTheta");
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Delta#theta#circ (Reconstructed - Thrown)");
				dHistDeque_DeltaThetaVsTheta[loc_i][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaThetaBins, dMinDeltaTheta, dMaxDeltaTheta);

				// DeltaPhi Vs P
				locHistName = string("DeltaPhiVsP");
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta#phi#circ (Reconstructed - Thrown)");
				dHistDeque_DeltaPhiVsP[loc_i][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNumDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

				// DeltaPhi Vs Theta
				locHistName = string("DeltaPhiVsTheta");
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Delta#phi#circ (Reconstructed - Thrown)");
				dHistDeque_DeltaPhiVsTheta[loc_i][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

				// DeltaT Vs Theta
				locHistName = string("DeltaTVsTheta");
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat (ns) (Reconstructed - Thrown)");
				dHistDeque_DeltaTVsTheta[loc_i][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

				// DeltaT Vs P
				locHistName = string("DeltaTVsP");
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat (ns) (Reconstructed - Thrown)");
				dHistDeque_DeltaTVsP[loc_i][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

				// DeltaVertexZ Vs Theta
				locHistName = string("DeltaVertexZVsTheta");
				locHistTitle = locParticleROOTName + string(";#theta#circ;#DeltaVertex-Z (cm) (Reconstructed - Thrown)");
				dHistDeque_DeltaVertexZVsTheta[loc_i][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaVertexZBins, dMinDeltaVertexZ, dMaxDeltaVertexZ);

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
					dHistDeque_Pulls[loc_i][locPID][dPullTypes[loc_j]] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, -10.0, 10.0);

					//Pull vs P
					locHistName = locPullNames[loc_j] + string("PullVsP");
					locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta") + locPullTitles[loc_j] + string("/#sigma_{") + locPullTitles[loc_j] + string("} (Reconstructed - Thrown)");
					dHistDeque_PullsVsP[loc_i][locPID][dPullTypes[loc_j]] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

					//Pull vs Theta
					locHistName = locPullNames[loc_j] + string("PullVsTheta");
					locHistTitle = locParticleROOTName + string(";#theta#circ;#Delta") + locPullTitles[loc_j] + string("/#sigma_{") + locPullTitles[loc_j] + string("} (Reconstructed - Thrown)");
					dHistDeque_PullsVsTheta[loc_i][locPID][dPullTypes[loc_j]] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);
				}

				//Delta-t Pulls - CDC & ST
				if(ParticleCharge(locPID) != 0)
				{
					//CDC
					locHistName = "TimePull_CDC";
					locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
					dHistDeque_TimePull_CDC[loc_i][locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, -10.0, 10.0);

					locHistName = "TimePullVsTheta_CDC";
					locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat/#sigma_{#Deltat}");
					dHistDeque_TimePullVsTheta_CDC[loc_i][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);

					locHistName = "TimePullVsP_CDC";
					locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
					dHistDeque_TimePullVsP_CDC[loc_i][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

					//ST
					locHistName = "TimePull_ST";
					locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
					dHistDeque_TimePull_ST[loc_i][locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, -10.0, 10.0);

					locHistName = "TimePullVsTheta_ST";
					locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat/#sigma_{#Deltat}");
					dHistDeque_TimePullVsTheta_ST[loc_i][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);

					locHistName = "TimePullVsP_ST";
					locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
					dHistDeque_TimePullVsP_ST[loc_i][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);
				}

				//Delta-t Pulls - BCAL
				locHistName = "TimePull_BCAL";
				locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
				dHistDeque_TimePull_BCAL[loc_i][locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, -10.0, 10.0);

				locHistName = "TimePullVsTheta_BCAL";
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat/#sigma_{#Deltat}");
				dHistDeque_TimePullVsTheta_BCAL[loc_i][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);

				locHistName = "TimePullVsP_BCAL";
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
				dHistDeque_TimePullVsP_BCAL[loc_i][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

				//Delta-t Pulls - TOF
				if(ParticleCharge(locPID) != 0) //TOF
				{
					locHistName = "TimePull_TOF";
					locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
					dHistDeque_TimePull_TOF[loc_i][locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, -10.0, 10.0);

					locHistName = "TimePullVsP_TOF";
					locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
					dHistDeque_TimePullVsP_TOF[loc_i][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);
				}

				//Delta-t Pulls - FCAL
				locHistName = "TimePull_FCAL";
				locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
				dHistDeque_TimePull_FCAL[loc_i][locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, -10.0, 10.0);

				locHistName = "TimePullVsP_FCAL";
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
				dHistDeque_TimePullVsP_FCAL[loc_i][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

				gDirectory->cd("..");

				gDirectory->cd("..");
			} //end of particle loop
			gDirectory->cd("..");
		} //end of step loop

		//Return to the base directory
		ChangeTo_BaseDirectory();
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
	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action();
	{
		dRFBeamBunchDeltaT_Hist->Fill(locRFDeltaT);
	}
	Unlock_Action();

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

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action();
	{
		dBeamParticleHist_DeltaPOverP->Fill(locDeltaPOverP);
		dBeamParticleHist_DeltaPOverPVsP->Fill(locThrownP, locDeltaPOverP);
		dBeamParticleHist_DeltaT->Fill(locDeltaT);
	}
	Unlock_Action();
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

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action();
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
	Unlock_Action();
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
	DVector3 locNeutralP3 = locNeutralParticleHypothesis->momentum();
	double locPMag = locNeutralP3.Mag();
	double locDeltaPOverP = (locPMag - locThrownP)/locThrownP;
	double locDeltaTheta = locNeutralP3.Theta()*180.0/TMath::Pi() - locThrownTheta;
	double locDeltaPhi = locNeutralP3.Phi()*180.0/TMath::Pi() - locMCThrown->momentum().Phi()*180.0/TMath::Pi();
	double locDeltaT = locNeutralParticleHypothesis->time() - locMCThrown->time(); //time comparison isn't fair if track comes from a detached vertex!!!
	double locDeltaVertexZ = locNeutralParticleHypothesis->position().Z() - locMCThrown->position().Z();
	const TMatrixDSym& locCovarianceMatrix = locNeutralParticleHypothesis->errorMatrix();

	double locStartTime = locThrownEventRFBunch->dTime + (locMCThrown->z() - dTargetZCenter)/29.9792458;
	double locTimePull = (locStartTime - locNeutralParticleHypothesis->time())/sqrt(locCovarianceMatrix(6, 6));

	double locEPull = 0.0;
	if(Get_UseKinFitResultsFlag())
	{
		DMatrix locJacobian(1, 7);
		locJacobian(0, 0) = locNeutralP3.Px()/locPMag;
		locJacobian(0, 1) = locNeutralP3.Py()/locPMag;
		locJacobian(0, 2) = locNeutralP3.Pz()/locPMag;
		for(unsigned int loc_i = 0; loc_i < 4; ++loc_i)
			locJacobian(0, 3 + loc_i) = 0.0;

		TMatrixDSym locCovCopy = locCovarianceMatrix;
		locCovCopy.Similarity(locJacobian);
		double locEUncertainty = sqrt(locCovCopy(0, 0));
		locEPull = (locNeutralParticleHypothesis->energy() - locMCThrown->energy())/locEUncertainty;
	}
	else
		locEPull = (locNeutralShower->dEnergy - locMCThrown->energy())/sqrt(locNeutralShower->dCovarianceMatrix(0, 0));

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action();
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
				locPull = locEPull;
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
	Unlock_Action();
}

void DHistogramAction_ThrownParticleKinematics::Initialize(JEventLoop* locEventLoop)
{
	string locHistName, locHistTitle, locParticleName, locParticleROOTName;
	Particle_t locPID;

	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
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
			dMCGENBeamParticle_P = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPBins, dMinP, dMaxP);

			locHistName = "Time";
			locHistTitle = string("MCGEN Thrown Beam ") + locParticleROOTName + string(";t (ns)");
			dMCGENBeamParticle_Time = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumTBins, dMinT, dMaxT);

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
			dAllBeamParticle_P = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPBins, dMinP, dMaxP);

			locHistName = "Time";
			locHistTitle = string("TRUTH Thrown Beam ") + locParticleROOTName + string(";t (ns)");
			dAllBeamParticle_Time = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumTBins, dMinT, dMaxT);

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
			dHistMap_P[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPBins, dMinP, dMaxP);

			// Theta
			locHistName = "Theta";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";#theta#circ");
			dHistMap_Theta[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumThetaBins, dMinTheta, dMaxTheta);

			// Phi
			locHistName = "Phi";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";#phi#circ");
			dHistMap_Phi[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPhiBins, dMinPhi, dMaxPhi);

			// P Vs Theta
			locHistName = "PVsTheta";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";#theta#circ;p (GeV/c)");
			dHistMap_PVsTheta[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			// Phi Vs Theta
			locHistName = "PhiVsTheta";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";#theta#circ;#phi#circ");
			dHistMap_PhiVsTheta[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

			// Vertex-Z
			locHistName = "VertexZ";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";Vertex-Z (cm)");
			dHistMap_VertexZ[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumVertexZBins, dMinVertexZ, dMaxVertexZ);

			// Vertex-Y Vs Vertex-X
			locHistName = "VertexYVsX";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";Vertex-X (cm);Vertex-Y (cm)");
			dHistMap_VertexYVsX[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY);

			// Vertex-T
			locHistName = "VertexT";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";Vertex-T (ns)");
			dHistMap_VertexT[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumTBins, dMinT, dMaxT);

			gDirectory->cd("..");
		}

		//Return to the base directory
		ChangeTo_BaseDirectory();
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

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action();
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
	Unlock_Action();

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

		//FILL HISTOGRAMS
		//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
		//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
		Lock_Action();
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
		Unlock_Action();
	}
	return true;
}

void DHistogramAction_ReconnedThrownKinematics::Initialize(JEventLoop* locEventLoop)
{
	string locHistName, locHistTitle, locParticleName, locParticleROOTName;
	Particle_t locPID;

	const DAnalysisUtilities* locAnalysisUtilities = NULL;
	locEventLoop->GetSingle(locAnalysisUtilities);

	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		dAnalysisUtilities = locAnalysisUtilities;

		CreateAndChangeTo_ActionDirectory();

		// Beam Particle
		locPID = Gamma;
		locParticleName = string("Beam_") + ParticleType(locPID);
		locParticleROOTName = ParticleName_ROOT(locPID);
		CreateAndChangeTo_Directory(locParticleName, locParticleName);

		locHistName = "Momentum";
		locHistTitle = string("Thrown Beam ") + locParticleROOTName + string(";p (GeV/c)");
		dBeamParticle_P = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPBins, dMinP, dMaxP);
		gDirectory->cd("..");

		//PID
		CreateAndChangeTo_Directory("PID", "PID");
		{
			//beta vs p
			locHistName = "BetaVsP_Q+";
			locHistTitle = "q^{+};p (GeV/c);#beta";
			dHistMap_QBetaVsP[1] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNumBetaBins, dMinBeta, dMaxBeta);

			locHistName = "BetaVsP_Q-";
			locHistTitle = "q^{-};p (GeV/c);#beta";
			dHistMap_QBetaVsP[-1] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNumBetaBins, dMinBeta, dMaxBeta);
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
			dHistMap_P[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPBins, dMinP, dMaxP);

			// Theta
			locHistName = "Theta";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";#theta#circ");
			dHistMap_Theta[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumThetaBins, dMinTheta, dMaxTheta);

			// Phi
			locHistName = "Phi";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";#phi#circ");
			dHistMap_Phi[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPhiBins, dMinPhi, dMaxPhi);

			// P Vs Theta
			locHistName = "PVsTheta";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";#theta#circ;p (GeV/c)");
			dHistMap_PVsTheta[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

			// Phi Vs Theta
			locHistName = "PhiVsTheta";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";#theta#circ;#phi#circ");
			dHistMap_PhiVsTheta[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

			// Vertex-Z
			locHistName = "VertexZ";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";Vertex-Z (cm)");
			dHistMap_VertexZ[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumVertexZBins, dMinVertexZ, dMaxVertexZ);

			// Vertex-Y Vs Vertex-X
			locHistName = "VertexYVsX";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";Vertex-X (cm);Vertex-Y (cm)");
			dHistMap_VertexYVsX[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY, dNumVertexXYBins, dMinVertexXY, dMaxVertexXY);

			// Vertex-T
			locHistName = "VertexT";
			locHistTitle = string("Thrown ") + locParticleROOTName + string(";Vertex-T (ns)");
			dHistMap_VertexT[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumTBins, dMinT, dMaxT);

			gDirectory->cd("..");
		}

		//Return to the base directory
		ChangeTo_BaseDirectory();
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

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action();
	{
		for(size_t loc_i = 0; loc_i < locBeamPhotons.size(); ++loc_i)
			dBeamParticle_P->Fill(locBeamPhotons[loc_i]->energy());
	}
	Unlock_Action();

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
			locBeta_Timing = locChargedTrackHypothesis->measuredBeta();
		}
		else
		{
			const DNeutralParticleHypothesis* locNeutralParticleHypothesis = locMCThrownMatching->Get_MatchingNeutralHypothesis(locMCThrown, locMatchFOM);
			if(locNeutralParticleHypothesis == NULL)
				continue; //not reconstructed
			locBeta_Timing = locNeutralParticleHypothesis->measuredBeta();
		}

		DVector3 locMomentum = locMCThrown->momentum();
		double locPhi = locMomentum.Phi()*180.0/TMath::Pi();
		double locTheta = locMomentum.Theta()*180.0/TMath::Pi();
		double locP = locMomentum.Mag();
		int locCharge = ParticleCharge(locPID);

		//FILL HISTOGRAMS
		//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
		//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
		Lock_Action();
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
		Unlock_Action();
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
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		CreateAndChangeTo_ActionDirectory();

		locGeometry->GetTargetZ(dTargetZCenter);

		locHistName = "DeltaT_RFBeamBunch";
		dRFBeamBunchDeltaT_Hist = GetOrCreate_Histogram<TH1I>(locHistName, ";RF #Deltat (Reconstructed - Thrown)", dNumRFDeltaTBins, dMinRFDeltaT, dMaxRFDeltaT);

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
			dHistMap_MatchFOM[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumMCMatchingFOMBins, 0.0, 1.0);

			// DeltaP/P
			locHistName = string("DeltaPOverP");
			locHistTitle = locParticleROOTName + string(";#Deltap/p (Reconstructed - Thrown)");
			dHistMap_DeltaPOverP[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

			// DeltaTheta
			locHistName = string("DeltaTheta");
			locHistTitle = locParticleROOTName + string(";#Delta#theta#circ (Reconstructed - Thrown)");
			dHistMap_DeltaTheta[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaThetaBins, dMinDeltaTheta, dMaxDeltaTheta);

			// DeltaPhi
			locHistName = string("DeltaPhi");
			locHistTitle = locParticleROOTName + string(";#Delta#phi#circ (Reconstructed - Thrown)");
			dHistMap_DeltaPhi[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

			// DeltaT
			locHistName = string("DeltaT");
			locHistTitle = locParticleROOTName + string(";#Deltat (ns) (Reconstructed - Thrown)");
			dHistMap_DeltaT[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

			// DeltaT - BCAL
			locHistName = string("DeltaT_BCAL");
			locHistTitle = locParticleROOTName + string(" in BCAL;#Deltat (ns) (Reconstructed - Thrown)");
			dHistMap_DeltaT_BCAL[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

			// DeltaT - TOF (charged only)
			if(ParticleCharge(locPID) != 0)
			{
				locHistName = string("DeltaT_TOF");
				locHistTitle = locParticleROOTName + string(" in TOF;#Deltat (ns) (Reconstructed - Thrown)");
				dHistMap_DeltaT_TOF[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
			}

			// DeltaT - FCAL (neutral only)
			if(ParticleCharge(locPID) == 0)
			{
				locHistName = string("DeltaT_FCAL");
				locHistTitle = locParticleROOTName + string(" in FCAL;#Deltat (ns) (Reconstructed - Thrown)");
				dHistMap_DeltaT_FCAL[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);
			}

			// DeltaVertexZ
			locHistName = string("DeltaVertexZ");
			locHistTitle = locParticleROOTName + string(";#DeltaVertex-Z (cm) (Reconstructed - Thrown)");
			dHistMap_DeltaVertexZ[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaVertexZBins, dMinDeltaVertexZ, dMaxDeltaVertexZ);

			// DeltaP/P Vs P
			locHistName = string("DeltaPOverPVsP");
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltap/p (Reconstructed - Thrown)");
			dHistMap_DeltaPOverPVsP[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

			// DeltaP/P Vs Theta
			locHistName = string("DeltaPOverPVsTheta");
			locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltap/p (Reconstructed - Thrown)");
			dHistMap_DeltaPOverPVsTheta[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

			// DeltaTheta Vs P
			locHistName = string("DeltaThetaVsP");
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta#theta#circ (Reconstructed - Thrown)");
			dHistMap_DeltaThetaVsP[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNumDeltaThetaBins, dMinDeltaTheta, dMaxDeltaTheta);

			// DeltaTheta Vs Theta
			locHistName = string("DeltaThetaVsTheta");
			locHistTitle = locParticleROOTName + string(";#theta#circ;#Delta#theta#circ (Reconstructed - Thrown)");
			dHistMap_DeltaThetaVsTheta[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaThetaBins, dMinDeltaTheta, dMaxDeltaTheta);

			// DeltaPhi Vs P
			locHistName = string("DeltaPhiVsP");
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta#phi#circ (Reconstructed - Thrown)");
			dHistMap_DeltaPhiVsP[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNumDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

			// DeltaPhi Vs Theta
			locHistName = string("DeltaPhiVsTheta");
			locHistTitle = locParticleROOTName + string(";#theta#circ;#Delta#phi#circ (Reconstructed - Thrown)");
			dHistMap_DeltaPhiVsTheta[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

			// DeltaT Vs Theta
			locHistName = string("DeltaTVsTheta");
			locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat (ns) (Reconstructed - Thrown)");
			dHistMap_DeltaTVsTheta[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

			// DeltaT Vs P
			locHistName = string("DeltaTVsP");
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat (ns) (Reconstructed - Thrown)");
			dHistMap_DeltaTVsP[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

			// DeltaVertexZ Vs Theta
			locHistName = string("DeltaVertexZVsTheta");
			locHistTitle = locParticleROOTName + string(";#theta#circ;#DeltaVertex-Z (cm) (Reconstructed - Thrown)");
			dHistMap_DeltaVertexZVsTheta[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaVertexZBins, dMinDeltaVertexZ, dMaxDeltaVertexZ);

			// P Vs Theta
			locHistName = "PVsTheta_LargeDeltaT";
			locHistTitle = locParticleROOTName + string(";#theta#circ;p (GeV/c)");
			dHistMap_PVsTheta_LargeDeltaT[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

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
				dHistMap_Pulls[locPID][dPullTypes[loc_j]] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, -10.0, 10.0);

				//Pull vs P
				locHistName = locPullNames[loc_j] + string("PullVsP");
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta") + locPullTitles[loc_j] + string("/#sigma_{") + locPullTitles[loc_j] + string("} (Reconstructed - Thrown)");
				dHistMap_PullsVsP[locPID][dPullTypes[loc_j]] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

				//Pull vs Theta
				locHistName = locPullNames[loc_j] + string("PullVsTheta");
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Delta") + locPullTitles[loc_j] + string("/#sigma_{") + locPullTitles[loc_j] + string("} (Reconstructed - Thrown)");
				dHistMap_PullsVsTheta[locPID][dPullTypes[loc_j]] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);
			}

			//Delta-t Pulls - CDC & ST
			if(ParticleCharge(locPID) != 0)
			{
				//CDC
				locHistName = "TimePull_CDC";
				locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
				dHistMap_TimePull_CDC[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, -10.0, 10.0);

				locHistName = "TimePullVsTheta_CDC";
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat/#sigma_{#Deltat}");
				dHistMap_TimePullVsTheta_CDC[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);

				locHistName = "TimePullVsP_CDC";
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
				dHistMap_TimePullVsP_CDC[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

				//ST
				locHistName = "TimePull_ST";
				locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
				dHistMap_TimePull_ST[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, -10.0, 10.0);

				locHistName = "TimePullVsTheta_ST";
				locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat/#sigma_{#Deltat}");
				dHistMap_TimePullVsTheta_ST[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);

				locHistName = "TimePullVsP_ST";
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
				dHistMap_TimePullVsP_ST[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);
			}

			//Delta-t Pulls - BCAL
			locHistName = "TimePull_BCAL";
			locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
			dHistMap_TimePull_BCAL[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, -10.0, 10.0);

			locHistName = "TimePullVsTheta_BCAL";
			locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltat/#sigma_{#Deltat}");
			dHistMap_TimePullVsTheta_BCAL[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPullBins, -10.0, 10.0);

			locHistName = "TimePullVsP_BCAL";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
			dHistMap_TimePullVsP_BCAL[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

			//Delta-t Pulls - TOF
			if(ParticleCharge(locPID) != 0) //TOF
			{
				locHistName = "TimePull_TOF";
				locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
				dHistMap_TimePull_TOF[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, -10.0, 10.0);

				locHistName = "TimePullVsP_TOF";
				locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
				dHistMap_TimePullVsP_TOF[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);
			}

			//Delta-t Pulls - FCAL
			locHistName = "TimePull_FCAL";
			locHistTitle = locParticleROOTName + string(";#Deltat/#sigma_{#Deltat}");
			dHistMap_TimePull_FCAL[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPullBins, -10.0, 10.0);

			locHistName = "TimePullVsP_FCAL";
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat/#sigma_{#Deltat}");
			dHistMap_TimePullVsP_FCAL[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DPullBins, -10.0, 10.0);

			gDirectory->cd("..");

			gDirectory->cd("..");
		}

		//Return to the base directory
		ChangeTo_BaseDirectory();
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

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action();
	{
		dRFBeamBunchDeltaT_Hist->Fill(locRFDeltaT);
	}
	Unlock_Action();

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

		//FILL HISTOGRAMS
		//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
		//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
		Lock_Action();
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
		Unlock_Action();
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

		//FILL HISTOGRAMS
		//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
		//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
		Lock_Action();
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
		Unlock_Action();
	}
	return true;
}

void DHistogramAction_TOFHitStudy::Initialize(JEventLoop* locEventLoop)
{
	string locHistName, locHistTitle, locParticleName, locParticleROOTName;
	Particle_t locPID;

	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
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
			dHistMap_DeltaT[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

			// DeltaX
			locHistName = string("DeltaX_") + locParticleName;
			locHistTitle = locParticleROOTName + string(";#Deltax (cm) (Reconstructed - Thrown)");
			dHistMap_DeltaX[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaXBins, dMinDeltaX, dMaxDeltaX);

			// DeltaY
			locHistName = string("DeltaY_") + locParticleName;
			locHistTitle = locParticleROOTName + string(";#Deltay (cm) (Reconstructed - Thrown)");
			dHistMap_DeltaY[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumDeltaXBins, dMinDeltaX, dMaxDeltaX);

			// dE
			locHistName = string("dE_") + locParticleName;
			locHistTitle = locParticleROOTName + string(";dE (MeV)");
			dHistMap_dE[locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumdEBins, dMindE, dMaxdE);

			// DeltaT Vs P
			locHistName = string("DeltaTVsP_") + locParticleName;
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltat (ns) (Reconstructed - Thrown)");
			dHistMap_DeltaTVsP[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNumDeltaTBins, dMinDeltaT, dMaxDeltaT);

			// DeltaX Vs P
			locHistName = string("DeltaXVsP_") + locParticleName;
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltax (cm) (Reconstructed - Thrown)");
			dHistMap_DeltaXVsP[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNumDeltaXBins, dMinDeltaX, dMaxDeltaX);

			// DeltaY Vs P
			locHistName = string("DeltaYVsP_") + locParticleName;
			locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltay (cm) (Reconstructed - Thrown)");
			dHistMap_DeltaYVsP[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNumDeltaXBins, dMinDeltaX, dMaxDeltaX);

			// dE Vs P
			locHistName = string("dEVsP_") + locParticleName;
			locHistTitle = locParticleROOTName + string(";p (GeV/c);dE (GeV)");
			dHistMap_dEVsP[locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNumdEBins, dMindE, dMaxdE);

			gDirectory->cd("..");
		}

		//Return to the base directory
		ChangeTo_BaseDirectory();
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

		//FILL HISTOGRAMS
		//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
		//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
		Lock_Action(); //ACQUIRE ROOT LOCK!!
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
		Unlock_Action(); //RELEASE ROOT LOCK!!
	}

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
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
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
				dHistDeque_P_CorrectID[loc_i][locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPBins, dMinP, dMaxP);

				//P of Incorrect ID
				locHistName = string("Momentum_IncorrectID_") + locParticleName;
				locHistTitle = string("Incorrect ") + locParticleROOTName + string(" ID, ") + locStepROOTName + string(";p (GeV/c)");
				dHistDeque_P_IncorrectID[loc_i][locPID] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumPBins, dMinP, dMaxP);

				//P Vs Theta of Correct ID
				locHistName = string("PVsTheta_CorrectID_") + locParticleName;
				locHistTitle = string("Correct ") + locParticleROOTName + string(" ID, ") + locStepROOTName + string(";#theta#circ;p (GeV/c)");
				dHistDeque_PVsTheta_CorrectID[loc_i][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNumThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

				//P Vs Theta of Incorrect ID
				locHistName = string("PVsTheta_IncorrectID_") + locParticleName;
				locHistTitle = string("Incorrect ") + locParticleROOTName + string(" ID, ") + locStepROOTName + string(";#theta#circ;p (GeV/c)");
				dHistDeque_PVsTheta_IncorrectID[loc_i][locPID] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNumThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);
			}
			gDirectory->cd("..");
		} //end of step loop

		//# Combos Pass/Fail All True PID
		locHistName = "Combo_TruePIDStatus";
		locHistTitle = Get_Reaction()->Get_ReactionName() + string(";# Combos;All Combo Particles True PID Status");
		dHist_TruePIDStatus = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, 2, -0.5, 1.5);

		//Return to the base directory
		ChangeTo_BaseDirectory();
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

			bool locCutResult = (locMCThrown == NULL) ? false : (((Particle_t)locMCThrown->type) == locPID);
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

			//FILL HISTOGRAMS
			//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
			//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
			Lock_Action();
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
			Unlock_Action();
		}
	}

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action();
	{
		dHist_TruePIDStatus->Fill(locComboTruePIDStatus);
	}
	Unlock_Action();

	return true;
}

