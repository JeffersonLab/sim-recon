#include "ANALYSIS/DHistogramActions.h"

void DHistogramAction_PID::Initialize(JEventLoop* locEventLoop)
{
	string locHistName, locHistTitle;
	string locParticleName, locParticleName2, locParticleROOTName, locParticleROOTName2;
	Particle_t locPID, locPID2;

	//setup pid deques: searched for and generated
	deque<Particle_t> locDesiredPIDs, locThrownPIDs;
	Get_Reaction()->Get_DetectedFinalPIDs(locDesiredPIDs);
	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);
	if(!locMCThrowns.empty()) //else real data, not MC!!
	{
		DMCThrownMatching_factory* locMCThrownMatchingFactory = static_cast<DMCThrownMatching_factory*>(locEventLoop->GetFactory("DMCThrownMatching"));
		locMCThrownMatchingFactory->Get_MCThrownComparisonPIDs(locThrownPIDs);
		locThrownPIDs.push_back(Unknown); //unmatched tracks
	}

	//CREATE THE HISTOGRAMS
	Get_Application()->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	CreateAndChangeTo_ActionDirectory();
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

		if(ParticleCharge(locPID) != 0) //no other sources of PID for neutrals
		{
			// TOF Confidence Level
			locHistName = "TOFConfidenceLevel";
			locHistTitle = locParticleROOTName + string(" PID;TOF Confidence Level");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_TOFFOM[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_TOFFOM[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumFOMBins, 0.0, 1.0);

			// DC dE/dx Confidence Level
			locHistName = "DCdEdxConfidenceLevel";
			locHistTitle = locParticleROOTName + string(" PID;DC #it{#frac{dE}{dx}} Confidence Level");
			if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
				dHistMap_DCdEdxFOM[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
			else
				dHistMap_DCdEdxFOM[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumFOMBins, 0.0, 1.0);
		}

		//one per thrown pid:
		for(size_t loc_j = 0; loc_j < locThrownPIDs.size(); ++loc_j)
		{
			locPID2 = locThrownPIDs[loc_j];
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
		}

		gDirectory->cd("..");
	} //end of PID loop

	Get_Application()->RootUnLock(); //RELEASE ROOT LOCK!!
} //end of Initialize() function

bool DHistogramAction_PID::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, const deque<pair<const DParticleCombo*, bool> >& locPreviousParticleCombos)
{
	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
	locEventLoop->Get(locMCThrownMatchingVector);
	const DMCThrownMatching* locMCThrownMatching = locMCThrownMatchingVector[0];

	deque<const DKinematicData*> locParticles;
	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);

		//charged tracks
		locParticleComboStep->Get_DetectedFinalChargedParticles_Measured(locParticles);
		for(size_t loc_j = 0; loc_j < locParticles.size(); ++loc_j)
		{
			if(!Get_AnalysisUtilities()->Find_SimilarCombos(locParticles[loc_j], locPreviousParticleCombos)) //else duplicate
				Fill_ChargedHists(static_cast<const DChargedTrackHypothesis*>(locParticles[loc_j]), locMCThrownMatching);
		}

		//neutral particles
		locParticleComboStep->Get_DetectedFinalNeutralParticles_Measured(locParticles);
		for(size_t loc_j = 0; loc_j < locParticles.size(); ++loc_j)
		{
			if(!Get_AnalysisUtilities()->Find_SimilarCombos(locParticles[loc_j], locPreviousParticleCombos)) //else duplicate
				Fill_NeutralHists(static_cast<const DNeutralParticleHypothesis*>(locParticles[loc_j]), locMCThrownMatching);
		}
	}

	return true;
}

void DHistogramAction_PID::Fill_ChargedHists(const DChargedTrackHypothesis* locChargedTrackHypothesis, const DMCThrownMatching* locMCThrownMatching)
{
	Particle_t locPID = locChargedTrackHypothesis->PID();
	double locDeltaBeta, locBeta_Timing, locStartTime = 0.0; //should be from RF!!

	locBeta_Timing = locChargedTrackHypothesis->pathLength()/(SPEED_OF_LIGHT*(locChargedTrackHypothesis->t1() - locStartTime));
	locDeltaBeta = locChargedTrackHypothesis->lorentzMomentum().Beta() - locBeta_Timing;
	double locFOM_Timing = (locChargedTrackHypothesis->dNDF_Timing > 0) ? TMath::Prob(locChargedTrackHypothesis->dChiSq_Timing, locChargedTrackHypothesis->dNDF_Timing) : NaN;
	double locFOM_DCdEdx = (locChargedTrackHypothesis->dNDF_DCdEdx > 0) ? TMath::Prob(locChargedTrackHypothesis->dChiSq_DCdEdx, locChargedTrackHypothesis->dNDF_DCdEdx) : NaN;

	double locP = locChargedTrackHypothesis->momentum().Mag();
	double locTheta = locChargedTrackHypothesis->momentum().Theta()*180.0/TMath::Pi();
	const DMCThrown* locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locChargedTrackHypothesis);

	Get_Application()->RootWriteLock();
	dHistMap_PIDFOM[locPID]->Fill(locChargedTrackHypothesis->dFOM);
	dHistMap_TOFFOM[locPID]->Fill(locFOM_Timing);
	dHistMap_DCdEdxFOM[locPID]->Fill(locFOM_DCdEdx);
	dHistMap_BetaVsP[locPID]->Fill(locP, locBeta_Timing);
	dHistMap_DeltaBetaVsP[locPID]->Fill(locP, locDeltaBeta);
	dHistMap_TOFFOMVsDeltaBeta[locPID]->Fill(locDeltaBeta, locFOM_Timing);
	pair<Particle_t, Particle_t> locPIDPair(locPID, Unknown); //default unless matched
	if(locMCThrown != NULL) //else bogus track (not matched to any thrown tracks)
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

	Get_Application()->RootUnLock();
}

void DHistogramAction_PID::Fill_NeutralHists(const DNeutralParticleHypothesis* locNeutralParticleHypothesis, const DMCThrownMatching* locMCThrownMatching)
{
	Particle_t locPID = locNeutralParticleHypothesis->PID();
	double locDeltaBeta, locBeta_Timing, locStartTime = 0.0; //should be from RF!!

	locBeta_Timing = locNeutralParticleHypothesis->pathLength()/(SPEED_OF_LIGHT*(locNeutralParticleHypothesis->t1() - locStartTime));
	locDeltaBeta = locNeutralParticleHypothesis->lorentzMomentum().Beta() - locBeta_Timing;

	double locP = locNeutralParticleHypothesis->momentum().Mag();
	double locTheta = locNeutralParticleHypothesis->momentum().Theta()*180.0/TMath::Pi();
	const DMCThrown* locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locNeutralParticleHypothesis);

	Get_Application()->RootWriteLock();
	dHistMap_PIDFOM[locPID]->Fill(locNeutralParticleHypothesis->dFOM);
	dHistMap_BetaVsP[locPID]->Fill(locP, locBeta_Timing);
	dHistMap_DeltaBetaVsP[locPID]->Fill(locP, locDeltaBeta);
	dHistMap_TOFFOMVsDeltaBeta[locPID]->Fill(locDeltaBeta, locNeutralParticleHypothesis->dFOM);
	pair<Particle_t, Particle_t> locPIDPair(locPID, Unknown); //default unless matched
	if(locMCThrown != NULL) //else bogus track (not matched to any thrown tracks)
		locPIDPair.second = (Particle_t)(locMCThrown->type); //matched
	if(dHistMap_PIDFOMForTruePID.find(locPIDPair) != dHistMap_PIDFOMForTruePID.end()) //else hist not created or PID is weird
		dHistMap_PIDFOMForTruePID[locPIDPair]->Fill(locNeutralParticleHypothesis->dFOM);

	if((locNeutralParticleHypothesis->dFOM < 0.01) && (locNeutralParticleHypothesis->dNDF > 0))
		dHistMap_PVsTheta_LowPIDFOM[locPID]->Fill(locTheta, locP);
	else if(locNeutralParticleHypothesis->dNDF == 0) //NaN
		dHistMap_PVsTheta_NaNPIDFOM[locPID]->Fill(locTheta, locP);

	Get_Application()->RootUnLock();
}

void DHistogramAction_TrackVertexComparison::Initialize(JEventLoop* locEventLoop)
{
	deque<deque<Particle_t> > locDetectedChargedPIDs;
	Get_Reaction()->Get_DetectedFinalChargedPIDs(locDetectedChargedPIDs);

	string locHistName, locHistTitle, locStepName, locStepROOTName, locParticleName, locParticleROOTName;
	Particle_t locPID, locHigherMassPID, locLowerMassPID;
	string locHigherMassParticleName, locLowerMassParticleName, locHigherMassParticleROOTName, locLowerMassParticleROOTName;

	size_t locNumSteps = Get_Reaction()->Get_NumReactionSteps();
	dHistDeque_TrackZToCommon.resize(locNumSteps);
	dHistDeque_TrackTToCommon.resize(locNumSteps);
	dHistDeque_TrackDOCAToCommon.resize(locNumSteps);
	dHistDeque_MaxTrackDeltaZ.resize(locNumSteps);
	dHistDeque_MaxTrackDeltaT.resize(locNumSteps);
	dHistDeque_MaxTrackDOCA.resize(locNumSteps);
	dHistDeque_TrackDeltaTVsP.resize(locNumSteps);

	//CREATE THE HISTOGRAMS
	Get_Application()->RootWriteLock(); //ACQUIRE ROOT LOCK!!
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
		}

		//delta-t vs p
		for(int loc_j = 0; loc_j < int(locDetectedChargedPIDs[loc_i].size()) - 1; ++loc_j)
		{
			locPID = locDetectedChargedPIDs[loc_i][loc_j];
			for(size_t loc_k = loc_j + 1; loc_k < locDetectedChargedPIDs[loc_i].size(); ++loc_k)
			{
				if(ParticleMass(locDetectedChargedPIDs[loc_i][loc_k]) > ParticleMass(locPID))
				{
					locHigherMassPID = locDetectedChargedPIDs[loc_i][loc_k];
					locLowerMassPID = locPID;
				}
				else
				{
					locHigherMassPID = locPID;
					locLowerMassPID = locDetectedChargedPIDs[loc_i][loc_k];
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

	Get_Application()->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_TrackVertexComparison::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, const deque<pair<const DParticleCombo*, bool> >& locPreviousParticleCombos)
{
	deque<const DKinematicData*> locParticles;
	DVector3 locVertex;
	double locVertexTime;
	Particle_t locPID;
	double locDOCA, locDeltaZ, locDeltaT, locMaxDOCA, locMaxDeltaZ, locMaxDeltaT;
	DKinFitType locKinFitType = Get_Reaction()->Get_KinFitType();
	const DKinFitResults* locKinFitResults = locParticleCombo->Get_KinFitResults();
	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);
		if(Get_AnalysisUtilities()->Find_SimilarCombos_FinalCharged(pair<const DParticleComboStep*, size_t>(locParticleComboStep, loc_i), locPreviousParticleCombos))
			continue; //duplicate, already histogrammed

		locParticleComboStep->Get_DetectedFinalChargedParticles_Measured(locParticles);
		if(locParticles.empty())
			continue;

		//Grab/Find common vertex & time
		if((locKinFitResults != NULL) && ((locKinFitType == d_VertexFit) || (locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit)))
		{
			//locKinFitResults NULL if failed kinfit
			//vertex kinfit
			deque<const DKinematicData*> locParticles_KinFit;
			locParticleComboStep->Get_DetectedFinalChargedParticles(locParticles_KinFit);
			locVertex = locParticles_KinFit[0]->position();
			if((locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndSpacetimeFit))
				locVertexTime = locParticles_KinFit[0]->time();
			else //do crude time determination
				locVertexTime = Get_AnalysisUtilities()->Calc_CrudeTime(locParticles_KinFit, locVertex);
		}
		else //do crude vertex & time determination
		{
			locVertex = Get_AnalysisUtilities()->Calc_CrudeVertex(locParticles);
			locVertexTime = Get_AnalysisUtilities()->Calc_CrudeTime(locParticles, locVertex);
		}

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

				locDOCA = Get_AnalysisUtilities()->Calc_DOCA(locParticles[loc_j], locParticles[loc_k]);
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
				if(Get_AnalysisUtilities()->Find_SimilarCombos(locParticlePairs, locPreviousParticleCombos))
					continue; //already histed!

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

			locDOCA = Get_AnalysisUtilities()->Calc_DOCAToVertex(locParticles[loc_j], locVertex);

			//HISTOGRAM //do all at once to reduce #locks & amount of time within the lock
			Get_Application()->RootWriteLock();
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
			Get_Application()->RootUnLock();
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

	string locHistName, locHistTitle, locStepName, locStepROOTName, locParticleName, locParticleROOTName;
	Particle_t locPID;

	size_t locNumSteps = Get_Reaction()->Get_NumReactionSteps();
	dHistDeque_PVsTheta.resize(locNumSteps);
	dHistDeque_PhiVsTheta.resize(locNumSteps);
	dHistDeque_P.resize(locNumSteps);
	dHistDeque_Theta.resize(locNumSteps);
	dHistDeque_Phi.resize(locNumSteps);
	dHistDeque_VertexZ.resize(locNumSteps);
	dHistDeque_VertexT.resize(locNumSteps);
	dHistDeque_VertexYVsX.resize(locNumSteps);

	deque<deque<Particle_t> > locDetectedFinalPIDs;
	Get_Reaction()->Get_DetectedFinalPIDs(locDetectedFinalPIDs);

	//CREATE THE HISTOGRAMS
	Get_Application()->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	CreateAndChangeTo_ActionDirectory();

	//beam particle
	locPID = Get_Reaction()->Get_ReactionStep(0)->Get_InitialParticleID();
	bool locBeamParticleUsed = (locPID == Gamma);
	if(locBeamParticleUsed)
	{
		locParticleName = string("Beam_") + ParticleType(locPID);
		CreateAndChangeTo_Directory(locParticleName, locParticleName);

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
	Get_Application()->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_ParticleComboKinematics::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, const deque<pair<const DParticleCombo*, bool> >& locPreviousParticleCombos)
{
	if(Get_UseKinFitResultsFlag() && (Get_Reaction()->Get_KinFitType() == d_NoFit))
	{
		cout << "WARNING: REQUESTED HISTOGRAM OF KINEMAITIC FIT RESULTS WHEN NO KINEMATIC FIT!!! Skipping histogram." << endl;
		return true; //no fit performed, but kinfit data requested!!
	}

	const DKinematicData* locKinematicData;
	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);

		//initial particle
		locKinematicData = locParticleComboStep->Get_InitialParticle();
		if(locKinematicData != NULL)
		{
			if(locParticleComboStep->Get_InitialParticleID() == Gamma)
			{
				if(Get_UseKinFitResultsFlag()) //kinfit, can be no duplicate
					Fill_BeamHists(locKinematicData);
				else if(!Get_AnalysisUtilities()->Find_SimilarCombos(pair<const DKinematicData*, size_t>(locKinematicData, loc_i), locPreviousParticleCombos)) //measured, check for dupe
					Fill_BeamHists(locKinematicData); //else duplicate
			}
			else if(Get_UseKinFitResultsFlag()) //decaying particle, but kinfit so can hist
				Fill_Hists(locKinematicData, loc_i);
		}

		//final particles
		for(size_t loc_j = 0; loc_j < locParticleComboStep->Get_NumFinalParticles(); ++loc_j)
		{
			if(Get_UseKinFitResultsFlag())
				locKinematicData = locParticleComboStep->Get_FinalParticle(loc_j);
			else
			{
				locKinematicData = locParticleComboStep->Get_FinalParticle_Measured(loc_j);
				if(locKinematicData == NULL)
					continue; //e.g. a decaying or missing particle whose params aren't set yet
				if(Get_AnalysisUtilities()->Find_SimilarCombos(pair<const DKinematicData*, size_t>(locKinematicData, loc_i), locPreviousParticleCombos))
					continue; //duplicate!
			}
			if(locKinematicData == NULL)
				continue;
			Fill_Hists(locKinematicData, loc_i);
		} //end of particle loop
	} //end of step loop
	return true;
}

void DHistogramAction_ParticleComboKinematics::Fill_Hists(const DKinematicData* locKinematicData, size_t locStepIndex)
{
	Particle_t locPID = locKinematicData->PID();
	DVector3 locMomentum = locKinematicData->momentum();
	double locPhi = locMomentum.Phi()*180.0/TMath::Pi();
	double locTheta = locMomentum.Theta()*180.0/TMath::Pi();
	double locP = locMomentum.Mag();

	Get_Application()->RootWriteLock();
	dHistDeque_P[locStepIndex][locPID]->Fill(locP);
	dHistDeque_Phi[locStepIndex][locPID]->Fill(locPhi);
	dHistDeque_Theta[locStepIndex][locPID]->Fill(locTheta);
	dHistDeque_PVsTheta[locStepIndex][locPID]->Fill(locTheta, locP);
	dHistDeque_PhiVsTheta[locStepIndex][locPID]->Fill(locTheta, locPhi);
	dHistDeque_VertexZ[locStepIndex][locPID]->Fill(locKinematicData->position().Z());
	dHistDeque_VertexYVsX[locStepIndex][locPID]->Fill(locKinematicData->position().X(), locKinematicData->position().Y());
	dHistDeque_VertexT[locStepIndex][locPID]->Fill(locKinematicData->time());
	Get_Application()->RootUnLock();
}

void DHistogramAction_ParticleComboKinematics::Fill_BeamHists(const DKinematicData* locKinematicData)
{
	DVector3 locMomentum = locKinematicData->momentum();
	double locPhi = locMomentum.Phi()*180.0/TMath::Pi();
	double locTheta = locMomentum.Theta()*180.0/TMath::Pi();
	double locP = locMomentum.Mag();

	Get_Application()->RootWriteLock();
	dBeamParticleHist_P->Fill(locP);
	dBeamParticleHist_Phi->Fill(locPhi);
	dBeamParticleHist_Theta->Fill(locTheta);
	dBeamParticleHist_PVsTheta->Fill(locTheta, locP);
	dBeamParticleHist_PhiVsTheta->Fill(locTheta, locPhi);
	dBeamParticleHist_VertexZ->Fill(locKinematicData->position().Z());
	dBeamParticleHist_VertexYVsX->Fill(locKinematicData->position().X(), locKinematicData->position().Y());
	dBeamParticleHist_VertexT->Fill(locKinematicData->time());
	Get_Application()->RootUnLock();
}

void DHistogramAction_ThrownParticleKinematics::Initialize(JEventLoop* locEventLoop)
{
	string locHistName, locHistTitle, locParticleName, locParticleROOTName;
	Particle_t locPID;

	dMCThrownMatchingFactory = static_cast<DMCThrownMatching_factory*>(locEventLoop->GetFactory("DMCThrownMatching"));

	//CREATE THE HISTOGRAMS
	Get_Application()->RootWriteLock(); //ACQUIRE ROOT LOCK!!
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
	Get_Application()->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_ThrownParticleKinematics::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, const deque<pair<const DParticleCombo*, bool> >& locPreviousParticleCombos)
{
	if(!locPreviousParticleCombos.empty())
		return true; //filled previously!

	Particle_t locPID;

	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons);
	for(size_t loc_i = 0; loc_i < locBeamPhotons.size(); ++loc_i)
		dBeamParticle_P->Fill(locBeamPhotons[loc_i]->energy());

	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);
	const DMCThrown* locMCThrown;

	for(size_t loc_i = 0; loc_i < locMCThrowns.size(); ++loc_i)
	{
		locMCThrown = locMCThrowns[loc_i];
		if(!dMCThrownMatchingFactory->Check_IsValidMCComparisonPID(locMCThrowns, locMCThrown))
			continue; //e.g. a muon from pion decay
		locPID = (Particle_t)locMCThrown->type;
		if(dHistMap_P.find(locPID) == dHistMap_P.end())
			continue; //not interested in histogramming

		DVector3 locMomentum = locMCThrown->momentum();
		double locPhi = locMomentum.Phi()*180.0/TMath::Pi();
		double locTheta = locMomentum.Theta()*180.0/TMath::Pi();
		double locP = locMomentum.Mag();
		Get_Application()->RootWriteLock();
		dHistMap_P[locPID]->Fill(locP);
		dHistMap_Phi[locPID]->Fill(locPhi);
		dHistMap_Theta[locPID]->Fill(locTheta);
		dHistMap_PVsTheta[locPID]->Fill(locTheta, locP);
		dHistMap_PhiVsTheta[locPID]->Fill(locTheta, locPhi);
		dHistMap_VertexZ[locPID]->Fill(locMCThrown->position().Z());
		dHistMap_VertexYVsX[locPID]->Fill(locMCThrown->position().X(), locMCThrown->position().Y());
		dHistMap_VertexT[locPID]->Fill(locMCThrown->time());
		Get_Application()->RootUnLock();
	}
	return true;
}

void DHistogramAction_DetectedParticleKinematics::Initialize(JEventLoop* locEventLoop)
{
	string locHistName, locHistTitle, locParticleName, locParticleROOTName;
	Particle_t locPID;

	//CREATE THE HISTOGRAMS
	Get_Application()->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	CreateAndChangeTo_ActionDirectory();

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

		// Vertex-Z
		locHistName = "VertexZ";
		locHistTitle = locParticleROOTName + string(";Vertex-Z (cm)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_VertexZ[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_VertexZ[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumVertexZBins, dMinVertexZ, dMaxVertexZ);

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
	Get_Application()->RootUnLock(); //RELEASE ROOT LOCK!!
}


bool DHistogramAction_DetectedParticleKinematics::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, const deque<pair<const DParticleCombo*, bool> >& locPreviousParticleCombos)
{
	if(!locPreviousParticleCombos.empty())
		return true; //filled previously!

	Particle_t locPID;

	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
	locEventLoop->Get(locMCThrownMatchingVector);
	const DMCThrownMatching* locMCThrownMatching = locMCThrownMatchingVector.empty() ? NULL : locMCThrownMatchingVector[0];

	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons);
	for(size_t loc_i = 0; loc_i < locBeamPhotons.size(); ++loc_i)
		dBeamParticle_P->Fill(locBeamPhotons[loc_i]->energy());

	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks);
	const DChargedTrackHypothesis* locChargedTrackHypothesis;

	for(size_t loc_i = 0; loc_i < locChargedTracks.size(); ++loc_i)
	{
		if(locMCThrownMatching != NULL)
		{
			const DMCThrown* locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locChargedTracks[loc_i]);
			if(locMCThrown == NULL)
			{
				locChargedTrackHypothesis = locChargedTracks[loc_i]->Get_BestFOM();
				locPID = locChargedTrackHypothesis->PID();
			}
			else
			{
				locPID = (Particle_t)locMCThrown->type;
				locChargedTrackHypothesis = locChargedTracks[loc_i]->Get_Hypothesis(locPID);
				if(locChargedTrackHypothesis == NULL)
					locChargedTrackHypothesis = locChargedTracks[loc_i]->Get_BestFOM();
			}
		}
		else
		{
			locChargedTrackHypothesis = locChargedTracks[loc_i]->Get_BestFOM();
			locPID = locChargedTrackHypothesis->PID();
		}

		if(dHistMap_P.find(locPID) == dHistMap_P.end())
			continue; //e.g. a decaying particle, or not interested in histogramming

		DVector3 locMomentum = locChargedTrackHypothesis->momentum();
		double locPhi = locMomentum.Phi()*180.0/TMath::Pi();
		double locTheta = locMomentum.Theta()*180.0/TMath::Pi();
		double locP = locMomentum.Mag();

		Get_Application()->RootWriteLock();
		dHistMap_P[locPID]->Fill(locP);
		dHistMap_Phi[locPID]->Fill(locPhi);
		dHistMap_Theta[locPID]->Fill(locTheta);
		dHistMap_PVsTheta[locPID]->Fill(locTheta, locP);
		dHistMap_PhiVsTheta[locPID]->Fill(locTheta, locPhi);
		dHistMap_VertexZ[locPID]->Fill(locChargedTrackHypothesis->position().Z());
		dHistMap_VertexYVsX[locPID]->Fill(locChargedTrackHypothesis->position().X(), locChargedTrackHypothesis->position().Y());
		dHistMap_VertexT[locPID]->Fill(locChargedTrackHypothesis->time());
		Get_Application()->RootUnLock();
	}

	vector<const DNeutralParticle*> locNeutralParticles;
	locEventLoop->Get(locNeutralParticles);
	const DNeutralParticleHypothesis* locNeutralParticleHypothesis;

	for(size_t loc_i = 0; loc_i < locNeutralParticles.size(); ++loc_i)
	{
		if(locMCThrownMatching != NULL)
		{
			const DMCThrown* locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locNeutralParticles[loc_i]);
			if(locMCThrown == NULL)
			{
				locNeutralParticleHypothesis = locNeutralParticles[loc_i]->Get_BestFOM();
				locPID = locNeutralParticleHypothesis->PID();
			}
			else
			{
				locPID = (Particle_t)locMCThrown->type;
				locNeutralParticleHypothesis = locNeutralParticles[loc_i]->Get_Hypothesis(locPID);
				if(locNeutralParticleHypothesis == NULL)
					locNeutralParticleHypothesis = locNeutralParticles[loc_i]->Get_BestFOM();
			}
		}
		else
		{
			locNeutralParticleHypothesis = locNeutralParticles[loc_i]->Get_BestFOM();
			locPID = locNeutralParticleHypothesis->PID();
		}

		if(dHistMap_P.find(locPID) == dHistMap_P.end())
			continue; //e.g. a decaying particle, or not interested in histogramming

		DVector3 locMomentum = locNeutralParticleHypothesis->momentum();
		double locPhi = locMomentum.Phi()*180.0/TMath::Pi();
		double locTheta = locMomentum.Theta()*180.0/TMath::Pi();
		double locP = locMomentum.Mag();

		Get_Application()->RootWriteLock();
		dHistMap_P[locPID]->Fill(locP);
		dHistMap_Phi[locPID]->Fill(locPhi);
		dHistMap_Theta[locPID]->Fill(locTheta);
		dHistMap_PVsTheta[locPID]->Fill(locTheta, locP);
		dHistMap_PhiVsTheta[locPID]->Fill(locTheta, locPhi);
		dHistMap_VertexZ[locPID]->Fill(locNeutralParticleHypothesis->position().Z());
		dHistMap_VertexYVsX[locPID]->Fill(locNeutralParticleHypothesis->position().X(), locNeutralParticleHypothesis->position().Y());
		dHistMap_VertexT[locPID]->Fill(locNeutralParticleHypothesis->time());
		Get_Application()->RootUnLock();
	}
	return true;
}

void DHistogramAction_GenReconTrackComparison::Initialize(JEventLoop* locEventLoop)
{
	string locHistName, locHistTitle, locParticleName, locParticleROOTName;
	Particle_t locPID;

	//CREATE THE HISTOGRAMS
	Get_Application()->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	CreateAndChangeTo_ActionDirectory();
	for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
	{
		locPID = dFinalStatePIDs[loc_i];
		locParticleName = ParticleType(locPID);
		locParticleROOTName = ParticleName_ROOT(locPID);
		CreateAndChangeTo_Directory(locParticleName, locParticleName);

		// DeltaP/P
		locHistName = string("DeltaPOverP_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";#Deltap/p (Reconstructed - Thrown)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_DeltaPOverP[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_DeltaPOverP[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

		// DeltaTheta
		locHistName = string("DeltaTheta_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";#Delta#theta#circ (Reconstructed - Thrown)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_DeltaTheta[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_DeltaTheta[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaThetaBins, dMinDeltaTheta, dMaxDeltaTheta);

		// DeltaPhi
		locHistName = string("DeltaPhi_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";#Delta#phi#circ (Reconstructed - Thrown)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_DeltaPhi[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_DeltaPhi[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

		// DeltaVertexZ
		locHistName = string("DeltaVertexZ_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";#DeltaVertex-Z (cm) (Reconstructed - Thrown)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_DeltaVertexZ[locPID] = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_DeltaVertexZ[locPID] = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumDeltaVertexZBins, dMinDeltaVertexZ, dMaxDeltaVertexZ);

		// DeltaP/P Vs P
		locHistName = string("DeltaPOverPVsP_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";p (GeV/c);#Deltap/p (Reconstructed - Thrown)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_DeltaPOverPVsP[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_DeltaPOverPVsP[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

		// DeltaP/P Vs Theta
		locHistName = string("DeltaPOverPVsTheta_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";#theta#circ;#Deltap/p (Reconstructed - Thrown)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_DeltaPOverPVsTheta[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_DeltaPOverPVsTheta[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

		// DeltaTheta Vs P
		locHistName = string("DeltaThetaVsP_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta#theta#circ (Reconstructed - Thrown)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_DeltaThetaVsP[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_DeltaThetaVsP[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaThetaBins, dMinDeltaTheta, dMaxDeltaTheta);

		// DeltaTheta Vs Theta
		locHistName = string("DeltaThetaVsTheta_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";#theta#circ;#Delta#theta#circ (Reconstructed - Thrown)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_DeltaThetaVsTheta[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_DeltaThetaVsTheta[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaThetaBins, dMinDeltaTheta, dMaxDeltaTheta);

		// DeltaPhi Vs P
		locHistName = string("DeltaPhiVsP_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";p (GeV/c);#Delta#phi#circ (Reconstructed - Thrown)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_DeltaPhiVsP[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_DeltaPhiVsP[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNumDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

		// DeltaPhi Vs Theta
		locHistName = string("DeltaPhiVsTheta_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";#theta#circ;#Delta#phi#circ (Reconstructed - Thrown)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_DeltaPhiVsTheta[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_DeltaPhiVsTheta[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

		// DeltaVertexZ Vs Theta
		locHistName = string("DeltaVertexZVsTheta_") + locParticleName;
		locHistTitle = locParticleROOTName + string(";#theta#circ;#DeltaVertex-Z (cm) (Reconstructed - Thrown)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistMap_DeltaVertexZVsTheta[locPID] = static_cast<TH2D*>(gDirectory->Get(locHistName.c_str()));
		else
			dHistMap_DeltaVertexZVsTheta[locPID] = new TH2D(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNumDeltaVertexZBins, dMinDeltaVertexZ, dMaxDeltaVertexZ);

		gDirectory->cd("..");
	}
	Get_Application()->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_GenReconTrackComparison::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, const deque<pair<const DParticleCombo*, bool> >& locPreviousParticleCombos)
{
	if(!locPreviousParticleCombos.empty())
		return true; //filled previously!

	Particle_t locPID;
	double locDeltaPOverP, locDeltaTheta, locDeltaPhi, locDeltaVertexZ;
	double locThrownP, locThrownTheta;

	const DMCThrown* locMCThrown;
	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
	locEventLoop->Get(locMCThrownMatchingVector);
	const DMCThrownMatching* locMCThrownMatching = locMCThrownMatchingVector[0];

	//charged particles
	const DChargedTrackHypothesis* locChargedTrackHypothesis;
	map<const DMCThrown*, const DChargedTrack*> locThrownToChargedMap;
	locMCThrownMatching->Get_ThrownToChargedMap(locThrownToChargedMap);
	for(map<const DMCThrown*, const DChargedTrack*>::iterator locIterator = locThrownToChargedMap.begin(); locIterator != locThrownToChargedMap.end(); ++locIterator)
	{
		locMCThrown = locIterator->first;
		locChargedTrackHypothesis = locIterator->second->Get_BestFOM();
		locPID = (Particle_t)locMCThrown->type;
		if(dHistMap_DeltaPOverP.find(locPID) == dHistMap_DeltaPOverP.end())
			continue; //e.g. not interested in histogramming

		locThrownP = locMCThrown->momentum().Mag();
		locThrownTheta = locMCThrown->momentum().Theta()*180.0/TMath::Pi();
		locDeltaPOverP = (locChargedTrackHypothesis->momentum().Mag() - locThrownP)/locThrownP;
		locDeltaTheta = locChargedTrackHypothesis->momentum().Theta()*180.0/TMath::Pi() - locThrownTheta;
		locDeltaPhi = locChargedTrackHypothesis->momentum().Phi()*180.0/TMath::Pi() - locMCThrown->momentum().Phi()*180.0/TMath::Pi();
		locDeltaVertexZ = locChargedTrackHypothesis->position().Z() - locMCThrown->position().Z();

		Get_Application()->RootWriteLock();
		dHistMap_DeltaPOverP[locPID]->Fill(locDeltaPOverP);
		dHistMap_DeltaTheta[locPID]->Fill(locDeltaTheta);
		dHistMap_DeltaPhi[locPID]->Fill(locDeltaPhi);
		dHistMap_DeltaVertexZ[locPID]->Fill(locDeltaVertexZ);
		dHistMap_DeltaPOverPVsP[locPID]->Fill(locThrownP, locDeltaPOverP);
		dHistMap_DeltaPOverPVsTheta[locPID]->Fill(locThrownTheta, locDeltaPOverP);
		dHistMap_DeltaThetaVsP[locPID]->Fill(locThrownP, locDeltaTheta);
		dHistMap_DeltaThetaVsTheta[locPID]->Fill(locThrownTheta, locDeltaTheta);
		dHistMap_DeltaPhiVsP[locPID]->Fill(locThrownP, locDeltaPhi);
		dHistMap_DeltaPhiVsTheta[locPID]->Fill(locThrownTheta, locDeltaPhi);
		dHistMap_DeltaVertexZVsTheta[locPID]->Fill(locThrownTheta, locDeltaVertexZ);
		Get_Application()->RootUnLock();
	}

	//neutral particles
	const DNeutralParticleHypothesis* locNeutralParticleHypothesis;
	map<const DMCThrown*, const DNeutralParticle*> locThrownToNeutralMap;
	locMCThrownMatching->Get_ThrownToNeutralMap(locThrownToNeutralMap);
	for(map<const DMCThrown*, const DNeutralParticle*>::iterator locIterator = locThrownToNeutralMap.begin(); locIterator != locThrownToNeutralMap.end(); ++locIterator)
	{
		locMCThrown = locIterator->first;
		locNeutralParticleHypothesis = locIterator->second->Get_BestFOM();
		locPID = (Particle_t)locMCThrown->type;
		if(dHistMap_DeltaPOverP.find(locPID) == dHistMap_DeltaPOverP.end())
			continue; //e.g. not interested in histogramming

		locThrownP = locMCThrown->momentum().Mag();
		locThrownTheta = locMCThrown->momentum().Theta()*180.0/TMath::Pi();
		locDeltaPOverP = (locNeutralParticleHypothesis->momentum().Mag() - locThrownP)/locThrownP;
		locDeltaTheta = locNeutralParticleHypothesis->momentum().Theta()*180.0/TMath::Pi() - locThrownTheta;
		locDeltaPhi = locNeutralParticleHypothesis->momentum().Phi()*180.0/TMath::Pi() - locMCThrown->momentum().Phi()*180.0/TMath::Pi();
		locDeltaVertexZ = locNeutralParticleHypothesis->position().Z() - locMCThrown->position().Z();

		Get_Application()->RootWriteLock();
		dHistMap_DeltaPOverP[locPID]->Fill(locDeltaPOverP);
		dHistMap_DeltaTheta[locPID]->Fill(locDeltaTheta);
		dHistMap_DeltaPhi[locPID]->Fill(locDeltaPhi);
		dHistMap_DeltaVertexZ[locPID]->Fill(locDeltaVertexZ);
		dHistMap_DeltaPOverPVsP[locPID]->Fill(locThrownP, locDeltaPOverP);
		dHistMap_DeltaPOverPVsTheta[locPID]->Fill(locThrownTheta, locDeltaPOverP);
		dHistMap_DeltaThetaVsP[locPID]->Fill(locThrownP, locDeltaTheta);
		dHistMap_DeltaThetaVsTheta[locPID]->Fill(locThrownTheta, locDeltaTheta);
		dHistMap_DeltaPhiVsP[locPID]->Fill(locThrownP, locDeltaPhi);
		dHistMap_DeltaPhiVsTheta[locPID]->Fill(locThrownTheta, locDeltaPhi);
		dHistMap_DeltaVertexZVsTheta[locPID]->Fill(locThrownTheta, locDeltaVertexZ);
		Get_Application()->RootUnLock();
	}
	return true;
}

void DHistogramAction_TrackMultiplicity::Initialize(JEventLoop* locEventLoop)
{
	//CREATE THE HISTOGRAMS
	Get_Application()->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	CreateAndChangeTo_ActionDirectory();

	string locLabelName;
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
			locLabelName = string("# ") + string(ParticleName_ROOT(dFinalStatePIDs[loc_i]));
			dHist_NumReconstructedTracks->GetXaxis()->SetBinLabel(6 + loc_i, locLabelName.c_str());
		}
	}

	Get_Application()->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_TrackMultiplicity::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, const deque<pair<const DParticleCombo*, bool> >& locPreviousParticleCombos)
{
	if(!locPreviousParticleCombos.empty())
		return true; //filled previously!

	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks);

	vector<const DNeutralParticle*> locNeutralParticles;
	locEventLoop->Get(locNeutralParticles);

	const DMCThrown* locMCThrown;
	Particle_t locPID;
	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
	locEventLoop->Get(locMCThrownMatchingVector);
	const DMCThrownMatching* locMCThrownMatching = locMCThrownMatchingVector[0];

	size_t locNumPositiveTracks = 0;
	size_t locNumNegativeTracks = 0;
	for(size_t loc_i = 0; loc_i < locChargedTracks.size(); ++loc_i)
	{
		if(locChargedTracks[loc_i]->Get_Charge() > 0.0)
			++locNumPositiveTracks;
		else
			++locNumNegativeTracks;
	}

	// get #tracks by pid type //USES MC PID IF EXISTS!!!
	map<Particle_t, size_t> locNumTracksByPID;
	for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
		locNumTracksByPID[dFinalStatePIDs[loc_i]] = 0;

	// charged by pid
	for(size_t loc_i = 0; loc_i < locChargedTracks.size(); ++loc_i)
	{
		if(locMCThrownMatching != NULL)
		{
			locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locChargedTracks[loc_i]);
			if(locMCThrown == NULL)
				locPID = locChargedTracks[loc_i]->Get_BestFOM()->PID();
			else
				locPID = (Particle_t)locMCThrown->type;
		}
		else
			locPID = locChargedTracks[loc_i]->Get_BestFOM()->PID();

		if(locNumTracksByPID.find(locPID) != locNumTracksByPID.end())
			++locNumTracksByPID[locPID];
	}

	// neutrals by pid
	for(size_t loc_i = 0; loc_i < locNeutralParticles.size(); ++loc_i)
	{
		if(locMCThrownMatching != NULL)
		{
			locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locNeutralParticles[loc_i]);
			if(locMCThrown == NULL)
				locPID = locNeutralParticles[loc_i]->Get_BestFOM()->PID();
			else
				locPID = (Particle_t)locMCThrown->type;
		}
		else
			locPID = locNeutralParticles[loc_i]->Get_BestFOM()->PID();

		if(locNumTracksByPID.find(locPID) != locNumTracksByPID.end())
			++locNumTracksByPID[locPID];
	}

	Get_Application()->RootWriteLock();
	dHist_NumReconstructedTracks->Fill(0.0, (Double_t)(locChargedTracks.size() + locNeutralParticles.size()));
	dHist_NumReconstructedTracks->Fill(1.0, (Double_t)locNumPositiveTracks);
	dHist_NumReconstructedTracks->Fill(2.0, (Double_t)locNumNegativeTracks);
	dHist_NumReconstructedTracks->Fill(3.0, (Double_t)locNeutralParticles.size());
	dHist_NumReconstructedTracks->Fill(4.0, (Double_t)locChargedTracks.size());
	for(size_t loc_i = 0; loc_i < dFinalStatePIDs.size(); ++loc_i)
		dHist_NumReconstructedTracks->Fill(5.0 + (Double_t)loc_i, (Double_t)locNumTracksByPID[dFinalStatePIDs[loc_i]]);
	Get_Application()->RootUnLock();

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

	//CREATE THE HISTOGRAMS
	Get_Application()->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	CreateAndChangeTo_ActionDirectory();
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
	Get_Application()->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_TruePID::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, const deque<pair<const DParticleCombo*, bool> >& locPreviousParticleCombos)
{
	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
	locEventLoop->Get(locMCThrownMatchingVector);
	const DMCThrownMatching* locMCThrownMatching = locMCThrownMatchingVector[0];
	double locP, locTheta;
	const DMCThrown* locMCThrown;
	Particle_t locPID;

	deque<const DKinematicData*> locParticles;
	int locComboTruePIDStatus = 1;
	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);
		locParticleComboStep->Get_DetectedFinalParticles_Measured(locParticles);

		for(size_t loc_j = 0; loc_j < locParticles.size(); ++loc_j)
		{
			locPID = locParticles[loc_j]->PID();

			if(ParticleCharge(locPID) == 0)
				locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(static_cast<const DNeutralParticleHypothesis*>(locParticles[loc_j]));
			else
				locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(static_cast<const DChargedTrackHypothesis*>(locParticles[loc_j]));

			bool locCutResult = (locMCThrown == NULL) ? false : (((Particle_t)locMCThrown->type) == locPID);
			if(!locCutResult)
				locComboTruePIDStatus = 0;

			if(Get_AnalysisUtilities()->Find_SimilarCombos(pair<const DKinematicData*, size_t>(locParticles[loc_j], loc_i), locPreviousParticleCombos))
				continue;

			locP = locParticles[loc_j]->momentum().Mag();
			locTheta = locParticles[loc_j]->momentum().Theta()*180.0/TMath::Pi();

			Get_Application()->RootWriteLock();
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
			Get_Application()->RootUnLock();
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

			DLorentzVector locFinalStateP4 = Get_AnalysisUtilities()->Calc_FinalStateP4(locParticleCombo, loc_i, Get_UseKinFitResultsFlag());
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

	//setup dFinalParticleNames: all possible pid combos that can give rise to the desired particle
	deque<string> locParticleNamesForHist;
	Get_Reaction()->Get_DecayChainFinalParticlesROOTNames(dInitialPID, dFinalParticleNames, locParticleNamesForHist, Get_UseKinFitResultsFlag());

	//CREATE THE HISTOGRAMS
	Get_Application()->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	CreateAndChangeTo_ActionDirectory();
	for(size_t loc_i = 0; loc_i < locParticleNamesForHist.size(); ++loc_i)
	{
		locHistName = "InvariantMass";
		if(locParticleNamesForHist.size() > 1)
			locHistName += string("_") + locParticleNamesForHist[loc_i];
		ostringstream locStream;
		locStream << locMassPerBin;
		locHistTitle = string(";") + locParticleNamesForHist[loc_i] + string(" Invariant Mass (GeV/c^{2});# Combos / ") + locStream.str() + string(" MeV/c^{2}");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
			dHistDeque_InvaraintMass.push_back(static_cast<TH1D*>(gDirectory->Get(locHistName.c_str())));
		else
			dHistDeque_InvaraintMass.push_back(new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumMassBins, dMinMass, dMaxMass));
	}
	Get_Application()->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_InvariantMass::Compare_ParticleNames(const deque<string>& locParticleNames1, const deque<string>& locParticleNames2) const
{
	deque<string>::const_iterator locIterator3;
	deque<string>::iterator locIterator4;
	if(locParticleNames1.size() != locParticleNames2.size())
		return false; //not same size, clearly can't be the same
	deque<string> locParticleNamesCopy = locParticleNames2;

	//loop over the lists of particles, see if they're identical
	for(locIterator3 = locParticleNames1.begin(); locIterator3 != locParticleNames1.end(); ++locIterator3)
	{
		for(locIterator4 = locParticleNamesCopy.begin(); locIterator4 != locParticleNamesCopy.end(); ++locIterator4)
		{
			if((*locIterator3) == (*locIterator4))
			{
				locParticleNamesCopy.erase(locIterator4); //particle name is identical, remove it from the list of remaining names
				break;
			}
		}
	}
	return locParticleNamesCopy.empty(); //all names removed means all names matched: duplicate
}

bool DHistogramAction_InvariantMass::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, const deque<pair<const DParticleCombo*, bool> >& locPreviousParticleCombos)
{
	deque<const DParticleComboStep*> locParticleComboSteps;
	deque<const DKinematicData*> locParticles;
	string locFinalParticleName;
	DLorentzVector locFinalStateP4;
	double locInvariantMass;

	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);
		if(locParticleComboStep->Get_InitialParticleID() != dInitialPID)
			continue;
		if(!Get_UseKinFitResultsFlag()) //measured
		{
			if((!dEnableDoubleCounting) && Get_AnalysisUtilities()->Find_SimilarCombos_AnyStep(locParticleCombo, loc_i, locPreviousParticleCombos))
				continue; //dupe: already histed!
		}

		locFinalStateP4 = Get_AnalysisUtilities()->Calc_FinalStateP4(locParticleCombo, loc_i, Get_UseKinFitResultsFlag());

		//get particle names so can select the correct histogram
		deque<string> locParticleNames;
		locParticleCombo->Get_DecayChainFinalParticlesROOTName(loc_i, locParticleNames, Get_UseKinFitResultsFlag());

		for(size_t loc_j = 0; loc_j < dFinalParticleNames.size(); ++loc_j)
		{
			if(Compare_ParticleNames(dFinalParticleNames[loc_j], locParticleNames))
			{
				locInvariantMass = locFinalStateP4.M();
				Get_Application()->RootWriteLock();
				dHistDeque_InvaraintMass[loc_j]->Fill(locInvariantMass);
				Get_Application()->RootUnLock();
				break;
			}
		}
	}
	return true;
}

void DHistogramAction_MissingMass::Initialize(JEventLoop* locEventLoop)
{
	string locHistName, locHistTitle;
	double locMassPerBin = 1000.0*(dMaxMass - dMinMass)/((double)dNumMassBins);
	string locFinalParticlesROOTName = Get_Reaction()->Get_DetectedParticlesROOTName();
	string locInitialParticlesROOTName = Get_Reaction()->Get_InitialParticlesROOTName();

	//CREATE THE HISTOGRAMS
	Get_Application()->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	CreateAndChangeTo_ActionDirectory();
	locHistName = "MissingMass";
	ostringstream locStream;
	locStream << locMassPerBin;
	locHistTitle = string(";") + locInitialParticlesROOTName + string("#rightarrow") + locFinalParticlesROOTName + string(" Missing Mass (GeV/c^{2});# Combos / ") + locStream.str() + string(" MeV/c^{2}");
	if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
		dHist_MissingMass = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
	else
		dHist_MissingMass = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumMassBins, dMinMass, dMaxMass);
	Get_Application()->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_MissingMass::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, const deque<pair<const DParticleCombo*, bool> >& locPreviousParticleCombos)
{
	//no need to check for duplicates: the missing mass for each DParticleCombo is guaranteed to be unique! (at least one different detected particle)
	double locMissingMass = (Get_AnalysisUtilities()->Calc_MissingP4(locParticleCombo, Get_UseKinFitResultsFlag())).M();
	Get_Application()->RootWriteLock();
	dHist_MissingMass->Fill(locMissingMass);
	Get_Application()->RootUnLock();
	return true;
}

void DHistogramAction_MissingMassSquared::Initialize(JEventLoop* locEventLoop)
{
	string locHistName, locHistTitle;
	double locMassSqPerBin = 1000.0*1000.0*(dMaxMassSq - dMinMassSq)/((double)dNumMassBins);
	string locFinalParticlesROOTName = Get_Reaction()->Get_DetectedParticlesROOTName();
	string locInitialParticlesROOTName = Get_Reaction()->Get_InitialParticlesROOTName();

	//CREATE THE HISTOGRAMS
	Get_Application()->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	CreateAndChangeTo_ActionDirectory();
	locHistName = "MissingMassSquared";
	ostringstream locStream;
	locStream << locMassSqPerBin;
	locHistTitle = string(";") + locInitialParticlesROOTName + string("#rightarrow") + locFinalParticlesROOTName + string(" Missing Mass Squared (GeV/c^{2})^{2};# Combos / ") + locStream.str() + string(" (MeV/c^{2})^{2}");
	if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
		dHist_MissingMassSquared = static_cast<TH1D*>(gDirectory->Get(locHistName.c_str()));
	else
		dHist_MissingMassSquared = new TH1D(locHistName.c_str(), locHistTitle.c_str(), dNumMassBins, dMinMassSq, dMaxMassSq);
	Get_Application()->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DHistogramAction_MissingMassSquared::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, const deque<pair<const DParticleCombo*, bool> >& locPreviousParticleCombos)
{
	//no need to check for duplicates: the missing mass for each DParticleCombo is guaranteed to be unique! (at least one different detected particle)
	double locMissingMassSq = (Get_AnalysisUtilities()->Calc_MissingP4(locParticleCombo, Get_UseKinFitResultsFlag())).M2();
	Get_Application()->RootWriteLock();
	dHist_MissingMassSquared->Fill(locMissingMassSq);
	Get_Application()->RootUnLock();
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

	//CREATE THE HISTOGRAMS
	Get_Application()->RootWriteLock(); //ACQUIRE ROOT LOCK!!
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
		Create_ParticlePulls("Beam", Gamma, dHistMap_BeamPulls, locKinFitTypeString);
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

			Create_ParticlePulls(locStepROOTName, locPID, locParticlePulls, locKinFitTypeString);
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

	Get_Application()->RootUnLock(); //RELEASE ROOT LOCK!!
}

void DHistogramAction_KinFitResults::Create_ParticlePulls(string locStepROOTName, Particle_t locPID, map<DKinFitPullType, TH1D*>& locParticlePulls, const string& locKinFitTypeString)
{
	string locHistName, locHistTitle, locParticleName, locParticleROOTName;
	locParticleName = ParticleType(locPID);
	locParticleROOTName = ParticleName_ROOT(locPID);

	DKinFitType locKinFitType = Get_Reaction()->Get_KinFitType();

	locParticlePulls.clear();
	if((ParticleCharge(locPID) == 0) && (locStepROOTName != "Beam"))
	{
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
	if((locKinFitType == d_VertexFit) || (locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit))
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

bool DHistogramAction_KinFitResults::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, const deque<pair<const DParticleCombo*, bool> >& locPreviousParticleCombos)
{
	//kinfit results are unique for each DParticleCombo: no need to check for duplicates
	const DKinFitResults* locKinFitResults = locParticleCombo->Get_KinFitResults();
	if(locKinFitResults == NULL)
		return true;

	// Confidence Level
	double locConfidenceLevel = locKinFitResults->Get_ConfidenceLevel();
	Get_Application()->RootWriteLock();
	dHist_ConfidenceLevel->Fill(locConfidenceLevel);
	Get_Application()->RootUnLock();

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
		Get_Application()->RootWriteLock();
		for(locIterator = locParticlePulls.begin(); locIterator != locParticlePulls.end(); ++locIterator)
			dHistMap_BeamPulls[locIterator->first]->Fill(locIterator->second);
		Get_Application()->RootUnLock();
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
			Get_Application()->RootWriteLock();
			for(locIterator = locParticlePulls.begin(); locIterator != locParticlePulls.end(); ++locIterator)
				(dHistMap_Pulls[locParticlePair])[locIterator->first]->Fill(locIterator->second);
			Get_Application()->RootUnLock();
		}
	}

	//rf time pull
	DKinFitType locKinFitType = Get_Reaction()->Get_KinFitType();
	if(locBeamFlag && ((locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndSpacetimeFit)))
	{
		locParticlePulls = locPulls[NULL];
		Get_Application()->RootWriteLock();
		dHist_RFTimePull->Fill(locParticlePulls[d_TPull]);
		Get_Application()->RootUnLock();
	}

	return true;
}

