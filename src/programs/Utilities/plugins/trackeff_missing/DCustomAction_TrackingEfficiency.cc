// $Id$
//
//    File: DCustomAction_TrackingEfficiency.cc
// Created: Wed Feb 25 09:38:06 EST 2015
// Creator: pmatt (on Linux pmattdesktop.jlab.org 2.6.32-504.8.1.el6.x86_64 x86_64)
//

#include "DCustomAction_TrackingEfficiency.h"

void DCustomAction_TrackingEfficiency::Initialize(JEventLoop* locEventLoop)
{
	//Optional: Create histograms and/or modify member variables.
	//Create any histograms/trees/etc. within a ROOT lock. 
		//This is so that when running multithreaded, only one thread is writing to the ROOT file at a time. 

	if(!Get_Reaction()->Get_MissingPID(dMissingPID))
		return; //invalid reaction setup
	if(dMissingPID == Unknown)
		return; //invalid reaction setup

	locEventLoop->GetSingle(dAnalysisUtilities);

	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry* locGeometry = locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());
	double locTargetZCenter = 0.0, locTargetLength = 1.0;
	locGeometry->GetTargetZ(locTargetZCenter);
	locGeometry->GetTargetLength(locTargetLength);
cout << "center, length = " << locTargetZCenter << ", " << locTargetLength << endl;
	dMinVertexZ = locTargetZCenter - 0.5*locTargetLength;
	dMaxVertexZ = locTargetZCenter + 0.5*locTargetLength;
	dVertexZBinSize = locTargetLength/double(dNumVertexZBins);

	bool locIsRESTEvent = (string(locEventLoop->GetJEvent().GetJEventSource()->className()) == string("DEventSourceREST"));

	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//Required: Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this action.
			//If another thread has already created the folder, it just changes to it. 
		CreateAndChangeTo_ActionDirectory();

		if(!locIsRESTEvent)
			Create_ResolutionHists(false);
		Create_ResolutionHists(true);

		if(!locIsRESTEvent)
			Create_EfficiencyHists(false);
		Create_EfficiencyHists(true);

		Create_PIDHists();
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

void DCustomAction_TrackingEfficiency::Create_ResolutionHists(bool locIsTimeBasedFlag)
{
	string locHistName, locHistTitle;

	string locDirectoryName = locIsTimeBasedFlag ? "Resolution_TimeBased" : "Resolution_WireBased";
	CreateAndChangeTo_Directory(locDirectoryName.c_str(), locDirectoryName.c_str());
	string locHistParticleName = locIsTimeBasedFlag ? "Time-Based " : "Wire-Based ";
	locHistParticleName += ParticleName_ROOT(dMissingPID);

	// DeltaP/P Vs P
	locHistName = string("DeltaPOverPVsP");
	locHistTitle = locHistParticleName + string(";p (GeV/c);#Deltap/p (Measured - Missing)");
	if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
		dHistMap_Resolution_DeltaPOverPVsP[locIsTimeBasedFlag] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
	else
		dHistMap_Resolution_DeltaPOverPVsP[locIsTimeBasedFlag] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

	// DeltaP/P Vs Theta
	locHistName = string("DeltaPOverPVsTheta");
	locHistTitle = locHistParticleName + string(";#theta#circ;#Deltap/p (Measured - Missing)");
	if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
		dHistMap_Resolution_DeltaPOverPVsTheta[locIsTimeBasedFlag] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
	else
		dHistMap_Resolution_DeltaPOverPVsTheta[locIsTimeBasedFlag] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

	// DeltaTheta Vs P
	locHistName = string("DeltaThetaVsP");
	locHistTitle = locHistParticleName + string(";p (GeV/c);#Delta#theta#circ (Measured - Missing)");
	if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
		dHistMap_Resolution_DeltaThetaVsP[locIsTimeBasedFlag] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
	else
		dHistMap_Resolution_DeltaThetaVsP[locIsTimeBasedFlag] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DDeltaThetaBins, dMinDeltaTheta, dMaxDeltaTheta);

	// DeltaTheta Vs Theta
	locHistName = string("DeltaThetaVsTheta");
	locHistTitle = locHistParticleName + string(";#theta#circ;#Delta#theta#circ (Measured - Missing)");
	if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
		dHistMap_Resolution_DeltaThetaVsTheta[locIsTimeBasedFlag] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
	else
		dHistMap_Resolution_DeltaThetaVsTheta[locIsTimeBasedFlag] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DDeltaThetaBins, dMinDeltaTheta, dMaxDeltaTheta);

	// DeltaPhi Vs P
	locHistName = string("DeltaPhiVsP");
	locHistTitle = locHistParticleName + string(";p (GeV/c);#Delta#phi#circ (Measured - Missing)");
	if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
		dHistMap_Resolution_DeltaPhiVsP[locIsTimeBasedFlag] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
	else
		dHistMap_Resolution_DeltaPhiVsP[locIsTimeBasedFlag] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

	// DeltaPhi Vs Theta
	locHistName = string("DeltaPhiVsTheta");
	locHistTitle = locHistParticleName + string(";#theta#circ;#Delta#phi#circ (Measured - Missing)");
	if(gDirectory->Get(locHistName.c_str()) != NULL) //already created by another thread, or directory name is duplicate (e.g. two identical steps)
		dHistMap_Resolution_DeltaPhiVsTheta[locIsTimeBasedFlag] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
	else
		dHistMap_Resolution_DeltaPhiVsTheta[locIsTimeBasedFlag] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

	gDirectory->cd(".."); //return to the action directory
}

void DCustomAction_TrackingEfficiency::Create_EfficiencyHists(bool locIsTimeBasedFlag)
{
	string locHistName, locHistTitle;

	string locDirectoryName = locIsTimeBasedFlag ? "Efficiency_TimeBased" : "Efficiency_WireBased";
	CreateAndChangeTo_Directory(locDirectoryName.c_str(), locDirectoryName.c_str());
	string locHistParticleName = locIsTimeBasedFlag ? "Time-Based " : "Wire-Based ";
	locHistParticleName += ParticleName_ROOT(dMissingPID);

	//Resize vectors
	dHistMap_TrackFound_PVsTheta[locIsTimeBasedFlag].resize(dNumVertexZBins);
	dHistMap_TrackMissing_PVsTheta[locIsTimeBasedFlag].resize(dNumVertexZBins);

	dHistMap_TrackFound_PhiVsTheta[locIsTimeBasedFlag].resize(dNumVertexZBins);
	dHistMap_TrackMissing_PhiVsTheta[locIsTimeBasedFlag].resize(dNumVertexZBins);

	dHistMap_FoundHasDetectorMatch_PVsTheta[locIsTimeBasedFlag].resize(dNumVertexZBins);
	dHistMap_FoundNoDetectorMatch_PVsTheta[locIsTimeBasedFlag].resize(dNumVertexZBins);

	dHistMap_FoundHasDetectorMatch_PhiVsTheta[locIsTimeBasedFlag].resize(dNumVertexZBins);
	dHistMap_FoundNoDetectorMatch_PhiVsTheta[locIsTimeBasedFlag].resize(dNumVertexZBins);

	locHistName = "MatchingFOM";
	locHistTitle = locHistParticleName + string(";Missing Match FOM");
	if(gDirectory->Get(locHistName.c_str()) != NULL) //check to see if already created by another thread
		dHistMap_MatchingFOM[locIsTimeBasedFlag] = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
	else //already created by another thread
		dHistMap_MatchingFOM[locIsTimeBasedFlag] = new TH1I(locHistName.c_str(), locHistTitle.c_str(), dNumMatchFOMBins, 0, 1.0);

	for(size_t locVertexZBin = 0; locVertexZBin < dNumVertexZBins; ++locVertexZBin)
	{
		double locBinMinVertZ = dMinVertexZ + double(locVertexZBin)*dVertexZBinSize;
		double locBinMaxVertZ = locBinMinVertZ + dVertexZBinSize;

		ostringstream locVertexZBinStream, locVertexZRangeStream;
		locVertexZBinStream << locVertexZBin;
		locVertexZRangeStream << locBinMinVertZ << " #leq Vertex-Z (cm) < " << locBinMaxVertZ;

		//Reconstruction
		locHistName = string("PVsTheta_TrackFound_VertexZBin") + locVertexZBinStream.str();
		locHistTitle = locHistParticleName + string(", ") + locVertexZRangeStream.str() + string(", Track Found;#theta#circ;p (GeV/c)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //check to see if already created by another thread
			dHistMap_TrackFound_PVsTheta[locIsTimeBasedFlag][locVertexZBin] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
		else //already created by another thread
			dHistMap_TrackFound_PVsTheta[locIsTimeBasedFlag][locVertexZBin] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

		locHistName = string("PVsTheta_TrackMissing_VertexZBin") + locVertexZBinStream.str();
		locHistTitle = locHistParticleName + string(", ") + locVertexZRangeStream.str() + string(", Track Missing;#theta#circ;p (GeV/c)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //check to see if already created by another thread
			dHistMap_TrackMissing_PVsTheta[locIsTimeBasedFlag][locVertexZBin] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
		else //already created by another thread
			dHistMap_TrackMissing_PVsTheta[locIsTimeBasedFlag][locVertexZBin] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

		locHistName = string("PhiVsTheta_TrackFound_VertexZBin") + locVertexZBinStream.str();
		locHistTitle = locHistParticleName + string(", ") + locVertexZRangeStream.str() + string(", Track Found;#theta#circ;#phi#circ");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //check to see if already created by another thread
			dHistMap_TrackFound_PhiVsTheta[locIsTimeBasedFlag][locVertexZBin] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
		else //already created by another thread
			dHistMap_TrackFound_PhiVsTheta[locIsTimeBasedFlag][locVertexZBin] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

		locHistName = string("PhiVsTheta_TrackMissing_VertexZBin") + locVertexZBinStream.str();
		locHistTitle = locHistParticleName + string(", ") + locVertexZRangeStream.str() + string(", Track Missing;#theta#circ;#phi#circ");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //check to see if already created by another thread
			dHistMap_TrackMissing_PhiVsTheta[locIsTimeBasedFlag][locVertexZBin] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
		else //already created by another thread
			dHistMap_TrackMissing_PhiVsTheta[locIsTimeBasedFlag][locVertexZBin] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

		//Detector Match
		locHistName = string("PVsTheta_FoundHasDetectorMatch_VertexZBin") + locVertexZBinStream.str();
		locHistTitle = locHistParticleName + string(", ") + locVertexZRangeStream.str() + string(", Has Detector Match;#theta#circ;p (GeV/c)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //check to see if already created by another thread
			dHistMap_FoundHasDetectorMatch_PVsTheta[locIsTimeBasedFlag][locVertexZBin] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
		else //already created by another thread
			dHistMap_FoundHasDetectorMatch_PVsTheta[locIsTimeBasedFlag][locVertexZBin] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

		locHistName = string("PVsTheta_FoundNoDetectorMatch_VertexZBin") + locVertexZBinStream.str();
		locHistTitle = locHistParticleName + string(", ") + locVertexZRangeStream.str() + string(", No Detector Match;#theta#circ;p (GeV/c)");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //check to see if already created by another thread
			dHistMap_FoundNoDetectorMatch_PVsTheta[locIsTimeBasedFlag][locVertexZBin] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
		else //already created by another thread
			dHistMap_FoundNoDetectorMatch_PVsTheta[locIsTimeBasedFlag][locVertexZBin] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

		locHistName = string("PhiVsTheta_FoundHasDetectorMatch_VertexZBin") + locVertexZBinStream.str();
		locHistTitle = locHistParticleName + string(", ") + locVertexZRangeStream.str() + string(", Has Detector Match;#theta#circ;#phi#circ");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //check to see if already created by another thread
			dHistMap_FoundHasDetectorMatch_PhiVsTheta[locIsTimeBasedFlag][locVertexZBin] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
		else //already created by another thread
			dHistMap_FoundHasDetectorMatch_PhiVsTheta[locIsTimeBasedFlag][locVertexZBin] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

		locHistName = string("PhiVsTheta_FoundNoDetectorMatch_VertexZBin") + locVertexZBinStream.str();
		locHistTitle = locHistParticleName + string(", ") + locVertexZRangeStream.str() + string(", No Detector Match;#theta#circ;#phi#circ");
		if(gDirectory->Get(locHistName.c_str()) != NULL) //check to see if already created by another thread
			dHistMap_FoundNoDetectorMatch_PhiVsTheta[locIsTimeBasedFlag][locVertexZBin] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
		else //already created by another thread
			dHistMap_FoundNoDetectorMatch_PhiVsTheta[locIsTimeBasedFlag][locVertexZBin] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);
	}

	gDirectory->cd(".."); //return to the action directory
}

void DCustomAction_TrackingEfficiency::Create_PIDHists(void)
{
	string locHistName, locHistTitle;

	CreateAndChangeTo_Directory("PID", "PID");
	string locHistParticleName = ParticleName_ROOT(dMissingPID);

	//dE/dx
	vector<DetectorSystem_t> locDetectors_dEdx;
	locDetectors_dEdx.push_back(SYS_CDC);  locDetectors_dEdx.push_back(SYS_FDC);  locDetectors_dEdx.push_back(SYS_TOF);  locDetectors_dEdx.push_back(SYS_START);
	for(size_t loc_i = 0; loc_i < locDetectors_dEdx.size(); ++loc_i)
	{
		DetectorSystem_t locSystem = locDetectors_dEdx[loc_i];
		string locUnitsString = ((locSystem == SYS_CDC) || (locSystem == SYS_FDC)) ? "(keV/cm)" : "(MeV/cm)";

		locHistName = string("dEdXVsP_") + SystemName(locSystem);
		locHistTitle = locHistParticleName + string(";p (GeV/c);") + SystemName(locSystem) + string(" dE/dX ") + locUnitsString;

		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHistMap_dEdXVsP[locSystem] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, dMaxP, dNum2DdEdxBins, dMindEdX, dMaxdEdX);
		else //already created by another thread
			dHistMap_dEdXVsP[locSystem] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
	}

	//beta vs p
	vector<DetectorSystem_t> locDetectors_BetaVsP;
	locDetectors_BetaVsP.push_back(SYS_BCAL);  locDetectors_BetaVsP.push_back(SYS_FCAL);  locDetectors_BetaVsP.push_back(SYS_TOF);
	for(size_t loc_i = 0; loc_i < locDetectors_BetaVsP.size(); ++loc_i)
	{
		DetectorSystem_t locSystem = locDetectors_BetaVsP[loc_i];
		double locMaxP = (locSystem == SYS_BCAL) ? dMaxPBCAL : dMaxP;

		locHistName = string("BetaVsP_") + SystemName(locSystem);
		locHistTitle = locHistParticleName + string(", ") + SystemName(locSystem) + string(";p (GeV/c);#beta");

		if(gDirectory->Get(locHistName.c_str()) == NULL) //check to see if already created by another thread
			dHistMap_BetaVsP[locSystem] = new TH2I(locHistName.c_str(), locHistTitle.c_str(), dNum2DPBins, dMinP, locMaxP, dNum2DBetaBins, dMinBeta, dMaxBeta);
		else //already created by another thread
			dHistMap_BetaVsP[locSystem] = static_cast<TH2I*>(gDirectory->Get(locHistName.c_str()));
	}
}

bool DCustomAction_TrackingEfficiency::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//Write custom code to perform an action on the INPUT DParticleCombo (DParticleCombo)
	//NEVER: Grab DParticleCombo or DAnalysisResults objects (of any tag!) from the JEventLoop within this function
	//NEVER: Grab objects that are created post-kinfit (e.g. DKinFitResults, etc.) from the JEventLoop if Get_UseKinFitResultsFlag() == false: CAN CAUSE INFINITE DEPENDENCY LOOP

	bool locUseKinFitResultsFlag = Get_UseKinFitResultsFlag();

	/*********************************************** MISSING PARTICLE INFO ***********************************************/

	if(dMissingPID == Unknown)
		return true; //invalid reaction setup
	if(ParticleCharge(dMissingPID) == 0)
		return true; //NOT SUPPORTED (YET?)
cout << "wega" << endl;
	const DKinematicData* locMissingParticle = locParticleCombo->Get_MissingParticle(); //is NULL if no kinfit!!

	// Get missing particle p4 & covariance
	DLorentzVector locMissingP4;
	DMatrixDSym locMissingCovarianceMatrix(3);
	if(locMissingParticle == NULL)
	{
		locMissingP4 = dAnalysisUtilities->Calc_MissingP4(locParticleCombo, locUseKinFitResultsFlag);
		locMissingCovarianceMatrix = dAnalysisUtilities->Calc_MissingP3Covariance(locParticleCombo);
	}
	else
	{
		locMissingP4 = locMissingParticle->lorentzMomentum();
		DMatrixDSym locKinFitCovarianceMatrix = locMissingParticle->errorMatrix();
		locKinFitCovarianceMatrix.ResizeTo(3, 3);
		locMissingCovarianceMatrix = locKinFitCovarianceMatrix;
	}

	double locVertexZ = locParticleCombo->Get_EventVertex().Z();
	int locVertexZBin = int((locVertexZ - dMinVertexZ)/dVertexZBinSize);
	cout << "wega, z, bin = " << locVertexZ << ", " << locVertexZBin << endl;
	if((locVertexZBin < 0) || (locVertexZBin >= int(dNumVertexZBins)))
		return true; //not within range
cout << "pass" << endl;
	/************************************************** TRACK CANDIDATES *************************************************/

	/*
		// VERY difficult to compare track candidates:
			// momentum is defined at some r > 0, so delta-p3 is not fair comparison (especially phi for slow tracks), also, no covariance matrix
		vector<const DTrackCandidate*> locUnusedTrackCandidates;
		dAnalysisUtilities->Get_UnusedTrackCandidates(locEventLoop, locParticleCombo, locUnusedTrackCandidates);
	*/

	/************************************************* WIRE-BASED TRACKS *************************************************/

	//Get unused tracks
	vector<const DTrackWireBased*> locUnusedWireBasedTracks;
	dAnalysisUtilities->Get_UnusedWireBasedTracks(locEventLoop, locParticleCombo, locUnusedWireBasedTracks);

	//find the best-matching wire-based track corresponding to the missing particle
	double locBestWireBasedMatchFOM = -1.0;
	const DTrackWireBased* locBestTrackWireBased = NULL;
	for(size_t loc_i = 0; loc_i < locUnusedWireBasedTracks.size(); ++loc_i)
	{
		if(locUnusedWireBasedTracks[loc_i]->PID() != dMissingPID)
			continue; //only use tracking results with correct PID

		DMatrixDSym locCovarianceMatrix = locUnusedWireBasedTracks[loc_i]->errorMatrix();
		locCovarianceMatrix.ResizeTo(3, 3);
		locCovarianceMatrix += locMissingCovarianceMatrix;

		//invert matrix
		TDecompLU locDecompLU(locCovarianceMatrix);
		//check to make sure that the matrix is decomposable and has a non-zero determinant
		if((!locDecompLU.Decompose()) || (fabs(locCovarianceMatrix.Determinant()) < 1.0E-300))
			continue; // matrix is not invertible
		locCovarianceMatrix.Invert();

		DVector3 locDeltaP3 = locUnusedWireBasedTracks[loc_i]->momentum() - locMissingP4.Vect();
		double locMatchFOM = Calc_MatchFOM(locDeltaP3, locCovarianceMatrix);

		if(locMatchFOM > locBestWireBasedMatchFOM)
		{
			locBestWireBasedMatchFOM = locMatchFOM;
			locBestTrackWireBased = locUnusedWireBasedTracks[loc_i];
		}
	}

	const DDetectorMatches* locDetectorMatches_WireBased = NULL;
	locEventLoop->GetSingle(locDetectorMatches_WireBased, "WireBased");

	bool locHasDetectorMatch_WireBased = (locBestTrackWireBased == NULL) ? false : locDetectorMatches_WireBased->Get_IsMatchedToHit(locBestTrackWireBased);
	Fill_NonPIDHistograms(locBestTrackWireBased, locMissingP4, locVertexZBin, locBestWireBasedMatchFOM, locHasDetectorMatch_WireBased, false);

	/************************************************* TIME-BASED TRACKS *************************************************/

	//Get unused tracks
	vector<const DTrackTimeBased*> locUnusedTimeBasedTracks;
	dAnalysisUtilities->Get_UnusedTimeBasedTracks(locEventLoop, locParticleCombo, locUnusedTimeBasedTracks);

	//find the best-matching time-based track corresponding to the missing particle
	double locBestTimeBasedMatchFOM = -1.0;
	const DTrackTimeBased* locBestTrackTimeBased = NULL;
	for(size_t loc_i = 0; loc_i < locUnusedTimeBasedTracks.size(); ++loc_i)
	{
		if(locUnusedTimeBasedTracks[loc_i]->PID() != dMissingPID)
			continue; //only use tracking results with correct PID

		DMatrixDSym locCovarianceMatrix = locUnusedTimeBasedTracks[loc_i]->errorMatrix();
		locCovarianceMatrix.ResizeTo(3, 3);
		locCovarianceMatrix += locMissingCovarianceMatrix;

		//invert matrix
		TDecompLU locDecompLU(locCovarianceMatrix);
		//check to make sure that the matrix is decomposable and has a non-zero determinant
		if((!locDecompLU.Decompose()) || (fabs(locCovarianceMatrix.Determinant()) < 1.0E-300))
			continue; // matrix is not invertible
		locCovarianceMatrix.Invert();

		DVector3 locDeltaP3 = locUnusedTimeBasedTracks[loc_i]->momentum() - locMissingP4.Vect();
		double locMatchFOM = Calc_MatchFOM(locDeltaP3, locCovarianceMatrix);

		if(locMatchFOM > locBestTimeBasedMatchFOM)
		{
			locBestTimeBasedMatchFOM = locMatchFOM;
			locBestTrackTimeBased = locUnusedTimeBasedTracks[loc_i];
		}
	}

	const DDetectorMatches* locDetectorMatches_TimeBased = NULL;
	locEventLoop->GetSingle(locDetectorMatches_TimeBased);

	bool locHasDetectorMatch_TimeBased = (locBestTrackTimeBased == NULL) ? false : locDetectorMatches_TimeBased->Get_IsMatchedToHit(locBestTrackTimeBased);
	Fill_NonPIDHistograms(locBestTrackTimeBased, locMissingP4, locVertexZBin, locBestTimeBasedMatchFOM, locHasDetectorMatch_TimeBased, true);

	/*********************************************** CHARGED HYPOTHESES: PID *********************************************/

	if((locBestTrackTimeBased == NULL) || (locBestTimeBasedMatchFOM < dMinTrackMatchFOM))
		return true; //No (good) charged hypothesis, don't bother

	//Get unused tracks
	vector<const DChargedTrack*> locUnusedChargedTracks;
	dAnalysisUtilities->Get_UnusedChargedTracks(locEventLoop, locParticleCombo, locUnusedChargedTracks);

	const DChargedTrackHypothesis* locBestChargedTrackHypothesis = NULL;
	for(size_t loc_i = 0; loc_i < locUnusedChargedTracks.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locUnusedChargedTracks[loc_i]->Get_Hypothesis(dMissingPID);
		if(locChargedTrackHypothesis == NULL)
			continue;

		//choose the one that matches the best time-based track
		const DTrackTimeBased* locTrackTimeBased = NULL;
		locChargedTrackHypothesis->GetSingle(locTrackTimeBased);
		if(locTrackTimeBased != locBestTrackTimeBased)
			continue;
		locBestChargedTrackHypothesis = locChargedTrackHypothesis;
		break;
	}
	if(locBestChargedTrackHypothesis == NULL)
		return true; //No charged hypothesis, don't bother //shouldn't be possible!

	//if RF time is indeterminate, start time will be NaN
	const DBCALShowerMatchParams& locBCALShowerMatchParams = locBestChargedTrackHypothesis->dBCALShowerMatchParams;
	const DFCALShowerMatchParams& locFCALShowerMatchParams = locBestChargedTrackHypothesis->dFCALShowerMatchParams;
	const DTOFHitMatchParams& locTOFHitMatchParams = locBestChargedTrackHypothesis->dTOFHitMatchParams;
	const DSCHitMatchParams& locSCHitMatchParams = locBestChargedTrackHypothesis->dSCHitMatchParams;

	double locP = locBestChargedTrackHypothesis->momentum().Mag();

	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		if(locSCHitMatchParams.dTrack != NULL)
			dHistMap_dEdXVsP[SYS_START]->Fill(locP, locSCHitMatchParams.dEdx*1.0E3);
		if(locTOFHitMatchParams.dTrack != NULL)
		{
			dHistMap_dEdXVsP[SYS_TOF]->Fill(locP, locTOFHitMatchParams.dEdx*1.0E3);
			double locBeta_Timing = locTOFHitMatchParams.dPathLength/(29.9792458*(locTOFHitMatchParams.dHitTime - locBestChargedTrackHypothesis->t0()));
			dHistMap_BetaVsP[SYS_TOF]->Fill(locP, locBeta_Timing);
		}
		if(locBCALShowerMatchParams.dTrack != NULL)
		{
			const DBCALShower* locBCALShower = locBCALShowerMatchParams.dBCALShower;
			double locBeta_Timing = locBCALShowerMatchParams.dPathLength/(29.9792458*(locBCALShower->t - locBestChargedTrackHypothesis->t0()));
			dHistMap_BetaVsP[SYS_BCAL]->Fill(locP, locBeta_Timing);
		}
		if(locFCALShowerMatchParams.dTrack != NULL)
		{
			const DFCALShower* locFCALShower = locFCALShowerMatchParams.dFCALShower;
			double locBeta_Timing = locFCALShowerMatchParams.dPathLength/(29.9792458*(locFCALShower->getTime() - locBestChargedTrackHypothesis->t0()));
			dHistMap_BetaVsP[SYS_FCAL]->Fill(locP, locBeta_Timing);
		}

		//Yes, these are time-based not hypothesis plots, but better to group here and minimize # of locks
		if(locBestTrackTimeBased->dNumHitsUsedFordEdx_CDC > 0)
			dHistMap_dEdXVsP[SYS_CDC]->Fill(locP, locBestTrackTimeBased->ddEdx_CDC*1.0E6);
		if(locBestTrackTimeBased->dNumHitsUsedFordEdx_FDC > 0)
			dHistMap_dEdXVsP[SYS_FDC]->Fill(locP, locBestTrackTimeBased->ddEdx_FDC*1.0E6);
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}

void DCustomAction_TrackingEfficiency::Fill_NonPIDHistograms(const DKinematicData* locTrack, DLorentzVector locMissingP4, size_t locVertexZBin, double locTrackMatchFOM, bool locHasDetectorMatch, bool locIsTimeBasedFlag)
{
	double locMeasuredP = (locTrack != NULL) ? locTrack->momentum().Mag() : numeric_limits<double>::quiet_NaN();
	double locMissingP = locMissingP4.P();
	double locMeasuredTheta = (locTrack != NULL) ? locTrack->momentum().Theta()*180.0/TMath::Pi() : numeric_limits<double>::quiet_NaN();
	double locMissingTheta = locMissingP4.Theta()*180.0/TMath::Pi();
	double locMeasuredPhi = (locTrack != NULL) ? locTrack->momentum().Phi()*180.0/TMath::Pi() : numeric_limits<double>::quiet_NaN();
	double locMissingPhi = locMissingP4.Phi()*180.0/TMath::Pi();

	double locDeltaPOverP = (locMeasuredP - locMissingP)/locMissingP;
	double locDeltaTheta = locMeasuredTheta - locMissingTheta;
	double locDeltaPhi = locMeasuredPhi - locMissingPhi;

	//Optional: Fill histograms
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//Resolution
		if(locTrack != NULL)
		{
			dHistMap_Resolution_DeltaPOverPVsP[locIsTimeBasedFlag]->Fill(locMissingP, locDeltaPOverP);
			dHistMap_Resolution_DeltaPOverPVsTheta[locIsTimeBasedFlag]->Fill(locMissingTheta, locDeltaPOverP);

			dHistMap_Resolution_DeltaThetaVsP[locIsTimeBasedFlag]->Fill(locMissingP, locDeltaTheta);
			dHistMap_Resolution_DeltaThetaVsTheta[locIsTimeBasedFlag]->Fill(locMissingTheta, locDeltaTheta);

			dHistMap_Resolution_DeltaPhiVsP[locIsTimeBasedFlag]->Fill(locMissingP, locDeltaPhi);
			dHistMap_Resolution_DeltaPhiVsTheta[locIsTimeBasedFlag]->Fill(locMissingTheta, locDeltaPhi);

			dHistMap_MatchingFOM[locIsTimeBasedFlag]->Fill(locTrackMatchFOM);
		}

		//Efficiency
		if((locTrack != NULL) && (locTrackMatchFOM >= dMinTrackMatchFOM))
		{
			//Found
			dHistMap_TrackFound_PVsTheta[locIsTimeBasedFlag][locVertexZBin]->Fill(locMissingTheta, locMissingP);
			dHistMap_TrackFound_PhiVsTheta[locIsTimeBasedFlag][locVertexZBin]->Fill(locMissingTheta, locMissingPhi);
			if(locHasDetectorMatch)
			{
				dHistMap_FoundHasDetectorMatch_PVsTheta[locIsTimeBasedFlag][locVertexZBin]->Fill(locMissingTheta, locMissingP);
				dHistMap_FoundHasDetectorMatch_PhiVsTheta[locIsTimeBasedFlag][locVertexZBin]->Fill(locMissingTheta, locMissingPhi);
			}
			else
			{
				dHistMap_FoundNoDetectorMatch_PVsTheta[locIsTimeBasedFlag][locVertexZBin]->Fill(locMissingTheta, locMissingP);
				dHistMap_FoundNoDetectorMatch_PhiVsTheta[locIsTimeBasedFlag][locVertexZBin]->Fill(locMissingTheta, locMissingPhi);
			}
		}
		else
		{
			//Missing
			dHistMap_TrackMissing_PVsTheta[locIsTimeBasedFlag][locVertexZBin]->Fill(locMissingTheta, locMissingP);
			dHistMap_TrackMissing_PhiVsTheta[locIsTimeBasedFlag][locVertexZBin]->Fill(locMissingTheta, locMissingPhi);
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

double DCustomAction_TrackingEfficiency::Calc_MatchFOM(const DVector3& locDeltaP3, DMatrixDSym locInverse3x3Matrix) const
{
	DMatrix locDeltas(3, 1);
	locDeltas(0, 0) = locDeltaP3.Px();
	locDeltas(1, 0) = locDeltaP3.Py();
	locDeltas(2, 0) = locDeltaP3.Pz();

	double locChiSq = (locInverse3x3Matrix.SimilarityT(locDeltas))(0, 0);
	return TMath::Prob(locChiSq, 3);
}
