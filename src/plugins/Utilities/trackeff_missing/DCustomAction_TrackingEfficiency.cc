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
	locEventLoop->GetSingle(dParticleID);

	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry* locGeometry = locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());
	double locTargetZCenter = 0.0, locTargetLength = 1.0;
	locGeometry->GetTargetZ(locTargetZCenter);
	locGeometry->GetTargetLength(locTargetLength);
	if(!(locTargetLength > 0.0)) //broken / invalid (e.g. commissioning data)
		locTargetLength = 30.0; //i dunno
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

		if(!locIsRESTEvent)
			Create_MatchingHists(false);
		Create_MatchingHists(true);

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
	locHistName = "DeltaPOverPVsP";
	locHistTitle = locHistParticleName + string(";p (GeV/c);#Deltap/p (Measured - Missing)");
	dHistMap_Resolution_DeltaPOverPVsP[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

	// DeltaP/P Vs Theta
	locHistName = string("DeltaPOverPVsTheta");
	locHistTitle = locHistParticleName + string(";#theta#circ;#Deltap/p (Measured - Missing)");
	dHistMap_Resolution_DeltaPOverPVsTheta[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DDeltaPOverPBins, dMinDeltaPOverP, dMaxDeltaPOverP);

	// DeltaTheta Vs P
	locHistName = string("DeltaThetaVsP");
	locHistTitle = locHistParticleName + string(";p (GeV/c);#Delta#theta#circ (Measured - Missing)");
	dHistMap_Resolution_DeltaThetaVsP[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltaThetaBins, dMinDeltaTheta, dMaxDeltaTheta);

	// DeltaTheta Vs Theta
	locHistName = string("DeltaThetaVsTheta");
	locHistTitle = locHistParticleName + string(";#theta#circ;#Delta#theta#circ (Measured - Missing)");
	dHistMap_Resolution_DeltaThetaVsTheta[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DDeltaThetaBins, dMinDeltaTheta, dMaxDeltaTheta);

	// DeltaPhi Vs P
	locHistName = string("DeltaPhiVsP");
	locHistTitle = locHistParticleName + string(";p (GeV/c);#Delta#phi#circ (Measured - Missing)");
	dHistMap_Resolution_DeltaPhiVsP[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

	// DeltaPhi Vs Theta
	locHistName = string("DeltaPhiVsTheta");
	locHistTitle = locHistParticleName + string(";#theta#circ;#Delta#phi#circ (Measured - Missing)");
	dHistMap_Resolution_DeltaPhiVsTheta[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

	gDirectory->cd(".."); //return to the action directory
}

void DCustomAction_TrackingEfficiency::Create_MatchingHists(bool locIsTimeBasedFlag)
{
	string locHistName, locHistTitle;

	string locDirectoryName = locIsTimeBasedFlag ? "Matching_TimeBased" : "Matching_WireBased";
	CreateAndChangeTo_Directory(locDirectoryName.c_str(), locDirectoryName.c_str());
	string locHistParticleName = locIsTimeBasedFlag ? "Time-Based " : "Wire-Based ";
	locHistParticleName += ParticleName_ROOT(dMissingPID);

	//Kinematics of has (no) hit
	vector<DetectorSystem_t> locDetectorSystems;
	locDetectorSystems.push_back(SYS_START);  locDetectorSystems.push_back(SYS_BCAL);
	locDetectorSystems.push_back(SYS_TOF);  locDetectorSystems.push_back(SYS_FCAL);
	for(size_t loc_i = 0; loc_i < locDetectorSystems.size(); ++loc_i)
	{
		DetectorSystem_t locSystem = locDetectorSystems[loc_i];

		double locMaxTheta = ((locSystem == SYS_FCAL) || (locSystem == SYS_TOF)) ? 12.0 : dMaxTheta;
		double locMaxP = (locSystem == SYS_BCAL) ? 3.0 : dMaxP;

		string locSystemName = SystemName(locSystem);
		if(locSystemName == "ST")
			locSystemName = "SC";
		string locDirName = locSystemName;
		if(locSystemName == "TOF")
			locDirName = "TOFPoint";

		CreateAndChangeTo_Directory(locDirName, locDirName);

		// PVsTheta Has Hit
		locHistName = "PVsTheta_HasHit";
		locHistTitle = locHistParticleName + string(", ") + locSystemName + string(" Has Hit;#theta#circ;p (GeV/c)");
		dHistMap_PVsTheta_HasHit[locSystem][locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, locMaxTheta, dNum2DPBins, dMinP, locMaxP);

		// PVsTheta Has No Hit
		locHistName = "PVsTheta_NoHit";
		locHistTitle = locHistParticleName + string(", ") + locSystemName + string(" No Hit;#theta#circ;p (GeV/c)");
		dHistMap_PVsTheta_NoHit[locSystem][locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, locMaxTheta, dNum2DPBins, dMinP, locMaxP);

		// PhiVsTheta Has Hit
		locHistName = "PhiVsTheta_HasHit";
		locHistTitle = locHistParticleName + string(", ") + locSystemName + string(" Has Hit;#theta#circ;#phi#circ");
		dHistMap_PhiVsTheta_HasHit[locSystem][locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, locMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

		// PhiVsTheta Has No Hit
		locHistName = "PhiVsTheta_NoHit";
		locHistTitle = locHistParticleName + string(", ") + locSystemName + string(" No Hit;#theta#circ;#phi#circ");
		dHistMap_PhiVsTheta_NoHit[locSystem][locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, locMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

		gDirectory->cd("..");
	}

	//SC
	CreateAndChangeTo_Directory("SC", "SC");
	locHistName = "SCPaddleVsTheta_HasHit";
	locHistTitle = locHistParticleName + string(", SC Has Hit;#theta#circ;Projected SC Paddle");
	dHistMap_SCPaddleVsTheta_HasHit[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, 30, 0.5, 30.5);

	locHistName = "SCPaddleVsTheta_NoHit";
	locHistTitle = locHistParticleName + string(", SC No Hit;#theta#circ;Projected SC Paddle");
	dHistMap_SCPaddleVsTheta_NoHit[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, 30, 0.5, 30.5);

	locHistName = "SCTrackDeltaPhiVsP";
	locHistTitle = locHistParticleName + string(";p (GeV/c);SC / Track #Delta#phi#circ");
	dHistMap_SCTrackDeltaPhiVsP[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltaPhiBins, dSCMatchMinDeltaPhi, dSCMatchMaxDeltaPhi);

	locHistName = "SCTrackDeltaPhiVsTheta";
	locHistTitle = locHistParticleName + string(";#theta#circ;SC / Track #Delta#phi#circ");
	dHistMap_SCTrackDeltaPhiVsTheta[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DDeltaPhiBins, dSCMatchMinDeltaPhi, dSCMatchMaxDeltaPhi);
	gDirectory->cd("..");

	//TOFPaddle
	CreateAndChangeTo_Directory("TOFPaddle", "TOFPaddle");
	locHistName = "VerticalPaddleTrackDeltaX";
	locHistTitle = locHistParticleName + string(";Vertical TOF Paddle / Track |#DeltaX| (cm)");
	dHistMap_TOFPaddleTrackDeltaX[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumTrackDOCABins/2, 0.0, dMaxTrackMatchDOCA);

	locHistName = "HorizontalPaddleTrackDeltaY";
	locHistTitle = locHistParticleName + string(";Horizontal TOF Paddle / Track |#DeltaY| (cm)");
	dHistMap_TOFPaddleTrackDeltaY[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumTrackDOCABins/2, 0.0, dMaxTrackMatchDOCA);

	locHistName = "TrackYVsVerticalPaddle_HasHit";
	locHistTitle = locHistParticleName + string(", TOF Paddle Has Hit;Projected Vertical Paddle;Projected TOF Hit Y (cm)");
	dHistMap_TOFPaddleTrackYVsVerticalPaddle_HasHit[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, 44, 0.5, 44.5, dNumFCALTOFXYBins, -130.0, 130.0);

	locHistName = "TrackYVsVerticalPaddle_NoHit";
	locHistTitle = locHistParticleName + string(", TOF Paddle No Hit;Projected Vertical Paddle;Projected TOF Hit Y (cm)");
	dHistMap_TOFPaddleTrackYVsVerticalPaddle_NoHit[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, 44, 0.5, 44.5, dNumFCALTOFXYBins, -130.0, 130.0);

	locHistName = "HorizontalPaddleVsTrackX_HasHit";
	locHistTitle = locHistParticleName + string(", TOF Paddle Has Hit;Projected TOF Hit X (cm);Projected Horizontal Paddle");
	dHistMap_TOFPaddleHorizontalPaddleVsTrackX_HasHit[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNumFCALTOFXYBins, -130.0, 130.0, 44, 0.5, 44.5);

	locHistName = "HorizontalPaddleVsTrackX_NoHit";
	locHistTitle = locHistParticleName + string(", TOF Paddle No Hit;Projected TOF Hit X (cm);Projected Horizontal Paddle");
	dHistMap_TOFPaddleHorizontalPaddleVsTrackX_NoHit[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNumFCALTOFXYBins, -130.0, 130.0, 44, 0.5, 44.5);
	gDirectory->cd("..");

	//TOFPoint
	CreateAndChangeTo_Directory("TOFPoint", "TOFPoint");
	locHistName = "TrackTOFYVsX_HasHit";
	locHistTitle = locHistParticleName + string(", TOF Has Hit;Projected TOF Hit X (cm);Projected TOF Hit Y (cm)");
	dHistMap_TrackTOFYVsX_HasHit[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNumFCALTOFXYBins, -130.0, 130.0, dNumFCALTOFXYBins, -130.0, 130.0);

	locHistName = "TrackTOFYVsX_NoHit";
	locHistTitle = locHistParticleName + string(", TOF No Hit;Projected TOF Hit X (cm);Projected TOF Hit Y (cm)");
	dHistMap_TrackTOFYVsX_NoHit[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNumFCALTOFXYBins, -130.0, 130.0, dNumFCALTOFXYBins, -130.0, 130.0);

	locHistName = "TrackTOF2DPaddles_HasHit";
	locHistTitle = locHistParticleName + string(", TOF Has Hit;Projected Vertical TOF Paddle;Projected Horizontal TOF Paddle");
	dHistMap_TrackTOF2DPaddles_HasHit[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, 44, 0.5, 44.5, 44, 0.5, 44.5);

	locHistName = "TrackTOF2DPaddles_NoHit";
	locHistTitle = locHistParticleName + string(", TOF No Hit;Projected Vertical TOF Paddle;Projected Horizontal TOF Paddle");
	dHistMap_TrackTOF2DPaddles_NoHit[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, 44, 0.5, 44.5, 44, 0.5, 44.5);

	locHistName = "TOFTrackDistanceVsP";
	locHistTitle = locHistParticleName + string(";p (GeV/c);TOF / Track Distance (cm)");
	dHistMap_TOFPointTrackDistanceVsP[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);

	locHistName = "TOFTrackDistanceVsTheta";
	locHistTitle = locHistParticleName + string(";#theta#circ;TOF / Track Distance (cm)");
	dHistMap_TOFPointTrackDistanceVsTheta[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, 20.0, dNum2DTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);

	locHistName = "TOFTrackDeltaXVsHorizontalPaddle";
	locHistTitle = locHistParticleName + string(";TOF Horizontal Paddle;TOF / Track #DeltaX (cm)");
	dHistMap_TOFPointTrackDeltaXVsHorizontalPaddle[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, 44, 0.5, 44.5, dNum2DTrackDOCABins, -1.0*dMaxTrackMatchDOCA, dMaxTrackMatchDOCA);

	locHistName = "TOFTrackDeltaXVsVerticalPaddle";
	locHistTitle = locHistParticleName + string(";TOF Vertical Paddle;TOF / Track #DeltaX (cm)");
	dHistMap_TOFPointTrackDeltaXVsVerticalPaddle[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, 44, 0.5, 44.5, dNum2DTrackDOCABins, -1.0*dMaxTrackMatchDOCA, dMaxTrackMatchDOCA);

	locHistName = "TOFTrackDeltaYVsHorizontalPaddle";
	locHistTitle = locHistParticleName + string(";TOF Horizontal Paddle;TOF / Track #DeltaY (cm)");
	dHistMap_TOFPointTrackDeltaYVsHorizontalPaddle[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, 44, 0.5, 44.5, dNum2DTrackDOCABins, -1.0*dMaxTrackMatchDOCA, dMaxTrackMatchDOCA);

	locHistName = "TOFTrackDeltaYVsVerticalPaddle";
	locHistTitle = locHistParticleName + string(";TOF Vertical Paddle;TOF / Track #DeltaY (cm)");
	dHistMap_TOFPointTrackDeltaYVsVerticalPaddle[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, 44, 0.5, 44.5, dNum2DTrackDOCABins, -1.0*dMaxTrackMatchDOCA, dMaxTrackMatchDOCA);

	locHistName = "TOFTrackDistance_BothPlanes";
	locHistTitle = locHistParticleName + string("TOF Hit in Both Planes;TOF / Track Distance (cm)");
	dHistMap_TOFPointTrackDistance_BothPlanes[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);

	locHistName = "TOFTrackDistance_OnePlane";
	locHistTitle = locHistParticleName + string("TOF Hit in One Plane;TOF / Track Distance (cm)");
	dHistMap_TOFPointTrackDistance_OnePlane[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);
	gDirectory->cd("..");

	//FCAL
	CreateAndChangeTo_Directory("FCAL", "FCAL");
	locHistName = "TrackFCALYVsX_HasHit";
	locHistTitle = locHistParticleName + string(", FCAL Has Hit;Projected FCAL Hit X (cm);Projected FCAL Hit Y (cm)");
	dHistMap_TrackFCALYVsX_HasHit[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNumFCALTOFXYBins, -130.0, 130.0, dNumFCALTOFXYBins, -130.0, 130.0);

	locHistName = "TrackFCALYVsX_NoHit";
	locHistTitle = locHistParticleName + string(", FCAL No Hit;Projected FCAL Hit X (cm);Projected FCAL Hit Y (cm)");
	dHistMap_TrackFCALYVsX_NoHit[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNumFCALTOFXYBins, -130.0, 130.0, dNumFCALTOFXYBins, -130.0, 130.0);

	locHistName = "TrackFCALRowVsColumn_HasHit";
	locHistTitle = locHistParticleName + string(", FCAL Has Hit;Projected FCAL Hit Column;Projected FCAL Hit Row");
	dHistMap_TrackFCALRowVsColumn_HasHit[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, 59, -0.5, 58.5, 59, -0.5, 58.5);

	locHistName = "TrackFCALRowVsColumn_NoHit";
	locHistTitle = locHistParticleName + string(", FCAL No Hit;Projected FCAL Hit Column;Projected FCAL Hit Row");
	dHistMap_TrackFCALRowVsColumn_NoHit[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, 59, -0.5, 58.5, 59, -0.5, 58.5);

	locHistName = "FCALTrackDistanceVsP";
	locHistTitle = locHistParticleName + string(";p (GeV/c);FCAL / Track Distance (cm)");
	dHistMap_FCALTrackDistanceVsP[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);

	locHistName = "FCALTrackDistanceVsTheta";
	locHistTitle = locHistParticleName + string(";#theta#circ;FCAL / Track Distance (cm)");
	dHistMap_FCALTrackDistanceVsTheta[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, 20.0, dNum2DTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);
	gDirectory->cd("..");

	//BCAL
	CreateAndChangeTo_Directory("BCAL", "BCAL");
	locHistName = "TrackBCALModuleVsZ_HasHit";
	locHistTitle = locHistParticleName + string(", BCAL Has Hit;Projected BCAL Hit Z (cm);Projected BCAL Hit Module");
	dHistMap_TrackBCALModuleVsZ_HasHit[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DBCALZBins, 0.0, 450.0, 48, 0.5, 48.5);

	locHistName = "TrackBCALModuleVsZ_NoHit";
	locHistTitle = locHistParticleName + string(", BCAL No Hit;Projected BCAL Hit Z (cm);Projected BCAL Hit Module");
	dHistMap_TrackBCALModuleVsZ_NoHit[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DBCALZBins, 0.0, 450.0, 48, 0.5, 48.5);

	locHistName = "TrackBCALPhiVsZ_HasHit";
	locHistTitle = locHistParticleName + string(", BCAL Has Hit;Projected BCAL Hit Z (cm);Projected BCAL Hit #phi#circ");
	dHistMap_TrackBCALPhiVsZ_HasHit[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DBCALZBins, 0.0, 450.0, dNum2DPhiBins, dMinPhi, dMaxPhi);

	locHistName = "TrackBCALPhiVsZ_NoHit";
	locHistTitle = locHistParticleName + string(", BCAL No Hit;Projected BCAL Hit Z (cm);Projected BCAL Hit #phi#circ");
	dHistMap_TrackBCALPhiVsZ_NoHit[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DBCALZBins, 0.0, 450.0, dNum2DPhiBins, dMinPhi, dMaxPhi);

	locHistName = "BCALDeltaPhiVsP";
	locHistTitle = locHistParticleName + string(";p (GeV/c);BCAL / Track #Delta#phi#circ");
	dHistMap_BCALDeltaPhiVsP[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, 4.0, dNum2DDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

	locHistName = "BCALDeltaZVsTheta";
	locHistTitle = locHistParticleName + string(";#theta#circ;BCAL / Track #Deltaz (cm)");
	dHistMap_BCALDeltaZVsTheta[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DDeltaZBins, dMinDeltaZ, dMaxDeltaZ);
	gDirectory->cd("..");

	gDirectory->cd("..");
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
	dHistMap_MatchingFOM[locIsTimeBasedFlag] = GetOrCreate_Histogram<TH1I>(locHistName, locHistTitle, dNumMatchFOMBins, 0, 1.0);

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
		dHistMap_TrackFound_PVsTheta[locIsTimeBasedFlag][locVertexZBin] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

		locHistName = string("PVsTheta_TrackMissing_VertexZBin") + locVertexZBinStream.str();
		locHistTitle = locHistParticleName + string(", ") + locVertexZRangeStream.str() + string(", Track Missing;#theta#circ;p (GeV/c)");
		dHistMap_TrackMissing_PVsTheta[locIsTimeBasedFlag][locVertexZBin] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

		locHistName = string("PhiVsTheta_TrackFound_VertexZBin") + locVertexZBinStream.str();
		locHistTitle = locHistParticleName + string(", ") + locVertexZRangeStream.str() + string(", Track Found;#theta#circ;#phi#circ");
		dHistMap_TrackFound_PhiVsTheta[locIsTimeBasedFlag][locVertexZBin] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

		locHistName = string("PhiVsTheta_TrackMissing_VertexZBin") + locVertexZBinStream.str();
		locHistTitle = locHistParticleName + string(", ") + locVertexZRangeStream.str() + string(", Track Missing;#theta#circ;#phi#circ");
		dHistMap_TrackMissing_PhiVsTheta[locIsTimeBasedFlag][locVertexZBin] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

		//Detector Match
		locHistName = string("PVsTheta_FoundHasDetectorMatch_VertexZBin") + locVertexZBinStream.str();
		locHistTitle = locHistParticleName + string(", ") + locVertexZRangeStream.str() + string(", Has Detector Match;#theta#circ;p (GeV/c)");
		dHistMap_FoundHasDetectorMatch_PVsTheta[locIsTimeBasedFlag][locVertexZBin] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

		locHistName = string("PVsTheta_FoundNoDetectorMatch_VertexZBin") + locVertexZBinStream.str();
		locHistTitle = locHistParticleName + string(", ") + locVertexZRangeStream.str() + string(", No Detector Match;#theta#circ;p (GeV/c)");
		dHistMap_FoundNoDetectorMatch_PVsTheta[locIsTimeBasedFlag][locVertexZBin] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPBins, dMinP, dMaxP);

		locHistName = string("PhiVsTheta_FoundHasDetectorMatch_VertexZBin") + locVertexZBinStream.str();
		locHistTitle = locHistParticleName + string(", ") + locVertexZRangeStream.str() + string(", Has Detector Match;#theta#circ;#phi#circ");
		dHistMap_FoundHasDetectorMatch_PhiVsTheta[locIsTimeBasedFlag][locVertexZBin] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);

		locHistName = string("PhiVsTheta_FoundNoDetectorMatch_VertexZBin") + locVertexZBinStream.str();
		locHistTitle = locHistParticleName + string(", ") + locVertexZRangeStream.str() + string(", No Detector Match;#theta#circ;#phi#circ");
		dHistMap_FoundNoDetectorMatch_PhiVsTheta[locIsTimeBasedFlag][locVertexZBin] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DPhiBins, dMinPhi, dMaxPhi);
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

		//dE/dx vs p
		locHistName = string("dEdXVsP_") + SystemName(locSystem);
		locHistTitle = locHistParticleName + string(";p (GeV/c);") + SystemName(locSystem) + string(" dE/dX ") + locUnitsString;
		dHistMap_dEdXVsP[locSystem] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DdEdxBins, dMindEdX, dMaxdEdX);

		//delta-dE/dx vs p
		locHistName = string("DeltadEdXVsP_") + SystemName(locSystem);
		locHistTitle = locHistParticleName + string(";p (GeV/c);") + SystemName(locSystem) + string(" #Delta(dE/dX) ") + locUnitsString;
		dHistMap_DeltadEdXVsP[locSystem] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltadEdxBins, dMinDeltadEdx, dMaxDeltadEdx);
	}

	//beta
	vector<DetectorSystem_t> locDetectors_BetaVsP;
	locDetectors_BetaVsP.push_back(SYS_BCAL);  locDetectors_BetaVsP.push_back(SYS_FCAL);  locDetectors_BetaVsP.push_back(SYS_TOF);
	for(size_t loc_i = 0; loc_i < locDetectors_BetaVsP.size(); ++loc_i)
	{
		DetectorSystem_t locSystem = locDetectors_BetaVsP[loc_i];
		double locMaxP = (locSystem == SYS_BCAL) ? dMaxPBCAL : dMaxP;

		//beta vs p
		locHistName = string("BetaVsP_") + SystemName(locSystem);
		locHistTitle = locHistParticleName + string(", ") + SystemName(locSystem) + string(";p (GeV/c);#beta");
		dHistMap_BetaVsP[locSystem] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, locMaxP, dNum2DBetaBins, dMinBeta, dMaxBeta);

		//delta-beta vs p
		locHistName = string("DeltaBetaVsP_") + SystemName(locSystem);
		locHistTitle = locHistParticleName + string(", ") + SystemName(locSystem) + string(";p (GeV/c);#Delta#beta");
		dHistMap_DeltaBetaVsP[locSystem] = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, locMaxP, dNum2DDeltaBetaBins, dMinDeltaBeta, dMaxDeltaBeta);
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
	if((locVertexZBin < 0) || (locVertexZBin >= int(dNumVertexZBins)))
		return true; //not within range

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

	if(locBestTrackWireBased != NULL)
	{
		vector<const DDetectorMatches*> locMatchVector_WireBased;
		locEventLoop->Get(locMatchVector_WireBased, "WireBased");

		const DDetectorMatches* locDetectorMatches_WireBased = locMatchVector_WireBased.empty() ? NULL : locMatchVector_WireBased[0];
		bool locHasDetectorMatch_WireBased = (locBestTrackWireBased == NULL) ? false : locDetectorMatches_WireBased->Get_IsMatchedToHit(locBestTrackWireBased);
		Fill_ResolutionAndTrackEff_Hists(locBestTrackWireBased, locMissingP4, locVertexZBin, locBestWireBasedMatchFOM, locHasDetectorMatch_WireBased, false);

		//Matching
		if(locBestWireBasedMatchFOM >= dMinTrackMatchFOM)
			Fill_MatchingHists(locEventLoop, locBestTrackWireBased, false);
	}

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
	Fill_ResolutionAndTrackEff_Hists(locBestTrackTimeBased, locMissingP4, locVertexZBin, locBestTimeBasedMatchFOM, locHasDetectorMatch_TimeBased, true);

	//Matching
	if((locBestTrackTimeBased != NULL) && (locBestTimeBasedMatchFOM >= dMinTrackMatchFOM))
		Fill_MatchingHists(locEventLoop, locBestTrackTimeBased, true);

	/*********************************************** CHARGED HYPOTHESES, PID *********************************************/

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
		return true; //No charged hypothesis, don't bother //e.g. reconned track didn't pass pre-select cuts

	//if RF time is indeterminate, start time will be NaN
	const DBCALShowerMatchParams* locBCALShowerMatchParams = locBestChargedTrackHypothesis->Get_BCALShowerMatchParams();
	const DFCALShowerMatchParams* locFCALShowerMatchParams = locBestChargedTrackHypothesis->Get_FCALShowerMatchParams();
	const DTOFHitMatchParams* locTOFHitMatchParams = locBestChargedTrackHypothesis->Get_TOFHitMatchParams();
	const DSCHitMatchParams* locSCHitMatchParams = locBestChargedTrackHypothesis->Get_SCHitMatchParams();

	double locP = locBestChargedTrackHypothesis->momentum().Mag();

	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//SC
		if(locSCHitMatchParams != NULL)
		{
			dHistMap_dEdXVsP[SYS_START]->Fill(locP, locSCHitMatchParams->dEdx*1.0E3);
			double locdx = locSCHitMatchParams->dHitEnergy/locSCHitMatchParams->dEdx;
			double locProbabledEdx = 0.0, locSigmadEdx = 0.0;
			dParticleID->GetScintMPdEandSigma(locP, locBestChargedTrackHypothesis->mass(), locdx, locProbabledEdx, locSigmadEdx);
			dHistMap_DeltadEdXVsP[SYS_START]->Fill(locP, (locSCHitMatchParams->dEdx - locProbabledEdx)*1.0E3);
		}

		//TOF
		if(locTOFHitMatchParams != NULL)
		{
			//energy
			dHistMap_dEdXVsP[SYS_TOF]->Fill(locP, locTOFHitMatchParams->dEdx*1.0E3);
			double locdx = locTOFHitMatchParams->dHitEnergy/locTOFHitMatchParams->dEdx;
			double locProbabledEdx = 0.0, locSigmadEdx = 0.0;
			dParticleID->GetScintMPdEandSigma(locP, locBestChargedTrackHypothesis->mass(), locdx, locProbabledEdx, locSigmadEdx);
			dHistMap_DeltadEdXVsP[SYS_TOF]->Fill(locP, (locTOFHitMatchParams->dEdx - locProbabledEdx)*1.0E3);

			//timing
			double locBeta_Timing = locTOFHitMatchParams->dPathLength/(29.9792458*(locTOFHitMatchParams->dHitTime - locBestChargedTrackHypothesis->t0()));
			dHistMap_BetaVsP[SYS_TOF]->Fill(locP, locBeta_Timing);
			double locDeltaBeta = locBestChargedTrackHypothesis->lorentzMomentum().Beta() - locBeta_Timing;
			dHistMap_DeltaBetaVsP[SYS_TOF]->Fill(locP, locDeltaBeta);
		}

		//BCAL
		if(locBCALShowerMatchParams != NULL)
		{
			const DBCALShower* locBCALShower = locBCALShowerMatchParams->dBCALShower;
			double locBeta_Timing = locBCALShowerMatchParams->dPathLength/(29.9792458*(locBCALShower->t - locBestChargedTrackHypothesis->t0()));
			dHistMap_BetaVsP[SYS_BCAL]->Fill(locP, locBeta_Timing);
			double locDeltaBeta = locBestChargedTrackHypothesis->lorentzMomentum().Beta() - locBeta_Timing;
			dHistMap_DeltaBetaVsP[SYS_BCAL]->Fill(locP, locDeltaBeta);
		}

		//FCAL
		if(locFCALShowerMatchParams != NULL)
		{
			const DFCALShower* locFCALShower = locFCALShowerMatchParams->dFCALShower;
			double locBeta_Timing = locFCALShowerMatchParams->dPathLength/(29.9792458*(locFCALShower->getTime() - locBestChargedTrackHypothesis->t0()));
			dHistMap_BetaVsP[SYS_FCAL]->Fill(locP, locBeta_Timing);
			double locDeltaBeta = locBestChargedTrackHypothesis->lorentzMomentum().Beta() - locBeta_Timing;
			dHistMap_DeltaBetaVsP[SYS_FCAL]->Fill(locP, locDeltaBeta);
		}

		//CDC
		if(locBestTrackTimeBased->dNumHitsUsedFordEdx_CDC > 0)
		{
			dHistMap_dEdXVsP[SYS_CDC]->Fill(locP, locBestTrackTimeBased->ddEdx_CDC*1.0E6);
			double locProbabledEdx = dParticleID->GetMostProbabledEdx_DC(locP, locBestChargedTrackHypothesis->mass(), locBestTrackTimeBased->ddx_CDC, true);
			dHistMap_DeltadEdXVsP[SYS_CDC]->Fill(locP, (locBestTrackTimeBased->ddEdx_CDC - locProbabledEdx)*1.0E6);
		}

		//FDC
		if(locBestTrackTimeBased->dNumHitsUsedFordEdx_FDC > 0)
		{
			dHistMap_dEdXVsP[SYS_FDC]->Fill(locP, locBestTrackTimeBased->ddEdx_FDC*1.0E6);
			double locProbabledEdx = dParticleID->GetMostProbabledEdx_DC(locP, locBestChargedTrackHypothesis->mass(), locBestTrackTimeBased->ddx_FDC, false);
			dHistMap_DeltadEdXVsP[SYS_FDC]->Fill(locP, (locBestTrackTimeBased->ddEdx_FDC - locProbabledEdx)*1.0E6);
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}

void DCustomAction_TrackingEfficiency::Fill_ResolutionAndTrackEff_Hists(const DKinematicData* locTrack, DLorentzVector locMissingP4, size_t locVertexZBin, double locTrackMatchFOM, bool locHasDetectorMatch, bool locIsTimeBasedFlag)
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
	while(locDeltaPhi > 180.0)
		locDeltaPhi -= 360.0;
	while(locDeltaPhi < -180.0)
		locDeltaPhi += 360.0;

	bool locTrackFoundFlag = false;
//	if((locTrack != NULL) && (locTrackMatchFOM >= dMinTrackMatchFOM))
//		locTrackFoundFlag = true;
	if((locTrack != NULL) && (locDeltaPOverP < 0.2) && (locDeltaTheta < 10.0) && (locDeltaPhi < 15.0))
		locTrackFoundFlag = true;

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
		if(locTrackFoundFlag)
		{
			//Found
			dHistMap_TrackFound_PVsTheta[locIsTimeBasedFlag][locVertexZBin]->Fill(locMissingTheta, locMissingP);
			dHistMap_TrackFound_PhiVsTheta[locIsTimeBasedFlag][locVertexZBin]->Fill(locMissingTheta, locMissingPhi);
			if(locHasDetectorMatch)
			{
				//Detector Match
				dHistMap_FoundHasDetectorMatch_PVsTheta[locIsTimeBasedFlag][locVertexZBin]->Fill(locMissingTheta, locMissingP);
				dHistMap_FoundHasDetectorMatch_PhiVsTheta[locIsTimeBasedFlag][locVertexZBin]->Fill(locMissingTheta, locMissingPhi);
			}
			else
			{
				//No Detector Match
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

void DCustomAction_TrackingEfficiency::Fill_MatchingHists(JEventLoop* locEventLoop, const DKinematicData* locTrack, bool locIsTimeBasedFlag)
{
	vector<const DBCALShower*> locBCALShowers;
	locEventLoop->Get(locBCALShowers);

	vector<const DFCALShower*> locFCALShowers;
	locEventLoop->Get(locFCALShowers);

	vector<const DTOFPoint*> locTOFPoints;
	locEventLoop->Get(locTOFPoints);

	vector<const DTOFPaddleHit*> locTOFPaddleHits;
	locEventLoop->Get(locTOFPaddleHits);

	vector<const DSCHit*> locSCHits;
	locEventLoop->Get(locSCHits);

	string locDetectorMatchesTag = locIsTimeBasedFlag ? "" : "WireBased";
	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches, locDetectorMatchesTag.c_str());

	//TRACK / BCAL CLOSEST MATCHES
	double locBestBCALMatchDeltaPhi = 999.9, locBestBCALMatchDeltaZ = 999.9;
	const DBCALShower* locClosestBCALShower = dParticleID->Get_ClosestToTrack_BCAL(locTrack, locBCALShowers, locBestBCALMatchDeltaPhi, locBestBCALMatchDeltaZ);

	//TRACK / FCAL CLOSEST MATCHES
	double locBestFCALDistance = 999.0;
	const DFCALShower* locClosestFCALShower = dParticleID->Get_ClosestToTrack_FCAL(locTrack, locFCALShowers, locBestFCALDistance);

	//TRACK / TOF PADDLE CLOSEST MATCHES
	double locBestTOFPaddleDeltaX = 999.9, locBestTOFPaddleDeltaY = 999.9;
	pair<const DTOFPaddleHit*, const DTOFPaddleHit*> locClosestTOFPaddleHits = dParticleID->Get_ClosestToTrack_TOFPaddles(locTrack, locTOFPaddleHits, locBestTOFPaddleDeltaX, locBestTOFPaddleDeltaY);

	//TRACK / TOF POINT CLOSEST MATCHES
	double locBestTOFPointDeltaX = 999.9, locBestTOFPointDeltaY = 999.9;
	const DTOFPoint* locClosestTOFPoint = dParticleID->Get_ClosestToTrack_TOFPoint(locTrack, locTOFPoints, locBestTOFPointDeltaX, locBestTOFPointDeltaY);

	//TRACK / SC CLOSEST MATCHES
	double locBestSCDeltaPhi = 999.0;
	const DSCHit* locClosestSCHit = dParticleID->Get_ClosestToTrack_SC(locTrack, locSCHits, locBestSCDeltaPhi);

	//PROJECTED HIT POSITIONS
	unsigned int locProjectedSCPaddle = 0;
	DVector3 locTOFIntersection;
	unsigned int locHorizontalTOFBar = 0, locVerticalTOFBar = 0;
	DVector3 locFCALIntersection;
	unsigned int locFCALRow = 0, locFCALColumn = 0;
	DVector3 locBCALIntersection;
	unsigned int locBCALModule = 0, locBCALSector = 0;

	//Get reference trajectory
	const DReferenceTrajectory* locReferenceTrajectory = NULL;
	if(locIsTimeBasedFlag)
		locReferenceTrajectory = (static_cast<const DTrackTimeBased*>(locTrack))->rt;
	else
		locReferenceTrajectory = (static_cast<const DTrackWireBased*>(locTrack))->rt;

	//Project
	if(locReferenceTrajectory != NULL)
	{
		locProjectedSCPaddle = dParticleID->PredictSCSector(locReferenceTrajectory, 999.9);
		dParticleID->PredictTOFPaddles(locReferenceTrajectory, locHorizontalTOFBar, locVerticalTOFBar, &locTOFIntersection);
		dParticleID->PredictFCALHit(locReferenceTrajectory, locFCALRow, locFCALColumn, &locFCALIntersection);
		dParticleID->PredictBCALWedge(locReferenceTrajectory, locBCALModule, locBCALSector, &locBCALIntersection);
	}

	//Fill Histograms
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		double locP = locTrack->momentum().Mag();
		double locTheta = locTrack->momentum().Theta()*180.0/TMath::Pi();
		double locPhi = locTrack->momentum().Phi()*180.0/TMath::Pi();

		/********************************************************** MATCHING DISTANCE **********************************************************/

		//BCAL
		if(locClosestBCALShower != NULL)
		{
			dHistMap_BCALDeltaPhiVsP[locIsTimeBasedFlag]->Fill(locP, locBestBCALMatchDeltaPhi*180.0/TMath::Pi());
			dHistMap_BCALDeltaZVsTheta[locIsTimeBasedFlag]->Fill(locTheta, locBestBCALMatchDeltaZ);
		}

		//FCAL
		if(locClosestFCALShower != NULL)
		{
			dHistMap_FCALTrackDistanceVsP[locIsTimeBasedFlag]->Fill(locP, locBestFCALDistance);
			dHistMap_FCALTrackDistanceVsTheta[locIsTimeBasedFlag]->Fill(locTheta, locBestFCALDistance);
		}

		//TOF Paddle
		if(locClosestTOFPaddleHits.second != NULL) //Horizontal
			dHistMap_TOFPaddleTrackDeltaY[locIsTimeBasedFlag]->Fill(locBestTOFPaddleDeltaY);
		if(locClosestTOFPaddleHits.first != NULL) //Vertical
			dHistMap_TOFPaddleTrackDeltaX[locIsTimeBasedFlag]->Fill(locBestTOFPaddleDeltaX);

		//TOF Point
		if(locClosestTOFPoint != NULL)
		{
			double locDistance = sqrt(locBestTOFPointDeltaX*locBestTOFPointDeltaX + locBestTOFPointDeltaY*locBestTOFPointDeltaY);
			if((locBestTOFPointDeltaX < 500.0) && (locBestTOFPointDeltaY < 500.0)) //else position not well-defined
			{
				dHistMap_TOFPointTrackDistanceVsP[locIsTimeBasedFlag]->Fill(locP, locDistance);
				dHistMap_TOFPointTrackDistanceVsTheta[locIsTimeBasedFlag]->Fill(locTheta, locDistance);
				if((locClosestTOFPoint->dHorizontalBar != 0) && (locClosestTOFPoint->dVerticalBar != 0))
					dHistMap_TOFPointTrackDistance_BothPlanes[locIsTimeBasedFlag]->Fill(locDistance);
				else
					dHistMap_TOFPointTrackDistance_OnePlane[locIsTimeBasedFlag]->Fill(locDistance);
			}

			dHistMap_TOFPointTrackDeltaXVsHorizontalPaddle[locIsTimeBasedFlag]->Fill(locClosestTOFPoint->dHorizontalBar, locBestTOFPointDeltaX);
			dHistMap_TOFPointTrackDeltaXVsVerticalPaddle[locIsTimeBasedFlag]->Fill(locClosestTOFPoint->dVerticalBar, locBestTOFPointDeltaX);

			dHistMap_TOFPointTrackDeltaYVsHorizontalPaddle[locIsTimeBasedFlag]->Fill(locClosestTOFPoint->dHorizontalBar, locBestTOFPointDeltaY);
			dHistMap_TOFPointTrackDeltaYVsVerticalPaddle[locIsTimeBasedFlag]->Fill(locClosestTOFPoint->dVerticalBar, locBestTOFPointDeltaY);
		}

		//SC
		if((locSCHits.size() <= 4) && (locClosestSCHit != NULL)) //don't fill if every paddle fired!
		{
			double locDeltaPhi = locBestSCDeltaPhi*180.0/TMath::Pi();
			dHistMap_SCTrackDeltaPhiVsP[locIsTimeBasedFlag]->Fill(locP, locDeltaPhi);
			dHistMap_SCTrackDeltaPhiVsTheta[locIsTimeBasedFlag]->Fill(locTheta, locDeltaPhi);
		}

		/********************************************************* MATCHING EFFICINECY *********************************************************/

		//BCAL
		if(locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_BCAL))
		{
			dHistMap_PVsTheta_HasHit[SYS_BCAL][locIsTimeBasedFlag]->Fill(locTheta, locP);
			dHistMap_PhiVsTheta_HasHit[SYS_BCAL][locIsTimeBasedFlag]->Fill(locTheta, locPhi);
			if(locBCALModule != 0)
			{
				dHistMap_TrackBCALPhiVsZ_HasHit[locIsTimeBasedFlag]->Fill(locBCALIntersection.Z(), locBCALIntersection.Phi()*180.0/TMath::Pi());
				dHistMap_TrackBCALModuleVsZ_HasHit[locIsTimeBasedFlag]->Fill(locBCALIntersection.Z(), locBCALModule);
			}
		}
		else
		{
			dHistMap_PVsTheta_NoHit[SYS_BCAL][locIsTimeBasedFlag]->Fill(locTheta, locP);
			dHistMap_PhiVsTheta_NoHit[SYS_BCAL][locIsTimeBasedFlag]->Fill(locTheta, locPhi);
			if(locBCALModule != 0)
			{
				dHistMap_TrackBCALPhiVsZ_NoHit[locIsTimeBasedFlag]->Fill(locBCALIntersection.Z(), locBCALIntersection.Phi()*180.0/TMath::Pi());
				dHistMap_TrackBCALModuleVsZ_NoHit[locIsTimeBasedFlag]->Fill(locBCALIntersection.Z(), locBCALModule);
			}
		}

		//FCAL
		if(locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_FCAL))
		{
			dHistMap_PVsTheta_HasHit[SYS_FCAL][locIsTimeBasedFlag]->Fill(locTheta, locP);
			dHistMap_PhiVsTheta_HasHit[SYS_FCAL][locIsTimeBasedFlag]->Fill(locTheta, locPhi);
			if(locFCALColumn != 0)
			{
				dHistMap_TrackFCALYVsX_HasHit[locIsTimeBasedFlag]->Fill(locFCALIntersection.X(), locFCALIntersection.Y());
				dHistMap_TrackFCALRowVsColumn_HasHit[locIsTimeBasedFlag]->Fill(locFCALColumn, locFCALRow);
			}
		}
		else
		{
			dHistMap_PVsTheta_NoHit[SYS_FCAL][locIsTimeBasedFlag]->Fill(locTheta, locP);
			dHistMap_PhiVsTheta_NoHit[SYS_FCAL][locIsTimeBasedFlag]->Fill(locTheta, locPhi);
			if(locFCALColumn != 0)
			{
				dHistMap_TrackFCALYVsX_NoHit[locIsTimeBasedFlag]->Fill(locFCALIntersection.X(), locFCALIntersection.Y());
				dHistMap_TrackFCALRowVsColumn_NoHit[locIsTimeBasedFlag]->Fill(locFCALColumn, locFCALRow);
			}
		}

		//TOF Paddle
		if(locHorizontalTOFBar != 0)
		{
			//Horizontal
			if(locClosestTOFPaddleHits.second != NULL)
			{
				if(locBestTOFPaddleDeltaY <= dMinTOFPaddleMatchDistance) //match
					dHistMap_TOFPaddleHorizontalPaddleVsTrackX_HasHit[locIsTimeBasedFlag]->Fill(locTOFIntersection.X(), locHorizontalTOFBar);
				else //no match
					dHistMap_TOFPaddleHorizontalPaddleVsTrackX_NoHit[locIsTimeBasedFlag]->Fill(locTOFIntersection.X(), locHorizontalTOFBar);
			}
			else // no match
				dHistMap_TOFPaddleHorizontalPaddleVsTrackX_NoHit[locIsTimeBasedFlag]->Fill(locTOFIntersection.X(), locHorizontalTOFBar);
		}
		if(locVerticalTOFBar != 0)
		{
			//Vertical
			if(locClosestTOFPaddleHits.first != NULL)
			{
				if(locBestTOFPaddleDeltaX <= dMinTOFPaddleMatchDistance) //match
					dHistMap_TOFPaddleTrackYVsVerticalPaddle_HasHit[locIsTimeBasedFlag]->Fill(locVerticalTOFBar, locTOFIntersection.Y());
				else //no match
					dHistMap_TOFPaddleTrackYVsVerticalPaddle_NoHit[locIsTimeBasedFlag]->Fill(locVerticalTOFBar, locTOFIntersection.Y());
			}
			else // no match
				dHistMap_TOFPaddleTrackYVsVerticalPaddle_NoHit[locIsTimeBasedFlag]->Fill(locVerticalTOFBar, locTOFIntersection.Y());
		}

		//TOF Point
		if(locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_TOF))
		{
			dHistMap_PVsTheta_HasHit[SYS_TOF][locIsTimeBasedFlag]->Fill(locTheta, locP);
			dHistMap_PhiVsTheta_HasHit[SYS_TOF][locIsTimeBasedFlag]->Fill(locTheta, locPhi);
			if(locHorizontalTOFBar != 0)
			{
				dHistMap_TrackTOFYVsX_HasHit[locIsTimeBasedFlag]->Fill(locTOFIntersection.X(), locTOFIntersection.Y());
				dHistMap_TrackTOF2DPaddles_HasHit[locIsTimeBasedFlag]->Fill(locVerticalTOFBar, locHorizontalTOFBar);
			}
		}
		else
		{
			dHistMap_PVsTheta_NoHit[SYS_TOF][locIsTimeBasedFlag]->Fill(locTheta, locP);
			dHistMap_PhiVsTheta_NoHit[SYS_TOF][locIsTimeBasedFlag]->Fill(locTheta, locPhi);
			if(locHorizontalTOFBar != 0)
			{
				dHistMap_TrackTOFYVsX_NoHit[locIsTimeBasedFlag]->Fill(locTOFIntersection.X(), locTOFIntersection.Y());
				dHistMap_TrackTOF2DPaddles_NoHit[locIsTimeBasedFlag]->Fill(locVerticalTOFBar, locHorizontalTOFBar);
			}
		}

		//SC
		if(locDetectorMatches->Get_IsMatchedToDetector(locTrack, SYS_START))
		{
			dHistMap_PVsTheta_HasHit[SYS_START][locIsTimeBasedFlag]->Fill(locTheta, locP);
			dHistMap_PhiVsTheta_HasHit[SYS_START][locIsTimeBasedFlag]->Fill(locTheta, locPhi);
			if(locProjectedSCPaddle != 0)
				dHistMap_SCPaddleVsTheta_HasHit[locIsTimeBasedFlag]->Fill(locTheta, locProjectedSCPaddle);
		}
		else
		{
			dHistMap_PVsTheta_NoHit[SYS_START][locIsTimeBasedFlag]->Fill(locTheta, locP);
			dHistMap_PhiVsTheta_NoHit[SYS_START][locIsTimeBasedFlag]->Fill(locTheta, locPhi);
			if(locProjectedSCPaddle != 0)
				dHistMap_SCPaddleVsTheta_NoHit[locIsTimeBasedFlag]->Fill(locTheta, locProjectedSCPaddle);
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

