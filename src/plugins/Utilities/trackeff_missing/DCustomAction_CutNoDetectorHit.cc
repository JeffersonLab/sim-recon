// $Id$
//
//    File: DCustomAction_CutNoDetectorHit.cc
// Created: Mon Feb 20 16:01:16 EST 2017
// Creator: pmatt (on Linux pmattdesktop.jlab.org 2.6.32-642.13.1.el6.x86_64 x86_64)
//

#include "DCustomAction_CutNoDetectorHit.h"

void DCustomAction_CutNoDetectorHit::Initialize(JEventLoop* locEventLoop)
{
	//Optional: Create histograms and/or modify member variables.
	//Create any histograms/trees/etc. within a ROOT lock. 
		//This is so that when running multithreaded, only one thread is writing to the ROOT file at a time. 
	//NEVER: Get anything from the JEventLoop while in a lock: May deadlock

	const DReaction* locReaction = Get_Reaction();
	if(!locReaction->Get_MissingPID(dMissingPID))
		return; //invalid reaction setup

	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	dMagneticFieldMap = locApplication->GetBfield(locEventLoop->GetJEvent().GetRunNumber());
	locEventLoop->GetSingle(dParticleID);

	string locHistName, locHistTitle;
	string locTrackString = string("Missing ") + ParticleName_ROOT(dMissingPID);

	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//Required: Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this action.
			//If another thread has already created the folder, it just changes to it. 
		CreateAndChangeTo_ActionDirectory();

		//SC MATCHING
		locHistName = "SCTrackDeltaPhiVsP";
		locHistTitle = locTrackString + string(";p (GeV/c);SC / Track #Delta#phi#circ");
		dHist_SCTrackDeltaPhiVsP = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DDeltaPhiBins, dSCMatchMinDeltaPhi, dSCMatchMaxDeltaPhi);

		locHistName = "SCTrackDeltaPhiVsTheta";
		locHistTitle = locTrackString + string(";#theta#circ;SC / Track #Delta#phi#circ");
		dHist_SCTrackDeltaPhiVsTheta = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DDeltaPhiBins, dSCMatchMinDeltaPhi, dSCMatchMaxDeltaPhi);

		locHistName = "SCTrackDeltaPhiVsZ";
		locHistTitle = locTrackString + string(";Projected SC Hit-Z (cm);SC / Track #Delta#phi#circ");
		dHistMap_SCTrackDeltaPhiVsZ = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DSCZBins, 0.0, 120.0, dNum2DDeltaPhiBins, dSCMatchMinDeltaPhi, dSCMatchMaxDeltaPhi);

		//TOF MATCHING
		locHistName = "TOFTrackDistanceVsP";
		locHistTitle = locTrackString + string(";p (GeV/c);TOF / Track Distance (cm)");
		dHist_TOFPointTrackDistanceVsP = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);

		locHistName = "TOFTrackDistanceVsTheta";
		locHistTitle = locTrackString + string(";#theta#circ;TOF / Track Distance (cm)");
		dHist_TOFPointTrackDistanceVsTheta = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, 20.0, dNum2DTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);

		//FCAL MATCHING
		locHistName = "FCALTrackDistanceVsP";
		locHistTitle = locTrackString + string(";p (GeV/c);FCAL / Track Distance (cm)");
		dHist_FCALTrackDistanceVsP = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, dMaxP, dNum2DTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);

		locHistName = "FCALTrackDistanceVsTheta";
		locHistTitle = locTrackString + string(";#theta#circ;FCAL / Track Distance (cm)");
		dHist_FCALTrackDistanceVsTheta = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, 20.0, dNum2DTrackDOCABins, dMinTrackDOCA, dMaxTrackMatchDOCA);

		//BCAL MATCHING
		locHistName = "BCALDeltaPhiVsP";
		locHistTitle = locTrackString + string(";p (GeV/c);BCAL / Track #Delta#phi#circ");
		dHist_BCALDeltaPhiVsP = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DPBins, dMinP, 4.0, dNum2DDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

		locHistName = "BCALDeltaPhiVsZ";
		locHistTitle = locTrackString + string(";Projected BCAL Hit-Z (cm);BCAL / Track #Delta#phi#circ");
		dHistMap_BCALDeltaPhiVsZ = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DBCALZBins, 0.0, 450.0, dNum2DDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

		locHistName = "BCALDeltaZVsTheta";
		locHistTitle = locTrackString + string(";#theta#circ;BCAL / Track #Deltaz (cm)");
		dHist_BCALDeltaZVsTheta = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DDeltaZBins, dMinDeltaZ, dMaxDeltaZ);

		locHistName = "BCALDeltaZVsZ";
		locHistTitle = locTrackString + string(";Projected BCAL Hit-Z (cm);BCAL / Track #Deltaz (cm)");
		dHistMap_BCALDeltaZVsZ = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DBCALZBins, 0.0, 450.0, dNum2DDeltaZBins, dMinDeltaZ, dMaxDeltaZ);
}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DCustomAction_CutNoDetectorHit::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//Write custom code to perform an action on the INPUT DParticleCombo (DParticleCombo)
	//NEVER: Grab DParticleCombo or DAnalysisResults objects (of any tag!) from the JEventLoop within this function
	//NEVER: Grab objects that are created post-kinfit (e.g. DKinFitResults, etc.) from the JEventLoop if Get_UseKinFitResultsFlag() == false: CAN CAUSE INFINITE DEPENDENCY LOOP
	//NEVER: Get anything from the JEventLoop while in a lock: May deadlock

	if(dMissingPID == Unknown)
		return false; //invalid reaction setup
	if(ParticleCharge(dMissingPID) == 0)
		return false; //NOT SUPPORTED

	const DKinematicData* locMissingParticle = locParticleCombo->Get_MissingParticle(); //is NULL if no kinfit!!
	if(locMissingParticle == nullptr)
		return false;

	double locP = locMissingParticle->momentum().Mag();
	double locTheta = locMissingParticle->momentum().Theta()*180.0/TMath::Pi();

	// Generate & swim reference trajectory for this
	DReferenceTrajectory rt(dMagneticFieldMap);
	rt.SetMass(ParticleMass(dMissingPID));
	rt.q = ParticleCharge(dMissingPID);
	rt.Swim(locMissingParticle->position(), locMissingParticle->momentum(), rt.q);

	vector<const DBCALShower*> locBCALShowers;
	locEventLoop->Get(locBCALShowers);

	vector<const DFCALShower*> locFCALShowers;
	locEventLoop->Get(locFCALShowers);

	vector<const DTOFPoint*> locTOFPoints;
	locEventLoop->Get(locTOFPoints);

	vector<const DSCHit*> locSCHits;
	locEventLoop->Get(locSCHits);

	// MATCHING: BCAL
	DBCALShowerMatchParams locBestBCALMatchParams;
	double locStartTimeVariance = 0.0;
	DVector3 locProjPos, locProjMom;
	double locStartTime = locParticleCombo->Get_ParticleComboStep(0)->Get_SpacetimeVertex().T();
	bool locBCALHitFoundFlag = dParticleID->Get_ClosestToTrack(&rt, locBCALShowers, false, locStartTime, locBestBCALMatchParams, &locStartTimeVariance, &locProjPos, &locProjMom);
	double locBCALProjectedZ = locProjPos.Z();
	double locBCALDeltaPhi = locBestBCALMatchParams.dDeltaPhiToShower*180.0/TMath::Pi();

	// MATCHING: FCAL
	DFCALShowerMatchParams locBestFCALMatchParams;
	locStartTime = locParticleCombo->Get_ParticleComboStep(0)->Get_SpacetimeVertex().T();
	bool locFCALHitFoundFlag = dParticleID->Get_ClosestToTrack(&rt, locFCALShowers, false, locStartTime, locBestFCALMatchParams);

	// MATCHING: SC
	DSCHitMatchParams locBestSCMatchParams;
	locStartTime = locParticleCombo->Get_ParticleComboStep(0)->Get_SpacetimeVertex().T();
	bool locSCHitFoundFlag = dParticleID->Get_ClosestToTrack(&rt, locSCHits, true, false, locStartTime, locBestSCMatchParams, &locStartTimeVariance, &locProjPos, &locProjMom);
	double locSCProjectedZ = locProjPos.Z();
	double locSCDeltaPhi = locBestSCMatchParams.dDeltaPhiToHit*180.0/TMath::Pi();

	// MATCHING: TOF
	DTOFHitMatchParams locBestTOFMatchParams;
	locStartTime = locParticleCombo->Get_ParticleComboStep(0)->Get_SpacetimeVertex().T();
	bool locTOFHitFoundFlag = dParticleID->Get_ClosestToTrack(&rt, locTOFPoints, false, locStartTime, locBestTOFMatchParams);
	double locDeltaX = locBestTOFMatchParams.dDeltaXToHit;
	double locDeltaY = locBestTOFMatchParams.dDeltaYToHit;
	double locTOFDistance = sqrt(locDeltaX*locDeltaX + locDeltaY*locDeltaY);

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action(); //ACQUIRE ROOT LOCK!!
	{
		/********************************************************** MATCHING DISTANCE **********************************************************/

		//BCAL
		if(locBCALHitFoundFlag)
		{
			dHist_BCALDeltaPhiVsP->Fill(locP, locBCALDeltaPhi);
			dHist_BCALDeltaZVsTheta->Fill(locTheta, locBestBCALMatchParams.dDeltaZToShower);
			dHistMap_BCALDeltaPhiVsZ->Fill(locBCALProjectedZ, locBCALDeltaPhi);
			dHistMap_BCALDeltaZVsZ->Fill(locBCALProjectedZ, locBestBCALMatchParams.dDeltaZToShower);
		}

		//FCAL
		if(locFCALHitFoundFlag)
		{
			dHist_FCALTrackDistanceVsP->Fill(locP, locBestFCALMatchParams.dDOCAToShower);
			dHist_FCALTrackDistanceVsTheta->Fill(locTheta, locBestFCALMatchParams.dDOCAToShower);
		}

		//TOF Point
		if(locTOFHitFoundFlag)
		{
			dHist_TOFPointTrackDistanceVsP->Fill(locP, locTOFDistance);
			dHist_TOFPointTrackDistanceVsTheta->Fill(locTheta, locTOFDistance);
		}

		//SC
		if(locSCHitFoundFlag)
		{
			dHist_SCTrackDeltaPhiVsP->Fill(locP, locSCDeltaPhi);
			dHist_SCTrackDeltaPhiVsTheta->Fill(locTheta, locSCDeltaPhi);
			dHistMap_SCTrackDeltaPhiVsZ->Fill(locSCProjectedZ, locSCDeltaPhi);
		}
	}
	Unlock_Action(); //RELEASE ROOT LOCK!!

	//REQUIRE: ST hit + EITHER BCAL hit OR TOF && FCAL hits
	if(!locSCHitFoundFlag)
		return false;
	if((locSCProjectedZ < 76.5) && (fabs(locSCDeltaPhi) > 10.0))
		return false; //HARD-CODED: BAD!
	if((locSCProjectedZ >= 76.5) && (fabs(locSCDeltaPhi) > (10.0 + locSCProjectedZ - 76.5)))
		return false;

	//BCAL Hit
	if(locBCALHitFoundFlag)
		return ((fabs(locBestBCALMatchParams.dDeltaZToShower) < 7.0) && (fabs(locBCALDeltaPhi) < 10.0));

	//MUST HAVE FCAL & TOF
	if(!locTOFHitFoundFlag || !locFCALHitFoundFlag)
		return false;
	return ((locTOFDistance < 10.0) && (locBestFCALMatchParams.dDOCAToShower < 10.0));
}
