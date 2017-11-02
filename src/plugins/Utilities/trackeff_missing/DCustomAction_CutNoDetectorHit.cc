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

	auto locMissingPIDs = Get_Reaction()->Get_MissingPIDs();
	dMissingPID = (locMissingPIDs.size() == 1) ? locMissingPIDs[0] : Unknown;
	if(locMissingPIDs.size() != 1)
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

		locHistName = "BCALDeltaPhiVsTheta";
		locHistTitle = locTrackString + string(";#theta#circ;BCAL / Track #Delta#phi#circ");
		dHistMap_BCALDeltaPhiVsTheta = GetOrCreate_Histogram<TH2I>(locHistName, locHistTitle, dNum2DThetaBins, dMinTheta, dMaxTheta, dNum2DDeltaPhiBins, dMinDeltaPhi, dMaxDeltaPhi);

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

	auto locMissingParticles = locParticleCombo->Get_MissingParticles(Get_Reaction());
	if(locMissingParticles.empty())
		return false; //kinfit failed to converge
	auto locMissingParticle = locMissingParticles[0];
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
	shared_ptr<const DBCALShowerMatchParams> locBestBCALMatchParams;
	double locStartTimeVariance = 0.0;
	DVector3 locProjPos, locProjMom;
	double locStartTime = locParticleCombo->Get_ParticleComboStep(0)->Get_SpacetimeVertex().T();
	bool locBCALHitFoundFlag = dParticleID->Get_ClosestToTrack(&rt, locBCALShowers, false, locStartTime, locBestBCALMatchParams, &locStartTimeVariance, &locProjPos, &locProjMom);
	double locBCALProjectedZ = locProjPos.Z();
	double locBCALDeltaPhi = locBCALHitFoundFlag ? locBestBCALMatchParams->dDeltaPhiToShower*180.0/TMath::Pi() : 9.9E9;

	// MATCHING: FCAL
	shared_ptr<const DFCALShowerMatchParams> locBestFCALMatchParams;
	locStartTime = locParticleCombo->Get_ParticleComboStep(0)->Get_SpacetimeVertex().T();
	DVector3 locProjPos_FCALMissing, locProjMom_FCALMissing;
	bool locFCALHitFoundFlag = dParticleID->Get_ClosestToTrack(&rt, locFCALShowers, false, locStartTime, locBestFCALMatchParams, &locStartTimeVariance, &locProjPos_FCALMissing, &locProjMom_FCALMissing);

	// MATCHING: SC
	shared_ptr<const DSCHitMatchParams> locBestSCMatchParams;
	locStartTime = locParticleCombo->Get_ParticleComboStep(0)->Get_SpacetimeVertex().T();
	bool locSCHitFoundFlag = dParticleID->Get_ClosestToTrack(&rt, locSCHits, true, false, locStartTime, locBestSCMatchParams, &locStartTimeVariance, &locProjPos, &locProjMom);
	double locSCProjectedZ = locProjPos.Z();
	double locSCDeltaPhi = locSCHitFoundFlag ? locBestSCMatchParams->dDeltaPhiToHit*180.0/TMath::Pi() : 9.9E9;

	// MATCHING: TOF
	shared_ptr<const DTOFHitMatchParams> locBestTOFMatchParams;
	DVector3 locProjPos_TOFMissing, locProjMom_TOFMissing;
	locStartTime = locParticleCombo->Get_ParticleComboStep(0)->Get_SpacetimeVertex().T();
	bool locTOFHitFoundFlag = dParticleID->Get_ClosestToTrack(&rt, locTOFPoints, false, locStartTime, locBestTOFMatchParams, &locStartTimeVariance, &locProjPos_TOFMissing, &locProjMom_TOFMissing);
	double locTOFDistance = 999.9;
	if(locTOFHitFoundFlag)
	{
		double locDeltaX = locBestTOFMatchParams->dDeltaXToHit;
		double locDeltaY = locBestTOFMatchParams->dDeltaYToHit;
		if(locBestTOFMatchParams->dTOFPoint->Is_XPositionWellDefined() == locBestTOFMatchParams->dTOFPoint->Is_YPositionWellDefined())
			locTOFDistance = sqrt(locDeltaX*locDeltaX + locDeltaY*locDeltaY);
		else
			locTOFDistance = locBestTOFMatchParams->dTOFPoint->Is_XPositionWellDefined() ? locDeltaX : locDeltaY;
	}


	/************************************************* TIME-BASED TRACKS *************************************************/

/*
	const DAnalysisUtilities* dAnalysisUtilities = nullptr;
	locEventLoop->GetSingle(dAnalysisUtilities);


	TMatrixDSym locMissingCovarianceMatrix(3);
	const TMatrixFSym& locKinFitCovarianceMatrix = *(locMissingParticle->errorMatrix());
	for(unsigned int loc_q = 0; loc_q < 3; ++loc_q)
	{
		for(unsigned int loc_r = 0; loc_r < 3; ++loc_r)
			locMissingCovarianceMatrix(loc_q, loc_r) = locKinFitCovarianceMatrix(loc_q, loc_r);
	}
	

	//Get unused tracks
	vector<const DTrackTimeBased*> locUnusedTimeBasedTracks;
	dAnalysisUtilities->Get_UnusedTimeBasedTracks(locEventLoop, locParticleCombo, locUnusedTimeBasedTracks);

	//loop over unused tracks
	double locBestMatchFOM = -1.0;

	const DTrackTimeBased* locBestTimeBasedTrack = nullptr;
	for(size_t loc_i = 0; loc_i < locUnusedTimeBasedTracks.size(); ++loc_i)
	{
		const DTrackTimeBased* locTimeBasedTrack = locUnusedTimeBasedTracks[loc_i];
		if(locTimeBasedTrack->PID() != dMissingPID)
			continue; //only use tracking results with correct PID

		DVector3 locTimeBasedDP3 = locTimeBasedTrack->momentum();
		TVector3 locTimeBasedP3(locTimeBasedDP3.X(), locTimeBasedDP3.Y(), locTimeBasedDP3.Z());

		const TMatrixFSym& locCovarianceMatrix = *(locTimeBasedTrack->errorMatrix());
		TMatrixDSym locDCovarianceMatrix(3);
		for(unsigned int loc_j = 0; loc_j < 3; ++loc_j)
		{
			for(unsigned int loc_k = 0; loc_k < 3; ++loc_k)
				locDCovarianceMatrix(loc_j, loc_k) = locCovarianceMatrix(loc_j, loc_k);
		}
		locDCovarianceMatrix += locMissingCovarianceMatrix;

		//invert matrix
		TDecompLU locDecompLU(locDCovarianceMatrix);
		//check to make sure that the matrix is decomposable and has a non-zero determinant
		if(!locDecompLU.Decompose() || (fabs(locDCovarianceMatrix.Determinant()) < 1.0E-300))
			continue;

		locDCovarianceMatrix.Invert();
		DVector3 locDeltaP3 = locTimeBasedTrack->momentum() - locMissingParticle->momentum();

		DMatrix locDeltas(3, 1);
		locDeltas(0, 0) = locDeltaP3.Px();
		locDeltas(1, 0) = locDeltaP3.Py();
		locDeltas(2, 0) = locDeltaP3.Pz();
		double locChiSq = (locDCovarianceMatrix.SimilarityT(locDeltas))(0, 0);
		double locMatchFOM = TMath::Prob(locChiSq, 3);

		if(locMatchFOM < locBestMatchFOM)
			continue;
		locBestMatchFOM = locMatchFOM;
		locBestTimeBasedTrack = locTimeBasedTrack;
	}

	if(locTOFHitFoundFlag && (locBestTimeBasedTrack != nullptr))
	{
		DVector3 locBestProjPos, locBestProjMom;
		double locStartTimeVariance = 0.0;
		shared_ptr<const DTOFHitMatchParams> locReconTOFMatchParams;
		locStartTime = locParticleCombo->Get_ParticleComboStep(0)->Get_SpacetimeVertex().T();
		if(dParticleID->Get_ClosestToTrack(locBestTimeBasedTrack->rt, locTOFPoints, false, locStartTime, locReconTOFMatchParams, &locStartTimeVariance, &locBestProjPos, &locBestProjMom))
		{
			//NOW CHECK FOR FCAL HITS
			locStartTime = locParticleCombo->Get_ParticleComboStep(0)->Get_SpacetimeVertex().T();
			DVector3 locBestProjPosFCAL, locBestProjMomFCAL;
			double locStartTimeVarianceFCAL = 0.0;
			shared_ptr<const DFCALShowerMatchParams> locReconFCALMatchParams;
			if(locFCALHitFoundFlag && dParticleID->Get_ClosestToTrack(locBestTimeBasedTrack->rt, locFCALShowers, false, locStartTime, locReconFCALMatchParams, &locStartTimeVarianceFCAL, &locBestProjPosFCAL, &locBestProjMomFCAL))
			{
				cout << "projected to hit tof AND fcal:" << endl;
				auto locReconP3 = locBestTimeBasedTrack->momentum();
				auto locMissingP3 = locMissingParticle->momentum();
				cout << "recon p/theta/phi: " << locReconP3.Mag() << ", " << locReconP3.Theta()*180.0/TMath::Pi() << ", " << locReconP3.Phi()*180.0/TMath::Pi() << endl;
				cout << "missing p/theta/phi: " << locMissingP3.Mag() << ", " << locMissingP3.Theta()*180.0/TMath::Pi() << ", " << locMissingP3.Phi()*180.0/TMath::Pi() << endl;
				cout << "TOF recon proj position = " << locBestProjPos.X() << ", " << locBestProjPos.Y() << ", " << locBestProjPos.Z() << endl;
				cout << "TOF missing proj position = " << locProjPos_TOFMissing.X() << ", " << locProjPos_TOFMissing.Y() << ", " << locProjPos_TOFMissing.Z() << endl;
				cout << "TOF: nearest hit: delta x/y, hor/vert bars: " << locReconTOFMatchParams->dDeltaXToHit << ", " << locReconTOFMatchParams->dDeltaYToHit << ", " << locReconTOFMatchParams->dTOFPoint->dHorizontalBar << ", " << locReconTOFMatchParams->dTOFPoint->dVerticalBar << endl;
				cout << "TOF: missing proj'd delta x/y, distance = " << locBestTOFMatchParams->dDeltaXToHit << ", " << locBestTOFMatchParams->dDeltaYToHit << ", " << locTOFDistance << endl;
				cout << "FCAL recon proj position = " << locBestProjPosFCAL.X() << ", " << locBestProjPosFCAL.Y() << ", " << locBestProjPosFCAL.Z() << endl;
				cout << "FCAL missing proj position = " << locProjPos_FCALMissing.X() << ", " << locProjPos_FCALMissing.Y() << ", " << locProjPos_FCALMissing.Z() << endl;
				cout << "FCAL: recon/missing doca to shower: " << locReconFCALMatchParams->dDOCAToShower << ", " << locBestFCALMatchParams->dDOCAToShower << endl;
			}
		}

//		const DDetectorMatches* locDetectorMatches = nullptr;
//		locEventLoop->GetSingle(locDetectorMatches);
	}
*/

	//FILL HISTOGRAMS
	auto locKinFitResults = locParticleCombo->Get_KinFitResults();
	if((locKinFitResults != nullptr) && (locKinFitResults->Get_ConfidenceLevel() > 0.0001))
	{
		//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
		//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
		Lock_Action(); //ACQUIRE ROOT LOCK!!
		{
			/********************************************************** MATCHING DISTANCE **********************************************************/

			//BCAL
			if(locBCALHitFoundFlag)
			{
				dHist_BCALDeltaPhiVsP->Fill(locP, locBCALDeltaPhi);
				dHistMap_BCALDeltaPhiVsZ->Fill(locBCALProjectedZ, locBCALDeltaPhi);
				dHistMap_BCALDeltaPhiVsTheta->Fill(locTheta, locBCALDeltaPhi);
				dHist_BCALDeltaZVsTheta->Fill(locTheta, locBestBCALMatchParams->dDeltaZToShower);
				dHistMap_BCALDeltaZVsZ->Fill(locBCALProjectedZ, locBestBCALMatchParams->dDeltaZToShower);
			}

			//FCAL
			if(locFCALHitFoundFlag)
			{
				dHist_FCALTrackDistanceVsP->Fill(locP, locBestFCALMatchParams->dDOCAToShower);
				dHist_FCALTrackDistanceVsTheta->Fill(locTheta, locBestFCALMatchParams->dDOCAToShower);
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
	}

	//Check for slow protons stopping in target: don't expect any hits
	bool locMassiveParticleFlag = (ParticleMass(dMissingPID) + 0.001 > ParticleMass(Proton));
	if((locP < 0.2) && locMassiveParticleFlag)
		return false; //don't even bother to try
	if((locP < 0.3) && locMassiveParticleFlag)
		return true;
	if((locP < 0.4) && (locTheta < 30.0) && locMassiveParticleFlag)
		return true;

	//REQUIRE: ST hit
	if(!locSCHitFoundFlag)
		return false;
	double locSCDeltaPhiCut = 25.0 + 0.12*exp(0.13*(locSCProjectedZ - 50.0));
	if(fabs(locSCDeltaPhi) > locSCDeltaPhiCut)
		return false; //could expand, would only effect 2pi missing pi forward-z

	//Check for slow protons stopping in the drift chambers: don't expect any further hits
	if((locP < 0.4) && locMassiveParticleFlag)
		return true;

	//REQUIRE EITHER BCAL hit OR TOF && FCAL hits

	//BCAL Hit
	if(locBCALHitFoundFlag)
		return ((fabs(locBCALDeltaPhi) < 2.0) && (fabs(locBestBCALMatchParams->dDeltaZToShower) < 15.0));

	//FCAL & TOF Hit:
	if(!locTOFHitFoundFlag || !locFCALHitFoundFlag)
		return false;

	if(locTOFHitFoundFlag)
		return (locTOFDistance < 30.0);

	return (locBestFCALMatchParams->dDOCAToShower < 50.0); //is always true: there's always another FCAL hit
}

