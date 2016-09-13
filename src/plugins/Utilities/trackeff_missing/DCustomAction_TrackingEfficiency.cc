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

	const DReaction* locReaction = Get_Reaction();
	if(!locReaction->Get_MissingPID(dMissingPID))
		return; //invalid reaction setup
	if(dMissingPID == Unknown)
		return; //invalid reaction setup

	locEventLoop->GetSingle(dAnalysisUtilities);
	locEventLoop->GetSingle(dParticleID);

	//CREATE TTREE, TFILE
	dTreeInterface = DTreeInterface::Create_DTreeInterface(locReaction->Get_ReactionName(), "tree_trackeff.root");

	//TTREE BRANCHES
	DTreeBranchRegister locBranchRegister;

	//USER INFO
	TList* locUserInfo = locBranchRegister.Get_UserInfo();
	TMap* locMiscInfoMap = new TMap(); //collection of pairs
	locMiscInfoMap->SetName("MiscInfoMap");
	locUserInfo->Add(locMiscInfoMap);
	//set pid
	ostringstream locPIDStream;
	locPIDStream << PDGtype(dMissingPID);
	locMiscInfoMap->Add(new TObjString("MissingPID_PDG"), new TObjString(locPIDStream.str().c_str()));

	//CHANNEL INFO
	locBranchRegister.Register_Single<Float_t>("BeamEnergy");
	locBranchRegister.Register_Single<Float_t>("BeamRFDeltaT");
	locBranchRegister.Register_Single<Float_t>("MissingMassSquared"); //is measured
	locBranchRegister.Register_Single<TVector3>("MissingP3"); //is kinfit if kinfit
	locBranchRegister.Register_Single<Float_t>("ComboVertexZ");
	locBranchRegister.Register_Single<UChar_t>("NumExtraTracks");

	//TRACKING INFO: //"Recon:" Time-based track
	locBranchRegister.Register_Single<Float_t>("ReconMatchFOM"); //FOM < 0 if nothing, no-match
	locBranchRegister.Register_Single<Float_t>("ReconTrackingFOM"); //FOM < 0 if nothing, no-match
	locBranchRegister.Register_Single<TVector3>("ReconP3");
	locBranchRegister.Register_Single<UInt_t>("TrackCDCRings"); //rings correspond to bits (1 -> 28)
	locBranchRegister.Register_Single<UInt_t>("TrackFDCPlanes"); //planes correspond to bits (1 -> 24)

	//HADRONIC BCAL SHOWER EFFICIENCY: TIMING, MATCHING //cannot get accurate PID without missing-track study
	locBranchRegister.Register_Single<UChar_t>("ProjectedBCALSector"); //4*(module - 1) + sector: 1 -> 192 //0 if proj to miss
	locBranchRegister.Register_Single<Float_t>("ProjectedBCALHitPhi"); //degrees
	locBranchRegister.Register_Single<Float_t>("ProjectedBCALHitZ");
	locBranchRegister.Register_Single<Float_t>("NearestShowerEnergy"); // < 0 if none in time: PID:OUT_OF_TIME_CUT
	locBranchRegister.Register_Single<Float_t>("TrackDeltaPhiToShower"); //is signed: BCAL - Track //degrees
	locBranchRegister.Register_Single<Float_t>("TrackDeltaZToShower"); //is signed: BCAL - Track
	locBranchRegister.Register_Single<Bool_t>("IsMatchedToBCALShower");
	locBranchRegister.Register_Single<Float_t>("BCALDeltaT");
	locBranchRegister.Register_Single<Float_t>("BCALTimeFOM");

	//REGISTER BRANCHES
	dTreeInterface->Create_Branches(locBranchRegister);
}

bool DCustomAction_TrackingEfficiency::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//Write custom code to perform an action on the INPUT DParticleCombo (DParticleCombo)
	//NEVER: Grab DParticleCombo or DAnalysisResults objects (of any tag!) from the JEventLoop within this function
	//NEVER: Grab objects that are created post-kinfit (e.g. DKinFitResults, etc.) from the JEventLoop if Get_UseKinFitResultsFlag() == false: CAN CAUSE INFINITE DEPENDENCY LOOP

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	vector<const DBCALShower*> locBCALShowers;
	locEventLoop->Get(locBCALShowers);

	vector<const DChargedTrack*> locUnusedChargedTracks;
	dAnalysisUtilities->Get_UnusedChargedTracks(locEventLoop, locParticleCombo, locUnusedChargedTracks);

	const DEventRFBunch* locEventRFBunch = locParticleCombo->Get_EventRFBunch();

	/*********************************************** MISSING PARTICLE INFO ***********************************************/

	if(dMissingPID == Unknown)
		return true; //invalid reaction setup
	if(ParticleCharge(dMissingPID) == 0)
		return true; //NOT SUPPORTED

	const DKinematicData* locMissingParticle = locParticleCombo->Get_MissingParticle(); //is NULL if no kinfit!!
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(0);

	// Get missing particle p4 & covariance
	DLorentzVector locMeasuredMissingP4 = dAnalysisUtilities->Calc_MissingP4(locParticleCombo, false);
	DVector3 locMissingP3 = locMeasuredMissingP4.Vect();
	DMatrixDSym locMissingCovarianceMatrix(3);
	const DKinematicData* locBeamParticle = NULL;
	if(locMissingParticle == NULL) //measured, or kinfit failed
	{
		locMissingCovarianceMatrix = dAnalysisUtilities->Calc_MissingP3Covariance(locParticleCombo);
		locBeamParticle = locParticleComboStep->Get_InitialParticle_Measured();
	}
	else //kinfit succeeded
	{
		locMissingP3 = locMissingParticle->momentum();
		DMatrixDSym locKinFitCovarianceMatrix = locMissingParticle->errorMatrix();
		locKinFitCovarianceMatrix.ResizeTo(3, 3);
		locMissingCovarianceMatrix = locKinFitCovarianceMatrix;
		locBeamParticle = locParticleComboStep->Get_InitialParticle();
	}
	
	double locVertexZ = locParticleCombo->Get_EventVertex().Z();
	double locBeamRFDeltaT = locBeamParticle->time() - locEventRFBunch->dTime;

	//get number of extra tracks that have at least 10 hits
	size_t locNumExtraTracks = 0;
	for(auto locChargedTrack : locUnusedChargedTracks)
	{
		auto locBestHypothesis = locChargedTrack->Get_BestFOM();
		auto locNumTrackHits = locBestHypothesis->dNDF_Track + 5;
		if(locNumTrackHits >= 10)
			++locNumExtraTracks;
	}

	//FILL CHANNEL INFO
	dTreeFillData.Fill_Single<Float_t>("BeamEnergy", locBeamParticle->energy());
	dTreeFillData.Fill_Single<Float_t>("BeamRFDeltaT", locBeamRFDeltaT);
	dTreeFillData.Fill_Single<UChar_t>("NumExtraTracks", (UChar_t)locNumExtraTracks);
	dTreeFillData.Fill_Single<Float_t>("MissingMassSquared", locMeasuredMissingP4.M2());
	TVector3 locTMissingP3(locMissingP3.Px(), locMissingP3.Py(), locMissingP3.Pz());
	dTreeFillData.Fill_Single<TVector3>("MissingP3", locTMissingP3); //is kinfit if kinfit
	dTreeFillData.Fill_Single<Float_t>("ComboVertexZ", locVertexZ);

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

		DVector3 locDeltaP3 = locUnusedTimeBasedTracks[loc_i]->momentum() - locMissingP3;
		double locMatchFOM = Calc_MatchFOM(locDeltaP3, locCovarianceMatrix);

		if(locMatchFOM > locBestTimeBasedMatchFOM)
		{
			locBestTimeBasedMatchFOM = locMatchFOM;
			locBestTrackTimeBased = locUnusedTimeBasedTracks[loc_i];
		}
	}

	//FILL TRACKING INFO: //"Recon:" Time-based track
	dTreeFillData.Fill_Single<Float_t>("ReconMatchFOM", locBestTimeBasedMatchFOM); //FOM < 0 if nothing, no-match
	if(locBestTrackTimeBased == nullptr)
	{
		dTreeFillData.Fill_Single<Float_t>("ReconTrackingFOM", -1.0);
		dTreeFillData.Fill_Single<TVector3>("ReconP3", TVector3());
		dTreeFillData.Fill_Single<UInt_t>("TrackCDCRings", 0);
		dTreeFillData.Fill_Single<UInt_t>("TrackFDCPlanes", 0);

		//HADRONIC BCAL SHOWER EFFICIENCY: TIMING, MATCHING
		dTreeFillData.Fill_Single<Float_t>("ProjectedBCALHitPhi", 0.0);
		dTreeFillData.Fill_Single<Float_t>("ProjectedBCALHitZ", 0.0);
		dTreeFillData.Fill_Single<UChar_t>("ProjectedBCALSector", 0);
		dTreeFillData.Fill_Single<Float_t>("NearestShowerEnergy", 0.0);
		dTreeFillData.Fill_Single<Float_t>("TrackDeltaPhiToShower", 0.0);
		dTreeFillData.Fill_Single<Float_t>("TrackDeltaZToShower", 0.0);
		dTreeFillData.Fill_Single<Bool_t>("IsMatchedToBCALShower", false);
		dTreeFillData.Fill_Single<Float_t>("BCALDeltaT", 0.0);
		dTreeFillData.Fill_Single<Float_t>("BCALTimeFOM", -1.0);

		//FILL TTREE
		dTreeInterface->Fill(dTreeFillData);
		return true;
	}

	//RESUME FILL TRACKING INFO: //"Recon:" Time-based track
	dTreeFillData.Fill_Single<Float_t>("ReconTrackingFOM", locBestTrackTimeBased->FOM);
	DVector3 locDP3 = locBestTrackTimeBased->momentum();
	TVector3 locP3(locDP3.X(), locDP3.Y(), locDP3.Z());
	dTreeFillData.Fill_Single<TVector3>("ReconP3", locP3);
	dTreeFillData.Fill_Single<UInt_t>("TrackCDCRings", locBestTrackTimeBased->dCDCRings);
	dTreeFillData.Fill_Single<UInt_t>("TrackFDCPlanes", locBestTrackTimeBased->dFDCPlanes);

	/********************************************* BCAL SHOWER RECONSTRUCTION ********************************************/

	//Predict BCAL Surface Hit Location
	unsigned int locPredictedSurfaceModule = 0, locPredictedSurfaceSector = 0;
	DVector3 locPredictedSurfacePosition;
	dParticleID->PredictBCALWedge(locBestTrackTimeBased->rt, locPredictedSurfaceModule, locPredictedSurfaceSector, &locPredictedSurfacePosition);
	unsigned int locTrackProjectedBCALSector = 4*(locPredictedSurfaceModule - 1) + locPredictedSurfaceSector; //0 if misses

	//FILL PROJECTED HIT POSITION
	dTreeFillData.Fill_Single<UChar_t>("ProjectedBCALSector", locTrackProjectedBCALSector);
	dTreeFillData.Fill_Single<Float_t>("ProjectedBCALHitPhi", locPredictedSurfacePosition.Phi()*180.0/TMath::Pi());
	dTreeFillData.Fill_Single<Float_t>("ProjectedBCALHitZ", locPredictedSurfacePosition.Z());

	//Get closest BCAL shower
	double locBestMatchDeltaPhi = 999.9, locBestMatchDeltaZ = 999.9;
	const DBCALShower* locBestBCALShower = dParticleID->Get_ClosestToTrack_BCAL(locBestTrackTimeBased, locBCALShowers, locBestMatchDeltaPhi, locBestMatchDeltaZ);

	//If none, fill and return
	if(locBestBCALShower == NULL)
	{
		//HADRONIC BCAL SHOWER EFFICIENCY
		dTreeFillData.Fill_Single<Float_t>("NearestShowerEnergy", -1.0);
		dTreeFillData.Fill_Single<Float_t>("TrackDeltaPhiToShower", 0.0);
		dTreeFillData.Fill_Single<Float_t>("TrackDeltaZToShower", 0.0);
		dTreeFillData.Fill_Single<Bool_t>("IsMatchedToBCALShower", false);
		dTreeFillData.Fill_Single<Float_t>("BCALDeltaT", 0.0);
		dTreeFillData.Fill_Single<Float_t>("BCALTimeFOM", -1.0);

		//FILL TTREE
		dTreeInterface->Fill(dTreeFillData);
		return true;
	}

	//HADRONIC BCAL SHOWER EFFICIENCY
	dTreeFillData.Fill_Single<Float_t>("NearestShowerEnergy", locBestBCALShower->E);
	dTreeFillData.Fill_Single<Float_t>("TrackDeltaPhiToShower", locBestMatchDeltaPhi*180.0/TMath::Pi());
	dTreeFillData.Fill_Single<Float_t>("TrackDeltaZToShower", locBestMatchDeltaZ);

	//See if shower is matched
	DBCALShowerMatchParams locBCALShowerMatchParams;
	if(!dParticleID->Get_BestBCALMatchParams(locBestTrackTimeBased, locDetectorMatches, locBCALShowerMatchParams))
	{
		//SHOWER NOT MATCHED
		dTreeFillData.Fill_Single<Bool_t>("IsMatchedToBCALShower", false);
		dTreeFillData.Fill_Single<Float_t>("BCALDeltaT", 0.0);
		dTreeFillData.Fill_Single<Float_t>("BCALTimeFOM", -1.0);

		//FILL TTREE
		dTreeInterface->Fill(dTreeFillData);
		return true;
	}

	double locStartTime = dParticleID->Calc_PropagatedRFTime(locBestTrackTimeBased, locEventRFBunch);
	double locDeltaT = locBestBCALShower->t - locBCALShowerMatchParams.dFlightTime - locStartTime;

	//FILL SHOWER MATCHED
	dTreeFillData.Fill_Single<Bool_t>("IsMatchedToBCALShower", true);
	dTreeFillData.Fill_Single<Float_t>("BCALDeltaT", locDeltaT);
	dTreeFillData.Fill_Single<Float_t>("BCALTimeFOM", -1.0);

	//FILL TTREE
	dTreeInterface->Fill(dTreeFillData);
	return true;
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

