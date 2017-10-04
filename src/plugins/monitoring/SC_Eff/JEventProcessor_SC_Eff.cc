// $Id$
//
//    File: JEventProcessor_SC_Eff.cc
//

#include "JEventProcessor_SC_Eff.h"

// Routine used to create our JEventProcessor
extern "C"
{
	void InitPlugin(JApplication *locApplication)
	{
		InitJANAPlugin(locApplication);
		locApplication->AddProcessor(new JEventProcessor_SC_Eff()); //register this plugin
	}
} // "C"

//define static local variable //declared in header file
thread_local DTreeFillData JEventProcessor_SC_Eff::dTreeFillData;

//------------------
// init
//------------------
jerror_t JEventProcessor_SC_Eff::init(void)
{
	//TRACK REQUIREMENTS
	dMinTrackingFOM = 5.73303E-7; // +/- 5 sigma
	dMinNumTrackHits = 14; //e.g. 6 in CDC, 8 in 
	dMinHitRingsPerCDCSuperlayer = 3;
	dMinHitPlanesPerFDCPackage = 4;
	dMaxVertexR = 1.0;
	dCutAction_TrackHitPattern = new DCutAction_TrackHitPattern(NULL, dMinHitRingsPerCDCSuperlayer, dMinHitPlanesPerFDCPackage);
	//action initialize not necessary: is empty
	dMaxPIDDeltaTMap[SYS_BCAL] = 1.0;
	dMaxPIDDeltaTMap[SYS_TOF] = 1.0;

	TDirectory* locOriginalDir = gDirectory;
	gDirectory->mkdir("SC_Eff")->cd();

	//Histograms
	dHist_HitFound = new TH2I("HitFound", "Hit Found;Projected Hit-Z (cm);Sector", 120, 40.0, 100.0, 30, 0.5, 30.5);
	dHist_HitTotal = new TH2I("HitTotal", "Hit Total;Projected Hit-Z (cm);Sector", 120, 40.0, 100.0, 30, 0.5, 30.5);
	
	// back to original dir
	locOriginalDir->cd();

	//TTREE INTERFACE
	//MUST DELETE WHEN FINISHED: OR ELSE DATA WON'T BE SAVED!!!
	dTreeInterface = DTreeInterface::Create_DTreeInterface("sc_eff", "tree_sc_eff.root");

	//TTREE BRANCHES
	DTreeBranchRegister locTreeBranchRegister;

	//TRACK
	locTreeBranchRegister.Register_Single<Int_t>("PID_PDG"); //gives charge, mass, beta
	locTreeBranchRegister.Register_Single<Float_t>("TrackVertexZ");
	locTreeBranchRegister.Register_Single<TVector3>("TrackP3");
	locTreeBranchRegister.Register_Single<UInt_t>("TrackCDCRings"); //rings correspond to bits (1 -> 28)
	locTreeBranchRegister.Register_Single<UInt_t>("TrackFDCPlanes"); //planes correspond to bits (1 -> 24)

	//PROJECTED HIT
	locTreeBranchRegister.Register_Single<Float_t>("ProjectedSCHitPhi"); //degrees
	locTreeBranchRegister.Register_Single<Float_t>("ProjectedSCHitZ");
	locTreeBranchRegister.Register_Single<UChar_t>("ProjectedSCHitSector");

	//SEARCH
	locTreeBranchRegister.Register_Single<UChar_t>("NumMatchingSCHits");
	locTreeBranchRegister.Register_FundamentalArray<UChar_t>("SCHitSector", "NumMatchingSCHits"); //0 if none in time: PID:OUT_OF_TIME_CUT
	locTreeBranchRegister.Register_FundamentalArray<Float_t>("TrackHitDeltaPhi", "NumMatchingSCHits"); //is signed: SC - Track
	locTreeBranchRegister.Register_FundamentalArray<Float_t>("TrackHitDeltaT", "NumMatchingSCHits"); //is signed: SC - Track
	locTreeBranchRegister.Register_FundamentalArray<Bool_t>("IsMatchedToTrack", "NumMatchingSCHits"); //false if not registered in DDetectorMatches

	//REGISTER BRANCHES
	dTreeInterface->Create_Branches(locTreeBranchRegister);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_SC_Eff::brun(jana::JEventLoop* locEventLoop, int locRunNumber)
{
	// This is called whenever the run number changes

	return NOERROR;
}

//------------------
// evnt
//------------------

jerror_t JEventProcessor_SC_Eff::evnt(jana::JEventLoop* locEventLoop, uint64_t locEventNumber)
{
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// locEventLoop->Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.

	// This plugin is used to determine the reconstruction efficiency of hits in the BCAL
		// Note, this is hit-level, not shower-level.  Hits: By sector/layer/module/end

	//CUT ON TRIGGER TYPE
	const DTrigger* locTrigger = NULL;
	locEventLoop->GetSingle(locTrigger);
	if(locTrigger->Get_L1FrontPanelTriggerBits() != 0)
		return NOERROR;

	//SEE IF SC REQUIRED TO TRIGGER
	uint32_t locTriggerBits = locTrigger->Get_L1TriggerBits();
	if(locTriggerBits >= 32)
		return NOERROR;

	const DEventRFBunch* locEventRFBunch = NULL;
	locEventLoop->GetSingle(locEventRFBunch);
	if(locEventRFBunch->dNumParticleVotes <= 1)
		return NOERROR; //don't trust PID: beta-dependence

	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks);

	vector<const DSCHit*> locSCHits;
	locEventLoop->Get(locSCHits);

	const DParticleID* locParticleID = NULL;
	locEventLoop->GetSingle(locParticleID);

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	//Try to select the most-pure sample of tracks possible
	set<const DChargedTrackHypothesis*> locBestTracks;
	for(auto& locChargedTrack : locChargedTracks)
	{
		//loop over PID hypotheses and find the best (if any good enough)
		const DChargedTrackHypothesis* locBestChargedTrackHypothesis = nullptr;
		double locBestTrackingFOM = -1.0;
		for(auto& locChargedTrackHypothesis : locChargedTrack->dChargedTrackHypotheses)
		{
			if((locChargedTrackHypothesis->PID() == Electron) || (locChargedTrackHypothesis->PID() == Positron))
				continue; //don't evaluate for these

			if(locChargedTrackHypothesis->position().Perp() > dMaxVertexR)
				continue; //don't trust reconstruction if not close to target

			//Need PID for beta-dependence
			if(!Cut_PIDDeltaT(locChargedTrackHypothesis))
				continue; //also requires match to BCAL or TOF: no need for separate check

			auto locTrackTimeBased = locChargedTrackHypothesis->Get_TrackTimeBased();
			if(locTrackTimeBased->FOM < dMinTrackingFOM)
				continue; //don't trust tracking results: bad tracking FOM

			if(!dCutAction_TrackHitPattern->Cut_TrackHitPattern(locParticleID, locTrackTimeBased))
				continue; //don't trust tracking results: not many grouped hits

			unsigned int locNumTrackHits = locTrackTimeBased->Ndof + 5;
			if(locNumTrackHits < dMinNumTrackHits)
				continue; //don't trust tracking results: not many hits

			if(locTrackTimeBased->FOM < locBestTrackingFOM)
				continue; //not the best mass hypothesis

			locBestTrackingFOM = locTrackTimeBased->FOM;
			locBestChargedTrackHypothesis = locChargedTrackHypothesis;
		}

		//if passed all cuts, save the best
		if(locBestChargedTrackHypothesis != nullptr)
			locBestTracks.insert(locBestChargedTrackHypothesis);
	}

	//for histograms, keep running total of what to fill
	//will fill at end: only one lock
	vector<pair<int, double> > locHitMap_HitFound, locHitMap_HitTotal; //int: sector, double: projected-z

	// Loop over the good tracks, using the best DTrackTimeBased object for each
	for(auto& locChargedTrackHypothesis : locBestTracks)
	{
		auto locTrackTimeBased = locChargedTrackHypothesis->Get_TrackTimeBased();

		//Predict ST Surface Hit Location
		DVector3 locPredictedSurfacePosition;
		bool locProjBarrelRegion = false;
		unsigned int locPredictedSCSector = locParticleID->PredictSCSector(locTrackTimeBased->rt, &locPredictedSurfacePosition, &locProjBarrelRegion);

		//Save for hist if is matched
		pair<int, double> locHitPair(locPredictedSCSector, locPredictedSurfacePosition.Z());
		locHitMap_HitTotal.push_back(locHitPair);
		if(locDetectorMatches->Get_IsMatchedToDetector(locTrackTimeBased, SYS_START))
			locHitMap_HitFound.push_back(locHitPair);

		//TRACK
		dTreeFillData.Fill_Single<Int_t>("PID_PDG", PDGtype(locChargedTrackHypothesis->PID()));
		dTreeFillData.Fill_Single<Float_t>("TrackVertexZ", locChargedTrackHypothesis->position().Z());
		dTreeFillData.Fill_Single<UInt_t>("TrackCDCRings", locTrackTimeBased->dCDCRings);
		dTreeFillData.Fill_Single<UInt_t>("TrackFDCPlanes", locTrackTimeBased->dFDCPlanes);
		DVector3 locDP3 = locChargedTrackHypothesis->momentum();
		TVector3 locP3(locDP3.X(), locDP3.Y(), locDP3.Z());
		dTreeFillData.Fill_Single<TVector3>("TrackP3", locP3);

		//PROJECTED
		dTreeFillData.Fill_Single<Float_t>("ProjectedSCHitPhi", locPredictedSurfacePosition.Phi()*180.0/TMath::Pi());
		dTreeFillData.Fill_Single<Float_t>("ProjectedSCHitZ", locPredictedSurfacePosition.Z());
		dTreeFillData.Fill_Single<UChar_t>("ProjectedSCHitSector", locPredictedSCSector);

		vector<shared_ptr<const DSCHitMatchParams>> locAllSCMatchParams;
		locDetectorMatches->Get_SCMatchParams(locTrackTimeBased, locAllSCMatchParams);

		//pre-sort
		set<const DSCHit*> locMatchedSCHits;
		for(auto locSCMatchParams : locAllSCMatchParams)
			locMatchedSCHits.insert(locSCMatchParams->dSCHit);

		//FILL ARRAYS
		size_t locArrayIndex = 0;
		for(size_t loc_i = 0; loc_i < locSCHits.size(); ++loc_i)
		{
			const DSCHit* locSCHit = locSCHits[loc_i];

			shared_ptr<DSCHitMatchParams> locSCHitMatchParams;
			if(!locParticleID->Distance_ToTrack(locTrackTimeBased->rt, locSCHit, locTrackTimeBased->t0(), locSCHitMatchParams))
				continue;

			double locDeltaT = locSCHitMatchParams->dHitTime - locSCHitMatchParams->dFlightTime - locChargedTrackHypothesis->time();
			bool locIsMatchedToTrack = (locMatchedSCHits.find(locSCHit) != locMatchedSCHits.end());

			dTreeFillData.Fill_Array<UChar_t>("SCHitSector", locSCHit->sector, locArrayIndex);
			dTreeFillData.Fill_Array<Float_t>("TrackHitDeltaPhi", locSCHitMatchParams->dDeltaPhiToHit, locArrayIndex); //is signed: SC - Track
			dTreeFillData.Fill_Array<Float_t>("TrackHitDeltaT", locDeltaT, locArrayIndex); //is signed: SC - Track
			dTreeFillData.Fill_Array<Bool_t>("IsMatchedToTrack", locIsMatchedToTrack, locArrayIndex);

			++locArrayIndex;
		}

		//FILL ARRAY SIZE
		dTreeFillData.Fill_Single<UChar_t>("NumMatchingSCHits", locArrayIndex);

		//FILL TTREE
		dTreeInterface->Fill(dTreeFillData);
	}

	// FILL HISTOGRAMS
	// Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
	{
		//Fill Found
		for(auto& locHitPair : locHitMap_HitFound)
			dHist_HitFound->Fill(locHitPair.second, locHitPair.first);

		//Fill Total
		for(auto& locHitPair : locHitMap_HitTotal)
			dHist_HitTotal->Fill(locHitPair.second, locHitPair.first);
	}
	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

	return NOERROR;
}

bool JEventProcessor_SC_Eff::Cut_PIDDeltaT(const DChargedTrackHypothesis* locChargedTrackHypothesis)
{
	DetectorSystem_t locSystem = locChargedTrackHypothesis->t1_detector();
	if(dMaxPIDDeltaTMap.find(locSystem) == dMaxPIDDeltaTMap.end())
		return false;

	double locDeltaT = locChargedTrackHypothesis->time() - locChargedTrackHypothesis->t0();
	return (fabs(locDeltaT) <= dMaxPIDDeltaTMap[locSystem]);
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_SC_Eff::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_SC_Eff::fini(void)
{
	// Called before program exit after event processing is finished.  

	delete dCutAction_TrackHitPattern;
	delete dTreeInterface; //saves trees to file, closes file

	return NOERROR;
}

