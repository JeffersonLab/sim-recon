// $Id$
//
//    File: JEventProcessor_FCAL_Hadronic_Eff.cc
//

#include "JEventProcessor_FCAL_Hadronic_Eff.h"

// Routine used to create our JEventProcessor
extern "C"
{
	void InitPlugin(JApplication *locApplication)
	{
		InitJANAPlugin(locApplication);
		locApplication->AddProcessor(new JEventProcessor_FCAL_Hadronic_Eff()); //register this plugin
	}
} // "C"

//define static local variable //declared in header file
thread_local DTreeFillData JEventProcessor_FCAL_Hadronic_Eff::dTreeFillData;

//------------------
// init
//------------------
jerror_t JEventProcessor_FCAL_Hadronic_Eff::init(void)
{
	//TRACK REQUIREMENTS
	dMaxTOFDeltaT = 1.0;
	dMinTrackingFOM = 5.73303E-7; // +/- 5 sigma
	dMinNumTrackHits = 16; //e.g. 4 in each FDC plane
	dMinHitRingsPerCDCSuperlayer = 0;
	dMinHitPlanesPerFDCPackage = 4;
	dMaxFCALThetaCut = 15.0;
	dMaxVertexR = 1.0;
	dCutAction_TrackHitPattern = new DCutAction_TrackHitPattern(NULL, dMinHitRingsPerCDCSuperlayer, dMinHitPlanesPerFDCPackage);
	//action initialize not necessary: is empty

	TDirectory* locOriginalDir = gDirectory;
	gDirectory->mkdir("FCAL_Hadronic_Eff")->cd();

	//Histograms
	string locHistName, locHistTitle;

	locHistName = "TrackFCALYVsX_HasHit";
	locHistTitle = "FCAL Has Hit;Projected FCAL Hit X (cm);Projected FCAL Hit Y (cm)";
	dHist_TrackFCALYVsX_HasHit = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 130, -130.0, 130.0, 130, -130.0, 130.0);

	locHistName = "TrackFCALYVsX_NoHit";
	locHistTitle = "FCAL Total Hit;Projected FCAL Hit X (cm);Projected FCAL Hit Y (cm)";
	dHist_TrackFCALYVsX_TotalHit = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 130, -130.0, 130.0, 130, -130.0, 130.0);

	locHistName = "TrackFCALRowVsColumn_HasHit";
	locHistTitle = "FCAL Has Hit;Projected FCAL Hit Column;Projected FCAL Hit Row";
	dHist_TrackFCALRowVsColumn_HasHit = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 59, -0.5, 58.5, 59, -0.5, 58.5);

	locHistName = "TrackFCALRowVsColumn_NoHit";
	locHistTitle = "FCAL Total Hit;Projected FCAL Hit Column;Projected FCAL Hit Row";
	dHist_TrackFCALRowVsColumn_TotalHit = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 59, -0.5, 58.5, 59, -0.5, 58.5);

	// back to original dir
	locOriginalDir->cd();

	//TTREE INTERFACE
	//MUST DELETE WHEN FINISHED: OR ELSE DATA WON'T BE SAVED!!!
	dTreeInterface = DTreeInterface::Create_DTreeInterface("fcal_hadronic_eff", "tree_fcal_hadronic_eff.root");

	//TTREE BRANCHES
	DTreeBranchRegister locTreeBranchRegister;

	//TRACK
	locTreeBranchRegister.Register_Single<Int_t>("PID_PDG"); //gives charge, mass, beta
	locTreeBranchRegister.Register_Single<Float_t>("TrackVertexZ");
	locTreeBranchRegister.Register_Single<TVector3>("TrackP3");
	locTreeBranchRegister.Register_Single<UInt_t>("TrackCDCRings"); //rings correspond to bits (1 -> 28)
	locTreeBranchRegister.Register_Single<UInt_t>("TrackFDCPlanes"); //planes correspond to bits (1 -> 24)

	//PROJECTED FCAL
	locTreeBranchRegister.Register_Single<Float_t>("ProjectedFCALX");
	locTreeBranchRegister.Register_Single<Float_t>("ProjectedFCALY");
	locTreeBranchRegister.Register_Single<UChar_t>("ProjectedFCALRow"); //0 if projected to miss
	locTreeBranchRegister.Register_Single<UChar_t>("ProjectedFCALColumn"); //0 if projected to miss

	//NEAREST FCAL SHOWER
	locTreeBranchRegister.Register_Single<Float_t>("NearestFCALEnergy"); //if < 0: nothing in time: PID:OUT_OF_TIME_CUT
	locTreeBranchRegister.Register_Single<Float_t>("NearestFCALX"); //ignore if E < 0
	locTreeBranchRegister.Register_Single<Float_t>("NearestFCALY"); //ignore if E < 0
	locTreeBranchRegister.Register_Single<Bool_t>("IsMatchedToTrack"); //false if not registered in DDetectorMatches

	//FOUND FCAL //only if matched: for evaluating PID quality
	locTreeBranchRegister.Register_Single<Float_t>("FCALDeltaT"); //FCAL - RF
	locTreeBranchRegister.Register_Single<Float_t>("FCALTimeFOM");

	//REGISTER BRANCHES
	dTreeInterface->Create_Branches(locTreeBranchRegister);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_FCAL_Hadronic_Eff::brun(jana::JEventLoop* locEventLoop, int locRunNumber)
{
	// This is called whenever the run number changes

	return NOERROR;
}

//------------------
// evnt
//------------------

jerror_t JEventProcessor_FCAL_Hadronic_Eff::evnt(jana::JEventLoop* locEventLoop, uint64_t locEventNumber)
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

	//SEE IF FCAL REQUIRED TO TRIGGER
	uint32_t locTriggerBits = locTrigger->Get_L1TriggerBits();
	if(((locTriggerBits & 4) != 4) && ((locTriggerBits & 64) != 64))
		return NOERROR;

	const DEventRFBunch* locEventRFBunch = NULL;
	locEventLoop->GetSingle(locEventRFBunch);
	if(locEventRFBunch->dNumParticleVotes <= 1)
		return NOERROR; //don't trust PID: beta-dependence

	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks);

	const DParticleID* locParticleID = NULL;
	locEventLoop->GetSingle(locParticleID);

	vector<const DFCALShower*> locFCALShowers;
	locEventLoop->Get(locFCALShowers);

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
			if(locChargedTrackHypothesis->position().Perp() > dMaxVertexR)
				continue; //don't trust reconstruction if not close to target

			//Need PID for beta-dependence, but cannot use FCAL info: would bias
			if(!Cut_TOFTiming(locChargedTrackHypothesis))
				continue; //also requires match to TOF: no need for separate check

			const DTrackTimeBased* locTrackTimeBased = NULL;
			locChargedTrackHypothesis->GetSingle(locTrackTimeBased);
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

	// Loop over the good tracks, using the best DTrackTimeBased object for each
	for(auto& locChargedTrackHypothesis : locBestTracks)
	{
		const DTrackTimeBased* locTrackTimeBased = NULL;
		locChargedTrackHypothesis->GetSingle(locTrackTimeBased);

		//Predict FCAL Surface Hit Location
		DVector3 locProjectedFCALIntersection;
		unsigned int locProjectedFCALRow = 0, locProjectedFCALColumn = 0;
		if(!locParticleID->PredictFCALHit(locTrackTimeBased->rt, locProjectedFCALRow, locProjectedFCALColumn, &locProjectedFCALIntersection))
		{
			if(locTrackTimeBased->momentum().Theta()*180.0/TMath::Pi() > dMaxFCALThetaCut)
				continue; //not predicted to hit FCAL
		}

		//Find closest FCAL Shower
		double locBestDistance = 999.0;
		const DFCALShower* locFCALShower = locParticleID->Get_ClosestToTrack_FCAL(locTrackTimeBased, locFCALShowers, locBestDistance);

		//Is match to FCAL shower?
		const DFCALShowerMatchParams* locFCALShowerMatchParams = locChargedTrackHypothesis->Get_FCALShowerMatchParams();
		bool locIsMatchedToTrack = (locFCALShowerMatchParams != nullptr);

		//If so, calc PID info: timing
		double locFCALDeltaT = 0.0;
		double locFCALTimeFOM = Calc_FCALTiming(locChargedTrackHypothesis, locParticleID, locEventRFBunch, locFCALDeltaT);

		// FILL HISTOGRAMS
		// Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
		japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
		{
			//TOFPoint
			dHist_TrackFCALYVsX_TotalHit->Fill(locProjectedFCALIntersection.X(), locProjectedFCALIntersection.Y());
			dHist_TrackFCALRowVsColumn_TotalHit->Fill(locProjectedFCALColumn, locProjectedFCALRow);

			//TOFPoint
			if(locIsMatchedToTrack)
			{
				dHist_TrackFCALYVsX_HasHit->Fill(locProjectedFCALIntersection.X(), locProjectedFCALIntersection.Y());
				dHist_TrackFCALRowVsColumn_HasHit->Fill(locProjectedFCALColumn, locProjectedFCALRow);
			}
		}
		japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

		//TRACK
		dTreeFillData.Fill_Single<Int_t>("PID_PDG", PDGtype(locChargedTrackHypothesis->PID()));
		dTreeFillData.Fill_Single<Float_t>("TrackVertexZ", locChargedTrackHypothesis->position().Z());
		dTreeFillData.Fill_Single<UInt_t>("TrackCDCRings", locTrackTimeBased->dCDCRings);
		dTreeFillData.Fill_Single<UInt_t>("TrackFDCPlanes", locTrackTimeBased->dFDCPlanes);
		DVector3 locDP3 = locChargedTrackHypothesis->momentum();
		TVector3 locP3(locDP3.X(), locDP3.Y(), locDP3.Z());
		dTreeFillData.Fill_Single<TVector3>("TrackP3", locP3);

		//PROJECTED FCAL
		dTreeFillData.Fill_Single<Float_t>("ProjectedFCALX", locProjectedFCALIntersection.X());
		dTreeFillData.Fill_Single<Float_t>("ProjectedFCALY", locProjectedFCALIntersection.Y());
		dTreeFillData.Fill_Single<UChar_t>("ProjectedFCALRow", locProjectedFCALRow);
		dTreeFillData.Fill_Single<UChar_t>("ProjectedFCALColumn", locProjectedFCALColumn);

		//NEAREST FCAL SHOWER
		if(locFCALShower != NULL)
		{
			dTreeFillData.Fill_Single<Float_t>("NearestFCALEnergy", locFCALShower->getEnergy());
			dTreeFillData.Fill_Single<Float_t>("NearestFCALX", locFCALShower->getPosition().X());
			dTreeFillData.Fill_Single<Float_t>("NearestFCALY", locFCALShower->getPosition().Y());
		}
		else
		{
			dTreeFillData.Fill_Single<Float_t>("NearestFCALEnergy", -1.0);
			dTreeFillData.Fill_Single<Float_t>("NearestFCALX", 0.0);
			dTreeFillData.Fill_Single<Float_t>("NearestFCALY", 0.0);
		}
		dTreeFillData.Fill_Single<Bool_t>("IsMatchedToTrack", locIsMatchedToTrack);

		//FOUND FCAL
		dTreeFillData.Fill_Single<Float_t>("FCALDeltaT", locFCALDeltaT);
		dTreeFillData.Fill_Single<Float_t>("FCALTimeFOM", locFCALTimeFOM);

		//FILL TTREE
		dTreeInterface->Fill(dTreeFillData);
	}

	return NOERROR;
}

double JEventProcessor_FCAL_Hadronic_Eff::Calc_FCALTiming(const DChargedTrackHypothesis* locChargedTrackHypothesis, const DParticleID* locParticleID, const DEventRFBunch* locEventRFBunch, double& locDeltaT)
{
	const DFCALShowerMatchParams* locFCALShowerMatchParams = locChargedTrackHypothesis->Get_FCALShowerMatchParams();
	if(locFCALShowerMatchParams == NULL)
		return false;

	double locStartTime = locParticleID->Calc_PropagatedRFTime(locChargedTrackHypothesis, locEventRFBunch);
	const DFCALShower* locFCALShower = locFCALShowerMatchParams->dFCALShower;
	locDeltaT = locFCALShower->getTime() - locFCALShowerMatchParams->dFlightTime - locStartTime;

	double locPIDFOM = -1.0; //not able to calc this correctly yet
	return locPIDFOM;
}

bool JEventProcessor_FCAL_Hadronic_Eff::Cut_TOFTiming(const DChargedTrackHypothesis* locChargedTrackHypothesis)
{
	if(locChargedTrackHypothesis->t1_detector() != SYS_TOF)
		return false;

	double locDeltaT = locChargedTrackHypothesis->time() - locChargedTrackHypothesis->t0();
	return (fabs(locDeltaT) <= dMaxTOFDeltaT);
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_FCAL_Hadronic_Eff::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_FCAL_Hadronic_Eff::fini(void)
{
	// Called before program exit after event processing is finished.  

	delete dCutAction_TrackHitPattern;
	delete dTreeInterface; //saves trees to file, closes file

	return NOERROR;
}

