// $Id$
//
//    File: JEventProcessor_TOF_Eff.cc
//

#include "JEventProcessor_TOF_Eff.h"

// Routine used to create our JEventProcessor
extern "C"
{
	void InitPlugin(JApplication *locApplication)
	{
		InitJANAPlugin(locApplication);
		locApplication->AddProcessor(new JEventProcessor_TOF_Eff()); //register this plugin
	}
} // "C"

//define static local variable //declared in header file
thread_local DTreeFillData JEventProcessor_TOF_Eff::dTreeFillData;

//------------------
// init
//------------------
jerror_t JEventProcessor_TOF_Eff::init(void)
{
	//TRACK REQUIREMENTS
	dMaxFCALDeltaT = 2.0;
	dMinTrackingFOM = 5.73303E-7; // +/- 5 sigma
	dMinNumTrackHits = 16; //e.g. 4 in each FDC plane
	dMinHitRingsPerCDCSuperlayer = 0;
	dMinHitPlanesPerFDCPackage = 4;
	dMaxVertexR = 1.0;
	dMaxTOFThetaCut = 15.0; //it could be projected to miss, but still hit it
	dCutAction_TrackHitPattern = new DCutAction_TrackHitPattern(NULL, dMinHitRingsPerCDCSuperlayer, dMinHitPlanesPerFDCPackage);
	//action initialize not necessary: is empty

	TDirectory* locOriginalDir = gDirectory;
	gDirectory->mkdir("TOF_Eff")->cd();

	//Histograms
	string locHistName, locHistTitle;
	dMinTOFPaddleMatchDistance = 9.0; //1.5*BARWIDTH //What DTOFPoint factory uses as of this writing

	//TOFPaddle
	gDirectory->mkdir("TOFPaddle")->cd();
	locHistName = "TrackYVsVerticalPaddle_HasHit_Top";
	locHistTitle = "TOF Paddle Has Top Hit;Projected Vertical Paddle;Projected TOF Hit Y (cm)";
	dHist_TOFPaddleTrackYVsVerticalPaddle_HasHit_Top = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 44, 0.5, 44.5, 130, -130.0, 130.0);

	locHistName = "TrackYVsVerticalPaddle_TotalHit_Top";
	locHistTitle = "TOF Paddle Total Top Hit;Projected Vertical Paddle;Projected TOF Hit Y (cm)";
	dHist_TOFPaddleTrackYVsVerticalPaddle_TotalHit_Top = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 44, 0.5, 44.5, 130, -130.0, 130.0);

	locHistName = "HorizontalPaddleVsTrackX_HasHit_North";
	locHistTitle = "TOF Paddle Has North Hit;Projected TOF Hit X (cm);Projected Horizontal Paddle";
	dHist_TOFPaddleHorizontalPaddleVsTrackX_HasHit_North = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 130, -130.0, 130.0, 44, 0.5, 44.5);

	locHistName = "HorizontalPaddleVsTrackX_TotalHit_North";
	locHistTitle = "TOF Paddle Total North Hit;Projected TOF Hit X (cm);Projected Horizontal Paddle";
	dHist_TOFPaddleHorizontalPaddleVsTrackX_TotalHit_North = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 130, -130.0, 130.0, 44, 0.5, 44.5);

	locHistName = "TrackYVsVerticalPaddle_HasHit_Bottom";
	locHistTitle = "TOF Paddle Has Bottom Hit;Projected Vertical Paddle;Projected TOF Hit Y (cm)";
	dHist_TOFPaddleTrackYVsVerticalPaddle_HasHit_Bottom = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 44, 0.5, 44.5, 130, -130.0, 130.0);

	locHistName = "TrackYVsVerticalPaddle_TotalHit_Bottom";
	locHistTitle = "TOF Paddle Total Bottom Hit;Projected Vertical Paddle;Projected TOF Hit Y (cm)";
	dHist_TOFPaddleTrackYVsVerticalPaddle_TotalHit_Bottom = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 44, 0.5, 44.5, 130, -130.0, 130.0);

	locHistName = "HorizontalPaddleVsTrackX_HasHit_South";
	locHistTitle = "TOF Paddle Has South Hit;Projected TOF Hit X (cm);Projected Horizontal Paddle";
	dHist_TOFPaddleHorizontalPaddleVsTrackX_HasHit_South = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 130, -130.0, 130.0, 44, 0.5, 44.5);

	locHistName = "HorizontalPaddleVsTrackX_TotalHit_South";
	locHistTitle = "TOF Paddle Total South Hit;Projected TOF Hit X (cm);Projected Horizontal Paddle";
	dHist_TOFPaddleHorizontalPaddleVsTrackX_TotalHit_South = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 130, -130.0, 130.0, 44, 0.5, 44.5);
	gDirectory->cd("..");

	//TOFPoint
	gDirectory->mkdir("TOFPoint")->cd();
	locHistName = "TrackTOFYVsX_HasHit";
	locHistTitle = "TOF Has Hit;Projected TOF Hit X (cm);Projected TOF Hit Y (cm)";
	dHist_TrackTOFYVsX_HasHit = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 130, -130.0, 130.0, 130, -130.0, 130.0);

	locHistName = "TrackTOFYVsX_TotalHit";
	locHistTitle = "TOF Total Hit;Projected TOF Hit X (cm);Projected TOF Hit Y (cm)";
	dHist_TrackTOFYVsX_TotalHit = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 130, -130.0, 130.0, 130, -130.0, 130.0);

	locHistName = "TrackTOF2DPaddles_HasHit";
	locHistTitle = "TOF Has Hit;Projected Vertical TOF Paddle;Projected Horizontal TOF Paddle";
	dHist_TrackTOF2DPaddles_HasHit = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 44, 0.5, 44.5, 44, 0.5, 44.5);

	locHistName = "TrackTOF2DPaddles_TotalHit";
	locHistTitle = "TOF Total Hit;Projected Vertical TOF Paddle;Projected Horizontal TOF Paddle";
	dHist_TrackTOF2DPaddles_TotalHit = new TH2I(locHistName.c_str(), locHistTitle.c_str(), 44, 0.5, 44.5, 44, 0.5, 44.5);
	gDirectory->cd("..");
	
	// back to original dir
	locOriginalDir->cd();

	//TTREE INTERFACE
	//MUST DELETE WHEN FINISHED: OR ELSE DATA WON'T BE SAVED!!!
	dTreeInterface = DTreeInterface::Create_DTreeInterface("tof_eff", "tree_tof_eff.root");

	//TTREE BRANCHES
	DTreeBranchRegister locTreeBranchRegister;

	//TRACK
	locTreeBranchRegister.Register_Single<Int_t>("PID_PDG"); //gives charge, mass, beta
	locTreeBranchRegister.Register_Single<Float_t>("TrackVertexZ");
	locTreeBranchRegister.Register_Single<TVector3>("TrackP3");
	locTreeBranchRegister.Register_Single<UInt_t>("TrackCDCRings"); //rings correspond to bits (1 -> 28)
	locTreeBranchRegister.Register_Single<UInt_t>("TrackFDCPlanes"); //planes correspond to bits (1 -> 24)

	//TOF
	locTreeBranchRegister.Register_Single<UChar_t>("NumTOFPoints");
	locTreeBranchRegister.Register_Single<Float_t>("ProjectedTOFX");
	locTreeBranchRegister.Register_Single<Float_t>("ProjectedTOFY");
	locTreeBranchRegister.Register_Single<UChar_t>("ProjectedTOFBar_Horizontal");
	locTreeBranchRegister.Register_Single<UChar_t>("ProjectedTOFBar_Vertical");

	//SEARCH TOF PADDLE
		//Nearest hit: //0 for none, 1 - 44 for both ends //101 - 144 for North/Top only, 201 - 244 for South/Bottom only (only = above threshold)
		//Nearest: No time cut (yet?)
	locTreeBranchRegister.Register_Single<UChar_t>("NearestTOFHit_Horizontal");
	locTreeBranchRegister.Register_Single<Float_t>("HorizontalTOFHitDeltaY");
	locTreeBranchRegister.Register_Single<UChar_t>("NearestTOFHit_Vertical");
	locTreeBranchRegister.Register_Single<Float_t>("VerticalTOFHitDeltaX");

	//SEARCH TOF POINT //nearest: must be in time: PID:OUT_OF_TIME_CUT
	locTreeBranchRegister.Register_Single<Float_t>("NearestTOFPointDeltaX");
	locTreeBranchRegister.Register_Single<Float_t>("NearestTOFPointDeltaY");
	locTreeBranchRegister.Register_Single<UShort_t>("NearestTOFPointStatus");
	/*
	NearestTOFPointStatus: horizontal_bar + 45*vertical_bar + 45*45*horizontal_status + 45*45*4*vertical_status
	*_bar = 0 -> 44 (0 for none (not matched to this plane))
	*_Status: 0 if no hit (not matched to this plane), 1 if only North/Top hit above threshold, 2 if only South/Bottom hit above threshold, 3 if both hits above threshold
	Note that if it's a single-ended bar, the status cannot be 3. 
	*/
	locTreeBranchRegister.Register_Single<Bool_t>("IsMatchedToTrack"); //false if not registered in DDetectorMatches

	//FOUND TOF POINT //only if matched: for evaluating PID quality
	locTreeBranchRegister.Register_Single<Float_t>("TOFPointDeltaT"); //TOF - RF
	locTreeBranchRegister.Register_Single<Float_t>("TOFPointTimeFOM");
	locTreeBranchRegister.Register_Single<Float_t>("TOFPointdEdX");

	//REGISTER BRANCHES
	dTreeInterface->Create_Branches(locTreeBranchRegister);

	return NOERROR;
}


//------------------
// brun
//------------------
jerror_t JEventProcessor_TOF_Eff::brun(jana::JEventLoop* locEventLoop, int locRunNumber)
{
	// This is called whenever the run number changes

	return NOERROR;
}

//------------------
// evnt
//------------------

jerror_t JEventProcessor_TOF_Eff::evnt(jana::JEventLoop* locEventLoop, uint64_t locEventNumber)
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

	const DEventRFBunch* locEventRFBunch = NULL;
	locEventLoop->GetSingle(locEventRFBunch);
	if(locEventRFBunch->dNumParticleVotes <= 1)
		return NOERROR; //don't trust PID: beta-dependence

	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks);

	const DParticleID* locParticleID = NULL;
	locEventLoop->GetSingle(locParticleID);

	vector<const DTOFPoint*> locTOFPoints;
	locEventLoop->Get(locTOFPoints);

	vector<const DTOFPaddleHit*> locTOFPaddleHits;
	locEventLoop->Get(locTOFPaddleHits);

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

			//Need PID for beta-dependence, but cannot use TOF info: would bias
			if(!Cut_FCALTiming(locChargedTrackHypothesis, locParticleID, locEventRFBunch))
				continue; //also requires match to FCAL: no need for separate check

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

		//Predict TOF Surface Hit Location
		DVector3 locProjectedTOFIntersection;
		unsigned int locProjectedHorizontalBar = 0, locProjectedVerticalBar = 0;
		if(!locParticleID->PredictTOFPaddles(locTrackTimeBased->rt, locProjectedHorizontalBar, locProjectedVerticalBar, &locProjectedTOFIntersection))
		{
			if(locTrackTimeBased->momentum().Theta()*180.0/TMath::Pi() > dMaxTOFThetaCut)
				continue; //not predicted to hit TOF
		}

		//Find closest TOF point
		double locBestPointDeltaX = 999.9, locBestPointDeltaY = 999.9;
		const DTOFPoint* locClosestTOFPoint = locParticleID->Get_ClosestToTrack_TOFPoint(locTrackTimeBased, locTOFPoints, locBestPointDeltaX, locBestPointDeltaY);

		//Find closest TOF paddles
		double locBestPaddleDeltaX = 999.9, locBestPaddleDeltaY = 999.9;
		//first in pair is vertical, second is horizontal // NULL If none / doesn't hit TOF
		double locStartTime = locParticleID->Calc_PropagatedRFTime(locChargedTrackHypothesis, locEventRFBunch);
		const DTOFPaddleHit* locClosestTOFPaddleHit_Vertical = locParticleID->Get_ClosestTOFPaddleHit_Vertical(locTrackTimeBased->rt, locTOFPaddleHits, locStartTime, locBestPaddleDeltaX);
		const DTOFPaddleHit* locClosestTOFPaddleHit_Horizontal = locParticleID->Get_ClosestTOFPaddleHit_Horizontal(locTrackTimeBased->rt, locTOFPaddleHits, locStartTime, locBestPaddleDeltaY);

		//Is match to TOF point?
		const DTOFHitMatchParams* locTOFHitMatchParams = locChargedTrackHypothesis->Get_TOFHitMatchParams();
		bool locIsMatchedToTrack = (locTOFHitMatchParams != nullptr);

		//If so, calc PID info: timing, dE/dx
		double locTOFPointDeltaT = 0.0;
		double locTOFTimeFOM = Calc_TOFTiming(locChargedTrackHypothesis, locParticleID, locEventRFBunch, locTOFPointDeltaT);
		double locTOFPointdEdX = locIsMatchedToTrack ? locTOFHitMatchParams->dEdx : 0.0;

		//calc nearest tof hit/point status, bars
		int locNearestTOFPointStatus = 0;
		if(locClosestTOFPoint != NULL)
		{
			locNearestTOFPointStatus = locClosestTOFPoint->dHorizontalBar + 45*locClosestTOFPoint->dVerticalBar;
			locNearestTOFPointStatus += 45*45*locClosestTOFPoint->dHorizontalBarStatus + 45*45*4*locClosestTOFPoint->dVerticalBarStatus;
		}

		int locNearestTOFHitHorizontal = Calc_NearestHit(locClosestTOFPaddleHit_Horizontal);
		int locNearestTOFHitVertical = Calc_NearestHit(locClosestTOFPaddleHit_Vertical);


		// FILL HISTOGRAMS
		// Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
		japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
		{
			//TOFPaddle
			dHist_TOFPaddleTrackYVsVerticalPaddle_TotalHit_Top->Fill(locProjectedVerticalBar, locProjectedTOFIntersection.Y());
			dHist_TOFPaddleHorizontalPaddleVsTrackX_TotalHit_North->Fill(locProjectedTOFIntersection.X(), locProjectedHorizontalBar);
			dHist_TOFPaddleTrackYVsVerticalPaddle_TotalHit_Bottom->Fill(locProjectedVerticalBar, locProjectedTOFIntersection.Y());
			dHist_TOFPaddleHorizontalPaddleVsTrackX_TotalHit_South->Fill(locProjectedTOFIntersection.X(), locProjectedHorizontalBar);

			//TOFPoint
			dHist_TrackTOFYVsX_TotalHit->Fill(locProjectedTOFIntersection.X(), locProjectedTOFIntersection.Y());
			dHist_TrackTOF2DPaddles_TotalHit->Fill(locProjectedVerticalBar, locProjectedHorizontalBar);

			//TOFPaddle: Horizontal
			if(locBestPaddleDeltaY <= dMinTOFPaddleMatchDistance) //horizontal match
			{
				dHist_TOFPaddleHorizontalPaddleVsTrackX_HasHit_North->Fill(locProjectedTOFIntersection.X(), locProjectedHorizontalBar);
				dHist_TOFPaddleHorizontalPaddleVsTrackX_HasHit_South->Fill(locProjectedTOFIntersection.X(), locProjectedHorizontalBar);
			}

			//TOFPaddle: Vertical
			if(locBestPaddleDeltaX <= dMinTOFPaddleMatchDistance) //vertical match
			{
				dHist_TOFPaddleTrackYVsVerticalPaddle_HasHit_Top->Fill(locProjectedVerticalBar, locProjectedTOFIntersection.Y());
				dHist_TOFPaddleTrackYVsVerticalPaddle_HasHit_Bottom->Fill(locProjectedVerticalBar, locProjectedTOFIntersection.Y());
			}

			//TOFPoint
			if(locIsMatchedToTrack)
			{
				dHist_TrackTOFYVsX_HasHit->Fill(locProjectedTOFIntersection.X(), locProjectedTOFIntersection.Y());
				dHist_TrackTOF2DPaddles_HasHit->Fill(locProjectedVerticalBar, locProjectedHorizontalBar);
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

		//TOF
		dTreeFillData.Fill_Single<UChar_t>("NumTOFPoints", locTOFPoints.size());
		dTreeFillData.Fill_Single<Float_t>("ProjectedTOFX", locProjectedTOFIntersection.X());
		dTreeFillData.Fill_Single<Float_t>("ProjectedTOFY", locProjectedTOFIntersection.Y());
		dTreeFillData.Fill_Single<UChar_t>("ProjectedTOFBar_Horizontal", locProjectedHorizontalBar);
		dTreeFillData.Fill_Single<UChar_t>("ProjectedTOFBar_Vertical", locProjectedVerticalBar);

		//SEARCH TOF PADDLE
		dTreeFillData.Fill_Single<UChar_t>("NearestTOFHit_Horizontal", locNearestTOFHitHorizontal);
		dTreeFillData.Fill_Single<Float_t>("HorizontalTOFHitDeltaY", locBestPaddleDeltaY);
		dTreeFillData.Fill_Single<UChar_t>("NearestTOFHit_Vertical", locNearestTOFHitVertical);
		dTreeFillData.Fill_Single<Float_t>("VerticalTOFHitDeltaX", locBestPaddleDeltaX);

		//SEARCH TOF POINT
		dTreeFillData.Fill_Single<Float_t>("NearestTOFPointDeltaX", locBestPointDeltaX);
		dTreeFillData.Fill_Single<Float_t>("NearestTOFPointDeltaY", locBestPointDeltaY);
		dTreeFillData.Fill_Single<UShort_t>("NearestTOFPointStatus", locNearestTOFPointStatus);
		dTreeFillData.Fill_Single<Bool_t>("IsMatchedToTrack", locIsMatchedToTrack);

		//FOUND TOF POINT
		dTreeFillData.Fill_Single<Float_t>("TOFPointDeltaT", locTOFPointDeltaT);
		dTreeFillData.Fill_Single<Float_t>("TOFPointTimeFOM", locTOFTimeFOM);
		dTreeFillData.Fill_Single<Float_t>("TOFPointdEdX", locTOFPointdEdX);

		//FILL TTREE
		dTreeInterface->Fill(dTreeFillData);
	}

	return NOERROR;
}

int JEventProcessor_TOF_Eff::Calc_NearestHit(const DTOFPaddleHit* locPaddleHit) const
{
	//Nearest hit: //0 for none, 1 - 44 for both ends //101 - 144 for North/Top only, 201 - 244 for South/Bottom only (only = above threshold)
	if(locPaddleHit == nullptr)
		return 0;

	int locNearestHit = locPaddleHit->bar;
	if((locPaddleHit->E_north > 0.00001) && !(locPaddleHit->E_south > 0.00001)) //0.00001: tolerance
		locNearestHit += 100.0;
	else if(!(locPaddleHit->E_north > 0.00001) && (locPaddleHit->E_south > 0.00001))
		locNearestHit += 200.0;

	return locNearestHit;
}

bool JEventProcessor_TOF_Eff::Cut_FCALTiming(const DChargedTrackHypothesis* locChargedTrackHypothesis, const DParticleID* locParticleID, const DEventRFBunch* locEventRFBunch)
{
	const DFCALShowerMatchParams* locFCALShowerMatchParams = locChargedTrackHypothesis->Get_FCALShowerMatchParams();
	if(locFCALShowerMatchParams == NULL)
		return false;

	double locStartTime = locParticleID->Calc_PropagatedRFTime(locChargedTrackHypothesis, locEventRFBunch);
	const DFCALShower* locFCALShower = locFCALShowerMatchParams->dFCALShower;

	double locDeltaT = locFCALShower->getTime() - locFCALShowerMatchParams->dFlightTime - locStartTime;
	return (fabs(locDeltaT) <= dMaxFCALDeltaT);
}

double JEventProcessor_TOF_Eff::Calc_TOFTiming(const DChargedTrackHypothesis* locChargedTrackHypothesis, const DParticleID* locParticleID, const DEventRFBunch* locEventRFBunch, double& locDeltaT)
{
	const DTOFHitMatchParams* locTOFHitMatchParams = locChargedTrackHypothesis->Get_TOFHitMatchParams();
	if(locTOFHitMatchParams == nullptr)
		return -1.0;

	double locStartTime = locParticleID->Calc_PropagatedRFTime(locChargedTrackHypothesis, locEventRFBunch);
	locDeltaT = locTOFHitMatchParams->dHitTime - locTOFHitMatchParams->dFlightTime - locStartTime;

	double locPIDFOM = -1.0; //not able to calc this correctly yet
	return locPIDFOM;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_TOF_Eff::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_TOF_Eff::fini(void)
{
	// Called before program exit after event processing is finished.  

	delete dCutAction_TrackHitPattern;
	delete dTreeInterface; //saves trees to file, closes file

	return NOERROR;
}

