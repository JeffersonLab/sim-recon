// $Id$
//
//    File: JEventProcessor_BCAL_Hadronic_Eff.cc
//

#include "JEventProcessor_BCAL_Hadronic_Eff.h"

// Routine used to create our JEventProcessor
extern "C"
{
	void InitPlugin(JApplication *locApplication)
	{
		InitJANAPlugin(locApplication);
		locApplication->AddProcessor(new JEventProcessor_BCAL_Hadronic_Eff()); //register this plugin
	}
} // "C"

//define static local variable //declared in header file
thread_local DTreeFillData JEventProcessor_BCAL_Hadronic_Eff::dTreeFillData;

//------------------
// init
//------------------
jerror_t JEventProcessor_BCAL_Hadronic_Eff::init(void)
{
	//TRACK REQUIREMENTS
	dMinPIDFOM = 5.73303E-7; //+/- 5 sigma
	dMinTrackingFOM = 5.73303E-7; // +/- 5 sigma
	dMinNumTrackHits = 14; //e.g. 6 in CDC, 8 in 
	dMinHitRingsPerCDCSuperlayer = 3;
	dMinHitPlanesPerFDCPackage = 4;
	dCutAction_TrackHitPattern = new DCutAction_TrackHitPattern(NULL, dMinHitRingsPerCDCSuperlayer, dMinHitPlanesPerFDCPackage);
	//action initialize not necessary: is empty

	TDirectory* locOriginalDir = gDirectory;
	gDirectory->mkdir("bcal_hadronic_eff")->cd();

	dHist_HadronicShowerMatched = new TH2I("HadronicShowerMatched", "Hadronic Shower Matched;Projected BCAL Hit Z (cm);Projected BCAL Hit #phi#circ;", 225, 0.0, 450.0, 180, -180.0, 180.0);
	dHist_HadronicShowerTotal = new TH2I("HadronicShowerTotal", "Hadronic Shower Total;Projected BCAL Hit Z (cm);Projected BCAL Hit #phi#circ;", 225, 0.0, 450.0, 180, -180.0, 180.0);

	dHistFoundDeltaSector = 4;
	for(int locLayer = 1; locLayer <= 4; ++locLayer)
	{
		ostringstream locHistName, locHistTitle;

		//Upstream, Found
		locHistName << "HitFound_Layer" << locLayer << "_Upstream";
		locHistTitle << "Hit Found, Layer " << locLayer << ", Upstream;Sector";
		dHistMap_HitFound[locLayer][true] = new TH1I(locHistName.str().c_str(), locHistTitle.str().c_str(), 192, 0.5, 192.5);

		//Upstream, Total
		locHistName.str("");
		locHistName << "HitTotal_Layer" << locLayer << "_Upstream";
		locHistTitle.str("");
		locHistTitle << "Hit Total, Layer " << locLayer << ", Upstream;Sector";
		dHistMap_HitTotal[locLayer][true] = new TH1I(locHistName.str().c_str(), locHistTitle.str().c_str(), 192, 0.5, 192.5);

		//Downstream, Found
		locHistName.str("");
		locHistName << "HitFound_Layer" << locLayer << "_Downstream";
		locHistTitle.str("");
		locHistTitle << "Hit Found, Layer " << locLayer << ", Downstream;Sector";
		dHistMap_HitFound[locLayer][false] = new TH1I(locHistName.str().c_str(), locHistTitle.str().c_str(), 192, 0.5, 192.5);

		//Downstream, Total
		locHistName.str("");
		locHistName << "HitTotal_Layer" << locLayer << "_Downstream";
		locHistTitle.str("");
		locHistTitle << "Hit Total, Layer " << locLayer << ", Downstream;Sector";
		dHistMap_HitTotal[locLayer][false] = new TH1I(locHistName.str().c_str(), locHistTitle.str().c_str(), 192, 0.5, 192.5);
	}
	
	// back to original dir
	locOriginalDir->cd();

	//TTREE INTERFACE
	//MUST DELETE WHEN FINISHED: OR ELSE DATA WON'T BE SAVED!!!
	dTreeInterface = DTreeInterface::Create_DTreeInterface("bcal_hadronic_eff", "tree_bcal_hadronic_eff.root");

	//TTREE BRANCHES
	DTreeBranchRegister locTreeBranchRegister;

	//TRACK
	locTreeBranchRegister.Register_Single<Int_t>("PID_PDG"); //gives charge, mass, beta
	locTreeBranchRegister.Register_Single<Float_t>("TimingBeta");
	locTreeBranchRegister.Register_Single<Float_t>("TrackVertexZ");
	locTreeBranchRegister.Register_Single<TVector3>("TrackP3");
	locTreeBranchRegister.Register_Single<UInt_t>("TrackCDCRings"); //rings correspond to bits (1 -> 28)
	locTreeBranchRegister.Register_Single<UInt_t>("TrackFDCPlanes"); //planes correspond to bits (1 -> 24)

	//SHOWER
	locTreeBranchRegister.Register_Single<Float_t>("NearestShowerEnergy"); //is zero if none
	locTreeBranchRegister.Register_Single<Float_t>("TrackDeltaPhiToShower"); //is signed: BCAL - Track
	locTreeBranchRegister.Register_Single<Float_t>("TrackDeltaZToShower"); //is signed: BCAL - Track
	locTreeBranchRegister.Register_Single<Float_t>("ProjectedBCALHitPhi"); //degrees
	locTreeBranchRegister.Register_Single<Float_t>("ProjectedBCALHitZ");
	locTreeBranchRegister.Register_Single<Bool_t>("IsMatchedToTrack"); //false if not registered in DDetectorMatches

	//HIT SEARCH
	//BCALClusterLayers: first 4 bits: point layers, next 4: unmatched-unified-hit layers
	locTreeBranchRegister.Register_Single<UChar_t>("BCALClusterLayers"); //is 0 if track not matched to shower: ignore hit efficiencies
	//IsHitInCluster: bits 1 -> 8 correspond to Upstream layer 1 -> 4, then Downstream layer 1 -> 4
	locTreeBranchRegister.Register_Single<UChar_t>("IsHitInCluster"); //1 if true, 0 if false or no hit
	//LAYER 1:
	//"Sector:" 4*(module - 1) + sector //sector: 1 -> 192
	locTreeBranchRegister.Register_Single<UChar_t>("ProjectedBCALSectors_Layer1"); //0 if biased or indeterminate
	locTreeBranchRegister.Register_Single<UChar_t>("NearestBCALSectors_Layer1_Downstream"); //0 if not found
	locTreeBranchRegister.Register_Single<UChar_t>("NearestBCALSectors_Layer1_Upstream"); //0 if not found
	//LAYER 2:
	locTreeBranchRegister.Register_Single<UChar_t>("ProjectedBCALSectors_Layer2"); //0 if biased or indeterminate
	locTreeBranchRegister.Register_Single<UChar_t>("NearestBCALSectors_Layer2_Downstream"); //0 if not found
	locTreeBranchRegister.Register_Single<UChar_t>("NearestBCALSectors_Layer2_Upstream"); //0 if not found
	//LAYER 3:
	locTreeBranchRegister.Register_Single<UChar_t>("ProjectedBCALSectors_Layer3"); //0 if biased or indeterminate
	locTreeBranchRegister.Register_Single<UChar_t>("NearestBCALSectors_Layer3_Downstream"); //0 if not found
	locTreeBranchRegister.Register_Single<UChar_t>("NearestBCALSectors_Layer3_Upstream"); //0 if not found
	//LAYER 4:
	locTreeBranchRegister.Register_Single<UChar_t>("ProjectedBCALSectors_Layer4"); //0 if biased or indeterminate
	locTreeBranchRegister.Register_Single<UChar_t>("NearestBCALSectors_Layer4_Downstream"); //0 if not found
	locTreeBranchRegister.Register_Single<UChar_t>("NearestBCALSectors_Layer4_Upstream"); //0 if not found

	//REGISTER BRANCHES
	dTreeInterface->Create_Branches(locTreeBranchRegister);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_BCAL_Hadronic_Eff::brun(jana::JEventLoop* locEventLoop, int locRunNumber)
{
	// This is called whenever the run number changes

	//Get Target Center Z
	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry* locGeometry = locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());
	locGeometry->GetTargetZ(dTargetCenterZ); //thread_local: each thread has own copy!

	//Effective velocities
	locEventLoop->GetCalib("/BCAL/effective_velocities", effective_velocities);

	//THE WORST THING EVER.  FIX THIS GEOMETRY CLASS.  NOTHING IN IT SHOULD BE STATIC.
	DBCALGeometry::Initialize(locRunNumber);

	return NOERROR;
}

//------------------
// evnt
//------------------

jerror_t JEventProcessor_BCAL_Hadronic_Eff::evnt(jana::JEventLoop* locEventLoop, uint64_t locEventNumber)
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

	//COMPUTE TOTAL FCAL ENERGY
	vector<const DFCALShower*> locFCALShowers;
	locEventLoop->Get(locFCALShowers);

	double locTotalFCALEnergy = 0.0;
	for(auto& locShower : locFCALShowers)
		locTotalFCALEnergy += locShower->getEnergy();

	//SEE IF BCAL REQUIRED TO TRIGGER
	uint32_t locTriggerBits = locTrigger->Get_L1TriggerBits();
	bool locBCALRequiredForTriggerFlag = true;
	if(((locTriggerBits & 2) == 2) || ((locTriggerBits & 32) == 32) || ((locTriggerBits & 64) == 64))
		locBCALRequiredForTriggerFlag = false;
	if(((locTriggerBits & 1) == 1) && (locTotalFCALEnergy > 1.0))
		locBCALRequiredForTriggerFlag = false;
	if(locBCALRequiredForTriggerFlag)
		return NOERROR;

	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks);

	const DParticleID* locParticleID = NULL;
	locEventLoop->GetSingle(locParticleID);

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	vector<const DBCALUnifiedHit*> locBCALUnifiedHits;
	locEventLoop->Get(locBCALUnifiedHits);

	vector<const DBCALPoint*> locBCALPoints;
	locEventLoop->Get(locBCALPoints);

	vector<const DBCALShower*> locBCALShowers;
	locEventLoop->Get(locBCALShowers);

	//sort bcal points by layer, total sector //first int: layer. second int: sector
	map<int, map<int, set<const DBCALPoint*> > > locSortedPoints;
	for(const auto& locPoint : locBCALPoints)
	{
		int locSector = (locPoint->module() - 1)*4 + locPoint->sector(); //1 -> 192
		locSortedPoints[locPoint->layer()][locSector].insert(locPoint);
	}

	//sort unified hits by layer, total sector //first int: layer. second int: sector
	map<int, map<int, set<const DBCALUnifiedHit*> > > locSortedHits_Upstream, locSortedHits_Downstream;
	for(const auto& locUnifiedHit : locBCALUnifiedHits)
	{
		int locSector = (locUnifiedHit->module - 1)*4 + locUnifiedHit->sector;
		if(locUnifiedHit->end == DBCALGeometry::kUpstream)
			locSortedHits_Upstream[locUnifiedHit->layer][locSector].insert(locUnifiedHit);
		else
			locSortedHits_Downstream[locUnifiedHit->layer][locSector].insert(locUnifiedHit);
	}

	//Try to select the most-pure sample of tracks possible
	//select the best DTrackTimeBased for each track: of tracks with good hit pattern, use best PID confidence level (need accurate beta)
		//also, must have min tracking FOM, min PID confidence level, min #-hits on track, and matching hit in SC
	set<const DChargedTrackHypothesis*> locBestTracks; //lowest tracking FOM for each candidate id
	for(size_t loc_i = 0; loc_i < locChargedTracks.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTracks[loc_i]->Get_BestFOM();
		if(locChargedTrackHypothesis->dFOM < dMinPIDFOM)
			continue; //don't trust PID

		const DTrackTimeBased* locTrackTimeBased = NULL;
		locChargedTrackHypothesis->GetSingle(locTrackTimeBased);
		if(locTrackTimeBased->FOM < dMinTrackingFOM)
			continue;

		if(!locDetectorMatches->Get_IsMatchedToDetector(locTrackTimeBased, SYS_START))
			continue; //not matched to SC

		if(!dCutAction_TrackHitPattern->Cut_TrackHitPattern(locParticleID, locTrackTimeBased))
			continue;

		unsigned int locNumTrackHits = locTrackTimeBased->Ndof + 5;
		if(locNumTrackHits < dMinNumTrackHits)
			continue;

		locBestTracks.insert(locChargedTrackHypothesis);
	}

	//for histograms, keep running total of what to fill //int = layer, bool = isUpstream
	//will fill at end: only one lock
	map<int, map<bool, vector<int> > > locHitMap_HitFound, locHitMap_HitTotal;
	vector<pair<double, double> > locShowerVector_ShowerFound, locShowerVector_ShowerTotal; //first is proj-z, second is proj-phi

	// Loop over the good tracks, using the best DTrackTimeBased object for each
	for(auto& locChargedTrackHypothesis : locBestTracks)
	{
		const DTrackTimeBased* locTrackTimeBased = NULL;
		locChargedTrackHypothesis->GetSingle(locTrackTimeBased);

		//Predict BCAL Surface Hit Location
		unsigned int locPredictedSurfaceModule = 0, locPredictedSurfaceSector = 0;
		DVector3 locPredictedSurfacePosition;
		locParticleID->PredictBCALWedge(locTrackTimeBased->rt, locPredictedSurfaceModule, locPredictedSurfaceSector, &locPredictedSurfacePosition);

		pair<double, double> locShowerPair(locPredictedSurfacePosition.Z(), locPredictedSurfacePosition.Phi()*180.0/TMath::Pi());
		locShowerVector_ShowerTotal.push_back(locShowerPair);

		/************************************************ CHECK SHOWER MATCH EFFICIENCY ************************************************/

		//get the best-matched DBCALShower for this track (if any)
		const DBCALShowerMatchParams* locBCALShowerMatchParams = locChargedTrackHypothesis->Get_BCALShowerMatchParams();
		if(locBCALShowerMatchParams == NULL)
		{
			Fill_NoClusterStudy(locChargedTrackHypothesis, locBCALShowers, locParticleID, false); //no match: fill for shower efficiency
			continue; //don't compute hit efficiencies
		}

		locShowerVector_ShowerFound.push_back(locShowerPair);

		/*************************************************** GET & SORT CLUSTER HITS ***************************************************/

		//use the DBCALPoints in this DBCALShower to help define where DBCALUnifiedHits are expected to be
		//The DBCALShower reconstruction does not care which BCAL layers fired
			//Therefore, given a DBCALShower hit layer, whether or not there are hits in the other layers is an unbiased check

		//get the clusters in the shower, reject if > 1
		vector<const DBCALCluster*> locClusters;
		locBCALShowerMatchParams->dBCALShower->Get(locClusters);
		if(locClusters.size() > 1)
		{
			//likely a messy shower, will probably mess up efficiency calculation. 
			Fill_NoClusterStudy(locChargedTrackHypothesis, locBCALShowers, locParticleID, true); //can't use for hits: fill for shower efficiency
			continue; //don't compute hit efficiencies
		}

		//get cluster, points, unmatched hits
		const DBCALCluster* locCluster = locClusters[0];
		vector<const DBCALPoint*> locClusterBCALPoints = locCluster->points();
		vector<pair<const DBCALUnifiedHit*, double> > locClusterBCALHits = locCluster->hits(); //double: corrected energy

		//sort cluster points by layer, total sector //first int: layer. second int: sector
		//also, build set of all cluster points
		set<const DBCALPoint*> locClusterPointSet;
		map<int, map<int, set<const DBCALPoint*> > > locSortedClusterPoints;
		for(const auto& locPoint : locClusterBCALPoints)
		{
			locClusterPointSet.insert(locPoint);
			int locSector = (locPoint->module() - 1)*4 + locPoint->sector(); //1 -> 192
			locSortedClusterPoints[locPoint->layer()][locSector].insert(locPoint);
		}

		//require points in at least 2 layers: make sure it's not a noise cluster
		//Note: This can introduce bias into the study. To avoid bias, only evaluate layer efficiency if at least 2 OTHER layers have hits
		if(locSortedClusterPoints.size() < 2)
		{
			Fill_NoClusterStudy(locChargedTrackHypothesis, locBCALShowers, locParticleID, true); //can't use for hits: fill for shower efficiency
			continue; //don't compute hit efficiencies
		}

		//sort cluster unified hits by layer, total sector //first int: layer. second int: sector
		//also, build set of all cluster unmatched unified hits
		set<const DBCALUnifiedHit*> locClusterUnifiedHitSet;
		map<int, map<int, set<const DBCALUnifiedHit*> > > locSortedClusterHits_Upstream, locSortedClusterHits_Downstream;
		for(const auto& locUnifiedHitPair : locClusterBCALHits)
		{
			const DBCALUnifiedHit* locUnifiedHit = locUnifiedHitPair.first;
			locClusterUnifiedHitSet.insert(locUnifiedHit);

			int locSector = (locUnifiedHit->module - 1)*4 + locUnifiedHit->sector;
			if(locUnifiedHit->end == DBCALGeometry::kUpstream)
				locSortedClusterHits_Upstream[locUnifiedHit->layer][locSector].insert(locUnifiedHit);
			else
				locSortedClusterHits_Downstream[locUnifiedHit->layer][locSector].insert(locUnifiedHit);
		}

		//register layers
		UChar_t locClusterLayers = 0; //which layers are in the cluster (both by points & by hits): 4bits each: 8total
		for(auto& locLayerPair : locSortedClusterPoints)
			locClusterLayers |= (1 << (locLayerPair.first - 1));
		for(auto& locLayerPair : locSortedClusterHits_Upstream)
			locClusterLayers |= (1 << (locLayerPair.first + 4 - 1));
		for(auto& locLayerPair : locSortedClusterHits_Downstream)
			locClusterLayers |= (1 << (locLayerPair.first + 4 - 1));

		/**************************************************** LOOP OVER BCAL LAYERS ****************************************************/

		//Tree-save variables
		map<int, UChar_t> locProjectedSectorsMap, locNearestSectorsMap_Downstream, locNearestSectorsMap_Upstream;
		//locProjectedSectors: (rounded from search): 4*(module - 1) + sector //sector: 1 -> 192, 0 if biased or indeterminate
		//locNearestSectors_Downstream & locNearestSectors_Upstream: 4*(module - 1) + sector //sector: 1 -> 192, 0 if not found
		UChar_t locIsHitInClusterBits = 0; //bits 1 -> 8 correspond to Upstream layer 1 -> 4, then Downstream layer 1 -> 4

		//loop over layers
		for(int locLayer = 1; locLayer <= 4; ++locLayer)
		{
			//initialize to 0's
			locProjectedSectorsMap[locLayer] = 0;
			locNearestSectorsMap_Downstream[locLayer] = 0;
			locNearestSectorsMap_Upstream[locLayer] = 0;

			//require at least 2 other layers in the cluster to have points
			//use the DBCALPoint's (double-ended hits) to define where hits are to be expected (not the single-ended ones!)
			//if no hit in this layer, then requirement is already met. if there is, must check
			auto locClusterPointsIterator = locSortedClusterPoints.find(locLayer);
			if((locClusterPointsIterator != locSortedClusterPoints.end()) && (locSortedClusterPoints.size() <= 2))
				continue; //not at least 2 hits in other layers

			//ideally, would require hits in adjacent layers, but cannot since one may be dead: won't be able to evaluate for this layer
				//can always constrain this later though, if desired

			//ok, can now evaluate whether or not there are hits in this layer in an unbiased manner

			/******************************************** PROJECT HIT LOCATION ********************************************/

			//in this layer, need to know where to search for a hit:
				//showers may contain hits in multiple sectors within a layer: 
				//cannot evaluate all found as "good," because when missing, don't know how many to mark as "bad"
				//so, pick the most-probable location: energy-weighted average of hit locations in the nearest layers

			//locProjectedSector is (module - 1)*4 + sector //double: average
			double locProjectedSector = Calc_ProjectedSector(locLayer, locSortedClusterPoints);
			int locProjectedSectorInt = int(locProjectedSector + 0.5);

			//register for tree & hists:
			locProjectedSectorsMap[locLayer] = UChar_t(locProjectedSectorInt); //tree
			locHitMap_HitTotal[locLayer][true].push_back(locProjectedSectorInt); //hist
			locHitMap_HitTotal[locLayer][false].push_back(locProjectedSectorInt); //hist

			/*************************************************** SEARCH ***************************************************/

			//now, find the closest hits in the layer, separately for in-cluster and not-in-cluster
				//time cuts copied from DBCALCluster_factory constructor

			//first, search ALL points
			auto locPointsIterator = locSortedPoints.find(locLayer);
			pair<const DBCALPoint*, double> locNearestPoint(NULL, 999.0);
			if(locPointsIterator != locSortedPoints.end())
				locNearestPoint = Find_NearestPoint(locProjectedSector, locPointsIterator->second, locCluster, 8.0);

			//next, search ALL hits: Upstream
			pair<const DBCALUnifiedHit*, double> locNearestHit_Upstream(NULL, 999.0);
			auto locHitsIterator = locSortedHits_Upstream.find(locLayer);
			if(locHitsIterator != locSortedHits_Upstream.end())
				locNearestHit_Upstream = Find_NearestHit(locProjectedSector, locHitsIterator->second, locCluster, 20.0);

			//finally, search ALL hits: Downstream
			pair<const DBCALUnifiedHit*, double> locNearestHit_Downstream(NULL, 999.0);
			locHitsIterator = locSortedHits_Downstream.find(locLayer);
			if(locHitsIterator != locSortedHits_Downstream.end())
				locNearestHit_Downstream = Find_NearestHit(locProjectedSector, locHitsIterator->second, locCluster, 20.0);

			/********************************************** REGISTER RESULTS **********************************************/

			//UPSTREAM
			if((locNearestPoint.first != nullptr) || (locNearestHit_Upstream.first != nullptr))
			{
				bool locUsePointFlag = (locNearestPoint.second <= locNearestHit_Upstream.second);

				int locFoundSector = 0;
				bool locIsInClusterFlag = false;
				if(locUsePointFlag)
				{
					locFoundSector = (locNearestPoint.first->module() - 1)*4 + locNearestPoint.first->sector();
					locIsInClusterFlag = (locClusterPointSet.find(locNearestPoint.first) != locClusterPointSet.end());
				}
				else
				{
					locFoundSector = (locNearestHit_Upstream.first->module - 1)*4 + locNearestHit_Upstream.first->sector;
					locIsInClusterFlag = (locClusterUnifiedHitSet.find(locNearestHit_Upstream.first) != locClusterUnifiedHitSet.end());
				}

				locNearestSectorsMap_Upstream[locLayer] = UChar_t(locFoundSector);
				if(locIsInClusterFlag)
					locIsHitInClusterBits |= (1 << (locLayer - 1)); //Upstream bits: 1 -> 4

				int locDeltaSector = Calc_DeltaSector<int>(locFoundSector, locProjectedSectorInt);
				if(abs(locDeltaSector) <= dHistFoundDeltaSector)
					locHitMap_HitFound[locLayer][true].push_back(locProjectedSectorInt); //true: upstream
			}

			//DOWNSTREAM
			if((locNearestPoint.first != nullptr) || (locNearestHit_Downstream.first != nullptr))
			{
				bool locUsePointFlag = (locNearestPoint.second <= locNearestHit_Downstream.second);

				int locFoundSector = 0;
				bool locIsInClusterFlag = false;
				if(locUsePointFlag)
				{
					locFoundSector = (locNearestPoint.first->module() - 1)*4 + locNearestPoint.first->sector();
					locIsInClusterFlag = (locClusterPointSet.find(locNearestPoint.first) != locClusterPointSet.end());
				}
				else
				{
					locFoundSector = (locNearestHit_Downstream.first->module - 1)*4 + locNearestHit_Downstream.first->sector;
					locIsInClusterFlag = (locClusterUnifiedHitSet.find(locNearestHit_Downstream.first) != locClusterUnifiedHitSet.end());
				}

				locNearestSectorsMap_Downstream[locLayer] = UChar_t(locFoundSector);
				if(locIsInClusterFlag)
					locIsHitInClusterBits |= (1 << (locLayer - 1 + 4)); //Downstream bits: 5 -> 8

				int locDeltaSector = Calc_DeltaSector<int>(locFoundSector, locProjectedSectorInt);
				if(abs(locDeltaSector) <= dHistFoundDeltaSector)
					locHitMap_HitFound[locLayer][false].push_back(locProjectedSectorInt); //false: downstream
			}
		}

		//STAGE DATA FOR TREE FILL

		//TRACK
		dTreeFillData.Fill_Single<Int_t>("PID_PDG", PDGtype(locChargedTrackHypothesis->PID()));
		dTreeFillData.Fill_Single<Float_t>("TimingBeta", locChargedTrackHypothesis->measuredBeta());
		dTreeFillData.Fill_Single<Float_t>("TrackVertexZ", locChargedTrackHypothesis->position().Z());
		DVector3 locDP3 = locChargedTrackHypothesis->momentum();
		TVector3 locP3(locDP3.X(), locDP3.Y(), locDP3.Z());
		dTreeFillData.Fill_Single<TVector3>("TrackP3", locP3);
		dTreeFillData.Fill_Single<UInt_t>("TrackCDCRings", locTrackTimeBased->dCDCRings);
		dTreeFillData.Fill_Single<UInt_t>("TrackFDCPlanes", locTrackTimeBased->dFDCPlanes);

		//SHOWER
		dTreeFillData.Fill_Single<Float_t>("NearestShowerEnergy", locBCALShowerMatchParams->dBCALShower->E);
		dTreeFillData.Fill_Single<Float_t>("TrackDeltaPhiToShower", locBCALShowerMatchParams->dDeltaPhiToShower); //is signed: BCAL - Track
		dTreeFillData.Fill_Single<Float_t>("TrackDeltaZToShower", locBCALShowerMatchParams->dDeltaZToShower); //is signed: BCAL - Track
		dTreeFillData.Fill_Single<Float_t>("ProjectedBCALHitPhi", locPredictedSurfacePosition.Phi()*180.0/TMath::Pi());
		dTreeFillData.Fill_Single<Float_t>("ProjectedBCALHitZ", locPredictedSurfacePosition.Z());
		dTreeFillData.Fill_Single<Bool_t>("IsMatchedToTrack", true);

		//HIT SEARCH
		//BCALClusterLayers: first 4 bits: point layers, next 4: unmatched-unified-hit layers
		dTreeFillData.Fill_Single<UChar_t>("BCALClusterLayers", locClusterLayers);
		dTreeFillData.Fill_Single<UChar_t>("IsHitInCluster", locIsHitInClusterBits);
		//Sector Branches: 4*(module - 1) + sector //sector: 1 -> 192
		//LAYER 1:
		dTreeFillData.Fill_Single<UInt_t>("ProjectedBCALSectors_Layer1", locProjectedSectorsMap[1]);
		dTreeFillData.Fill_Single<UInt_t>("NearestBCALSectors_Layer1_Downstream", locNearestSectorsMap_Downstream[1]);
		dTreeFillData.Fill_Single<UInt_t>("NearestBCALSectors_Layer1_Upstream", locNearestSectorsMap_Upstream[1]);
		//LAYER 2:
		dTreeFillData.Fill_Single<UInt_t>("ProjectedBCALSectors_Layer2", locProjectedSectorsMap[2]);
		dTreeFillData.Fill_Single<UInt_t>("NearestBCALSectors_Layer2_Downstream", locNearestSectorsMap_Downstream[2]);
		dTreeFillData.Fill_Single<UInt_t>("NearestBCALSectors_Layer2_Upstream", locNearestSectorsMap_Upstream[2]);
		//LAYER 3:
		dTreeFillData.Fill_Single<UInt_t>("ProjectedBCALSectors_Layer3", locProjectedSectorsMap[3]);
		dTreeFillData.Fill_Single<UInt_t>("NearestBCALSectors_Layer3_Downstream", locNearestSectorsMap_Downstream[3]);
		dTreeFillData.Fill_Single<UInt_t>("NearestBCALSectors_Layer3_Upstream", locNearestSectorsMap_Upstream[3]);
		//LAYER 4:
		dTreeFillData.Fill_Single<UInt_t>("ProjectedBCALSectors_Layer4", locProjectedSectorsMap[4]);
		dTreeFillData.Fill_Single<UInt_t>("NearestBCALSectors_Layer4_Downstream", locNearestSectorsMap_Downstream[4]);
		dTreeFillData.Fill_Single<UInt_t>("NearestBCALSectors_Layer4_Upstream", locNearestSectorsMap_Upstream[4]);

		//FILL TTREE
		dTreeInterface->Fill(dTreeFillData);
	}

	// FILL HISTOGRAMS
	// Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
	{
		//Fill Shower
		for(auto& locShowerPair : locShowerVector_ShowerFound)
			dHist_HadronicShowerMatched->Fill(locShowerPair.first, locShowerPair.second);
		for(auto& locShowerPair : locShowerVector_ShowerTotal)
			dHist_HadronicShowerTotal->Fill(locShowerPair.first, locShowerPair.second);

		//Fill Hit Found
		for(auto& locLayerPair : locHitMap_HitFound)
		{
			for(auto& locEndPair : locLayerPair.second)
			{
				for(int& locSector : locEndPair.second)
					dHistMap_HitFound[locLayerPair.first][locEndPair.first]->Fill(locSector);
			}
		}

		//Fill Hit Total
		for(auto& locLayerPair : locHitMap_HitTotal)
		{
			for(auto& locEndPair : locLayerPair.second)
			{
				for(int& locSector : locEndPair.second)
					dHistMap_HitTotal[locLayerPair.first][locEndPair.first]->Fill(locSector);
			}
		}
	}
	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

	return NOERROR;
}

void JEventProcessor_BCAL_Hadronic_Eff::Fill_NoClusterStudy(const DChargedTrackHypothesis* locChargedTrackHypothesis, vector<const DBCALShower*>& locBCALShowers, const DParticleID* locParticleID, bool locIsMatchedToTrackFlag)
{
	const DTrackTimeBased* locTrackTimeBased = NULL;
	locChargedTrackHypothesis->GetSingle(locTrackTimeBased);

	//Predict BCAL Surface Hit Location
	unsigned int locPredictedSurfaceModule = 0, locPredictedSurfaceSector = 0;
	DVector3 locPredictedSurfacePosition;
	locParticleID->PredictBCALWedge(locTrackTimeBased->rt, locPredictedSurfaceModule, locPredictedSurfaceSector, &locPredictedSurfacePosition);

	//Find closest shower match for BCAL
	double locBestMatchDeltaPhi = 0.0, locBestMatchDeltaZ = 0.0;
	const DBCALShower* locClosestBCALShower = locParticleID->Get_ClosestToTrack_BCAL(locTrackTimeBased, locBCALShowers, locBestMatchDeltaPhi, locBestMatchDeltaZ);

	//TRACK
	dTreeFillData.Fill_Single<Int_t>("PID_PDG", PDGtype(locChargedTrackHypothesis->PID()));
	dTreeFillData.Fill_Single<Float_t>("TimingBeta", locChargedTrackHypothesis->measuredBeta());
	dTreeFillData.Fill_Single<Float_t>("TrackVertexZ", locChargedTrackHypothesis->position().Z());
	DVector3 locDP3 = locChargedTrackHypothesis->momentum();
	TVector3 locP3(locDP3.X(), locDP3.Y(), locDP3.Z());
	dTreeFillData.Fill_Single<TVector3>("TrackP3", locP3);
	dTreeFillData.Fill_Single<UInt_t>("TrackCDCRings", locTrackTimeBased->dCDCRings);
	dTreeFillData.Fill_Single<UInt_t>("TrackFDCPlanes", locTrackTimeBased->dFDCPlanes);

	//SHOWER
	double locShowerEnergy = (locClosestBCALShower != nullptr) ? locClosestBCALShower->E : 0.0;
	dTreeFillData.Fill_Single<Float_t>("NearestShowerEnergy", locShowerEnergy);
	dTreeFillData.Fill_Single<Float_t>("TrackDeltaPhiToShower", locBestMatchDeltaPhi); //is signed: BCAL - Track
	dTreeFillData.Fill_Single<Float_t>("TrackDeltaZToShower", locBestMatchDeltaZ); //is signed: BCAL - Track
	dTreeFillData.Fill_Single<Float_t>("ProjectedBCALHitPhi", locPredictedSurfacePosition.Phi()*180.0/TMath::Pi());
	dTreeFillData.Fill_Single<Float_t>("ProjectedBCALHitZ", locPredictedSurfacePosition.Z());
	dTreeFillData.Fill_Single<Bool_t>("IsMatchedToTrack", locIsMatchedToTrackFlag);

	//HIT SEARCH
	dTreeFillData.Fill_Single<UChar_t>("BCALClusterLayers", 0);
	dTreeFillData.Fill_Single<UChar_t>("IsHitInCluster", 0);
	dTreeFillData.Fill_Single<UInt_t>("ProjectedBCALSectors_Layer1", 0);
	dTreeFillData.Fill_Single<UInt_t>("NearestBCALSectors_Layer1_Downstream", 0);
	dTreeFillData.Fill_Single<UInt_t>("NearestBCALSectors_Layer1_Upstream", 0);
	//LAYER 2:
	dTreeFillData.Fill_Single<UInt_t>("ProjectedBCALSectors_Layer2", 0);
	dTreeFillData.Fill_Single<UInt_t>("NearestBCALSectors_Layer2_Downstream", 0);
	dTreeFillData.Fill_Single<UInt_t>("NearestBCALSectors_Layer2_Upstream", 0);
	//LAYER 3:
	dTreeFillData.Fill_Single<UInt_t>("ProjectedBCALSectors_Layer3", 0);
	dTreeFillData.Fill_Single<UInt_t>("NearestBCALSectors_Layer3_Downstream", 0);
	dTreeFillData.Fill_Single<UInt_t>("NearestBCALSectors_Layer3_Upstream", 0);
	//LAYER 4:
	dTreeFillData.Fill_Single<UInt_t>("ProjectedBCALSectors_Layer4", 0);
	dTreeFillData.Fill_Single<UInt_t>("NearestBCALSectors_Layer4_Downstream", 0);
	dTreeFillData.Fill_Single<UInt_t>("NearestBCALSectors_Layer4_Upstream", 0);

	//FILL TTREE
	dTreeInterface->Fill(dTreeFillData);
}

double JEventProcessor_BCAL_Hadronic_Eff::Calc_ProjectedSector(int locLayer, const map<int, map<int, set<const DBCALPoint*> > >& locSortedPoints)
{
	//in the nearest layer of the cluster, find the average hit-sector
		//if two layers are adjacent: average them

	//ProjectedSector is (module - 1)*4 + sector //double: can be average
	if((locLayer == 2) || (locLayer == 3))
	{
		//if layer 2 or 3, there must be an adjacent layer present (must be 2 others)
		vector<double> locSectors;
		auto locLayerIterator = locSortedPoints.find(locLayer - 1);
		if(locLayerIterator != locSortedPoints.end())
			locSectors.push_back(Calc_AverageSector(locLayerIterator->second));

		locLayerIterator = locSortedPoints.find(locLayer + 1);
		if(locLayerIterator != locSortedPoints.end())
			locSectors.push_back(Calc_AverageSector(locLayerIterator->second));

		//average results if hits in 2 layers, but beware 0/2pi wrap-around
		if(locSectors.size() == 1)
			return locSectors[0];

		double locProjectedSector = locSectors[0] + locSectors[1];
		if(abs(locSectors[0] - locSectors[1]) > 96.0) //wrap-around!
			locProjectedSector += 192.0;
		locProjectedSector /= 2.0;
		while(locProjectedSector >= 193.0)
			locProjectedSector -= 192.0;
		return locProjectedSector;
	}
	else if(locLayer == 1)
	{
		auto locLayerIterator = locSortedPoints.find(2);
		if(locLayerIterator != locSortedPoints.end())
			return Calc_AverageSector(locLayerIterator->second);
		else
			return Calc_AverageSector(locSortedPoints.at(3));
	}
	else //layer 4
	{
		auto locLayerIterator = locSortedPoints.find(3);
		if(locLayerIterator != locSortedPoints.end())
			return Calc_AverageSector(locLayerIterator->second);
		else
			return Calc_AverageSector(locSortedPoints.at(2));
	}
}

double JEventProcessor_BCAL_Hadronic_Eff::Calc_AverageSector(const map<int, set<const DBCALPoint*> >& locBCALPoints)
{
	//int: full sector
	if(locBCALPoints.empty())
		return -1.0;

	//Has low/high: Beware 0/2pi wrap-around showers!
	vector<pair<double, double> > locSectors; //double: energy
	bool locHasLowPhiHits = false, locHasHighPhiHits = false;
	for(auto& locPointPair : locBCALPoints)
	{
		//sum energy of hits in the sector //if somehow > 1
		double locEnergySum = 0.0;
		for(auto& locPoint : locPointPair.second)
			locEnergySum += locPoint->E();

		locSectors.push_back(pair<double, double>(double(locPointPair.first), locEnergySum));
		if(locPointPair.first <= 12.0) //first 3 modules
			locHasLowPhiHits = true;
		else if(locPointPair.first >= 179) //last 3 modules
			locHasHighPhiHits = true;
	}

	//compute weighted average: numerator & denominator
	double locNumerator = 0.0, locDenominator = 0.0;
	for(auto& locSectorPair : locSectors)
	{
		double locWeight = 1.0/(locSectorPair.second*locSectorPair.second);

		double locSector = locSectorPair.first;
		if(locHasLowPhiHits && locHasHighPhiHits && (locSector < 96.0)) //96: half-way point
			locSector += 192.0; //add 2pi
		locNumerator += locSector*locWeight;

		locDenominator += locWeight;
	}

	//calc average
	double locAverageSector = locNumerator/locDenominator;
	while(locAverageSector >= 193.0)
		locAverageSector -= 192.0;

	//return
	return locAverageSector;
}

pair<const DBCALPoint*, double> JEventProcessor_BCAL_Hadronic_Eff::Find_NearestPoint(double locProjectedSector, const map<int, set<const DBCALPoint*> >& locLayerBCALPoints, const DBCALCluster* locBCALCluster, double locTimeCut)
{
	//input map int: full_sector
	//returned double: delta-sector
	const DBCALPoint* locBestBCALPoint = NULL;
	double locBestDeltaSector = 999.0;

	for(auto& locPointPair : locLayerBCALPoints)
	{
		//do time cut!
		const DBCALPoint* locBCALPoint = Find_ClosestTimePoint(locPointPair.second, locBCALCluster, locTimeCut);
		if(locBCALPoint == NULL)
			continue; //no points, or failed time cut //not applied if cut <= 0 //cluster

		double locDeltaSector = Calc_DeltaSector<double>(locPointPair.first, locProjectedSector);
		if(fabs(locDeltaSector) >= fabs(locBestDeltaSector))
			continue;
		locBestDeltaSector = locDeltaSector;
		locBestBCALPoint = locBCALPoint;
	}

	return pair<const DBCALPoint*, double>(locBestBCALPoint, locBestDeltaSector);
}

pair<const DBCALUnifiedHit*, double> JEventProcessor_BCAL_Hadronic_Eff::Find_NearestHit(double locProjectedSector, const map<int, set<const DBCALUnifiedHit*> >& locLayerUnifiedHits, const DBCALCluster* locBCALCluster, double locTimeCut)
{
	//input map int: full_sector
	//returned double: delta-sector
	const DBCALUnifiedHit* locBestBCALHit = NULL;
	double locBestDeltaSector = 999.0;

	for(auto& locHitPair : locLayerUnifiedHits)
	{
		const DBCALUnifiedHit* locHit = Find_ClosestTimeHit(locHitPair.second, locBCALCluster, locTimeCut);
		if(locHit == NULL)
			continue; //no hits, or failed time cut //not applied if cut <= 0 //cluster

		double locDeltaSector = Calc_DeltaSector<double>(locHitPair.first, locProjectedSector);
		if(fabs(locDeltaSector) >= fabs(locBestDeltaSector))
			continue;
		locBestDeltaSector = locDeltaSector;
		locBestBCALHit = locHit;
	}

	return pair<const DBCALUnifiedHit*, double>(locBestBCALHit, locBestDeltaSector);
}

const DBCALPoint* JEventProcessor_BCAL_Hadronic_Eff::Find_ClosestTimePoint(const set<const DBCALPoint*>& locPoints, const DBCALCluster* locBCALCluster, double locTimeCut)
{
	if(locPoints.empty())
		return nullptr;
	if(locTimeCut <= 0.0)
		return *(locPoints.begin()); //no cut, doesn't matter which one

	const DBCALPoint* locBestPoint = NULL;
	for(auto& locPoint : locPoints)
	{
		double locDeltaT = fabs(locBCALCluster->t() - locPoint->t());
		if(locDeltaT > locTimeCut)
			continue;

		//if best time, register new delta-t
		locTimeCut = locDeltaT; //cut becomes tighter to improve match
		locBestPoint = locPoint;
	}

	return locBestPoint;
}

const DBCALUnifiedHit* JEventProcessor_BCAL_Hadronic_Eff::Find_ClosestTimeHit(const set<const DBCALUnifiedHit*>& locHits, const DBCALCluster* locBCALCluster, double locTimeCut)
{
	if(locHits.empty())
		return nullptr;
	if(locTimeCut <= 0.0)
		return *(locHits.begin()); //no cut, doesn't matter which one

	//THIS CODE HAS BEEN COPIED FROM DBCALCluster_factory.cc
	const DBCALUnifiedHit* locBestHit = NULL;
	for(auto& locHit : locHits)
	{
		// given the location of the cluster, we need the best guess for z with respect to target at this radius
		double locClusterZ = locBCALCluster->rho()*cos(locBCALCluster->theta()) + dTargetCenterZ;

		// calc the distance to upstream or downstream end of BCAL depending on where the hit was with respect to the cluster z position.
		double locDistance = DBCALGeometry::GetBCAL_length()/2.0;
		if(locHit->end == 0)
			locDistance += locClusterZ - DBCALGeometry::GetBCAL_center();
		else
			locDistance += DBCALGeometry::GetBCAL_center() - locClusterZ;

		//correct time, calc delta-t
		int channel_calib = 16*(locHit->module - 1) + 4*(locHit->layer - 1) + locHit->sector - 1; //for CCDB
		double locCorrectedHitTime = locHit->t - locDistance/effective_velocities[channel_calib];  // to the interaction point in the bar.
		double locDeltaT = fabs(locCorrectedHitTime - locBCALCluster->t());

		//if best time, register new delta-t
		if(locDeltaT > locTimeCut)
			continue;
		locTimeCut = locDeltaT; //cut becomes tighter to improve match
		locBestHit = locHit;
	}

	return locBestHit;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_BCAL_Hadronic_Eff::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_BCAL_Hadronic_Eff::fini(void)
{
	// Called before program exit after event processing is finished.  

	delete dCutAction_TrackHitPattern;
	delete dTreeInterface; //saves trees to file, closes file

	return NOERROR;
}
