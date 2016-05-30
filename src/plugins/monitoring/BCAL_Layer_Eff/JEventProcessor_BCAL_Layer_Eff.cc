// $Id$
//
//    File: JEventProcessor_BCAL_Layer_Eff.cc
//

#include "JEventProcessor_BCAL_Layer_Eff.h"

// Routine used to create our JEventProcessor
extern "C"
{
	void InitPlugin(JApplication *locApplication)
	{
		InitJANAPlugin(locApplication);
		locApplication->AddProcessor(new JEventProcessor_BCAL_Layer_Eff()); //register this plugin
	}
} // "C"

//define static local variable //declared in header file
thread_local DTreeFillData JEventProcessor_BCAL_Layer_Eff::dTreeFillData;

//------------------
// init
//------------------
jerror_t JEventProcessor_BCAL_Layer_Eff::init(void)
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
	gDirectory->mkdir("bcal_layer_eff")->cd();

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
	dTreeInterface = DTreeInterface::Create_DTreeInterface("bcal_layer_eff", "tree_bcal_layer_eff.root");

	//TTREE BRANCHES
	DTreeBranchRegister locTreeBranchRegister;

	//TRACK
	locTreeBranchRegister.Register_Branch_Single<Int_t>("PID_PDG"); //gives charge, mass, beta
	locTreeBranchRegister.Register_Branch_Single<Float_t>("TrackVertexZ");
	locTreeBranchRegister.Register_Branch_Single<TVector3>("TrackP3");
	locTreeBranchRegister.Register_Branch_Single<Float_t>("TrackDeltaPhiToShower"); //is signed: BCAL - Track
	locTreeBranchRegister.Register_Branch_Single<Float_t>("TrackDeltaZToShower"); //is signed: BCAL - Track
	locTreeBranchRegister.Register_Branch_Single<Float_t>("ProjectedBCALHitPhi"); //degrees
	locTreeBranchRegister.Register_Branch_Single<Float_t>("ProjectedBCALHitZ");

	//HIT SEARCH
	//BCALClusterLayers: first 4 bits: point layers, next 4: unmatched-unified-hit layers
	locTreeBranchRegister.Register_Branch_Single<UChar_t>("BCALClusterLayers");
	//LAYER 1:
	//"Sector:" 4*(module - 1) + sector //sector: 1 -> 192
	locTreeBranchRegister.Register_Branch_Single<UChar_t>("ProjectedBCALSectors_Layer1"); //0 if biased or indeterminate
	locTreeBranchRegister.Register_Branch_Single<UChar_t>("NearestBCALSectors_Layer1_Downstream"); //0 if not found
	locTreeBranchRegister.Register_Branch_Single<UChar_t>("NearestBCALSectors_Layer1_Upstream"); //0 if not found
	//LAYER 2:
	locTreeBranchRegister.Register_Branch_Single<UChar_t>("ProjectedBCALSectors_Layer2"); //0 if biased or indeterminate
	locTreeBranchRegister.Register_Branch_Single<UChar_t>("NearestBCALSectors_Layer2_Downstream"); //0 if not found
	locTreeBranchRegister.Register_Branch_Single<UChar_t>("NearestBCALSectors_Layer2_Upstream"); //0 if not found
	//LAYER 3:
	locTreeBranchRegister.Register_Branch_Single<UChar_t>("ProjectedBCALSectors_Layer3"); //0 if biased or indeterminate
	locTreeBranchRegister.Register_Branch_Single<UChar_t>("NearestBCALSectors_Layer3_Downstream"); //0 if not found
	locTreeBranchRegister.Register_Branch_Single<UChar_t>("NearestBCALSectors_Layer3_Upstream"); //0 if not found
	//LAYER 4:
	locTreeBranchRegister.Register_Branch_Single<UChar_t>("ProjectedBCALSectors_Layer4"); //0 if biased or indeterminate
	locTreeBranchRegister.Register_Branch_Single<UChar_t>("NearestBCALSectors_Layer4_Downstream"); //0 if not found
	locTreeBranchRegister.Register_Branch_Single<UChar_t>("NearestBCALSectors_Layer4_Upstream"); //0 if not found

	//REGISTER BRANCHES
	dTreeInterface->Create_Branches(locTreeBranchRegister);

	return NOERROR;
}


//------------------
// brun
//------------------
jerror_t JEventProcessor_BCAL_Layer_Eff::brun(jana::JEventLoop* locEventLoop, int locRunNumber)
{
	// This is called whenever the run number changes

	return NOERROR;
}

//------------------
// evnt
//------------------

jerror_t JEventProcessor_BCAL_Layer_Eff::evnt(jana::JEventLoop* locEventLoop, uint64_t locEventNumber)
{
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// locEventLoop->Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.

	// This plugin is used to determine the reconstruction efficiency of hits in the BCAL
		// Note, this is hit-level, not shower-level.  Hits: By sector/layer/module/end

/*
	//CUT ON TRIGGER TYPE
	const DTrigger* locTrigger = NULL;
	locEventLoop->GetSingle(locTrigger);
*/
	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks);

	const DParticleID* locParticleID = NULL;
	locEventLoop->GetSingle(locParticleID);

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

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

	// Loop over the good tracks, using the best DTrackTimeBased object for each
	for(auto& locChargedTrackHypothesis : locBestTracks)
	{
		//get the best-matched DBCALShower for this track (if any)
		const DBCALShowerMatchParams* locBCALShowerMatchParams = locChargedTrackHypothesis->Get_BCALShowerMatchParams();
		if(locBCALShowerMatchParams == NULL)
			continue; //nope

		//use the DBCALPoints in this DBCALShower to help define where DBCALUnifiedHits are expected to be
		//The DBCALShower reconstruction does not care which BCAL layers fired
			//Therefore, given a DBCALShower hit layer, whether or not there are hits in the other layers is an unbiased check

		const DTrackTimeBased* locTrackTimeBased = NULL;
		locChargedTrackHypothesis->GetSingle(locTrackTimeBased);

		//get the clusters in the shower, reject if > 1
		vector<const DBCALCluster*> locClusters;
		locBCALShowerMatchParams->dBCALShower->Get(locClusters);
		if(locClusters.size() > 1)
			continue; //likely a messy shower, will probably mess up efficiency calculation. 

		//get cluster, points, unmatched hits
		const DBCALCluster* locCluster = locClusters[0];
		vector<const DBCALPoint*> locBCALPoints = locCluster->points();
		vector<pair<const DBCALUnifiedHit*, double> > locBCALHits = locCluster->hits(); //double: corrected energy

		//sort points by layer
		map<int, set<const DBCALPoint*> > locSortedPoints; //int = layer
		for(auto& locBCALPoint : locBCALPoints)
			locSortedPoints[locBCALPoint->layer()].insert(locBCALPoint);

		//require hits in at least 2 layers: make sure it's not a noise cluster
		//Note: This can introduce bias into the study. To avoid bias, only evaluate layer efficiency if at least 2 OTHER layers have hits
		if(locSortedPoints.size() < 2)
			continue;

		//sort hits by layer
		map<int, set<const DBCALUnifiedHit*> > locSortedHits; //int = layer
		for(auto& locHitPair : locBCALHits)
			locSortedHits[locHitPair.first->layer].insert(locHitPair.first);

		//Tree-save variables
		map<int, UChar_t> locProjectedSectorsMap, locNearestSectorsMap_Downstream, locNearestSectorsMap_Upstream;
		//locProjectedSectors: (rounded from search): 4*(module - 1) + sector //sector: 1 -> 192, 0 if biased or indeterminate
		//locNearestSectors_Downstream & locNearestSectors_Upstream: 4*(module - 1) + sector //sector: 1 -> 192, 0 if not found

		//register layers
		UChar_t locClusterLayers = 0; //which layers are in the cluster (both by points & by hits): 4bits each: 8total
		for(auto& locLayerPair : locSortedPoints)
			locClusterLayers |= (1 << (locLayerPair.first - 1));
		for(auto& locLayerPair : locSortedHits)
			locClusterLayers |= (1 << (locLayerPair.first + 4 - 1));

		//loop over layers
		for(int locLayer = 1; locLayer <= 4; ++locLayer)
		{
			//initialize to 0's
			locProjectedSectorsMap[locLayer] = 0;
			locNearestSectorsMap_Downstream[locLayer] = 0;
			locNearestSectorsMap_Upstream[locLayer] = 0;

			//require at least 2 other layers to have points
			//use the DBCALPoint's (double-ended hits) to define where hits are to be expected (not the single-ended ones!)
			//if no hit in this layer, then requirement is already met. if there is, must check
			auto locPointsIterator = locSortedPoints.find(locLayer);
			if((locPointsIterator != locSortedPoints.end()) && (locSortedPoints.size() <= 2))
				continue; //not at least 2 hits in other layers

			//ideally, would require hits in adjacent layers, but cannot since one may be dead: won't be able to evaluate for this layer
				//can always constrain this later though, if desired

			//ok, can now evaluate whether or not there are hits in this layer in an unbiased manner

			//in this layer, need to know where to search for a hit:
				//showers may contain hits in multiple sectors within a layer: 
				//cannot evaluate all found as "good," because when missing, don't know how many to mark as "bad"
				//so, pick the most-probable location: energy-weighted average of hit locations in the nearest layers

			//locProjectedSector is (module - 1)*4 + sector //double: average
			double locProjectedSector = Calc_ProjectedSector(locLayer, locSortedPoints);
			int locProjectedSectorInt = int(locProjectedSector + 0.5);

			//register for tree & hists:
			locProjectedSectorsMap[locLayer] = UChar_t(locProjectedSectorInt); //tree
			locHitMap_HitTotal[locLayer][true].push_back(locProjectedSectorInt); //hist
			locHitMap_HitTotal[locLayer][false].push_back(locProjectedSectorInt); //hist

			//the clustering algorithm is a far more sophisticated & accurate method of finding hits than anything we can come up with here
				//in fact, the algorithm for determining whether a hit is a "good" match for efficiency purposes may as well be the clustering algorithm
				//if a hit was reconstructed and usable, it will be in the cluster

			//now, see if can find a matching hit in the layer
			//first, search the cluster points //returned double: signed distance (point sector - search sector)
			pair<const DBCALPoint*, double> locNearestClusterPoint(NULL, 999.0);
			if(locPointsIterator != locSortedPoints.end())
				locNearestClusterPoint = Find_NearestClusterPoint(locProjectedSector, locPointsIterator->second);

			//next, search the unmatched cluster hits
			pair<const DBCALUnifiedHit*, double> locNearestClusterHit(NULL, 999.0);
			auto locHitsIterator = locSortedHits.find(locLayer);
			if(locHitsIterator != locSortedHits.end())
				locNearestClusterHit = Find_NearestClusterHit(locProjectedSector, locHitsIterator->second);

			if((locNearestClusterPoint.first == nullptr) && (locNearestClusterHit.first == nullptr))
				continue; //nothing found for this layer

			//determine whether to register the results for the nearest point or the nearest hit (the closest)
			bool locRegisterPointFlag = (locNearestClusterPoint.first != nullptr);
			if(locRegisterPointFlag && (locNearestClusterHit.first != nullptr)) //both present, see which is closer
				locRegisterPointFlag = (locNearestClusterPoint.second < locNearestClusterHit.second);

			//ok, now register results, preparing to write to tree
			if(locRegisterPointFlag)
			{
				//register point: both ends
				int locSector = (locNearestClusterPoint.first->module() - 1)*4 + locNearestClusterPoint.first->sector();
				//tree
				locNearestSectorsMap_Upstream[locLayer] = UChar_t(locSector);
				locNearestSectorsMap_Downstream[locLayer] = UChar_t(locSector);
				//hists
				locHitMap_HitFound[locLayer][true].push_back(locProjectedSectorInt); //hist
				locHitMap_HitFound[locLayer][false].push_back(locProjectedSectorInt); //hist
			}
			else
			{
				//register hit: one end
				int locSector = (locNearestClusterHit.first->module - 1)*4 + locNearestClusterHit.first->sector;
				if(locNearestClusterHit.first->end == DBCALGeometry::kUpstream)
				{
					locNearestSectorsMap_Upstream[locLayer] = UChar_t(locSector);
					locHitMap_HitFound[locLayer][true].push_back(locProjectedSectorInt); //hist
				}
				else
				{
					locNearestSectorsMap_Downstream[locLayer] = UChar_t(locSector);
					locHitMap_HitFound[locLayer][false].push_back(locProjectedSectorInt); //hist
				}
			}
		}

		//Predict BCAL Surface Hit Location
		unsigned int locPredictedSurfaceModule = 0, locPredictedSurfaceSector = 0;
		DVector3 locPredictedSurfacePosition;
		locParticleID->PredictBCALWedge(locTrackTimeBased->rt, locPredictedSurfaceModule, locPredictedSurfaceSector, &locPredictedSurfacePosition);

		//STAGE DATA FOR TREE FILL

		//TRACK
		dTreeFillData.Fill_Single<Int_t>("PID_PDG", PDGtype(locChargedTrackHypothesis->PID()));
		dTreeFillData.Fill_Single<Float_t>("TrackVertexZ", locChargedTrackHypothesis->position().Z());
		DVector3 locDP3 = locChargedTrackHypothesis->momentum();
		TVector3 locP3(locDP3.X(), locDP3.Y(), locDP3.Z());
		dTreeFillData.Fill_Single<TVector3>("TrackP3", locP3);
		dTreeFillData.Fill_Single<Float_t>("TrackDeltaPhiToShower", locBCALShowerMatchParams->dDeltaPhiToShower); //is signed: BCAL - Track
		dTreeFillData.Fill_Single<Float_t>("TrackDeltaZToShower", locBCALShowerMatchParams->dDeltaZToShower); //is signed: BCAL - Track
		dTreeFillData.Fill_Single<Float_t>("ProjectedBCALHitPhi", locPredictedSurfacePosition.Phi()*180.0/TMath::Pi());
		dTreeFillData.Fill_Single<Float_t>("ProjectedBCALHitZ", locPredictedSurfacePosition.Z());

		//HIT SEARCH
		//BCALClusterLayers: first 4 bits: point layers, next 4: unmatched-unified-hit layers
		dTreeFillData.Fill_Single<UChar_t>("BCALClusterLayers", locClusterLayers);
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
		//Fill Found
		for(auto& locLayerPair : locHitMap_HitFound)
		{
			for(auto& locEndPair : locLayerPair.second)
			{
				for(int& locSector : locEndPair.second)
					dHistMap_HitFound[locLayerPair.first][locEndPair.first]->Fill(locSector);
			}
		}

		//Fill Total
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

double JEventProcessor_BCAL_Layer_Eff::Calc_ProjectedSector(int locLayer, const map<int, set<const DBCALPoint*> >& locSortedPoints)
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

double JEventProcessor_BCAL_Layer_Eff::Calc_AverageSector(const set<const DBCALPoint*>& locBCALPoints)
{
	if(locBCALPoints.empty())
		return -1.0;

	//Beware 0/2pi wrap-around showers!
	vector<pair<double, double> > locSectors; //double: energy
	bool locHasLowPhiHits = false, locHasHighPhiHits = false;
	for(auto& locBCALPoint : locBCALPoints)
	{
		int locSector = (locBCALPoint->module() - 1)*4 + locBCALPoint->sector(); //1 -> 192
		locSectors.push_back(pair<double, double>(double(locSector), locBCALPoint->E()));
		if(locSector <= 8) //first 2 modules
			locHasLowPhiHits = true;
		else if(locSector >= 183) //last 2 modules
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

pair<const DBCALUnifiedHit*, double> JEventProcessor_BCAL_Layer_Eff::Find_NearestClusterHit(double locProjectedSector, const set<const DBCALUnifiedHit*>& locBCALUnifiedHits)
{
	//report all hits within 0.99 sectors of the best delta-sector
	//input map int: full_sector
	const DBCALUnifiedHit* locBestBCALHit = NULL;
	double locBestDeltaSector = 16.0; //if you are >= 4 modules away, don't bother reporting
	for(auto& locBCALHit : locBCALUnifiedHits)
	{
		double locSector = (locBCALHit->module - 1)*4 + locBCALHit->sector;
		double locDeltaSector = locSector - locProjectedSector;
		if(fabs(locDeltaSector) >= fabs(locBestDeltaSector))
			continue;
		locBestDeltaSector = locDeltaSector;
		locBestBCALHit = locBCALHit;
	}

	return pair<const DBCALUnifiedHit*, double>(locBestBCALHit, locBestDeltaSector);
}

pair<const DBCALPoint*, double> JEventProcessor_BCAL_Layer_Eff::Find_NearestClusterPoint(double locProjectedSector, const set<const DBCALPoint*>& locClusterLayerBCALPoints)
{
	//returned double: delta-sector
	//no time cut: was already placed when forming the cluster
	const DBCALPoint* locBestBCALPoint = NULL;
	double locBestDeltaSector = 16.0; //if you are >= 4 modules away, don't bother reporting

	for(auto& locBCALPoint : locClusterLayerBCALPoints)
	{
		double locSector = double((locBCALPoint->module() - 1)*4 + locBCALPoint->sector());
		double locDeltaSector = locSector - locProjectedSector;
		if(fabs(locDeltaSector) >= fabs(locBestDeltaSector))
			continue;
		locBestDeltaSector = locDeltaSector;
		locBestBCALPoint = locBCALPoint;
	}

	return pair<const DBCALPoint*, double>(locBestBCALPoint, locBestDeltaSector);
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_BCAL_Layer_Eff::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_BCAL_Layer_Eff::fini(void)
{
	// Called before program exit after event processing is finished.  

	delete dCutAction_TrackHitPattern;
	delete dTreeInterface; //saves trees to file, closes file

	return NOERROR;
}

