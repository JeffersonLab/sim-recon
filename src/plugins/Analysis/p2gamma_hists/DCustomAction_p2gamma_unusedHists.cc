// $Id$
//
//    File: DCustomAction_p2gamma_unusedHists.cc
// Created: Thu Jan 22 08:06:18 EST 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "DCustomAction_p2gamma_unusedHists.h"

void DCustomAction_p2gamma_unusedHists::Initialize(JEventLoop* locEventLoop)
{

	// get PID algos
	const DParticleID* locParticleID = NULL;
        locEventLoop->GetSingle(locParticleID);
	dParticleID = locParticleID;

	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//Required: Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this action.
			//If another thread has already created the folder, it just changes to it. 
		CreateAndChangeTo_ActionDirectory();

		dNShowerBCAL_FCAL = GetOrCreate_Histogram<TH2I>("NShowerBCAL_FCAL", "Number of non-matched showers; FCAL; BCAL", 10, 0, 10, 10, 0, 10);

		TString locHistName;
		TString locHistTitle;
		TString match[2] = {"Match",""};
		map<int, TString> charge;
		charge[-1] = "Neg"; charge[+1] = "Pos";

		for(int locDummy = 0; locDummy < 2; locDummy++) {
                        bool locMatch = (locDummy == 0);
			
			for(int locCharge = -1; locCharge < 2; locCharge+=2) {
				locHistName = "TrackNhits_Theta" + match[locDummy] + charge[locCharge];
				locHistTitle = "Number of hits on track vs #theta: " + match[locDummy] + charge[locCharge];
				locHistTitle += "; #theta; # DOF";
				dHistMap_TrackNhits_Theta[locMatch][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 150, 0., 150., 50, 0, 50);

				locHistName = "TrackChiSq_Theta" + match[locDummy] + charge[locCharge];
				locHistTitle = "FOM for track #chi^{2} vs #theta: " + match[locDummy] + charge[locCharge];
				locHistTitle += "; #theta; #chi^{2}";
				dHistMap_TrackChiSq_Theta[locMatch][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 150, 0., 150., 500, 0., 100.);

				locHistName = "TrackFOM_Theta" + match[locDummy] + charge[locCharge];
				locHistTitle = "FOM for track fit vs #theta: " + match[locDummy] + charge[locCharge];
				locHistTitle += "; #theta; FOM";
				dHistMap_TrackFOM_Theta[locMatch][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 150, 0., 150., 1000, 0., 1.);

				locHistName = "TrackP_Theta" + match[locDummy] + charge[locCharge];
				locHistTitle = "P vs #theta: " + match[locDummy] + charge[locCharge];
				locHistTitle += "; #theta; momentum";
				dHistMap_TrackP_Theta[locMatch][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 150, 0., 150., 120, 0., 6.);

				locHistName = "TrackPOCAXY" + match[locDummy] + charge[locCharge];
				locHistTitle = "POCA to beamline Y vs X: " + match[locDummy] + charge[locCharge];
				locHistTitle += "; POCA X; POCAY";
				dHistMap_TrackPOCAXY[locMatch][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 200, -20., 20., 200, -20., 20.);

				locHistName = "TrackPOCAZ" + match[locDummy] + charge[locCharge];
				locHistTitle = "POCA to beamline Z: " + match[locDummy] + charge[locCharge];
				locHistTitle += "; POCA Z";
				dHistMap_TrackPOCAZ[locMatch][locCharge] = GetOrCreate_Histogram<TH1I>(locHistName.Data(), locHistTitle.Data(), 500, 0., 250.);

				/*
				  locHistName = "TrackNposs" + match[locDummy] + charge[locCharge];
				  locHistTitle = "Number possible hits on track vs #theta: " + match[locDummy] + charge[locCharge];
				  locHistTitle += "; #theta; # possible hits";
				  dHistMap_TrackNposs_Theta[locMatch][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 150, 0., 150., 50, 0, 50);
				  
				  locHistName = "TrackHitFrac" + match[locDummy] + charge[locCharge];
				  locHistTitle = "Fraction of possible hits used in track fit vs #theta: " + match[locDummy] + charge[locCharge];
				  locHistTitle += "; #theta; # hits fit / # hits possible";
				  dHistMap_TrackHitFrac_Theta[locMatch][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 150, 0., 150., 100, 0., 1.);
				*/
			}
			
			vector<DetectorSystem_t> locDetectorSystems;
			locDetectorSystems.push_back(SYS_FCAL); locDetectorSystems.push_back(SYS_BCAL); 
			for(size_t loc_i = 0; loc_i < locDetectorSystems.size(); ++loc_i)
                        {
				DetectorSystem_t locSystem = locDetectorSystems[loc_i];
				TString locSystemName = (locSystem == SYS_FCAL) ? "FCAL" : "BCAL";
				double thetaMax = (locSystem == SYS_FCAL) ? 15. : 150.;

				locHistName = "ShowerEnergy_Theta" + match[locDummy] + locSystemName;
				locHistTitle = "Shower Energy vs #theta: " + match[locDummy] + locSystemName;
				locHistTitle += ";#theta; Energy";
				dHistMap_ShowerEnergy_Theta[locMatch][locSystem] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 150, 0., thetaMax, 100., 0., 5.);			

				locHistName = "ShowerPhi_Theta" + match[locDummy] + locSystemName;
				locHistTitle = "Shower #phi vs #theta: " + match[locDummy] + locSystemName;
				locHistTitle += ";#theta; #phi";
				dHistMap_ShowerPhi_Theta[locMatch][locSystem] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 150, 0., thetaMax, 360., -180., 180.);			

				locHistName = "ShowerNClusters" + match[locDummy] + locSystemName;
				locHistTitle = "Number of clusters in shower: " + match[locDummy] + locSystemName;
				dHistMap_ShowerNclusters[locMatch][locSystem] = GetOrCreate_Histogram<TH1I>(locHistName.Data(), locHistTitle.Data(), 15, 0, 15);


				locHistName = "ShowerNHits" + match[locDummy] + locSystemName;
				locHistTitle = "Number of hits in shower: " + match[locDummy] + locSystemName;
				dHistMap_ShowerNhits[locMatch][locSystem] = GetOrCreate_Histogram<TH1I>(locHistName.Data(), locHistTitle.Data(), 15, 0, 15);

				locHistName = "ShowerMaxEnergy_NHits" + match[locDummy] + locSystemName;
				locHistTitle = "E_{max}/E_{tot} vs number of hits in shower: " + match[locDummy] + locSystemName;
				dHistMap_ShowerMaxEnergy_Nhits[locMatch][locSystem] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 15, 0, 15, 500, 0., 1.);

				locHistName = "ShowerDeltaT_NHits" + match[locDummy] + locSystemName;
				locHistTitle = "t_{shower} - t_{#gamma} vs number of hits in shower: " + match[locDummy] + locSystemName;
				dHistMap_ShowerDeltaT_Nhits[locMatch][locSystem] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 15, 0, 15, 200, -100., 100.);

				locHistName = "ShowerDeltaT_E" + match[locDummy] + locSystemName;
				locHistTitle = "t_{shower} - t_{#gamma} vs shower energy: " + match[locDummy] + locSystemName;
				dHistMap_ShowerDeltaT_E[locMatch][locSystem] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 100., 0., 5., 200, -100., 100.);

				locHistName = "ShowerE_Theta" + match[locDummy] + locSystemName;
				locHistTitle = "t_{shower} - t_{#gamma} vs shower #theta: " + match[locDummy] + locSystemName;
				dHistMap_ShowerE_Theta[locMatch][locSystem] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 150, 0., thetaMax, 120., 0., 6.);

				if(locSystem == SYS_BCAL){
					locHistName = "ShowerNcell" + match[locDummy] + locSystemName;
					locHistTitle = "Number of points in shower: " + match[locDummy] + locSystemName;
					dHistMap_BCALShowerNcell[locMatch] = GetOrCreate_Histogram<TH1I>(locHistName.Data(), locHistTitle.Data(), 15, 0, 15);			

					locHistName = "Layer1Energy_Theta" + match[locDummy] + locSystemName;
					locHistTitle = "Shower Layer 1 Energy vs #theta: " + match[locDummy] + locSystemName;
					locHistTitle += ";#theta; Energy Layer 1";
					dHistMap_Layer1Energy_Theta[locMatch] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 150, 0., thetaMax, 500., 0., 1.);

					locHistName = "Layer2Energy_Theta" + match[locDummy] + locSystemName;
					locHistTitle = "Shower Layer 2 Energy vs #theta: " + match[locDummy] + locSystemName;
					locHistTitle += ";#theta; Energy Layer 2";
					dHistMap_Layer2Energy_Theta[locMatch] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 150, 0., thetaMax, 500., 0., 1.);

					locHistName = "Layer3Energy_Theta" + match[locDummy] + locSystemName;
					locHistTitle = "Shower Layer 3 Energy vs #theta: " + match[locDummy] + locSystemName;
					locHistTitle += ";#theta; Energy Layer 3";
					dHistMap_Layer3Energy_Theta[locMatch] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 150, 0., thetaMax, 500., 0., 1.);

					locHistName = "Layer4Energy_Theta" + match[locDummy] + locSystemName;
					locHistTitle = "Shower Layer 4 Energy vs #theta: " + match[locDummy] + locSystemName;
					locHistTitle += ";#theta; Energy Layer 4";
					dHistMap_Layer4Energy_Theta[locMatch] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 150, 0., thetaMax, 500., 0., 1.);
				}
			}
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DCustomAction_p2gamma_unusedHists::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	
	// should only have one reaction step
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(0);

	// get beam photon energy and final state particles
        const DKinematicData* locBeamPhoton = NULL;
        deque<const DKinematicData*> locParticles;
        if(!Get_UseKinFitResultsFlag()) { //measured
		locBeamPhoton = locParticleComboStep->Get_InitialParticle_Measured();
                locParticleComboStep->Get_FinalParticles_Measured(locParticles);
	}
	else {
		locBeamPhoton = locParticleComboStep->Get_InitialParticle();
		locParticleComboStep->Get_FinalParticles(locParticles);
	}
	double locBeamPhotonTime = locBeamPhoton->time();

	// detector matches for charged track -to- shower matching
	const DDetectorMatches* locDetectorMatches = NULL;
        locEventLoop->GetSingle(locDetectorMatches);

	// get all charged tracks and showers for comparison to combo (no "pre-select" cut here)
        vector<const DChargedTrack*> locChargedTracks;
        locEventLoop->Get(locChargedTracks); //, "PreSelect"); 
	vector<const DNeutralShower*> locNeutralShowers;
        locEventLoop->Get(locNeutralShowers); //, "PreSelect");  

	// loop over charged tracks to make list of unused
        for(size_t loc_j = 0; loc_j < locChargedTracks.size(); ++loc_j) {
                int nMatched = 0;

                // add tracks not in combo to vector of unused tracks
                for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i) {
			if(locParticles[loc_i]->PID() != Proton) continue;

                        const DChargedTrack* locChargedTrack = static_cast<const DChargedTrack*>(locParticleComboStep->Get_FinalParticle_SourceObject(loc_i));
                        if(locChargedTrack == NULL) continue; // should never happen

                        const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_BestFOM();
                        if(locChargedTracks[loc_j]->Get_BestFOM()->candidateid == locChargedTrackHypothesis->candidateid) {
                                nMatched++;
				FillTrack(locChargedTracks[loc_j], true);
			}
                }
		
		// plot properties of unused charged tracks
                if(nMatched==0) {
			FillTrack(locChargedTracks[loc_j], false);
		}
        }

	int nUnmatchedFCAL = 0;
	int nUnmatchedBCAL = 0;

	/////////////////////////////////////////////////////////////////////////////////////////
	// loop over neutral showers which should all be unused (unless bad matching criteria) //
	/////////////////////////////////////////////////////////////////////////////////////////
	for(size_t loc_j = 0; loc_j < locNeutralShowers.size(); ++loc_j) {
		int nMatched = 0;

		// compare shower to pi+pi- in DReaction
		for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i) {
			if(locParticles[loc_i]->PID() == Proton) continue;

			const DNeutralShower* locNeutralShower = static_cast<const DNeutralShower*>(locParticleComboStep->Get_FinalParticle_SourceObject(loc_i));
			if(locNeutralShower == locNeutralShowers[loc_j]){
				nMatched++;
				double locFlightTime = locNeutralShowers[loc_j]->dSpacetimeVertex.Vect().Mag()/SPEED_OF_LIGHT;
				FillShower(locNeutralShowers[loc_j], true, locBeamPhotonTime, locFlightTime);
			}
		}

		if(nMatched==0) {
			double locFlightTime = locNeutralShowers[loc_j]->dSpacetimeVertex.Vect().Mag()/SPEED_OF_LIGHT;
			FillShower(locNeutralShowers[loc_j], false, locBeamPhotonTime, locFlightTime);
			if(locNeutralShowers[loc_j]->dDetectorSystem == SYS_FCAL) nUnmatchedFCAL++;
			if(locNeutralShowers[loc_j]->dDetectorSystem == SYS_BCAL) nUnmatchedBCAL++;
		}	
	}
	
	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action(); //ACQUIRE ROOT LOCK!!
	{
		dNShowerBCAL_FCAL->Fill(nUnmatchedFCAL, nUnmatchedBCAL);
	}
	Unlock_Action(); //RELEASE ROOT LOCK!!

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}


void DCustomAction_p2gamma_unusedHists::FillTrack(const DChargedTrack* locChargedTrack, bool locMatch)
{

	const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_BestFOM();
	int locCharge = locChargedTrackHypothesis->charge();

	const DTrackTimeBased* locTrackTimeBased = NULL;
	locChargedTrackHypothesis->GetSingleT(locTrackTimeBased);
	
	float nHits = locTrackTimeBased->Ndof + 5.;
	
	/* Number of possible hits not implemented yet?  Always returns 0
	  float nPoss = locTrackTimeBased->cdc_hit_usage.total_hits + locTrackTimeBased->fdc_hit_usage.total_hits;
	  float fitFrac = 0;
	  if(nPoss > 0.) fitFrac = nHits/nPoss;
	  else cout<<nPoss<<" "<<nHits<<endl;
	*/

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action(); //ACQUIRE ROOT LOCK!!
	{
		dHistMap_TrackNhits_Theta[locMatch][locCharge]->Fill(locTrackTimeBased->momentum().Theta()*180/TMath::Pi(), nHits);
		dHistMap_TrackChiSq_Theta[locMatch][locCharge]->Fill(locTrackTimeBased->momentum().Theta()*180/TMath::Pi(), locTrackTimeBased->chisq);
		dHistMap_TrackFOM_Theta[locMatch][locCharge]->Fill(locTrackTimeBased->momentum().Theta()*180/TMath::Pi(), locTrackTimeBased->FOM);
		dHistMap_TrackP_Theta[locMatch][locCharge]->Fill(locTrackTimeBased->momentum().Theta()*180/TMath::Pi(), locTrackTimeBased->momentum().Mag());
		dHistMap_TrackPOCAXY[locMatch][locCharge]->Fill(locTrackTimeBased->position().X(), locTrackTimeBased->position().Y());
		dHistMap_TrackPOCAZ[locMatch][locCharge]->Fill(locTrackTimeBased->position().Z());

		//dHistMap_TrackNposs_Theta[locMatch][locCharge]->Fill(locTrackTimeBased->momentum().Theta()*180/TMath::Pi(), nPoss);
		//dHistMap_TrackHitFrac_Theta[locMatch][locCharge]->Fill(locTrackTimeBased->momentum().Theta()*180/TMath::Pi(), fitFrac);
	}
	Unlock_Action(); //RELEASE ROOT LOCK!!

	return;
}


void DCustomAction_p2gamma_unusedHists::FillShower(const DNeutralShower* locNeutralShower, bool locMatch, double locBeamPhotonTime, double locFlightTime)
{
	
	int nClusters = 0;
	int nHits = 0;

	double layerE[4] = {0., 0., 0., 0.};

	double locEnergy = locNeutralShower->dEnergy;
	DVector3 locPosition = locNeutralShower->dSpacetimeVertex.Vect();
	double locDeltaT = locNeutralShower->dSpacetimeVertex.T() - locFlightTime - locBeamPhotonTime;
	
	double locEnergyCluster = 0.;
	double locMaxEnergyCluster = 0.;

	DetectorSystem_t locSystem = locNeutralShower->dDetectorSystem;
	if(locSystem == SYS_FCAL) {
			
		const DFCALShower* locFCALShower = NULL;
		locNeutralShower->GetSingleT(locFCALShower);

		vector<const DFCALCluster*> locFCALClusters;
		locFCALShower->Get(locFCALClusters);
		nClusters = locFCALClusters.size();

		// get hits in FCAL shower
		for(unsigned int i=0; i<locFCALClusters.size(); i++){
			const vector<DFCALCluster::DFCALClusterHit_t> locFCALHits = locFCALClusters[i]->GetHits();

			locEnergyCluster = locFCALClusters[i]->getEnergy();
			locMaxEnergyCluster = locFCALClusters[i]->getEmax();

			for(unsigned int j=0; j<locFCALHits.size(); j++){
//				const DFCALCluster::DFCALClusterHit_t hit = locFCALHits[j];
				nHits++;
			}
		}
	}
	if(locSystem == SYS_BCAL) {

		const DBCALShower* locBCALShower = NULL;
		locNeutralShower->GetSingleT(locBCALShower);

		//FILL HISTOGRAMS
		//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
		//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
		Lock_Action(); //ACQUIRE ROOT LOCK!!
		{
			dHistMap_BCALShowerNcell[locMatch]->Fill(locBCALShower->N_cell);
		}
		Unlock_Action(); //RELEASE ROOT LOCK!!

		vector<const DBCALCluster*> locBCALClusters;
		locBCALShower->Get(locBCALClusters);
		nClusters = locBCALClusters.size();
	
		// get points in BCAL shower
		for(unsigned int i=0; i<locBCALClusters.size(); i++){
			vector<const DBCALPoint*> locBCALPoints;
			locBCALClusters[i]->Get(locBCALPoints);

			locEnergyCluster = locBCALClusters[i]->E();
			if(locBCALPoints.size() == 1) locMaxEnergyCluster = locBCALClusters[i]->E();

			for(unsigned int j=0; j<locBCALPoints.size(); j++){
				const DBCALPoint *point = locBCALPoints[j];
				nHits++;

				if(point->E() > locMaxEnergyCluster)
					locMaxEnergyCluster = point->E();

				layerE[point->layer()-1] += point->E();
			}
		}
	}

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action(); //ACQUIRE ROOT LOCK!!
	{
		dHistMap_ShowerEnergy_Theta[locMatch][locSystem]->Fill(locPosition.Theta()*180./TMath::Pi(), locEnergy);
		dHistMap_ShowerPhi_Theta[locMatch][locSystem]->Fill(locPosition.Theta()*180./TMath::Pi(), locPosition.Phi()*180./TMath::Pi());
		dHistMap_ShowerNclusters[locMatch][locSystem]->Fill(nClusters);
		dHistMap_ShowerNhits[locMatch][locSystem]->Fill(nHits);
		
		dHistMap_ShowerMaxEnergy_Nhits[locMatch][locSystem]->Fill(nHits, locMaxEnergyCluster/locEnergyCluster);
		dHistMap_ShowerDeltaT_Nhits[locMatch][locSystem]->Fill(nHits, locDeltaT);
		dHistMap_ShowerDeltaT_E[locMatch][locSystem]->Fill(locEnergy, locDeltaT);
		if(locDeltaT > -5. && locDeltaT < 5. && locSystem == SYS_FCAL)
			dHistMap_ShowerE_Theta[locMatch][locSystem]->Fill(locPosition.Theta()*180./TMath::Pi(), locEnergy);
		if(locDeltaT > -5. && locDeltaT < 5. && locSystem == SYS_BCAL)
			dHistMap_ShowerE_Theta[locMatch][locSystem]->Fill(locPosition.Theta()*180./TMath::Pi(), locEnergy);

		if(locSystem == SYS_BCAL) {
			dHistMap_Layer1Energy_Theta[locMatch]->Fill(locPosition.Theta()*180./TMath::Pi(), layerE[0]/locEnergy);
			dHistMap_Layer2Energy_Theta[locMatch]->Fill(locPosition.Theta()*180./TMath::Pi(), layerE[1]/locEnergy);
			dHistMap_Layer3Energy_Theta[locMatch]->Fill(locPosition.Theta()*180./TMath::Pi(), layerE[2]/locEnergy);
			dHistMap_Layer4Energy_Theta[locMatch]->Fill(locPosition.Theta()*180./TMath::Pi(), layerE[3]/locEnergy);
		}

	}
	Unlock_Action(); //RELEASE ROOT LOCK!!

	return;
}


/*

        // loop over neutral showers to make list of unused showers
        for(size_t loc_j = 0; loc_j < locNeutralShowers.size(); ++loc_j) {
                int nMatched = 0;

                // add showers not matched to tracks in combo to vector of unused showers
                for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i) {

                        const DChargedTrack* locChargedTrack = static_cast<const DChargedTrack*>(locParticleComboStep->Get_FinalParticle_SourceObject(loc_i));
                        if(locChargedTrack == NULL) continue; // should never happen

                        const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_BestFOM();
			const DTrackTimeBased* locTrackTimeBased = NULL;
			locChargedTrackHypothesis->GetSingleT(locTrackTimeBased);
			
			if(locNeutralShowers[loc_j]->dDetectorSystem == SYS_FCAL) {
				DFCALShowerMatchParams locFCALShowerMatchParams;
				bool foundFCAL = dParticleID->Get_BestFCALMatchParams(locTrackTimeBased, locDetectorMatches, locFCALShowerMatchParams);
				if(foundFCAL){
					const DFCALShower *locFCALShower = NULL;
					locNeutralShowers[loc_j]->GetSingleT(locFCALShower);

					if(locFCALShowerMatchParams.dFCALShower == locFCALShower) {
						nMatched++;
						FillShower(locNeutralShowers[loc_j], true, locBeamPhotonTime);
					}
				}
			}
			if(locNeutralShowers[loc_j]->dDetectorSystem == SYS_BCAL) {
				DBCALShowerMatchParams locBCALShowerMatchParams;
				bool foundBCAL = dParticleID->Get_BestBCALMatchParams(locTrackTimeBased, locDetectorMatches, locBCALShowerMatchParams);
				if(foundBCAL){
					const DBCALShower *locBCALShower = NULL;
					locNeutralShowers[loc_j]->GetSingleT(locBCALShower);

					if(locBCALShowerMatchParams.dBCALShower == locBCALShower) {
						nMatched++;
						FillShower(locNeutralShowers[loc_j], true, locBeamPhotonTime);
					}
				}
			}
                }
		
		
		// plot properties of unused showers
                if(nMatched==0) {
		
		}
*/
