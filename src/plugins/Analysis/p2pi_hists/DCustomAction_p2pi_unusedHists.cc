// $Id$
//
//    File: DCustomAction_p2pi_unusedHists.cc
// Created: Thu Jan 22 08:06:18 EST 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "DCustomAction_p2pi_unusedHists.h"

void DCustomAction_p2pi_unusedHists::Initialize(JEventLoop* locEventLoop)
{

	// get PID algos
	const DParticleID* locParticleID = NULL;
        locEventLoop->GetSingle(locParticleID);
	dParticleID = locParticleID;

	vector<const DFCALGeometry*> locFCALGeometry;
	locEventLoop->Get(locFCALGeometry);
	dFCALGeometry = locFCALGeometry[0];
	
	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//Required: Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this action.
			//If another thread has already created the folder, it just changes to it. 
		CreateAndChangeTo_ActionDirectory();

		dNShowerBCAL_FCAL = GetOrCreate_Histogram<TH2I>("NShowerBCAL_FCAL", "Number of non-matched showers; FCAL; BCAL", 10, 0, 10, 10, 0, 10);
		dHistMap_FCALShowerDeltaR_P = GetOrCreate_Histogram<TH2I>("FCALShowerDeltaR_P", "Unused FCAL showers #DeltaR vs momentum; momentum; #DeltaR", 120, 0., 6., 200, 0, 50); 
		dHistMap_FCALShowerDeltaR_Theta = GetOrCreate_Histogram<TH2I>("FCALShowerDeltaR_Theta", "Unused FCAL showers #DeltaR vs #theta; #theta; #DeltaR", 150, 0, 15, 200, 0, 50);
		dHistMap_FCALShowerDeltaD_P = GetOrCreate_Histogram<TH2I>("FCALShowerDeltaD_P", "Unused FCAL showers #DeltaD vs momentum; momentum; #DeltaR", 120, 0., 6., 200, 0, 50); 
		dHistMap_FCALShowerDeltaD_Theta = GetOrCreate_Histogram<TH2I>("FCALShowerDeltaD_Theta", "Unused FCAL showers #DeltaD vs #theta; #theta; #DeltaD", 150, 0, 15, 200, 0, 50);
		dHistMap_FCALShowerDeltaD_DeltaT = GetOrCreate_Histogram<TH2I>("FCALShowerDeltaD_DeltaT", "Unused FCAL showers #DeltaD vs #Deltat; #Deltat; #DeltaD", 200, -100., 100., 200, 0, 50);

		dHistMap_BCALShowerDeltaPhi_DeltaZ = GetOrCreate_Histogram<TH2I>("BCALShowerDeltaPhi_DeltaZ", "Unused BCAL showers #Delta#phi vs #DeltaZ; #DeltaZ; #Delta#phi", 200, -50, 50, 360, -180, 180);
		dHistMap_BCALShowerDeltaPhi_P = GetOrCreate_Histogram<TH2I>("BCALShowerDeltaPhi_P", "Unused BCAL showers #Delta#phi vs momentum; momentum; #Delta#phi", 120, 0., 6., 360, -180, 180); 
		dHistMap_BCALShowerDeltaPhi_Theta = GetOrCreate_Histogram<TH2I>("BCALShowerDeltaPhi_Theta", "Unused BCAL showers #Delta#phi vs #theta; #theta; #Delta#phi", 150, 0, 150, 360, -180, 180);
		dHistMap_BCALShowerDeltaPhi_DeltaT = GetOrCreate_Histogram<TH2I>("BCALShowerDeltaPhi_DeltaT", "Unused BCAL showers #Delta#phi vs #Deltat; #Deltat; #Delta$phi", 200, -100., 100., 360, -180, 180);
		dHistMap_BCALShowerDeltaZ_DeltaT = GetOrCreate_Histogram<TH2I>("BCALShowerDeltaZ_DeltaT", "Unused FCAL showers #DeltaZ vs #Deltat; #Deltat; #DeltaZ", 200, -100., 100., 200, -50, 50);

		int nbins = 100;
		h1_deltaX = GetOrCreate_Histogram<TH1I>("h1_deltaX", "Fcal X - Track X", nbins,-25,25);
		h1_deltaY = GetOrCreate_Histogram<TH1I>("h1_deltaY", "Fcal Y - Track Y", nbins,-25,25);
		h1_Efcal = GetOrCreate_Histogram<TH1I>("h1_Efcal", "Hit energy", nbins,0,4);
		h1_tfcal = GetOrCreate_Histogram<TH1I>("h1_tfcal", "Hit time", 250,-25,225);
	       
		h1_N9 = GetOrCreate_Histogram<TH1I>("h1_N9", "Hit N9", 25, 0, 25);
		h1_E9 = GetOrCreate_Histogram<TH1I>("h1_E9", "Energy E9 (GeV)", nbins, 0, 2);
		h1_t9 = GetOrCreate_Histogram<TH1I>("h1_t9", "Time t9 (ns)", 250, -25, 225);
		h1_t9sigma = GetOrCreate_Histogram<TH1I>("h1_t9sigma", "Time t9sigma",nbins,0,10);
		
		h1_N1 = GetOrCreate_Histogram<TH1I>("h1_N1", "Hit N1", 25, 0, 25);
		h1_E1 = GetOrCreate_Histogram<TH1I>("h1_E1", "Energy E1 (GeV)", nbins, 0, 2);
		h1_t1 = GetOrCreate_Histogram<TH1I>("h1_t1", "Time t1 (ns)", 250, -25, 225);
		h1_t1sigma = GetOrCreate_Histogram<TH1I>("h1_t1sigma", "Time t1sigma",nbins,0,10);

		h2_YvsX9 = GetOrCreate_Histogram<TH2I>("h2_YvsX9", "Hit position Y vs X, E9",240,-120,120,240,-120,120);
		h2_dEdx9_vsp= GetOrCreate_Histogram<TH2I>("h2_dEdx9_vsp", "Track dEdx9 vs p",nbins,0,4,nbins,0,10);
		h2_E9vsp= GetOrCreate_Histogram<TH2I>("h2_E9vsp", "E9 vs p",nbins,0,4,nbins,0,4);
		h2_dEdxvsE9= GetOrCreate_Histogram<TH2I>("h2_dEdxvsE9", "dEdx vs E9 energy",nbins,0,4,nbins,0,4);
		h2_E9_vsTheta= GetOrCreate_Histogram<TH2I>("h2_E9_vsTheta", "E9 vs Theta",90,0,30,nbins,0,4);
		h2_E9_vsPhi= GetOrCreate_Histogram<TH2I>("h2_E9_vsPhi", "E9 vs Phi",90,-180,180,nbins,0,4);
		
		h2_YvsX1 = GetOrCreate_Histogram<TH2I>("h2_YvsX1", "Hit position Y vs X, E1",240,-120,120,240,-120,120);
		h2_dEdx1_vsp = GetOrCreate_Histogram<TH2I>("h2_dEdx1_vsp", "Track dEdx1 vs p",nbins,0,4,nbins,0,10);
		h2_E1vsp = GetOrCreate_Histogram<TH2I>("h2_E1vsp", "E1 vs p",nbins,0,4,nbins,0,4);
		h2_dEdxvsE1 = GetOrCreate_Histogram<TH2I>("h2_dEdxvsE1", "dEdx vs E1 energy",nbins,0,4,nbins,0,4);
		h2_E1_vsTheta = GetOrCreate_Histogram<TH2I>("h2_E1_vsTheta", "E1 vs Theta",90,0,30,nbins,0,4);
		h2_E1_vsPhi = GetOrCreate_Histogram<TH2I>("h2_E1_vsPhi", "E1 vs Phi",90,-180,180,nbins,0,4);
		h2_E1vsRlocal = GetOrCreate_Histogram<TH2I>("h2_E1vsRlocal", "E1 vs Rtrk rel to Block center (cm)",nbins,0,5,nbins,0,4);
		h2_YvsXcheck = GetOrCreate_Histogram<TH2I>("h2_E1vsRlocal", "E1 vs Rtrk rel to Block center (cm)",nbins,0,5,nbins,0,4);

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
				dHistMap_TrackNhits_Theta[locMatch][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 90, 0., 90., 35, 0, 35);

				locHistName = "TrackNhitsCDC_Theta" + match[locDummy] + charge[locCharge];
				locHistTitle = "Number of CDC ring hits on track vs #theta: " + match[locDummy] + charge[locCharge];
				locHistTitle += "; #theta; # CDC ring hits";
				dHistMap_TrackNhitsCDC_Theta[locMatch][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 90, 0., 90., 35, 0, 35);

				locHistName = "TrackNhitsFDC_Theta" + match[locDummy] + charge[locCharge];
				locHistTitle = "Number of FDC plane hits on track vs #theta: " + match[locDummy] + charge[locCharge];
				locHistTitle += "; #theta; # FDC plane hits";
				dHistMap_TrackNhitsFDC_Theta[locMatch][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 90, 0., 90., 35, 0, 35);

				locHistName = "TrackNhitsFDCVsCDC_Theta13_16" + match[locDummy] + charge[locCharge];
				locHistTitle = "Number of FDC plane hits vs CDC ring hits (13<#theta<16): " + match[locDummy] + charge[locCharge];
				locHistTitle += "; # CDC ring hits; # FDC plane hits";
				dHistMap_TrackNhitsFDCVsCDC_Theta13_16[locMatch][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 35, 0, 35, 35, 0, 35);

				locHistName = "TrackNhitsFDCVsCDC_Theta16_20" + match[locDummy] + charge[locCharge];
				locHistTitle = "Number of FDC plane hits vs CDC ring hits (16<#theta<20): " + match[locDummy] + charge[locCharge];
				locHistTitle += "; # CDC ring hits; # FDC plane hits";
				dHistMap_TrackNhitsFDCVsCDC_Theta16_20[locMatch][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 35, 0, 35, 35, 0, 35);

				locHistName = "TrackChiSq_Theta" + match[locDummy] + charge[locCharge];
				locHistTitle = "FOM for track #chi^{2} vs #theta: " + match[locDummy] + charge[locCharge];
				locHistTitle += "; #theta; #chi^{2}";
				dHistMap_TrackChiSq_Theta[locMatch][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 90, 0., 90., 500, 0., 100.);

				locHistName = "TrackFOM_Theta" + match[locDummy] + charge[locCharge];
				locHistTitle = "FOM for track fit vs #theta: " + match[locDummy] + charge[locCharge];
				locHistTitle += "; #theta; FOM";
				dHistMap_TrackFOM_Theta[locMatch][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 90, 0., 90., 1000, 0., 1.);

				locHistName = "TrackP_Theta" + match[locDummy] + charge[locCharge];
				locHistTitle = "P vs #theta: " + match[locDummy] + charge[locCharge];
				locHistTitle += "; #theta; momentum";
				dHistMap_TrackP_Theta[locMatch][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 90, 0., 90., 120, 0., 3.);

				locHistName = "TrackPOCAXY" + match[locDummy] + charge[locCharge];
				locHistTitle = "POCA to beamline Y vs X: " + match[locDummy] + charge[locCharge];
				locHistTitle += "; POCA X; POCAY";
				dHistMap_TrackPOCAXY[locMatch][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 200, -20., 20., 200, -20., 20.);

				locHistName = "TrackPOCAZ" + match[locDummy] + charge[locCharge];
				locHistTitle = "POCA to beamline Z: " + match[locDummy] + charge[locCharge];
				locHistTitle += "; POCA Z";
				dHistMap_TrackPOCAZ[locMatch][locCharge] = GetOrCreate_Histogram<TH1I>(locHistName.Data(), locHistTitle.Data(), 500, 0., 250.);

				locHistName = "TrackCDCHitRadius_NHits" + match[locDummy] + charge[locCharge];
				locHistTitle = "CDC Hit radius vs # hits: " + match[locDummy] + charge[locCharge];
				locHistTitle += "; # hits; CDC Hit Radius";
				dHistMap_TrackCDCHitRadius_Nhits[locMatch][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 50, 0, 50, 50, 0, 50);
				
				locHistName = "HighNDFTrackCDCHitRadius_PT" + match[locDummy] + charge[locCharge];
				locHistTitle = "High NDF tracks CDC Hit radius vs P_{T}: " + match[locDummy] + charge[locCharge];
				locHistTitle += "; P_{T} (GeV); CDC Hit Radius";
				dHistMap_HighNDFTrackCDCHitRadius_PT[locMatch][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 120, 0., 3., 50, 0, 50);

				locHistName = "LowNDFTrackCDCHitRadius_PT" + match[locDummy] + charge[locCharge];
				locHistTitle = "Low NDF tracks CDC Hit radius vs P_{T}: " + match[locDummy] + charge[locCharge];
				locHistTitle += "; P_{T} (GeV); CDC Hit Radius";
				dHistMap_LowNDFTrackCDCHitRadius_PT[locMatch][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 120, 0., 3., 50, 0, 50);

				locHistName = "LowNDFTrackP_VertZ" + match[locDummy] + charge[locCharge];
				locHistTitle = "Low NDF tracks P vs vertex Z: " + match[locDummy] + charge[locCharge];
				locHistTitle += "; Vertex Z (cm); momentum";
				dHistMap_LowNDFTrackP_VertZ[locMatch][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 500, 0., 250., 120, 0., 3.);

				locHistName = "LowNDFTrackPT_Phi" + match[locDummy] + charge[locCharge];
                                locHistTitle = "Low NDF tracks P_{T} vs #phi: " + match[locDummy] + charge[locCharge];
                                locHistTitle += "; #phi; transverse momentum";			
				dHistMap_LowNDFTrackPT_Phi[locMatch][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 360., -180., 180., 120, 0., 3.);

				locHistName = "TrackMCHitFraction_Theta" + match[locDummy] + charge[locCharge];
                                locHistTitle = "Track MC Hit Fraction vs #theta: " + match[locDummy] + charge[locCharge];
                                locHistTitle += "; #theta; MC Hit Fraction";
				dHistMap_TrackMCHitFraction_Theta[locMatch][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 90., 0., 90., 100, 0., 1.);
				
				locHistName = "TrackMCMomRes_Theta" + match[locDummy] + charge[locCharge];
                                locHistTitle = "Track MC Momentum Resolution vs #theta: " + match[locDummy] + charge[locCharge];
                                locHistTitle += "; #theta; MC Momentum Resolution #Delta p/p";
				dHistMap_TrackMCMomRes_Theta[locMatch][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 90., 0., 90., 100, -0.4, 0.4);

				locHistName = "TrackMCThetaRes_Theta" + match[locDummy] + charge[locCharge];
                                locHistTitle = "Track MC Theta Resolution vs #theta: " + match[locDummy] + charge[locCharge];
                                locHistTitle += "; #theta; MC Theta Resolution #Delta#theta (degrees)";
				dHistMap_TrackMCThetaRes_Theta[locMatch][locCharge] = GetOrCreate_Histogram<TH2I>(locHistName.Data(), locHistTitle.Data(), 90., 0., 90., 100, -5, 5);
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

bool DCustomAction_p2pi_unusedHists::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	
	// should only have one reaction step
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(0);

	// get beam photon energy and final state particles
	auto locBeamPhoton = Get_UseKinFitResultsFlag() ? locParticleComboStep->Get_InitialParticle() : locParticleComboStep->Get_InitialParticle_Measured();
	auto locParticles = Get_UseKinFitResultsFlag() ? locParticleComboStep->Get_FinalParticles() : locParticleComboStep->Get_FinalParticles_Measured();
	double locBeamPhotonTime = locBeamPhoton->time();

	// detector matches for charged track -to- shower matching
	const DDetectorMatches* locDetectorMatches = NULL;
        locEventLoop->GetSingle(locDetectorMatches);

	// thrown particles for track resolution
	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
        locEventLoop->Get(locMCThrownMatchingVector);
	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	// get all charged tracks and showers for comparison to combo (no "pre-select" cut here)
        vector<const DChargedTrack*> locChargedTracks;
        locEventLoop->Get(locChargedTracks); //, "PreSelect"); 
	vector<const DNeutralShower*> locNeutralShowers;
        locEventLoop->Get(locNeutralShowers); //, "PreSelect"); 
	vector<const DFCALShower*> locFCALShowers;
        locEventLoop->Get(locFCALShowers); 
	vector<const DBCALShower*> locBCALShowers;
        locEventLoop->Get(locBCALShowers); 

        vector<const DChargedTrack*> locUnusedChargedTracks;
	vector<const DNeutralShower*> locUnusedNeutralShowers;

	// loop over charged tracks to make list of unused
        for(size_t loc_j = 0; loc_j < locChargedTracks.size(); ++loc_j) {
                int nMatched = 0;

		double locMatchFOM;
		const DMCThrown* locMCThrown = NULL; 
		if(!locMCThrowns.empty()) 
			locMCThrown = locMCThrownMatchingVector[0]->Get_MatchingMCThrown(locChargedTracks[loc_j], locMatchFOM);

                // add tracks not in combo to vector of unused tracks
                for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i) {
			
                        const DChargedTrack* locChargedTrack = static_cast<const DChargedTrack*>(locParticleComboStep->Get_FinalParticle_SourceObject(loc_i));
                        if(locChargedTrack == NULL) continue; // should never happen

			if(locChargedTracks[loc_j]->candidateid == locChargedTrack->candidateid) {
                                nMatched++;
				FillTrack(locEventLoop, locChargedTracks[loc_j], true, locMCThrown);
			}
                }
		
		// plot properties of unused charged tracks
                if(nMatched==0) {
			locUnusedChargedTracks.push_back(locChargedTracks[loc_j]);
			FillTrack(locEventLoop, locChargedTracks[loc_j], false, locMCThrown);
		}
        }

	// loop over BCAL and FCAL showers to fill histograms for matched tracks
	for(size_t loc_j = 0; loc_j < locBCALShowers.size(); ++loc_j) {

		// add showers matched to tracks in combo
                for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i) {

                        const DChargedTrack* locChargedTrack = static_cast<const DChargedTrack*>(locParticleComboStep->Get_FinalParticle_SourceObject(loc_i));
                        if(locChargedTrack == NULL) continue; // should never happen

                        const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_BestFOM();
	auto locTrackTimeBased = locChargedTrackHypothesis->Get_TrackTimeBased();
			const DReferenceTrajectory* rt = locTrackTimeBased->rt;
			if(rt == NULL) continue;

			const DBCALShower *locBCALShower = locBCALShowers[loc_j];
			DVector3 bcal_pos(locBCALShower->x, locBCALShower->y, locBCALShower->z);
				
			double locFlightTime = 9.9E9, locPathLength = 9.9E9, locFlightTimeVariance = 9.9E9;
			rt->DistToRTwithTime(bcal_pos, &locPathLength, &locFlightTime, &locFlightTimeVariance, SYS_BCAL);
			
			DVector3 proj_pos = rt->GetLastDOCAPoint();
			if(proj_pos.Perp() < 65.0)
				continue;  // not inside BCAL!
			
			double dz = bcal_pos.z() - proj_pos.z();
			double dphi = bcal_pos.Phi() - proj_pos.Phi();
			while(dphi >    TMath::Pi())
				dphi -= 2*TMath::Pi();
			while(dphi < -1.*TMath::Pi())
				dphi += 2*TMath::Pi();
			
			double locDeltaT = locBCALShowers[loc_j]->t - locBeamPhotonTime;

			//FILL HISTOGRAMS
			//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
			//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
			Lock_Action(); //ACQUIRE ROOT LOCK!!
			{
				dHistMap_BCALShowerDeltaPhi_DeltaZ->Fill(dz, dphi*180./TMath::Pi());
				dHistMap_BCALShowerDeltaPhi_P->Fill(locTrackTimeBased->momentum().Mag(), dphi*180./TMath::Pi());
				dHistMap_BCALShowerDeltaPhi_Theta->Fill(locTrackTimeBased->momentum().Theta()*180./TMath::Pi(), dphi*180./TMath::Pi());
				dHistMap_BCALShowerDeltaPhi_DeltaT->Fill(locDeltaT, dphi*180./TMath::Pi());
				dHistMap_BCALShowerDeltaZ_DeltaT->Fill(locDeltaT, dz);
			}
			Unlock_Action(); //RELEASE ROOT LOCK!!

			shared_ptr<const DBCALShowerMatchParams> locBCALShowerMatchParams;
			bool foundBCAL = dParticleID->Get_BestBCALMatchParams(locTrackTimeBased, locDetectorMatches, locBCALShowerMatchParams);
			if(foundBCAL){
				if(locBCALShowerMatchParams->dBCALShower == locBCALShower) {
					DNeutralShower* locNeutralShower = new DNeutralShower();
					locNeutralShower->dDetectorSystem = SYS_BCAL;
					locNeutralShower->dEnergy = locBCALShowers[loc_j]->E;
					locNeutralShower->dSpacetimeVertex.SetXYZT(locBCALShowers[loc_j]->x, locBCALShowers[loc_j]->y, locBCALShowers[loc_j]->z, locBCALShowers[loc_j]->t);
					locNeutralShower->AddAssociatedObject(locBCALShowers[loc_j]);
					double locFlightTime = locBCALShowerMatchParams->dFlightTime;
					FillShower(locNeutralShower, true, locBeamPhotonTime, locFlightTime);
					delete locNeutralShower;
				}
			}
		}
	}

	for(size_t loc_j = 0; loc_j < locFCALShowers.size(); ++loc_j) {

		// add showers matched to tracks in combo
                for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i) {
			
                        const DChargedTrack* locChargedTrack = static_cast<const DChargedTrack*>(locParticleComboStep->Get_FinalParticle_SourceObject(loc_i));
                        if(locChargedTrack == NULL) continue; // should never happen
			
                        const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_BestFOM();
	auto locTrackTimeBased = locChargedTrackHypothesis->Get_TrackTimeBased();
			const DReferenceTrajectory* rt = locTrackTimeBased->rt;
			if(rt == NULL) continue;

			const DFCALShower *locFCALShower = locFCALShowers[loc_j];
			DVector3 fcal_pos = locFCALShower->getPosition();
			DVector3 norm(0.0, 0.0, 1.0); //normal vector for the FCAL plane
			DVector3 proj_pos, proj_mom;
			double locPathLength = 9.9E9, locFlightTime = 9.9E9, locFlightTimeVariance = 9.9e9;
			if(rt->GetIntersectionWithPlane(fcal_pos, norm, proj_pos, proj_mom, &locPathLength, &locFlightTime, &locFlightTimeVariance, SYS_FCAL) != NOERROR)
				continue;
			
			double dd = (fcal_pos - proj_pos).Mag();
			double dr = (fcal_pos - proj_pos).Perp();
			double dphi = fcal_pos.Phi() - proj_pos.Phi();
			while(dphi >    TMath::Pi())
				dphi -= 2*TMath::Pi();
			while(dphi < -1.*TMath::Pi())
				dphi += 2*TMath::Pi();
			
			double locDeltaT = locFCALShowers[loc_j]->getTime() - locBeamPhotonTime;

			//FILL HISTOGRAMS
			//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
			//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
			Lock_Action(); //ACQUIRE ROOT LOCK!!
			{
				dHistMap_FCALShowerDeltaR_P->Fill(locTrackTimeBased->momentum().Mag(), dr);
				dHistMap_FCALShowerDeltaR_Theta->Fill(locTrackTimeBased->momentum().Theta()*180./TMath::Pi(), dr);
				dHistMap_FCALShowerDeltaD_P->Fill(locTrackTimeBased->momentum().Mag(), dd);
				dHistMap_FCALShowerDeltaD_Theta->Fill(locTrackTimeBased->momentum().Theta()*180./TMath::Pi(), dd);
				dHistMap_FCALShowerDeltaD_DeltaT->Fill(locDeltaT, dd);
			}
			Unlock_Action(); //RELEASE ROOT LOCK!!

			shared_ptr<const DFCALShowerMatchParams> locFCALShowerMatchParams;
			bool foundFCAL = dParticleID->Get_BestFCALMatchParams(locTrackTimeBased, locDetectorMatches, locFCALShowerMatchParams);
			if(foundFCAL){
				if(locFCALShowerMatchParams->dFCALShower == locFCALShower) {
					DNeutralShower* locNeutralShower = new DNeutralShower();
					locNeutralShower->dDetectorSystem = SYS_FCAL;
					locNeutralShower->dEnergy = locFCALShowers[loc_j]->getEnergy();
					locNeutralShower->dSpacetimeVertex.SetVect(locFCALShowers[loc_j]->getPosition());
					locNeutralShower->dSpacetimeVertex.SetT(locFCALShowers[loc_j]->getTime());
					locNeutralShower->AddAssociatedObject(locFCALShowers[loc_j]);
					double locFlightTime = locFCALShowerMatchParams->dFlightTime;
					FillShower(locNeutralShower, true, locBeamPhotonTime, locFlightTime);
					delete locNeutralShower;
				}
			}
		}
	}

	int nUnmatchedFCAL = 0;
	int nUnmatchedBCAL = 0;

	/////////////////////////////////////////////////////////////////////////////////////////
	// loop over neutral showers which should all be unused (unless bad matching criteria) //
	/////////////////////////////////////////////////////////////////////////////////////////
	for(size_t loc_j = 0; loc_j < locNeutralShowers.size(); ++loc_j) {
		locUnusedNeutralShowers.push_back(locNeutralShowers[loc_j]);
		double locFlightTime = locNeutralShowers[loc_j]->dSpacetimeVertex.Vect().Mag()/SPEED_OF_LIGHT;
		FillShower(locNeutralShowers[loc_j], false, locBeamPhotonTime, locFlightTime);
		
		if(locNeutralShowers[loc_j]->dDetectorSystem == SYS_FCAL) {
			nUnmatchedFCAL++;
		}
		if(locNeutralShowers[loc_j]->dDetectorSystem == SYS_BCAL) {
			nUnmatchedBCAL++;
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


void DCustomAction_p2pi_unusedHists::FillTrack(JEventLoop* locEventLoop, const DChargedTrack* locChargedTrack, bool locMatch, const DMCThrown* locMCThrown)
{

	const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_BestFOM();
	int locCharge = locChargedTrackHypothesis->charge();

	auto locTrackTimeBased = locChargedTrackHypothesis->Get_TrackTimeBased();
	
	double nHits = locTrackTimeBased->Ndof + 5.;
	double locTheta = locTrackTimeBased->momentum().Theta()*180/TMath::Pi();

	set<int> locCDCRings;
	dParticleID->Get_CDCRings(locTrackTimeBased->dCDCRings, locCDCRings);
	
	set<int> locFDCPlanes;
	dParticleID->Get_FDCPlanes(locTrackTimeBased->dFDCPlanes, locFDCPlanes);

	double locHitFraction = locTrackTimeBased->dNumHitsMatchedToThrown/nHits;
	double locMomentumRes = -999;
	double locThetaRes = -999;
	if(locMCThrown) {
	  locMomentumRes = (locTrackTimeBased->momentum().Mag() - locMCThrown->momentum().Mag())/locMCThrown->momentum().Mag();
	  locThetaRes = (locTrackTimeBased->momentum().Theta() - locMCThrown->momentum().Theta())*180./TMath::Pi();
	}

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action(); //ACQUIRE ROOT LOCK!!
	{
		dHistMap_TrackNhits_Theta[locMatch][locCharge]->Fill(locTheta, nHits);
		dHistMap_TrackNhitsCDC_Theta[locMatch][locCharge]->Fill(locTheta, locCDCRings.size());
		dHistMap_TrackNhitsFDC_Theta[locMatch][locCharge]->Fill(locTheta, locFDCPlanes.size());
		if(locTheta>13. && locTheta<16.)
		  dHistMap_TrackNhitsFDCVsCDC_Theta13_16[locMatch][locCharge]->Fill(locCDCRings.size(), locFDCPlanes.size());
		if(locTheta>16. && locTheta<20.)
		  dHistMap_TrackNhitsFDCVsCDC_Theta16_20[locMatch][locCharge]->Fill(locCDCRings.size(), locFDCPlanes.size());
		dHistMap_TrackChiSq_Theta[locMatch][locCharge]->Fill(locTheta, locTrackTimeBased->chisq);
		dHistMap_TrackFOM_Theta[locMatch][locCharge]->Fill(locTheta, locTrackTimeBased->FOM);
		dHistMap_TrackP_Theta[locMatch][locCharge]->Fill(locTheta, locTrackTimeBased->momentum().Mag());
		dHistMap_TrackPOCAXY[locMatch][locCharge]->Fill(locTrackTimeBased->position().X(), locTrackTimeBased->position().Y());
		dHistMap_TrackPOCAZ[locMatch][locCharge]->Fill(locTrackTimeBased->position().Z());

		dHistMap_TrackMCHitFraction_Theta[locMatch][locCharge]->Fill(locTrackTimeBased->momentum().Theta()*180.0/TMath::Pi(), locHitFraction);
		if(locMCThrown) {
		  dHistMap_TrackMCMomRes_Theta[locMatch][locCharge]->Fill(locTheta, locMomentumRes);
		  dHistMap_TrackMCThetaRes_Theta[locMatch][locCharge]->Fill(locTheta, locThetaRes);
		}

		for(set<int>::iterator locIterator = locFDCPlanes.begin(); locIterator != locFDCPlanes.end(); ++locIterator) {
		  //dHistMap_TrackFDCPlaneVs_Nhits->Fill(nHits, *locIterator);
		}

		// understand low Nhits tracks
		if(locTrackTimeBased->momentum().Theta()*180./TMath::Pi() > 30.){
			for(set<int>::iterator locIterator = locCDCRings.begin(); locIterator != locCDCRings.end(); ++locIterator) {
				dHistMap_TrackCDCHitRadius_Nhits[locMatch][locCharge]->Fill(nHits, *locIterator);
				if(nHits < 20) dHistMap_LowNDFTrackCDCHitRadius_PT[locMatch][locCharge]->Fill(locTrackTimeBased->momentum().Perp(), *locIterator);
				else dHistMap_HighNDFTrackCDCHitRadius_PT[locMatch][locCharge]->Fill(locTrackTimeBased->momentum().Perp(), *locIterator);
			}

			if(nHits < 20){
				dHistMap_LowNDFTrackP_VertZ[locMatch][locCharge]->Fill(locTrackTimeBased->position().Z(), locTrackTimeBased->momentum().Mag());
				dHistMap_LowNDFTrackPT_Phi[locMatch][locCharge]->Fill(locTrackTimeBased->momentum().Phi()*180/TMath::Pi(), locTrackTimeBased->momentum().Perp());
			}
		}
	}
	Unlock_Action(); //RELEASE ROOT LOCK!!

	if(!locMatch) return;

	// do FCAL energy sums to compare with Elton's E9 and E25 in fcal_charged
	Double_t zfcal=625;
        DVector3 fcal_origin(0.0,0.0,zfcal);
        DVector3 fcal_normal(0.0,0.0,1.0);
	DVector3 trkpos(0.0,0.0,0.0);
        DVector3 proj_mom(0.0,0.0,0.0);
	double theta = locTrackTimeBased->momentum().Theta()*180./TMath::Pi();
	double phi = locTrackTimeBased->momentum().Phi()*180./TMath::Pi();
	double p = locTrackTimeBased->momentum().Mag();
        double dEdx = locTrackTimeBased->dEdx()*1e6;

	vector<const DFCALHit*> fcalhits;
        locEventLoop->Get(fcalhits);
	if(fcalhits.empty()) return;

	if (locTrackTimeBased->rt->GetIntersectionWithPlane(fcal_origin,fcal_normal,trkpos,proj_mom,NULL,NULL,NULL,SYS_FCAL)==NOERROR){
		double trkposX = trkpos.X();
		double trkposY = trkpos.Y();
		int trkrow = dFCALGeometry->row((float)trkposY);
		int trkcol = dFCALGeometry->column((float)trkposX);
		
		// loop over fcal hits

		double E9=0;  // energy, x, y of 9 blocks surrounding track
		double E9peak=0;
		double x9=0;
		double y9=0;
		double t9=0;
		double t9sq=0;
		double t9sigma=0;
		int N9=0;
		int Delta_block=1;   // =1 for 3x3, =2 for 5x5
		double dX_E1=-1000;
		double dY_E1=-1000;
		
		int row_offset = 0;   // offset actual row to look for randoms
		int col_offset = 0;
		trkrow += row_offset;
		trkcol += col_offset;
		
		for (unsigned int j=0; j < fcalhits.size(); ++j) {
			int row = fcalhits[j]->row;
			int col = fcalhits[j]->column;
			double x = fcalhits[j]->x;
			double y = fcalhits[j]->y;
			double Efcal = fcalhits[j]->E;
			double tfcal= fcalhits[j]->t;
			double intOverPeak = fcalhits[j]->intOverPeak;
			
			// fill histograms
			int drow = abs(row - trkrow);
			int dcol = abs(col - trkcol);
			
			h1_deltaX->Fill(x - trkposX);
			h1_deltaY->Fill(y - trkposY);
			h1_Efcal->Fill(Efcal);
			h1_tfcal->Fill(tfcal);
		
			// select hits 
			if (drow<=Delta_block && dcol<=Delta_block && (tfcal>-15 && tfcal<50) && (intOverPeak>2.5 && intOverPeak<9)) {
				E9 += Efcal;
				E9peak += Efcal*6/intOverPeak;   // factor of 6 so that E9peak ~ E9
				x9 += x;
				y9 += y;
				t9 += tfcal;
				t9sq += tfcal*tfcal;
				N9 += 1;

				dX_E1 = x - trkposX;
				dY_E1 = y - trkposY;
			}
			
		} // end loop over fcal hits

//		x9 = N9>0? x9/N9 : 0;
//		y9 = N9>0? y9/N9 : 0;
		t9 = N9>0? t9/N9 : 0;
		t9sigma = N9>0? sqrt(t9sq/N9 - t9*t9): 0;
		
		//FILL HISTOGRAMS
		//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
		//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
		Lock_Action(); //ACQUIRE ROOT LOCK!!
		{
			if (N9>0) {
				h1_N9->Fill(N9);
				h1_E9->Fill(E9);
				h1_t9->Fill(t9);
				h1_t9sigma->Fill(t9sigma);
				h2_YvsX9->Fill(trkposX,trkposY);
				h2_dEdx9_vsp->Fill(p,dEdx);
				h2_E9vsp->Fill(p,E9);
				h2_dEdxvsE9->Fill(E9,dEdx);
				h2_E9_vsTheta->Fill(theta,E9);
				h2_E9_vsPhi->Fill(phi,E9);
			}
			if (N9==1) {
				h1_N1->Fill(N9);
				h1_E1->Fill(E9);
				h1_t1->Fill(t9);
				h1_t1sigma->Fill(t9sigma);
				h2_YvsX1->Fill(trkposX,trkposY);
				h2_dEdx1_vsp->Fill(p,dEdx);
				h2_E1vsp->Fill(p,E9);
				h2_dEdxvsE1->Fill(E9,dEdx);
				h2_E1_vsTheta->Fill(theta,E9);
				h2_E1_vsPhi->Fill(phi,E9);
				double Rlocal = sqrt(dX_E1*dX_E1 + dY_E1*dY_E1);
				h2_E1vsRlocal->Fill(Rlocal,E9);
				h2_YvsXcheck->Fill(dX_E1,dY_E1);
			  }
		}
		Unlock_Action(); //RELEASE ROOT LOCK!!
	}

	return;
}


void DCustomAction_p2pi_unusedHists::FillShower(const DNeutralShower* locNeutralShower, bool locMatch, double locBeamPhotonTime, double locFlightTime)
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
