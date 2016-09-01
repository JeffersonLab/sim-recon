// $Id$
//
//    File: DCustomAction_p2k_hists.cc
// Created: Wed Jan 21 16:53:41 EST 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "DCustomAction_p2k_hists.h"

void DCustomAction_p2k_hists::Initialize(JEventLoop* locEventLoop)
{
	DApplication* dapp=dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
        JCalibration *jcalib = dapp->GetJCalibration((locEventLoop->GetJEvent()).GetRunNumber());

        // Parameters for event selection to fill histograms
        endpoint_energy = 12.;
        map<string, double> photon_endpoint_energy;
        if(jcalib->Get("/PHOTON_BEAM/endpoint_energy", photon_endpoint_energy) == false) {
                endpoint_energy = photon_endpoint_energy["PHOTON_BEAM_ENDPOINT_ENERGY"];
        }
        endpoint_energy_bins = (int)(20*endpoint_energy);

        cohmin_energy = 0.;
        cohedge_energy = 12.;
        map<string, double> photon_beam_param;
        if(jcalib->Get("/ANALYSIS/beam_asymmetry/coherent_energy", photon_beam_param) == false) {
                cohmin_energy = photon_beam_param["cohmin_energy"];
                cohedge_energy = photon_beam_param["cohedge_energy"];
        }

        dEdxCut = 2.2;
        minMM2Cut = -0.01;
        maxMM2Cut = 0.01;
        maxPhiMassCut = 1.05;

	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//Required: Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this action.
			//If another thread has already created the folder, it just changes to it. 
		CreateAndChangeTo_ActionDirectory();
		
		dEgamma = GetOrCreate_Histogram<TH1I>("Egamma", "TAGGER photon energy; E_{#gamma}", endpoint_energy_bins, 0., endpoint_energy);
		dInvariantMass = GetOrCreate_Histogram<TH1I>("InvariantMass", "Invariant Mass of K^{+}K^{-};Invariant Mass of K^{+}K^{-} (GeV/c^{2}",270,0.9,1.8);
		dMissingMassSq = GetOrCreate_Histogram<TH1I>("MissingMassSq", "Missing Mass Sq. of K^{+}K^{-};MM^{2}_{K^{+}K^{-}} (GeV/c^{2}",600,-0.6,0.6);
		dKplus_deltaInvBeta_P = GetOrCreate_Histogram<TH2I>("Kplus_deltaInvBeta_P","K+ Inverse #Delta#beta vs p; p; 1/#beta_{exp} - 1/#beta_{meas}",200,0,10,500,-0.5,0.5);
		dKminus_deltaInvBeta_P = GetOrCreate_Histogram<TH2I>("Kminus_deltaInvBeta_P","K- Inverse #Delta#beta vs p; p; 1/#beta_{exp} - 1/#beta_{meas}",200,0,10,500,-0.5,0.5);

		dMM2_M2k = GetOrCreate_Histogram<TH2I>("MM2_M2k", "MM^{2} off k^{+}k^{-} vs M_{k^{+}k^{-}}; M_{k^{+}k^{-}}; MM^{2}", 800, 0.0, 2.0, 200, -1., 1.);
		dProton_dEdx_P = GetOrCreate_Histogram<TH2I>("Proton_dEdx_P","dE/dx vs p; p; dE/dx",200,0,2,500,0,5);
		dProton_P_Theta = GetOrCreate_Histogram<TH2I>("Proton_P_Theta","p vs #theta; #theta; p (GeV/c)",180,0,180,120,0,12);
		dDeltaE_M2k = GetOrCreate_Histogram<TH2I>("DeltaE_M2k", "#DeltaE vs M_{k^{+}k^{-}}; M_{k^{+}k^{-}}; #DeltaE (tagger - tracks)", 800, 0.0, 2.0, 200, -5., 5.);
		dM2pi_M2k = GetOrCreate_Histogram<TH2I>("M2pi_M2k", "M_{#pi^{+}#pi^{-}} vs M_{k^{+}k^{-}}; M_{k^{+}k^{-}}; M_{#pi^{+}#pi^{-}}", 800, 0.0, 2.0, 800, 0.0, 2.0);
		
		dMM2_M2k_ProtonTag = GetOrCreate_Histogram<TH2I>("MM2_M2k_ProtonTag", "MM^{2} off k^{+}k^{-} vs M_{k^{+}k^{-}}; M_{k^{+}k^{-}}; MM^{2}", 800, 0.0, 2.0, 200, -1., 1.);
		dDeltaE_M2k_ProtonTag = GetOrCreate_Histogram<TH2I>("DeltaE_M2k_ProtonTag", "#DeltaE vs M_{k^{+}k^{-}}; M_{k^{+}k^{-}}; #DeltaE (tagger - tracks)", 800, 0.0, 2.0, 200, -5., 5.);
		dM2pi_M2k_ProtonTag = GetOrCreate_Histogram<TH2I>("M2pi_M2k_ProtonTag", "M_{#pi^{+}#pi^{-}} vs M_{k^{+}k^{-}}; M_{k^{+}k^{-}}; M_{#pi^{+}#pi^{-}}", 800, 0.0, 2.0, 800, 0.0, 2.0);
		
		dEgamma_M2k_ProtonTag = GetOrCreate_Histogram<TH2I>("Egamma_M2k", "TAGGER photon energy vs M_{k^{+}k^{-}}; M_{k^{+}k^{-}}; E_{#gamma}", 800, 0.0, 2.0, 400, 2., 12.);
		dKplus_deltaInvBeta_P_ProtonTag = GetOrCreate_Histogram<TH2I>("Kplus_deltaInvBeta_P_ProtonTag","K+ Inverse #Delta#beta vs p; p; 1/#beta_{exp} - 1/#beta_{meas}",200,0,10,500,-0.5,0.5);
		dKminus_deltaInvBeta_P_ProtonTag = GetOrCreate_Histogram<TH2I>("Kminus_deltaInvBeta_P_ProtonTag","K- Inverse #Delta#beta vs p; p; 1/#beta_{exp} - 1/#beta_{meas}",200,0,10,500,-0.5,0.5);
		
		dKplus_deltaInvBeta_P_PhiTag = GetOrCreate_Histogram<TH2I>("Kplus_deltaInvBeta_P_PhiTag","K+ Inverse #Delta#beta vs p; p; 1/#beta_{exp} - 1/#beta_{meas}",200,0,10,500,-0.5,0.5);
		dKminus_deltaInvBeta_P_PhiTag = GetOrCreate_Histogram<TH2I>("Kminus_deltaInvBeta_P_PhiTag","K- Inverse #Delta#beta vs p; p; 1/#beta_{exp} - 1/#beta_{meas}",200,0,10,500,-0.5,0.5);
		dKplus_Beta_P_PhiTag = GetOrCreate_Histogram<TH2I>("Kplus_Beta_P_PhiTag","K+ #beta vs p; p; #beta",200,0,10,500,-0.5,1.5);
		dKminus_Beta_P_PhiTag = GetOrCreate_Histogram<TH2I>("Kminus_Beta_P_PhiTag","K- #beta vs p; p; #beta",200,0,10,500,-0.5,1.5);
		
		dKplus_deltaInvBeta_P_RhoTag = GetOrCreate_Histogram<TH2I>("Kplus_deltaInvBeta_P_RhoTag","K+ Inverse #Delta#beta vs p; p; 1/#beta_{exp} - 1/#beta_{meas}",200,0,10,500,-0.5,0.5);
		dKminus_deltaInvBeta_P_RhoTag = GetOrCreate_Histogram<TH2I>("Kminus_deltaInvBeta_P_RhoTag","K- Inverse #Delta#beta vs p; p; 1/#beta_{exp} - 1/#beta_{meas}",200,0,10,500,-0.5,0.5);
		dKplus_Beta_P_RhoTag = GetOrCreate_Histogram<TH2I>("Kplus_Beta_P_RhoTag","K+ #beta vs p; p; #beta",200,0,10,500,-0.5,1.5);
		dKminus_Beta_P_RhoTag = GetOrCreate_Histogram<TH2I>("Kminus_Beta_P_RhoTag","K- #beta vs p; p; #beta",200,0,10,500,-0.5,1.5);
		
		dKplus_P_Theta_PhiTag = GetOrCreate_Histogram<TH2I>("Kplus_P_Theta_PhiTag","K+ p vs #theta; #theta; p (GeV)",180,0,180,200,0.,10.);
		dKminus_P_Theta_PhiTag = GetOrCreate_Histogram<TH2I>("Kminus_P_Theta_PhiTag","K- p vs #theta; #theta; p (GeV)",180,0,180,200,0.,10.);
		dKplus_P_Theta_RhoTag = GetOrCreate_Histogram<TH2I>("Kplus_P_Theta_RhoTag","K+ p vs #theta; #theta; p (GeV)",180,0,180,200,0.,10.);
		dKminus_P_Theta_RhoTag = GetOrCreate_Histogram<TH2I>("Kminus_P_Theta_RhoTag","K- p vs #theta; #theta; p (GeV)",180,0,180,200,0.,10.);

	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DCustomAction_p2k_hists::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
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
        double locBeamPhotonEnergy = locBeamPhoton->energy();

	// get k+ and k- tracks
	double kplus_beta = -1.;
	double kminus_beta = -1.;
	double kplus_deltaInvBeta = -1.;
	double kminus_deltaInvBeta = -1.;
	DLorentzVector kplus_PionHypothesis, kminus_PionHypothesis;
	for(size_t loc_i = 0; loc_i < 2; ++loc_i) {
		const DChargedTrack* locChargedTrack = static_cast<const DChargedTrack*>(locParticleComboStep->Get_FinalParticle_SourceObject(loc_i));
		const DChargedTrackHypothesis* locChargedTrackHypothesis;
		if(loc_i == 0) {
			locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(KPlus);
			if(locChargedTrackHypothesis == NULL) { cout<<"Missing K+"<<endl; continue; }
			kplus_beta = locChargedTrackHypothesis->measuredBeta();
			kplus_deltaInvBeta = locChargedTrackHypothesis->deltaInvBeta();

			const DChargedTrackHypothesis* locChargedTrackPion = locChargedTrack->Get_Hypothesis(PiPlus);
			if(locChargedTrackPion != NULL)
				kplus_PionHypothesis = locChargedTrackPion->lorentzMomentum();
			//else 
				
			kplus_PionHypothesis.SetXYZM(locParticles[0]->lorentzMomentum().X(), locParticles[0]->lorentzMomentum().Y(), locParticles[0]->lorentzMomentum().Z(), 0.13957);
		}
		else if(loc_i == 1){
			locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(KMinus);
			if(locChargedTrackHypothesis == NULL) { cout<<"Missing K-"<<endl; continue; }
			kminus_beta = locChargedTrackHypothesis->measuredBeta();
			kminus_deltaInvBeta = locChargedTrackHypothesis->deltaInvBeta();

			const DChargedTrackHypothesis* locChargedTrackPion = locChargedTrack->Get_Hypothesis(PiMinus);
			if(locChargedTrackPion != NULL)
				kminus_PionHypothesis = locChargedTrackPion->lorentzMomentum();
			//else 
				
			kminus_PionHypothesis.SetXYZM(locParticles[1]->lorentzMomentum().X(), locParticles[1]->lorentzMomentum().Y(), locParticles[1]->lorentzMomentum().Z(), 0.13957);
		}
	}

	// calculate 2k P4
	DLorentzVector locP4_2k = locParticles[0]->lorentzMomentum() + locParticles[1]->lorentzMomentum();
	DLorentzVector locP4_2pi = kplus_PionHypothesis + kminus_PionHypothesis;

	// calculate missing P4
	DLorentzVector locMissingP4 = dAnalysisUtilities->Calc_MissingP4(locParticleCombo,Get_UseKinFitResultsFlag());

	// reconstructed proton variables
	double dEdx = 0.;
	DLorentzVector locProtonP4;
	const DChargedTrack* locChargedTrack = static_cast<const DChargedTrack*>(locParticleComboStep->Get_FinalParticle_SourceObject(2));
	const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(Proton);
	dEdx = locChargedTrackHypothesis->dEdx()*1e6;
	locProtonP4 = locChargedTrackHypothesis->lorentzMomentum();	

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action(); //ACQUIRE ROOT LOCK!!
	{
		// Fill histograms here
		dEgamma->Fill(locBeamPhotonEnergy);

		dKplus_deltaInvBeta_P->Fill(locParticles[0]->lorentzMomentum().Vect().Mag(), kplus_deltaInvBeta);
		dKminus_deltaInvBeta_P->Fill(locParticles[1]->lorentzMomentum().Vect().Mag(), kminus_deltaInvBeta);

		dMM2_M2k->Fill(locP4_2k.M(), locMissingP4.M2());
		dDeltaE_M2k->Fill(locP4_2k.M(),locMissingP4.E());
		dM2pi_M2k->Fill(locP4_2k.M(), locP4_2pi.M());

		dMissingMassSq->Fill(locMissingP4.M2());
		
		if(locMissingP4.M2() > minMM2Cut && locMissingP4.M2() < maxMM2Cut) {
		
			dProton_dEdx_P->Fill(locProtonP4.Vect().Mag(), dEdx);
			dProton_P_Theta->Fill(locProtonP4.Vect().Theta()*180/TMath::Pi(), locProtonP4.Vect().Mag());
			
			// plots for proton tagged events
			if(dEdx > dEdxCut) {
				
				dMM2_M2k_ProtonTag->Fill(locP4_2k.M(), locMissingP4.M2());
				dDeltaE_M2k_ProtonTag->Fill(locP4_2k.M(),locMissingP4.E());
				dM2pi_M2k_ProtonTag->Fill(locP4_2k.M(), locP4_2pi.M());
				
				dKplus_deltaInvBeta_P_ProtonTag->Fill(locParticles[0]->lorentzMomentum().Vect().Mag(), kplus_deltaInvBeta);
				dKminus_deltaInvBeta_P_ProtonTag->Fill(locParticles[1]->lorentzMomentum().Vect().Mag(), kminus_deltaInvBeta);
				
				// plots vs phi mass
				dEgamma_M2k_ProtonTag->Fill(locP4_2k.M(),locBeamPhotonEnergy);

				dInvariantMass->Fill(locP4_2k.M());
				
				// correlations for rho vs phi mass regions
				if(locP4_2k.M() < maxPhiMassCut){
					dKplus_deltaInvBeta_P_PhiTag->Fill(locParticles[0]->lorentzMomentum().Vect().Mag(), kplus_deltaInvBeta);
					dKminus_deltaInvBeta_P_PhiTag->Fill(locParticles[1]->lorentzMomentum().Vect().Mag(), kminus_deltaInvBeta);
					
					dKplus_P_Theta_PhiTag->Fill(locParticles[0]->lorentzMomentum().Theta()*180./TMath::Pi(),locParticles[0]->lorentzMomentum().Vect().Mag());
					dKminus_P_Theta_PhiTag->Fill(locParticles[1]->lorentzMomentum().Theta()*180./TMath::Pi(),locParticles[1]->lorentzMomentum().Vect().Mag());
					
					dKplus_Beta_P_PhiTag->Fill(locParticles[0]->lorentzMomentum().Vect().Mag(), kplus_beta);
					dKminus_Beta_P_PhiTag->Fill(locParticles[1]->lorentzMomentum().Vect().Mag(), kminus_beta);
				}
				else {
					dKplus_deltaInvBeta_P_RhoTag->Fill(locParticles[0]->lorentzMomentum().Vect().Mag(), kplus_deltaInvBeta);
					dKminus_deltaInvBeta_P_RhoTag->Fill(locParticles[1]->lorentzMomentum().Vect().Mag(), kminus_deltaInvBeta);
					
					dKplus_P_Theta_RhoTag->Fill(locParticles[0]->lorentzMomentum().Theta()*180./TMath::Pi(),locParticles[0]->lorentzMomentum().Vect().Mag());
					dKminus_P_Theta_RhoTag->Fill(locParticles[1]->lorentzMomentum().Theta()*180./TMath::Pi(),locParticles[1]->lorentzMomentum().Vect().Mag());
					
					dKplus_Beta_P_RhoTag->Fill(locParticles[0]->lorentzMomentum().Vect().Mag(), kplus_beta);
					dKminus_Beta_P_RhoTag->Fill(locParticles[1]->lorentzMomentum().Vect().Mag(), kminus_beta);
				}
			}
		}
	}
	Unlock_Action(); //RELEASE ROOT LOCK!!

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}

