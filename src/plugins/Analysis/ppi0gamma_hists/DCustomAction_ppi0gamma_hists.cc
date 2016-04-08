// $Id$
//
//    File: DCustomAction_ppi0gamma_hists.cc
// Created: Wed Jan 21 16:53:41 EST 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "DCustomAction_ppi0gamma_hists.h"

void DCustomAction_ppi0gamma_hists::Initialize(JEventLoop* locEventLoop)
{
	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//Required: Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this action.
			//If another thread has already created the folder, it just changes to it. 
		CreateAndChangeTo_ActionDirectory();

		dEgamma = GetOrCreate_Histogram<TH1I>("Egamma", "TAGGER photon energy; E_{#gamma}", 400, 0., 12.);
		
		
		dMM2_MPi0 = GetOrCreate_Histogram<TH2I>("MM2_MPi0", "MM^{2} off p #pi^{0}#gamma vs M_{#gamma#gamma}; M_{#gamma#gamma}; MM^{2}", 400, 0.0, 2.0, 200, -1., 1.);
		dMM2_MOmega = GetOrCreate_Histogram<TH2I>("MM2_MOmega", "MM^{2} off p #pi^{0}#gamma vs M_{#pi^{0}#gamma}; M_{#pi^{0}#gamma}; MM^{2}", 500, 0.0, 1.0, 200, -1., 1.);
		dProton_dEdx_P = GetOrCreate_Histogram<TH2I>("Proton_dEdx_P","dE/dx vs p; p; dE/dx",200,0,2,500,0,5);
		dProton_P_Theta = GetOrCreate_Histogram<TH2I>("Proton_P_Theta","p vs #theta; #theta; p (GeV/c)",180,0,180,120,0,12);

		dDeltaE_MOmega = GetOrCreate_Histogram<TH2I>("dDeltaE_MOmega", "#Delta E vs M_{#pi^{0}#gamma}; M_{#pi^{0}#gamma}; #Delta E (tagger - p#pi^{0}#gamma)", 400, 0.0, 2.0, 200, -5., 5.);
		dPhiOmega_PhiP = GetOrCreate_Histogram<TH2I>("PhiOmega_PhiP", "#phi_{#pi^{0}#gamma} vs. #phi_{p}; #phi_{p}; #phi_{#pi^{0}#gamma}", 360, -180.0, 180.0, 360, -180., 180.);
		dDeltaPhi_MOmega = GetOrCreate_Histogram<TH2I>("DeltaPhi_MOmega", "#Delta#phi vs M_{#pi^{0}#gamma}; M_{#pi^{0}#gamma}; #Delta#phi", 400, 0.0, 2.0, 360, 0.0, 360.0);

		dEgamma_MOmegaProtonTag = GetOrCreate_Histogram<TH2I>("Egamma_MOmegaProtonTag", "E_{#gamma} vs M_{#pi^{0}#gamma}; M_{#pi^{0}#gamma}; E_{#gamma}", 400, 0.0, 2.0, 240, 0., 6.);
		dEgamma_MOmegaKinTag = GetOrCreate_Histogram<TH2I>("Egamma_MOmegaKinTag", "E_{#gamma} vs M_{#pi^{0}#gamma}; M_{#pi^{0}#gamma}; E_{#gamma}", 400, 0.0, 2.0, 240, 0., 6.);

		dMM2_MOmegaCoplanarTag = GetOrCreate_Histogram<TH2I>("MM2_MOmegaCoplanarTag", "MM^{2} off p #pi^{0}#gamma vs M_{#pi^{0}#gamma}: Coplanar Tag; M_{#gamma#gamma}; MM^{2}", 400, 0.0, 2.0, 200, -1., 1.);
                dDeltaE_MOmegaCoplanarTag = GetOrCreate_Histogram<TH2I>("dDeltaE_MOmegaCoplanarTag", "#Delta E vs M_{#pi^{0}#gamma}: Coplanar Tag; M_{#pi^{0}#gamma}; #Delta E (tagger - p#pi^{0}#gamma)", 400, 0.0, 2.0, 200, -5., 5.);
		dMM2_DeltaE_CoplanarTag = GetOrCreate_Histogram<TH2I>("dMM2_DeltaE_CoplanarTag", "MM^{2} vs #Delta E: Coplanar Tag; #Delta E (tagger - p#pi^{0}#gamma); MM^{2}", 200, -5., 5., 200, -1., 1.);

		dMM2_MOmegaProtonTag = GetOrCreate_Histogram<TH2I>("MM2_MOmegaProtonTag", "MM^{2} off p #pi^{0}#gamma vs M_{#pi^{0}#gamma}: Proton Tag; M_{#pi^{0}#gamma}; MM^{2}", 400, 0.0, 2.0, 200, -1., 1.);
		dDeltaE_MOmegaProtonTag = GetOrCreate_Histogram<TH2I>("dDeltaE_MOmegaProtonTag", "#Delta E vs M_{#pi^{0}#gamma}: Proton Tag; M_{#pi^{0}#gamma}; #Delta E (tagger - p#gamma#gamma)", 400, 0.0, 2.0, 200, -5., 5.);
		dMM2_DeltaE_ProtonTag = GetOrCreate_Histogram<TH2I>("dMM2_DeltaE_ProtonTag", "MM^{2} vs #Delta E: Proton Tag; #Delta E (tagger - p#pi^{0}#gamma); MM^{2}", 200, -5., 5., 200, -1., 1.);
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DCustomAction_ppi0gamma_hists::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
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
	
	// calculate missing mass
        DLorentzVector locMissingP4 = dAnalysisUtilities->Calc_MissingP4(locParticleCombo, Get_UseKinFitResultsFlag());

        // calculate masses
        DLorentzVector locOmegaP4 = dAnalysisUtilities->Calc_FinalStateP4(locParticleCombo, 1, Get_UseKinFitResultsFlag());
	DLorentzVector locPi0P4 = dAnalysisUtilities->Calc_FinalStateP4(locParticleCombo, 2, Get_UseKinFitResultsFlag());

	// reconstructed proton variables
	const DChargedTrack* locChargedTrack = static_cast<const DChargedTrack*>(locParticleComboStep->Get_FinalParticle_SourceObject(1));
	const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(Proton);
	double dEdx = locChargedTrackHypothesis->dEdx()*1e6;
	DLorentzVector locProtonP4 = locChargedTrackHypothesis->lorentzMomentum();	

	double dEdxCut = 2.2;
	
	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action(); //ACQUIRE ROOT LOCK!!
	{
		// Fill histograms here
		dEgamma->Fill(locBeamPhotonEnergy);
		dMM2_MPi0->Fill(locPi0P4.M(), locMissingP4.M2());
		
		// cut on pi0 mass
		if(locPi0P4.M() < 0.12 || locPi0P4.M() > 0.15) return true;
		
		dMM2_MOmega->Fill(locOmegaP4.M(), locMissingP4.M2());
		dDeltaE_MOmega->Fill(locOmegaP4.M(),locMissingP4.E());

		double locDeltaPhi = (locProtonP4.Phi() - locOmegaP4.Phi())*180./TMath::Pi();
		if(locDeltaPhi > 360.) locDeltaPhi -= 360.;
		if(locDeltaPhi < 0.) locDeltaPhi += 360.;
		dPhiOmega_PhiP->Fill(locProtonP4.Phi()*180./TMath::Pi(), locOmegaP4.Phi()*180./TMath::Pi());
		dDeltaPhi_MOmega->Fill(locOmegaP4.M(), locDeltaPhi);

		// require proton and omega are back-to-back
		if(locDeltaPhi < 170. || locDeltaPhi > 190.) return true;
		dMM2_MOmegaCoplanarTag->Fill(locOmegaP4.M(), locMissingP4.M2());
		dDeltaE_MOmegaCoplanarTag->Fill(locOmegaP4.M(), locMissingP4.E());
		if(locOmegaP4.M() > 0.7 && locOmegaP4.M() < 0.9)
			dMM2_DeltaE_CoplanarTag->Fill(locMissingP4.E(), locMissingP4.M2());

		// tag proton with dE/dx
		if(dEdx > dEdxCut) {
			dMM2_MOmegaProtonTag->Fill(locOmegaP4.M(), locMissingP4.M2());
			dDeltaE_MOmegaProtonTag->Fill(locOmegaP4.M(), locMissingP4.E());
			if(locOmegaP4.M() > 0.7 && locOmegaP4.M() < 0.9)
				dMM2_DeltaE_ProtonTag->Fill(locMissingP4.E(), locMissingP4.M2());

			// invariant mass for exclusive g+p -> p + omega
			if(fabs(locMissingP4.M2()) < 0.01 && fabs(locMissingP4.E()) < 0.1) {
				dEgamma_MOmegaProtonTag->Fill(locOmegaP4.M(),locBeamPhotonEnergy);
			}
		}

		// for pi0gamma candidates require recoil proton
		if(fabs(locMissingP4.M2()) < 0.01 && fabs(locMissingP4.E()) < 0.1) {
			dProton_dEdx_P->Fill(locProtonP4.Vect().Mag(), dEdx);
			dProton_P_Theta->Fill(locProtonP4.Vect().Theta()*180/TMath::Pi(), locProtonP4.Vect().Mag());			
			dEgamma_MOmegaKinTag->Fill(locOmegaP4.M(),locBeamPhotonEnergy);

			if(locBeamPhotonEnergy > 2.5 && locBeamPhotonEnergy < 3.0){ // coherent peak
			
			}
			
			Unlock_Action(); //RELEASE ROOT LOCK!!
			return true;
		}
	}
	Unlock_Action(); //RELEASE ROOT LOCK!!

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}
