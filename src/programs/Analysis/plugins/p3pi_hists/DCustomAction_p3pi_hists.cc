// $Id$
//
//    File: DCustomAction_p3pi_hists.cc
// Created: Wed Jan 21 16:53:41 EST 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "DCustomAction_p3pi_hists.h"

void DCustomAction_p3pi_hists::Initialize(JEventLoop* locEventLoop)
{

	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//Required: Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this action.
			//If another thread has already created the folder, it just changes to it. 
		CreateAndChangeTo_ActionDirectory();

		dEgamma = GetOrCreate_Histogram<TH1I>("Egamma", "TAGGER photon energy; E_{#gamma}", 400, 0., 12.);
		
		dMM2_M3pi = GetOrCreate_Histogram<TH2I>("MM2_M3pi", "MM^{2} off #pi^{+}#pi^{-}#pi^{0} vs M_{#pi^{+}#pi^{-}#pi^{0}}; M_{#pi^{+}#pi^{-}#pi^{0}}; MM^{2}", 200, 0.0, 2.0, 200, -1., 1.);
		dProton_dEdx_P = GetOrCreate_Histogram<TH2I>("Proton_dEdx_P","dE/dx vs p; p; dE/dx",200,0,2,500,0,5);
		dProton_P_Theta = GetOrCreate_Histogram<TH2I>("Proton_P_Theta","p vs #theta; #theta; p (GeV/c)",180,0,180,120,0,12);
		
		dDeltaE_M3pi = GetOrCreate_Histogram<TH2I>("dDeltaE_M3pi", "#Delta E vs M_{#pi^{+}#pi^{-}#pi^{0}}; M_{#pi^{+}#pi^{-}#pi^{0}}; #Delta E (tagger - tracks)", 200, 0.0, 2.0, 200, -5., 5.);
		dPhi3pi_PhiP = GetOrCreate_Histogram<TH2I>("Phi3pi_PhiP", "#phi_{#pi^{+}#pi^{-}#pi^{0}} vs. #phi_{p}; #phi_{p}; #phi_{#pi^{+}#pi^{-}#pi^{0}}", 360, -180.0, 180.0, 360, -180., 180.);
		dDeltaPhi_M3pi = GetOrCreate_Histogram<TH2I>("DeltaPhi_M3pi", "#Delta#phi vs M_{#pi^{+}#pi^{-}#pi^{0}}; M_{#pi^{+}#pi^{-}#pi^{0}}; #Delta#phi", 200, 0.0, 2.0, 360, 0.0, 360.0);
		
		dEgamma_M3pi_ProtonTag = GetOrCreate_Histogram<TH2I>("dEgamma_M3pi_ProtonTag", "E_{#gamma} vs M_{#pi^{+}#pi^{-}#pi^{0}}; M_{#pi^{+}#pi^{-}#pi^{0}}; E_{#gamma}", 200, 0.0, 2.0, 240, 0., 6.);
		dMppizero_M2pi = GetOrCreate_Histogram<TH2I>("Mppizero_M2pi", "Psuedo-Dalitz p#pi^{+}#pi^{-}#pi^{0}; M_{#pi^{+}#pi^{-}}; M_{p#pi^{0}}", 200, 0.0, 1.0, 200, 1.0, 4.0);
		dMppiplus_M2pi = GetOrCreate_Histogram<TH2I>("Mppiplus_M2pi", "Psuedo-Dalitz p#pi^{+}#pi^{-}#pi^{0}; M_{#pi^{0}#pi^{-}}; M_{p#pi^{+}}", 200, 0.0, 1.0, 200, 1.0, 4.0);
		dMppiminus_M2pi = GetOrCreate_Histogram<TH2I>("Mppiminus_M2pi", "Psuedo-Dalitz p#pi^{+}#pi^{-}#pi^{0}; M_{#pi^{0}#pi^{+}}; M_{p#pi^{-}}", 200, 0.0, 1.0, 200, 1.0, 4.0);
		
		dMM2_M3pi_CoplanarTag = GetOrCreate_Histogram<TH2I>("MM2_M3pi_CoplanarTag", "MM^{2} off p #pi^{+}#pi^{-}#pi^{0} vs M_{#pi^{+}#pi^{-}#pi^{0}}: Coplanar Tag; M_{#pi^{+}#pi^{-}#pi^{0}}; MM^{2}", 200, 0.0, 2.0, 200, -1., 1.);
		dDeltaE_M3pi_CoplanarTag = GetOrCreate_Histogram<TH2I>("dDeltaE_M3pi_CoplanarTag", "#Delta E vs M_{#pi^{+}#pi^{-}#pi^{0}}: Coplanar Tag; M_{#pi^{+}#pi^{-}#pi^{0}}; #Delta E (tagger - p#pi^{+}#pi^{-}#pi^{0})", 200, 0.0, 2.0, 200, -5., 5.);
		dMM2_DeltaE_CoplanarTag = GetOrCreate_Histogram<TH2I>("dMM2_DeltaE_CoplanarTag", "MM^{2} vs #Delta E: Coplanar Tag; #Delta E (tagger - p#pi^{+}#pi^{-}#pi^{0}); MM^{2}", 200, -5., 5., 200, -1., 1.);
		
		dMM2_M3pi_ProtonTag = GetOrCreate_Histogram<TH2I>("MM2_M3pi_ProtonTag", "MM^{2} off p #pi^{+}#pi^{-}#pi^{0} vs M_{#pi^{+}#pi^{-}#pi^{0}}: Proton Tag; M_{#pi^{+}#pi^{-}#pi^{0}}; MM^{2}", 200, 0.0, 2.0, 200, -1., 1.);
		dDeltaE_M3pi_ProtonTag = GetOrCreate_Histogram<TH2I>("dDeltaE_M3pi_ProtonTag", "#Delta E vs M_{#pi^{+}#pi^{-}#pi^{0}}: Proton Tag; M_{#pi^{+}#pi^{-}#pi^{0}}; #Delta E (tagger - p#pi^{+}#pi^{-}#pi^{0})", 200, 0.0, 2.0, 200, -5., 5.);
		dMM2_DeltaE_ProtonTag = GetOrCreate_Histogram<TH2I>("dMM2_DeltaE_ProtonTag", "MM^{2} vs #Delta E: Proton Tag; #Delta E (tagger - p#pi^{+}#pi^{-}#pi^{0}); MM^{2}", 200, -5., 5., 200, -1., 1.);
		
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DCustomAction_p3pi_hists::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
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

	DLorentzVector locSumInitP4;
	locSumInitP4.SetXYZM(0, 0, 0, 0.938);
        locSumInitP4 += locBeamPhoton->lorentzMomentum();

	// cut on tagger energy
	if(locBeamPhotonEnergy < 1.5) 
		return true;

	// calculate missing mass
	DLorentzVector locMissingP4 = dAnalysisUtilities->Calc_MissingP4(locParticleCombo, Get_UseKinFitResultsFlag());
	
	// calculate 3pi mass
	DLorentzVector locOmegaP4 = dAnalysisUtilities->Calc_FinalStateP4(locParticleCombo, 1, Get_UseKinFitResultsFlag());
	DLorentzVector locPi0P4 = dAnalysisUtilities->Calc_FinalStateP4(locParticleCombo, 2, Get_UseKinFitResultsFlag());
	DLorentzVector locPiPlusP4 = locParticleCombo->Get_ParticleComboStep(1)->Get_FinalParticle_Measured(0)->lorentzMomentum();
	DLorentzVector locPiMinusP4 = locParticleCombo->Get_ParticleComboStep(1)->Get_FinalParticle_Measured(1)->lorentzMomentum();

	// reconstructed proton variables
	double dEdx = 0.;
	DLorentzVector locProtonP4;
	const DChargedTrack* locChargedTrack = static_cast<const DChargedTrack*>(locParticleComboStep->Get_FinalParticle_SourceObject(1));
	const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(Proton);
	dEdx = locChargedTrackHypothesis->dEdx()*1e6;
	locProtonP4 = locChargedTrackHypothesis->lorentzMomentum();	

	double dEdxCut = 2.2;

	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		// Fill histograms here
		dEgamma->Fill(locBeamPhotonEnergy);

		dMM2_M3pi->Fill(locOmegaP4.M(), locMissingP4.M2());
		dDeltaE_M3pi->Fill(locOmegaP4.M(),locMissingP4.E());

		double locDeltaPhi = (locProtonP4.Phi() - locOmegaP4.Phi())*180./TMath::Pi();
		if(locDeltaPhi > 360.) locDeltaPhi -= 360.;
		if(locDeltaPhi < 0.) locDeltaPhi += 360.;
		dPhi3pi_PhiP->Fill(locProtonP4.Phi()*180./TMath::Pi(), locOmegaP4.Phi()*180./TMath::Pi());
		dDeltaPhi_M3pi->Fill(locOmegaP4.M(), locDeltaPhi);
		
		// require proton and omega are back-to-back
		if(locDeltaPhi < 175. || locDeltaPhi > 185.) return true;
		dMM2_M3pi_CoplanarTag->Fill(locOmegaP4.M(), locMissingP4.M2());
		dDeltaE_M3pi_CoplanarTag->Fill(locOmegaP4.M(),locMissingP4.E());
		if(locOmegaP4.M() > 0.7 && locOmegaP4.M() < 0.9)
			dMM2_DeltaE_CoplanarTag->Fill(locMissingP4.E(), locMissingP4.M2());
		
		// tag proton with dE/dx
		if(dEdx > dEdxCut) {
			dMM2_M3pi_ProtonTag->Fill(locOmegaP4.M(), locMissingP4.M2());
			dDeltaE_M3pi_ProtonTag->Fill(locOmegaP4.M(),locMissingP4.E());
			if(locOmegaP4.M() > 0.7 && locOmegaP4.M() < 0.9)
				dMM2_DeltaE_ProtonTag->Fill(locMissingP4.E(), locMissingP4.M2());
			
			// 3pi invariant mass for exclusive
			if(fabs(locMissingP4.M2()) < 0.05 && fabs(locMissingP4.E()) < 0.5) {
				dEgamma_M3pi_ProtonTag->Fill(locOmegaP4.M(),locBeamPhotonEnergy);
				
				if(locBeamPhotonEnergy > 2.5 && locBeamPhotonEnergy < 3.0){
					dMppizero_M2pi->Fill((locPiPlusP4+locPiMinusP4).M(), (locProtonP4+locPi0P4).M());
					dMppiplus_M2pi->Fill((locPi0P4+locPiMinusP4).M(), (locProtonP4+locPiPlusP4).M());
					dMppiminus_M2pi->Fill((locPiPlusP4+locPi0P4).M(), (locProtonP4+locPiMinusP4).M());
				}
			}
		}
		
		if(fabs(locMissingP4.M2()) < 0.05 && fabs(locMissingP4.E()) < 0.5) {
			
			dProton_dEdx_P->Fill(locProtonP4.Vect().Mag(), dEdx);
			dProton_P_Theta->Fill(locProtonP4.Vect().Theta()*180/TMath::Pi(), locProtonP4.Vect().Mag());
			
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}
