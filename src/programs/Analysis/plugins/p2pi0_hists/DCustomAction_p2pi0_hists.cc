// $Id$
//
//    File: DCustomAction_p2pi0_hists.cc
// Created: Wed Jan 21 16:53:41 EST 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "DCustomAction_p2pi0_hists.h"

void DCustomAction_p2pi0_hists::Initialize(JEventLoop* locEventLoop)
{
	
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//Required: Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this action.
			//If another thread has already created the folder, it just changes to it. 
		CreateAndChangeTo_ActionDirectory();

		dEgamma = GetOrCreate_Histogram<TH1I>("Egamma", "TAGGER photon energy; E_{#gamma}", 400, 0., 12.);
		
		dMM2_M2pi0 = GetOrCreate_Histogram<TH2I>("MM2_M2pi0", "MM^{2} off p#pi^{0}#pi^{0} vs M_{#pi^{0}#pi^{0}}; M_{#pi^{0}#pi^{0}}; MM^{2}", 200, 0.0, 2.0, 200, -1., 1.);
		dDeltaE_M2pi0 = GetOrCreate_Histogram<TH2I>("dDeltaE_M2pi0", "#Delta E vs M_{#pi^{0}#pi^{0}}; M_{#pi^{0}}; #Delta E (tagger - tracks)", 200, 0.0, 2.0, 200, -5., 5.);
		dMpi0_corr = GetOrCreate_Histogram<TH2I>("Mpi0_corr", "M_{#gamma#gamma}(1) vs M_{#gamma#gamma}(2); M_{#gamma#gamma}(1); M_{#gamma#gamma}(2)", 200, 0.0, 1.0, 200, 0., 1.);
	
		dPhi2pi0_PhiP = GetOrCreate_Histogram<TH2I>("Phi2pi0_PhiP", "#phi_{#pi^{0}#pi^{0}} vs. #phi_{p}; #phi_{p}; #phi_{#pi^{0}#pi^{0}}", 360, -180.0, 180.0, 360, -180., 180.);
		dDeltaPhi_M2pi0 = GetOrCreate_Histogram<TH2I>("DeltaPhi_M2pi0", "#Delta#phi vs M_{#pi^{0}#pi^{0}}; M_{#pi^{0}#pi^{0}}; #Delta#phi", 200, 0.0, 2.0, 360, 0.0, 360.0);
		
		dMM2_M2pi0_CoplanarTag = GetOrCreate_Histogram<TH2I>("MM2_M2pi0_CoplanarTag", "MM^{2} off p #pi^{0}#pi^{0} vs M_{#pi^{0}#pi^{0}}: Coplanar Tag; M_{#pi^{0}#pi^{0}}; MM^{2}", 200, 0.0, 2.0, 200, -1., 1.);
		dDeltaE_M2pi0_CoplanarTag = GetOrCreate_Histogram<TH2I>("dDeltaE_M2pi0_CoplanarTag", "#Delta E vs M_{#pi^{0}#pi^{0}}: Coplanar Tag; M_{#pi^{0}#pi^{0}}; #Delta E (tagger - p#pi^{0}#pi^{0})", 200, 0.0, 2.0, 200, -5., 5.);
		dMpi0_corr_CoplanarTag = GetOrCreate_Histogram<TH2I>("Mpi0_Corr_CoplanarTag", "M_{#gamma#gamma}(1) vs M_{#gamma#gamma}(2); M_{#gamma#gamma}(1); M_{#gamma#gamma}(2)", 200, 0.0, 1.0, 200, 0., 1.);

		dMM2_M2pi0_ProtonTag = GetOrCreate_Histogram<TH2I>("MM2_M2pi0_ProtonTag", "MM^{2} off p #pi^{0}#pi^{0} vs M_{#pi^{0}#pi^{0}}: Proton Tag; M_{#pi^{0}#pi^{0}}; MM^{2}", 200, 0.0, 2.0, 200, -1., 1.);
		dDeltaE_M2pi0_ProtonTag = GetOrCreate_Histogram<TH2I>("dDeltaE_M2pi0_ProtonTag", "#Delta E vs M_{#pi^{0}#pi^{0}}: Proton Tag; M_{#pi^{0}#pi^{0}}; #Delta E (tagger - p#pi^{0}#pi^{0})", 200, 0.0, 2.0, 200, -5., 5.);
		dMpi0_corr_ProtonTag = GetOrCreate_Histogram<TH2I>("Mpi0_corr_ProtonTag", "M_{#gamma#gamma}(1) vs M_{#gamma#gamma}(2); M_{#gamma#gamma}(1); M_{#gamma#gamma}(2)", 200, 0.0, 1.0, 200, 0., 1.);

		dMpi0_corr_MMTag = GetOrCreate_Histogram<TH2I>("Mpi0_corr_MMTag", "M_{#gamma#gamma}(1) vs M_{#gamma#gamma}(2); M_{#gamma#gamma}(1); M_{#gamma#gamma}(2)", 200, 0.0, 1.0, 200, 0., 1.);
		
		dMM2_M2pi0_Pi0Tag = GetOrCreate_Histogram<TH2I>("MM2_M2pi0_Pi0Tag", "MM^{2} off p #pi^{0}#pi^{0} vs M_{#pi^{0}#pi^{0}}: Pi0 Tag; M_{#pi^{0}#pi^{0}}; MM^{2}", 200, 0.0, 2.0, 200, -1., 1.);
		dDeltaE_M2pi0_Pi0Tag = GetOrCreate_Histogram<TH2I>("dDeltaE_M2pi0_Pi0Tag", "#Delta E vs M_{#pi^{0}#pi^{0}}: Pi0 Tag; M_{#pi^{0}#pi^{0}}; #Delta E (tagger - p#pi^{0}#pi^{0})", 200, 0.0, 2.0, 200, -5., 5.);
		dMM2_DeltaE_Pi0Tag = GetOrCreate_Histogram<TH2I>("dMM2_DeltaE_Pi0Tag", "MM^{2} vs #Delta E: Coplanar Tag; #Delta E (tagger - p#gamma#gamma); MM^{2}", 200, -5., 5., 200, -1., 1.);

		dEgamma_M2pi0_Pi0Tag = GetOrCreate_Histogram<TH2I>("Egamma_M2pi0_Pi0Tag", "E_{#gamma} vs M_{#pi^{0}#pi^{0}}: Pi0 Tag; M_{#pi^{0}#pi^{0}}; E_{#gamma}", 200, 0.0, 2.0, 400, 0., 12.);		
		dDalitz1_Pi0Tag = GetOrCreate_Histogram<TH2I>("Dalitz1_Pi0Tag", "Dalitz Pi0 Tag; M^{2}_{p#pi^{0}}; M^{2}_{p#pi^{0}}", 200, 1.0, 11.0, 200, 1., 11.);
		dDalitz2_Pi0Tag = GetOrCreate_Histogram<TH2I>("Dalitz2_Pi0Tag", "Dalitz Pi0 Tag; M^{2}_{p#pi^{0}}; M^{2}_{p#pi^{0}}", 200, 1.0, 11.0, 200, 1., 11.);
		dDalitz3_Pi0Tag = GetOrCreate_Histogram<TH2I>("Dalitz3_Pi0Tag", "Dalitz Pi0 Tag; M^{2}_{p#pi^{0}}; M^{2}_{p#pi^{0}}", 200, 1.0, 11.0, 200, 1., 11.);
		dDalitz4_Pi0Tag = GetOrCreate_Histogram<TH2I>("Dalitz4_Pi0Tag", "Dalitz Pi0 Tag; M^{2}_{p#pi^{0}}; M^{2}_{p#pi^{0}}", 200, 1.0, 11.0, 200, 1., 11.);
		

		dProton_dEdx_P = GetOrCreate_Histogram<TH2I>("Proton_dEdx_P","dE/dx vs p; p; dE/dx",200,0,2,500,0,5);
		dProton_P_Theta = GetOrCreate_Histogram<TH2I>("Proton_P_Theta","p vs #theta; #theta; p (GeV/c)",180,0,180,120,0,12);
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DCustomAction_p2pi0_hists::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
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

	// cut on tagger energy
	if(locBeamPhotonEnergy < 1.5) 
		return true;

	// calculate missing mass
	DLorentzVector locMissingP4 = dAnalysisUtilities->Calc_MissingP4(locParticleCombo, Get_UseKinFitResultsFlag());
	
	// calculate Pi0 masses
	DLorentzVector locPi01_P4 = dAnalysisUtilities->Calc_FinalStateP4(locParticleCombo, 1, Get_UseKinFitResultsFlag());
	DLorentzVector locPi02_P4 = dAnalysisUtilities->Calc_FinalStateP4(locParticleCombo, 2, Get_UseKinFitResultsFlag());
	DLorentzVector loc2pi0_P4 = locPi01_P4; loc2pi0_P4 += locPi02_P4;

	// reconstructed proton variables
	double dEdx = 0.;
	DLorentzVector locProtonP4;
	const DChargedTrack* locChargedTrack = static_cast<const DChargedTrack*>(locParticleComboStep->Get_FinalParticle_SourceObject(2));
	const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(Proton);
	dEdx = locChargedTrackHypothesis->dEdx()*1e6;
	locProtonP4 = locChargedTrackHypothesis->lorentzMomentum();	

	DLorentzVector locPi01P_P4 = locProtonP4; locPi01P_P4 += locPi01_P4;
	DLorentzVector locPi02P_P4 = locProtonP4; locPi02P_P4 += locPi02_P4;

	double dEdxCut = 2.2;

	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		// Fill histograms here
		dEgamma->Fill(locBeamPhotonEnergy);

		dMM2_M2pi0->Fill(loc2pi0_P4.M(), locMissingP4.M2());
		dDeltaE_M2pi0->Fill(loc2pi0_P4.M(),locMissingP4.E());
		dMpi0_corr->Fill(locPi01_P4.M(), locPi02_P4.M());

		double locDeltaPhi = (locProtonP4.Phi() - loc2pi0_P4.Phi())*180./TMath::Pi();
		if(locDeltaPhi > 360.) locDeltaPhi -= 360.;
		if(locDeltaPhi < 0.) locDeltaPhi += 360.;
		dPhi2pi0_PhiP->Fill(locProtonP4.Phi()*180./TMath::Pi(), loc2pi0_P4.Phi()*180./TMath::Pi());
		dDeltaPhi_M2pi0->Fill(loc2pi0_P4.M(), locDeltaPhi);

		// require proton and omega are back-to-back
		if(locDeltaPhi < 170. || locDeltaPhi > 190.) return true;
		dMM2_M2pi0_CoplanarTag->Fill(loc2pi0_P4.M(), locMissingP4.M2());
		dDeltaE_M2pi0_CoplanarTag->Fill(loc2pi0_P4.M(),locMissingP4.E());
		dMpi0_corr_CoplanarTag->Fill(locPi01_P4.M(), locPi02_P4.M());

		// tag proton with dE/dx
		if(dEdx > dEdxCut) {
			dMM2_M2pi0_ProtonTag->Fill(loc2pi0_P4.M(), locMissingP4.M2());
			dDeltaE_M2pi0_ProtonTag->Fill(loc2pi0_P4.M(),locMissingP4.E());
			dMpi0_corr_ProtonTag->Fill(locPi01_P4.M(), locPi02_P4.M());
		}

		// cut on pi0 masses
		if(locPi01_P4.M() > 0.10 && locPi01_P4.M() < 0.16 && locPi02_P4.M() > 0.10 && locPi02_P4.M() < 0.16){
			dMM2_M2pi0_Pi0Tag->Fill(loc2pi0_P4.M(), locMissingP4.M2());
			dDeltaE_M2pi0_Pi0Tag->Fill(loc2pi0_P4.M(),locMissingP4.E());
			dMM2_DeltaE_Pi0Tag->Fill(locMissingP4.E(), locMissingP4.M2());
		}
		
		// 2pi invariant mass for exclusive
		if(fabs(locMissingP4.M2()) < 0.05 && fabs(locMissingP4.E()) < 0.5) {
			dMpi0_corr_MMTag->Fill(locPi01_P4.M(), locPi02_P4.M());
			
			if(locPi01_P4.M() > 0.10 && locPi01_P4.M() < 0.16 && locPi02_P4.M() > 0.10 && locPi02_P4.M() < 0.16){
				dEgamma_M2pi0_Pi0Tag->Fill(loc2pi0_P4.M(), locBeamPhotonEnergy);
				if(locBeamPhotonEnergy > 1.5 && locBeamPhotonEnergy < 2.0)
					dDalitz1_Pi0Tag->Fill(locPi01P_P4.M2(),locPi02P_P4.M2());
				if(locBeamPhotonEnergy > 2.0 && locBeamPhotonEnergy < 2.5)
					dDalitz2_Pi0Tag->Fill(locPi01P_P4.M2(),locPi02P_P4.M2());
				if(locBeamPhotonEnergy > 2.5 && locBeamPhotonEnergy < 3.0)
					dDalitz3_Pi0Tag->Fill(locPi01P_P4.M2(),locPi02P_P4.M2());
				if(locBeamPhotonEnergy > 3.0 && locBeamPhotonEnergy < 5.5)
					dDalitz4_Pi0Tag->Fill(locPi01P_P4.M2(),locPi02P_P4.M2());
			}
		}
		
		if(fabs(locMissingP4.M2()) < 0.05 && fabs(locMissingP4.E()) < 0.5) {
			
			dProton_dEdx_P->Fill(locProtonP4.Vect().Mag(), dEdx);
			dProton_P_Theta->Fill(locProtonP4.Vect().Theta()*180/TMath::Pi(), locProtonP4.Vect().Mag());
			
			if(dEdx < dEdxCut) return true;
			
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}
