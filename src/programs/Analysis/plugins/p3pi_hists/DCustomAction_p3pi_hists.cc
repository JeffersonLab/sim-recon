// $Id$
//
//    File: DCustomAction_p3pi_hists.cc
// Created: Wed Jan 21 16:53:41 EST 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "DCustomAction_p3pi_hists.h"

void DCustomAction_p3pi_hists::Initialize(JEventLoop* locEventLoop)
{
	// check if a particle is missing
        Get_Reaction()->Get_MissingPID(dMissingPID);

	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//Required: Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this action.
			//If another thread has already created the folder, it just changes to it. 
		CreateAndChangeTo_ActionDirectory();

		dEgamma = GetOrCreate_Histogram<TH1I>("Egamma", "TAGGER photon energy; E_{gamma}", 400, 2., 12.);
		
		if(dMissingPID == Proton)
			dMM_M3pi = GetOrCreate_Histogram<TH2I>("MM_M3pi", "MM off #pi^{+}#pi^{-}#pi^{0} vs M_{#pi^{+}#pi^{-}#pi^{0}}; M_{#pi^{+}#pi^{-}#pi^{0}}; MM", 200, 0.0, 2.0, 200, 0., 4.);
		else {
			dMM2_M3pi = GetOrCreate_Histogram<TH2I>("MM2_M3pi", "MM^{2} off #pi^{+}#pi^{-}#pi^{0} vs M_{#pi^{+}#pi^{-}#pi^{0}}; M_{#pi^{+}#pi^{-}#pi^{0}}; MM^2", 200, 0.0, 2.0, 200, -1., 1.);
                        dProton_dEdx_P = GetOrCreate_Histogram<TH2I>("Proton_dEdx_P","dE/dx vs p; p; dE/dx",200,0,2,500,0,5);
                        dProton_P_Theta = GetOrCreate_Histogram<TH2I>("Proton_P_Theta","p vs #theta; #theta; p (GeV/c)",180,0,180,120,0,12);
                        dDeltaE_M3pi = GetOrCreate_Histogram<TH2I>("dDeltaE_M3pi", "#Delta E vs M_{#pi^{+}#pi^{-}#pi^{0}}; M_{#pi^{+}#pi^{-}#pi^{0}}; #Delta E (tagger - tracks)", 200, 0.0, 2.0, 200, -5., 5.);
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DCustomAction_p3pi_hists::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(0);

        // get beam photon energy
        const DKinematicData* locBeamPhoton = locParticleComboStep->Get_InitialParticle(); 
        double locBeamPhotonEnergy = locBeamPhoton->energy();

	// calculate missing mass
	DLorentzVector locMissingP4 = dAnalysisUtilities->Calc_MissingP4(locParticleCombo, Get_UseKinFitResultsFlag());
	
	// calculate 3pi mass
	DLorentzVector locOmegaP4 = dAnalysisUtilities->Calc_FinalStateP4(locParticleCombo, 1, Get_UseKinFitResultsFlag());

	// reconstructed proton variables
	double dEdx = 0.;
	DLorentzVector locProtonP4;
	if(dMissingPID != Proton) {
		const DChargedTrack* locChargedTrack = static_cast<const DChargedTrack*>(locParticleComboStep->Get_FinalParticle_SourceObject(1));
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(Proton);
		dEdx = locChargedTrackHypothesis->dEdx()*1e3;
		locProtonP4 = locChargedTrackHypothesis->lorentzMomentum();	
	}

	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		// Fill histograms here
		dEgamma->Fill(locBeamPhotonEnergy);

		if(dMissingPID == Proton) {
			dMM_M3pi->Fill(locOmegaP4.M(), locMissingP4.M());
		}
		else {
			dMM2_M3pi->Fill(locOmegaP4.M(), locMissingP4.M2());
                        dDeltaE_M3pi->Fill(locOmegaP4.M(),locMissingP4.E());

			// for omega candidates look at recoil proton
			if(locOmegaP4.M() > 0.75 && locOmegaP4.M() < 0.875){

				dProton_dEdx_P->Fill(locProtonP4.Vect().Mag(), dEdx);
				dProton_P_Theta->Fill(locProtonP4.Vect().Theta()*180/TMath::Pi(), locProtonP4.Vect().Mag());
			}
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}
