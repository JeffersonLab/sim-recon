// $Id$
//
//    File: DCustomAction_p2pi_hists.cc
// Created: Wed Jan 21 16:53:41 EST 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "DCustomAction_p2pi_hists.h"

void DCustomAction_p2pi_hists::Initialize(JEventLoop* locEventLoop)
{
	// get PID algos
        const DParticleID* locParticleID = NULL;
        locEventLoop->GetSingle(locParticleID);
        dParticleID = locParticleID;

	// check if a particle is missing
	Get_Reaction()->Get_MissingPID(dMissingPID);

	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//Required: Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this action.
			//If another thread has already created the folder, it just changes to it. 
		CreateAndChangeTo_ActionDirectory();
		
		dEgamma = GetOrCreate_Histogram<TH1I>("Egamma", "TAGGER photon energy; E_{#gamma}", 400, 0., 12.);
		
		if(dMissingPID == Proton) {
			dMM_M2pi = GetOrCreate_Histogram<TH2I>("MM_M2pi", "MM off #pi^{+}#pi^{-} vs M_{#pi^{+}#pi^{-}}; M_{#pi^{+}#pi^{-}}; MM", 200, 0.0, 2.0, 200, 0., 4.);
		}
		else {
			dMM2_M2pi = GetOrCreate_Histogram<TH2I>("MM2_M2pi", "MM^{2} off #pi^{+}#pi^{-} vs M_{#pi^{+}#pi^{-}}; M_{#pi^{+}#pi^{-}}; MM^{2}", 200, 0.0, 2.0, 200, -1., 1.);
			dProton_dEdx_P = GetOrCreate_Histogram<TH2I>("Proton_dEdx_P","dE/dx vs p; p; dE/dx",200,0,2,500,0,5);
			dProton_P_Theta = GetOrCreate_Histogram<TH2I>("Proton_P_Theta","p vs #theta; #theta; p (GeV/c)",180,0,180,120,0,12);
			dDeltaE_M2pi = GetOrCreate_Histogram<TH2I>("DeltaE_M2pi", "#DeltaE vs M_{#pi^{+}#pi^{-}}; M_{#pi^{+}#pi^{-}}; #DeltaE (tagger - tracks)", 200, 0.0, 2.0, 200, -5., 5.);
			dDeltaE_M2pi_ProtonTag = GetOrCreate_Histogram<TH2I>("DeltaE_M2pi_ProtonTag", "#DeltaE vs M_{#pi^{+}#pi^{-}}; M_{#pi^{+}#pi^{-}}; #DeltaE (tagger - tracks)", 200, 0.0, 2.0, 200, -5., 5.);

			dDalitz_p2pi = GetOrCreate_Histogram<TH2I>("Dalitz_p2pi", "Dalitz p#pi^{+}#pi^{-}; M_{#pi^{+}#pi^{-}}^{2}; M_{p#pi^{+}}^{2}", 200, 0.0, 5.0, 200, 0., 5.);
			dMppiplus_M2pi = GetOrCreate_Histogram<TH2I>("Mppiplus_M2pi", "Psuedo-Dalitz p#pi^{+}#pi^{-}; M_{#pi^{+}#pi^{-}}; M_{p#pi^{+}}", 200, 0.0, 2.0, 200, 1.0, 3.0);
			dMppiminus_M2pi = GetOrCreate_Histogram<TH2I>("Mppiminus_M2pi", "Psuedo-Dalitz p#pi^{+}#pi^{-}; M_{#pi^{+}#pi^{-}}; M_{p#pi^{-}}", 200, 0.0, 2.0, 200, 1.0, 3.0);

			dProtonPhi_Egamma = GetOrCreate_Histogram<TH2I>("ProtonPhi_Egamma","#phi vs E_{#gamma}; E_{#gamma}; proton #phi",240,0,6,360,-180,180);	
			dPiPlusPsi_Egamma = GetOrCreate_Histogram<TH2I>("PiPlusPsi_Egamma","#psi vs E_{#gamma}; E_{#gamma}; #pi^{+} #phi",240,0,6,360,-180,180);
			dPiPlusPsi_t = GetOrCreate_Histogram<TH2I>("PiPlusPsi_t","#psi vs |t|; |t|; #pi^{+} #psi",1000,0.,5.,360,-180,180);

			dEgamma_M2pi = GetOrCreate_Histogram<TH2I>("Egamma_M2pi", "E_{#gamma} vs M_{#pi^{+}#pi^{-}}; M_{#pi^{+}#pi^{-}}; E_{#gamma}", 200, 0.0, 2.0, 240, 0., 6.);		
			
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DCustomAction_p2pi_hists::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{

	const DDetectorMatches* locDetectorMatches = NULL;
        locEventLoop->GetSingle(locDetectorMatches);

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

	// cut on tagger energy
	if(locBeamPhotonEnergy < 1.5) 
		return false;

	// calculate 2pi P4
	DLorentzVector locP4_2pi = locParticles[0]->lorentzMomentum() + locParticles[1]->lorentzMomentum();

	// calculate missing P4
	DLorentzVector locMissingP4 = dAnalysisUtilities->Calc_MissingP4(locParticleCombo,Get_UseKinFitResultsFlag());

	// reconstructed proton variables
	double dEdx = 0.;
	DLorentzVector locProtonP4;
	if(dMissingPID != Proton) {
		const DChargedTrack* locChargedTrack = static_cast<const DChargedTrack*>(locParticleComboStep->Get_FinalParticle_SourceObject(2));
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(Proton);
		dEdx = locChargedTrackHypothesis->dEdx()*1e6; // convert to keV
		locProtonP4 = locChargedTrackHypothesis->lorentzMomentum();	
	}

	// production kinematics
	DLorentzVector locSumInitP4; 
	DLorentzVector locProtonP4Init(0,0,0,0.938);
	locSumInitP4 += locProtonP4Init;
	locSumInitP4 += locBeamPhoton->lorentzMomentum();
	TLorentzVector locDelta = (locProtonP4 - locProtonP4Init);
	double t = locDelta.M2();

	TVector3 locBoostVector = -1.*locP4_2pi.BoostVector();
	TLorentzVector locPiPlus_P4 = locParticles[0]->lorentzMomentum();
	locPiPlus_P4.Boost(locBoostVector);
	double locPsi = locPiPlus_P4.Phi();

	double dEdxCut = 0.8;
	if(locEventLoop->GetJEvent().GetRunNumber() == 2931)
		dEdxCut = 2.0;
	if(locEventLoop->GetJEvent().GetRunNumber() == 3079)
		dEdxCut = 0.9;

	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		// Fill histograms here
		dEgamma->Fill(locBeamPhotonEnergy);

		if(dMissingPID == Proton) {
			dMM_M2pi->Fill(locP4_2pi.M(), locMissingP4.M());

		}
		else {
			dMM2_M2pi->Fill(locP4_2pi.M(), locMissingP4.M2());
			dDeltaE_M2pi->Fill(locP4_2pi.M(),locMissingP4.E());

			if(fabs(locMissingP4.M2()) < 0.02){
				dProton_dEdx_P->Fill(locProtonP4.Vect().Mag(), dEdx);
				dProton_P_Theta->Fill(locProtonP4.Vect().Theta()*180/TMath::Pi(), locProtonP4.Vect().Mag());

				if(dEdx < dEdxCut) return false;

				dEgamma_M2pi->Fill(locP4_2pi.M(), locBeamPhotonEnergy);
				dDeltaE_M2pi_ProtonTag->Fill(locP4_2pi.M(),locMissingP4.E());

				DLorentzVector locP4_ppiplus = locProtonP4;
				DLorentzVector locP4_ppiminus = locProtonP4;
				locP4_ppiplus += locParticles[0]->lorentzMomentum();
				locP4_ppiminus += locParticles[1]->lorentzMomentum();
				if(locBeamPhotonEnergy > 2.5 && locBeamPhotonEnergy < 3.0){
					dDalitz_p2pi->Fill(locP4_2pi.M2(), locP4_ppiplus.M2());
					dMppiplus_M2pi->Fill(locP4_2pi.M(), locP4_ppiplus.M());
					dMppiminus_M2pi->Fill(locP4_2pi.M(), locP4_ppiminus.M());
				}

				if(locP4_2pi.M() > 0.6 && locP4_2pi.M() < 0.9){
					dProtonPhi_Egamma->Fill(locBeamPhotonEnergy, locProtonP4.Phi()*180/TMath::Pi());
					dPiPlusPsi_Egamma->Fill(locBeamPhotonEnergy, locPsi*180/TMath::Pi());
					
					if(locBeamPhotonEnergy > 2.5 && locBeamPhotonEnergy < 3.0){
						dPiPlusPsi_t->Fill(fabs(t), locPsi*180/TMath::Pi());
					}
					
					japp->RootUnLock(); //RELEASE ROOT LOCK!!
					return true;
				}			
			}
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!

	return false; //return false if you want to use this action to apply a cut (and it fails the cut!)
}
