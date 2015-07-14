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

			dBaryonM_CosTheta_Egamma1 = GetOrCreate_Histogram<TH2I>("BaryonM_CosTheta_Egamma1", "Baryon M_{p#pi^{+}} vs cos#theta: E_{#gamma} < 2.5; cos#theta; M_{p#pi^{+}}", 100, -1, 1, 100, 1.0, 2.5);
			dBaryonM_CosTheta_Egamma2 = GetOrCreate_Histogram<TH2I>("BaryonM_CosTheta_Egamma2", "Baryon M_{p#pi^{+}} vs cos#theta: 2.5 < E_{#gamma} < 3.0; cos#theta; M_{p#pi^{+}}", 100, -1, 1, 100, 1.0, 2.5);
			dBaryonM_CosTheta_Egamma3 = GetOrCreate_Histogram<TH2I>("BaryonM_CosTheta_Egamma3", "Baryon M_{p#pi^{+}} vs cos#theta: E_{#gamma} > 3.0; cos#theta; M_{p#pi^{+}}", 100, -1, 1, 100, 1.0, 2.5);
			
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

	// calculate 2pi P4
	DLorentzVector locP4_2pi;
	for(size_t loc_i = 0; loc_i < 3; ++loc_i) {
		if(locParticles[loc_i] == NULL) continue; // missing proton
                if(locParticles[loc_i]->PID() == PiPlus || locParticles[loc_i]->PID() == PiMinus)
                        locP4_2pi += locParticles[loc_i]->lorentzMomentum();
        }

	// calculate missing P4
	DLorentzVector locMissingP4 = dAnalysisUtilities->Calc_MissingP4(locParticleCombo,Get_UseKinFitResultsFlag());

	// reconstructed proton variables
	double dEdx = 0.;
	DLorentzVector locProtonP4;
	if(dMissingPID != Proton) {
		const DChargedTrack* locChargedTrack = static_cast<const DChargedTrack*>(locParticleComboStep->Get_FinalParticle_SourceObject(0));
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(Proton);
		dEdx = locChargedTrackHypothesis->dEdx()*1e6; // convert to keV
		locProtonP4 = locChargedTrackHypothesis->lorentzMomentum();	
	}

	double dEdxCut = 2.2;

	// production kinematics
	DLorentzVector locSumInitP4; 
	DLorentzVector locProtonP4Init(0,0,0,0.938);
	locSumInitP4 += locProtonP4Init;
	locSumInitP4 += locBeamPhoton->lorentzMomentum();
	TLorentzVector locDelta = (locProtonP4 - locProtonP4Init);
	double t = locDelta.M2();

	TVector3 locBoostVector = -1.*locP4_2pi.BoostVector();
	TLorentzVector locPiPlus_P4 = locParticles[1]->lorentzMomentum();
	locPiPlus_P4.Boost(locBoostVector);
	double locPsi = locPiPlus_P4.Phi();

	// copied from TwoPiAngles
	TLorentzRotation resonanceBoost( -locP4_2pi.BoostVector() );
	
	TLorentzVector beam_res = resonanceBoost * locBeamPhoton->lorentzMomentum();
	TLorentzVector recoil_res = resonanceBoost * locProtonP4;
	TLorentzVector p1_res = resonanceBoost * locParticles[1]->lorentzMomentum();
	
	TVector3 z = -recoil_res.Vect().Unit();
	TVector3 y = beam_res.Vect().Cross(z).Unit();
	TVector3 x = y.Cross(z).Unit();
	
	TVector3 angles(   (p1_res.Vect()).Dot(x),
			   (p1_res.Vect()).Dot(y),
			   (p1_res.Vect()).Dot(z) );

	double cosTheta = angles.CosTheta();

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

			if(fabs(locMissingP4.M2()) < 0.02 && fabs(locMissingP4.E()) < 0.2){
				dProton_dEdx_P->Fill(locProtonP4.Vect().Mag(), dEdx);
				dProton_P_Theta->Fill(locProtonP4.Vect().Theta()*180/TMath::Pi(), locProtonP4.Vect().Mag());

				if(dEdx > dEdxCut) {
					
					dEgamma_M2pi->Fill(locP4_2pi.M(), locBeamPhotonEnergy);
					dDeltaE_M2pi_ProtonTag->Fill(locP4_2pi.M(),locMissingP4.E());
					
					DLorentzVector locP4_ppiplus = locProtonP4;
					DLorentzVector locP4_ppiminus = locProtonP4;
					locP4_ppiplus += locParticles[1]->lorentzMomentum();
					locP4_ppiminus += locParticles[2]->lorentzMomentum();
					if(locBeamPhotonEnergy > 2.5 && locBeamPhotonEnergy < 3.0){
						dDalitz_p2pi->Fill(locP4_2pi.M2(), locP4_ppiplus.M2());
						dMppiplus_M2pi->Fill(locP4_2pi.M(), locP4_ppiplus.M());
						dMppiminus_M2pi->Fill(locP4_2pi.M(), locP4_ppiminus.M());
					}
					
					// some p pi+ mass distributions (hunt for Delta++...)
					if(locBeamPhotonEnergy < 2.5)
						dBaryonM_CosTheta_Egamma1->Fill(cosTheta,locP4_ppiplus.M());
					else if(locBeamPhotonEnergy < 3.0)
						dBaryonM_CosTheta_Egamma2->Fill(cosTheta,locP4_ppiplus.M());
					else 
						dBaryonM_CosTheta_Egamma3->Fill(cosTheta,locP4_ppiplus.M());


					if(locP4_2pi.M() > 0.6 && locP4_2pi.M() < 0.9){
						dProtonPhi_Egamma->Fill(locBeamPhotonEnergy, locProtonP4.Phi()*180/TMath::Pi());
						dPiPlusPsi_Egamma->Fill(locBeamPhotonEnergy, locPsi*180/TMath::Pi());
						
						if(locBeamPhotonEnergy > 2.5 && locBeamPhotonEnergy < 3.0){
							dPiPlusPsi_t->Fill(fabs(t), locPsi*180/TMath::Pi());	
						}
					}	
				}	
			}
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}
