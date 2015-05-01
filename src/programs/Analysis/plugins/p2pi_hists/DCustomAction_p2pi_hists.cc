// $Id$
//
//    File: DCustomAction_p2pi_hists.cc
// Created: Wed Jan 21 16:53:41 EST 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "DCustomAction_p2pi_hists.h"

void DCustomAction_p2pi_hists::Initialize(JEventLoop* locEventLoop)
{
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
			dVertexDeltaZ_M2pi = GetOrCreate_Histogram<TH2I>("DeltaZ_M2pi", "Vertex #DeltaZ of #pi^{+}#pi^{-} vs M_{#pi^{+}#pi^{-}}; M_{#pi^{+}#pi^{-}}; #pi^{+}#pi^{-} vertex #DeltaZ (cm)", 200, 0.0, 2.0, 200, -50, 50.);
			dVertexDeltaMag_M2pi = GetOrCreate_Histogram<TH2I>("DeltaMag_M2pi", "Vertex #DeltaMag of #pi^{+}#pi^{-} vs M_{#pi^{+}#pi^{-}}; M_{#pi^{+}#pi^{-}}; #pi^{+}#pi^{-} vertex #Delta Magnitude (cm)", 200, 0.0, 2.0, 200, 0, 50.);
		}
		else {
			dMM2_M2pi = GetOrCreate_Histogram<TH2I>("MM2_M2pi", "MM^{2} off #pi^{+}#pi^{-} vs M_{#pi^{+}#pi^{-}}; M_{#pi^{+}#pi^{-}}; MM^{2}", 200, 0.0, 2.0, 200, -1., 1.);
			dProton_dEdx_P = GetOrCreate_Histogram<TH2I>("Proton_dEdx_P","dE/dx vs p; p; dE/dx",200,0,2,500,0,5);
			dProton_P_Theta = GetOrCreate_Histogram<TH2I>("Proton_P_Theta","p vs #theta; #theta; p (GeV/c)",180,0,180,120,0,12);
			dDeltaE_M2pi = GetOrCreate_Histogram<TH2I>("DeltaE_M2pi", "#DeltaE vs M_{#pi^{+}#pi^{-}}; M_{#pi^{+}#pi^{-}}; #DeltaE (tagger - tracks)", 200, 0.0, 2.0, 200, -5., 5.);
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DCustomAction_p2pi_hists::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
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
	
	// cut on tagger energy
	if(locBeamPhotonEnergy < 2.5) 
		return true;

	// calculate 2pi P4
	DLorentzVector locP4_2pi = locParticles[0]->lorentzMomentum() + locParticles[1]->lorentzMomentum();

	// calculate 2pi vertex
        DLorentzVector locX4_2pi = locParticleComboStep->Get_SpacetimeVertex();
        DVector3 locDiplacedVertex = locX4_2pi.Vect() - DVector3(0, 0, 65);

	// calculate missing P4
	DLorentzVector locMissingP4 = dAnalysisUtilities->Calc_MissingP4(locParticleCombo,Get_UseKinFitResultsFlag());

	// reconstructed proton variables
	double dEdx = 0.;
	DLorentzVector locProtonP4;
	if(dMissingPID != Proton) {
		const DChargedTrack* locChargedTrack = static_cast<const DChargedTrack*>(locParticleComboStep->Get_FinalParticle_SourceObject(2));
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(Proton);
		dEdx = locChargedTrackHypothesis->dEdx()*1e6;
		locProtonP4 = locChargedTrackHypothesis->lorentzMomentum();	
	}

	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		// Fill histograms here
		dEgamma->Fill(locBeamPhotonEnergy);

		if(dMissingPID == Proton) {
			dMM_M2pi->Fill(locP4_2pi.M(), locMissingP4.M());

			dVertexDeltaZ_M2pi->Fill(locP4_2pi.M(), locDiplacedVertex.Z());
			dVertexDeltaMag_M2pi->Fill(locP4_2pi.M(), locDiplacedVertex.Mag());
		}
		else {
			dMM2_M2pi->Fill(locP4_2pi.M(), locMissingP4.M2());
			dDeltaE_M2pi->Fill(locP4_2pi.M(),locMissingP4.E());

			if(fabs(locMissingP4.M2()) < 0.05){
				dProton_dEdx_P->Fill(locProtonP4.Vect().Mag(), dEdx);
				dProton_P_Theta->Fill(locProtonP4.Vect().Theta()*180/TMath::Pi(), locProtonP4.Vect().Mag());
			}
		}
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}
