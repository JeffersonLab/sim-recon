// $Id$
//
//  DCustomAction_Z2pi_trees.cc. Use for Primakoff events where the heavy Z target is never detected
//    custom action modeled after Justin's 
//
//    File: DCustomAction_p2pi_hists.cc
// Created: Wed Jan 21 16:53:41 EST 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "DCustomAction_Z2pi_trees.h"

void DCustomAction_Z2pi_trees::Initialize(JEventLoop* locEventLoop)
{
	DApplication* dapp=dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	JCalibration *jcalib = dapp->GetJCalibration((locEventLoop->GetJEvent()).GetRunNumber());

	// Parameters for event selection to fill histograms
	endpoint_energy = 12.;
	map<string, double> photon_endpoint_energy;
	if(jcalib->Get("/PHOTON_BEAM/endpoint_energy", photon_endpoint_energy) == false) {
		endpoint_energy = photon_endpoint_energy["PHOTON_BEAM_ENDPOINT_ENERGY"];
	}
	else {
		jout<<"No /PHOTON_BEAM/endpoint_energy for this run number: using default of 12 GeV"<<endl;
	}
	endpoint_energy_bins = (int)(20*endpoint_energy);

	cohmin_energy = 0.;
	cohedge_energy = 12.;
	map<string, double> photon_beam_param;
	if(jcalib->Get("/ANALYSIS/beam_asymmetry/coherent_energy", photon_beam_param) == false) {
		cohmin_energy = photon_beam_param["cohmin_energy"];
		cohedge_energy = photon_beam_param["cohedge_energy"];
	}
	else {
		jout<<"No /ANALYSIS/beam_asymmetry/coherent_energy for this run number: using default range of 0-12 GeV"<<endl;
	}

	dEdxCut = 2.2;
	minMMCut = 0.8;
	maxMMCut = 1.05;
	minMM2Cut = 30000;  // recal target particle is missing
	maxMM2Cut = 50000;
	missingEnergyCut = 20.;
	min2piMassCut = 0.2;
	max2piMassCut = 0.6;
	Mtarg = 208.;

	// get PID algos
        const DParticleID* locParticleID = NULL;
        locEventLoop->GetSingle(locParticleID);
        dParticleID = locParticleID;

	locEventLoop->GetSingle(dAnalysisUtilities);

	// check if a particle is missing
	Get_Reaction()->Get_MissingPID(dMissingPID);

	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!

	//Required: Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this action.
	//If another thread has already created the folder, it just changes to it. 
	CreateAndChangeTo_ActionDirectory();
		
	dEgamma = GetOrCreate_Histogram<TH1I>("Egamma", "TAGGER photon energy; E_{#gamma}", endpoint_energy_bins, 0., endpoint_energy);

	dMM2_M2pi = GetOrCreate_Histogram<TH2I>("MM2_M2pi", "MM^{2} off #pi^{+}#pi^{-} vs M_{#pi^{+}#pi^{-}}; M_{#pi^{+}#pi^{-}}; MM^{2}", 200, 0.0, 2.0, 200, 30000,50000);
	dDeltaE_M2pi = GetOrCreate_Histogram<TH2I>("DeltaE_M2pi", "#DeltaE vs M_{#pi^{+}#pi^{-}}; M_{#pi^{+}#pi^{-}}; #DeltaE (tagger - tracks -Mtarg)", 200, 0.0, 2.0, 200, -20,20);

	dPiPlusPhi_Egamma = GetOrCreate_Histogram<TH2I>("PiPlusPhi_Egamma","#phi vs E_{#gamma}; E_{#gamma}; #pi^{+} #phi",endpoint_energy_bins,0,endpoint_energy,360,-180,180);	
	dPiPlusPsi_Egamma = GetOrCreate_Histogram<TH2I>("PiPlusPsi_Egamma","#psi vs E_{#gamma}; E_{#gamma}; #pi^{+} #phi",endpoint_energy_bins,0,endpoint_energy,360,-180,180);
	dPiPlusPsi_t = GetOrCreate_Histogram<TH2I>("PiPlusPsi_t","#psi vs |t|; |t|; #pi^{+} #psi",1000,0.,5.,360,-180,180);

	dEgamma_M2pi = GetOrCreate_Histogram<TH2I>("Egamma_M2pi", "E_{#gamma} vs M_{#pi^{+}#pi^{-}}; M_{#pi^{+}#pi^{-}}; E_{#gamma}", 200, 0.0, 1.0, endpoint_energy_bins,0,endpoint_energy);		

	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DCustomAction_Z2pi_trees::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
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
	DLorentzVector locP4_2pi, locP4_p1;
	for(size_t loc_i = 0; loc_i < 3; ++loc_i) {
		if(locParticles[loc_i] == NULL) continue; // missing particle (i.e. target)
                if(locParticles[loc_i]->PID() == PiPlus || locParticles[loc_i]->PID() == PiMinus)
                        locP4_2pi += locParticles[loc_i]->lorentzMomentum();
                if(locParticles[loc_i]->PID() == PiPlus)
                        locP4_p1 = locParticles[loc_i]->lorentzMomentum();
        }

	// calculate missing P4
	DLorentzVector locMissingP4 = dAnalysisUtilities->Calc_MissingP4(locParticleCombo,Get_UseKinFitResultsFlag());

	// cout << "Missing P4 "; locMissingP4.Print();

	// reconstructed target  variables
	DLorentzVector locZP4;
	locZP4 = locMissingP4;

	// production kinematics
	DLorentzVector locSumInitP4; 
	DLorentzVector locZP4Init(0,0,0,Mtarg);   // assume Pb target
	locSumInitP4 += locZP4Init;
	locSumInitP4 += locBeamPhoton->lorentzMomentum();
	TLorentzVector locDelta = (locP4_2pi - locBeamPhoton->lorentzMomentum()); // calculate t from 2pi measured momenta
	double t = locDelta.M2();

        double phipol = 0;                           // *** Note assumes horizontal polarization plane.
        TVector3 eps(cos(phipol), sin(phipol), 0.0); // beam polarization vector in lab

	// boost to resonance frame for angular distributions 	
	TLorentzRotation resonanceBoost( -locP4_2pi.BoostVector() );
	TLorentzVector recoil_res = resonanceBoost * locZP4;     // this quantity is poorly defined from missing momenta but calculable
	TLorentzVector p1_res = resonanceBoost * locP4_p1;

        // choose helicity frame: z-axis opposite recoil target in rho rest frame. Note that for Primakoff recoil is missin P4
	TVector3 y = (locBeamPhoton->lorentzMomentum().Vect().Unit().Cross(-locZP4.Vect().Unit())).Unit();
	
	// choose helicity frame: z-axis opposite recoil proton in rho rest frame
	TVector3 z = -1. * recoil_res.Vect().Unit();
	TVector3 x = y.Cross(z).Unit();
	TVector3 angles( (p1_res.Vect()).Dot(x),
			 (p1_res.Vect()).Dot(y),
			 (p1_res.Vect()).Dot(z) );
	
	// double CosTheta = angles.CosTheta();

	double phi = angles.Phi();

	double Phi = atan2(y.Dot(eps), locBeamPhoton->lorentzMomentum().Vect().Unit().Dot(eps.Cross(y)));

        double locPsi = Phi - phi;               // angle difference
	if(locPsi < -1*TMath::Pi()) locPsi += 2*TMath::Pi();
        if(locPsi > TMath::Pi()) locPsi -= 2*TMath::Pi();

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action(); //ACQUIRE ROOT LOCK!!
	{
		// Fill histograms here
		dEgamma->Fill(locBeamPhotonEnergy);

		dMM2_M2pi->Fill(locP4_2pi.M(), locMissingP4.M2());
		dDeltaE_M2pi->Fill(locP4_2pi.M(),locMissingP4.E()-Mtarg);

		if(locMissingP4.M2() > minMM2Cut && locMissingP4.M2() < maxMM2Cut && fabs(locMissingP4.E()) -Mtarg  < missingEnergyCut){  // subtract target mass

		  // fill histograms if missing mass2/energy are within limits

			    dEgamma_M2pi->Fill(locP4_2pi.M(), locBeamPhotonEnergy);
					
			    if(locP4_2pi.M() > min2piMassCut && locP4_2pi.M() < max2piMassCut){
			      dPiPlusPhi_Egamma->Fill(locBeamPhotonEnergy, locP4_p1.Vect().Phi()*180/TMath::Pi());
				 dPiPlusPsi_Egamma->Fill(locBeamPhotonEnergy, locPsi*180/TMath::Pi());
						
			         if(locBeamPhotonEnergy > cohmin_energy && locBeamPhotonEnergy < cohedge_energy){
				      dPiPlusPsi_t->Fill(fabs(t), locPsi*180/TMath::Pi());	
			         }
			    }		
		}
	}
	Unlock_Action(); //RELEASE ROOT LOCK!!

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}
