// $Id$
//
//    File: DEventProcessor_bcal_hists.cc
// Created: Mon Apr  3 11:38:03 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.4.0 powerpc)
//

#include "DEventProcessor_bcal_hists.h"

#include <TLorentzVector.h>

#include "DANA/DApplication.h"
#include "BCAL/DBCALShower.h"
#include "BCAL/DBCALTruthShower.h"
#include "BCAL/DHDDMBCALHit.h"
#include "FCAL/DFCALCluster.h"
#include "TRACKING/DMCThrown.h"

// The executable should define the ROOTfile global variable. It will
// be automatically linked when dlopen is called.
extern TFile *ROOTfile;

// Routine used to create our DEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_bcal_hists());
}
} // "C"


#define FCAL_Z_OFFSET 640.0-65.0 // I don't know what this value is ???

//------------------
// init
//------------------
jerror_t DEventProcessor_bcal_hists::init(void)
{
	// Create THROWN directory
	TDirectory *dir = new TDirectoryFile("BCAL","BCAL");
	dir->cd();
	
	two_gamma_mass = new TH1F("two_gamma_mass","two_gamma_mass",100, 0.0, 0.300);
	two_gamma_mass_corr = new TH1F("two_gamma_mass_corr","two_gamma_mass_corr",100, 0.0, 0.300);
	two_gamma_mass_cut = new TH1F("two_gamma_mass_cut","two_gamma_mass_cut",100, 0.0, 0.300);
	bcal_fcal_two_gamma_mass = new TH1F("fcal_two_gamma_mass","bcal_fcal_two_gamma_mass",100, 0.0, 0.300);
	bcal_fcal_two_gamma_mass_cut = new TH1F("fcal_two_gamma_mass_cut","bcal_fcal_two_gamma_mass_cut",100, 0.0, 0.300);
	xy_shower = new TH2F("xy_shower","xy_shower",100, -100.0, 100., 100 , -100.0, 100.0);
	z_shower = new TH1F("z_shower","z_shower",450, -50.0, 400);
	E_shower = new TH1F("E_shower","E_shower", 200, 0.0, 6.0);
	
	Erec_over_Ethrown_vs_z = new TH2F("Erec_over_Ethrown_vs_z","Erec_over_Ethrown_vs_z", 200, -50.0, 600.0, 200, 0.0, 2.0);
	Ecorr_over_Erec_vs_z = new TH2F("Ecorr_over_Erec_vs_z","Ecorr_over_Erec_vs_z", 200, -50.0, 600.0, 200, 0.0, 4.0);
	Ereconstructed_vs_Ethrown = new TH2F("Ereconstructed_vs_Ethrown","BCAL total reconstructed E to total thrown E", 200, 0.0, 6.0, 200, 0.0, 6.0);

	Etot_truth = new TH1F("Etot_truth", "Sum of all truth showers (GeV)", 200, 0.0, 6.0);
	Etot_hits = new TH1F("Etot_hits", "Sum of all hits (GeV)", 200, 0.0, 6.0);
	Etruth_over_Ethrown_vs_z = new TH2F("Etruth_over_Ethrown_vs_z","Etruth_over_Ethrown_vs_z", 200, -50.0, 600.0, 200, 0.0, 2.0);

	Edeposited_over_Ethrown_vs_z = new TH2F("Edeposited_over_Ethrown_vs_z","Edeposited_over_Ethrown_vs_z", 80, -50.0, 450.0, 200, 0.0, 2.0);

	// Go back up to the parent directory
	dir->cd("../");

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_bcal_hists::brun(JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_bcal_hists::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DBCALShower*> showers;
	vector<const DFCALCluster*> fcal_showers;
	vector<const DBCALTruthShower*> truthshowers;	
	vector<const DMCThrown*> mcthrowns;
	vector<const DHDDMBCALHit*> bcalhits;
	loop->Get(showers, "KLOE" );
	//loop->Get(fcal_showers);
	loop->Get(truthshowers);
	loop->Get(mcthrowns);
	loop->Get(bcalhits);
	
	LockState();
	
	// Single shower params
	double Etot_reconstructed = 0.0;
	for(unsigned int i=0; i<showers.size(); i++){
		const DBCALShower *s = showers[i];
		xy_shower->Fill(s->x, s->y);
		z_shower->Fill(s->z);
		E_shower->Fill(s->Ecorr);
		Etot_reconstructed += s->Ecorr;
	}
	
	// 2-gamma inv. mass
	for(unsigned int i=0; i<showers.size(); i++){
		const DBCALShower *s1 = showers[i];
		double dx = s1->x;
		double dy = s1->y;
		double dz = s1->z - 65.0;
		double R = sqrt(dx*dx + dy*dy + dz*dz);
		double E = s1->Ecorr;
		double Edave = s1->Ecorr*(1.106+(dz+65.0-208.4)*(dz+65.0-208.4)*6.851E-6);
		TLorentzVector p1(E*dx/R, E*dy/R, E*dz/R, E);		
		TLorentzVector p1dave(Edave*dx/R, Edave*dy/R, Edave*dz/R, Edave);		
		
		for(unsigned int j=i+1; j<showers.size(); j++){
			const DBCALShower *s2 = showers[j];
			dx = s2->x;
			dy = s2->y;
			dz = s2->z - 65.0; // shift to coordinate relative to center of target
			R = sqrt(dx*dx + dy*dy + dz*dz);
			double E = s2->Ecorr;
			double Edave = s2->Ecorr*(1.106+(dz+65.0-208.4)*(dz+65.0-208.4)*6.851E-6);
			TLorentzVector p2(E*dx/R, E*dy/R, E*dz/R, E);		
			TLorentzVector p2dave(Edave*dx/R, Edave*dy/R, Edave*dz/R, Edave);		
			
			TLorentzVector ptot = p1+p2;
			two_gamma_mass->Fill(ptot.M());
			TLorentzVector ptotdave = p1dave+p2dave;
			two_gamma_mass_corr->Fill(ptotdave.M());
			
			if(showers.size()==2)two_gamma_mass_cut->Fill(ptotdave.M());
		}
		
		for(unsigned int j=0; j<fcal_showers.size(); j++){
			const DFCALCluster *s2 = fcal_showers[j];
			dx = s2->getCentroid().X();
			dy = s2->getCentroid().Y();
			dz = FCAL_Z_OFFSET; // shift to coordinate relative to center of target
			R = sqrt(dx*dx + dy*dy + dz*dz);
			double E = s2->getEnergy();
			TLorentzVector p2(E*dx/R, E*dy/R, E*dz/R, E);		
			
			TLorentzVector ptot = p1dave+p2;
			bcal_fcal_two_gamma_mass->Fill(ptot.M());

			if(showers.size()==1 && fcal_showers.size()==1){
				bcal_fcal_two_gamma_mass_cut->Fill(ptot.M());
			}
		}
		
	}
	
	// Total energy in hits
	double Ehit_tot = 0.0;
	for(unsigned int i=0; i<bcalhits.size(); i++){
		Ehit_tot += bcalhits[i]->E;
	}
	Etot_hits->Fill(Ehit_tot);
	
	// Truth values
	double Etruth_tot = 0.0;
	double z_truth = 0.0;
	for(unsigned int i=0; i<truthshowers.size(); i++){
		Etruth_tot += truthshowers[i]->E;
		z_truth += truthshowers[i]->E*truthshowers[i]->z;
	}
	z_truth/=Etruth_tot;
	Etot_truth->Fill(Etruth_tot);

	// Compare to thrown values
	double Etot_thrown=0.0;
	for(unsigned int i=0; i<mcthrowns.size(); i++){
		Etot_thrown += mcthrowns[i]->energy();
		for(unsigned int j=0; j<showers.size(); j++){
			double z = showers[j]->z;
			Erec_over_Ethrown_vs_z->Fill(z, showers[j]->Ecorr/mcthrowns[i]->energy());

			double Ecorr = showers[j]->Ecorr*(1.106+(z-208.4)*(z-208.4)*6.851E-6);
			Ecorr_over_Erec_vs_z->Fill(z, mcthrowns[i]->energy()/Ecorr);
		}
	}
	
	Ereconstructed_vs_Ethrown->Fill(Etot_thrown, Etot_reconstructed);
	Etruth_over_Ethrown_vs_z->Fill(z_truth, Etruth_tot/Etot_thrown);

	// Single thrown particle
	if(mcthrowns.size()==1){
		const DMCThrown* mcthrown = mcthrowns[0];
		if(mcthrown->momentum().Theta()>0.0001){
			double z = mcthrown->position().Z() + 65.0/tan(mcthrown->momentum().Theta());
			double Ethrown = 1.0; // for some reason, mcthrown->E is zero right now.
			// I fudge this for now since I know all of the events thrw 1.0GeV
			Edeposited_over_Ethrown_vs_z->Fill(z, Ehit_tot/Ethrown);
		}
	}

	UnlockState();	

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_bcal_hists::erun(void)
{
	// Any final calculations on histograms (like dividing them)
	// should be done here. This may get called more than once.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_bcal_hists::fini(void)
{
	return NOERROR;
}

