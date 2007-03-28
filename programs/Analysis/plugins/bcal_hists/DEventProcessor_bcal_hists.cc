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
#include "FCAL/DFCALShower.h"
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


#define BCAL_Z_OFFSET 26.03+123.4+65.0 //convert regina's coordinated to HDGeant's
#define FCAL_Z_OFFSET 640.0-65.0 // I don't know what this value is ???

//------------------
// init
//------------------
jerror_t DEventProcessor_bcal_hists::init(void)
{
	// open ROOT file
	//if(ROOTfile != NULL) ROOTfile->cd();
	
	two_gamma_mass = new TH1F("bcal_two_gamma_mass","two_gamma_mass",100, 0.0, 0.300);
	two_gamma_mass_corr = new TH1F("bcal_two_gamma_mass_corr","two_gamma_mass_corr",100, 0.0, 0.300);
	two_gamma_mass_cut = new TH1F("bcal_two_gamma_mass_cut","two_gamma_mass_cut",100, 0.0, 0.300);
	bcal_fcal_two_gamma_mass = new TH1F("bcal_fcal_two_gamma_mass","bcal_fcal_two_gamma_mass",100, 0.0, 0.300);
	bcal_fcal_two_gamma_mass_cut = new TH1F("bcal_fcal_two_gamma_mass_cut","bcal_fcal_two_gamma_mass_cut",100, 0.0, 0.300);
	xy_shower = new TH2F("bcal_xy_shower","xy_shower",100, -100.0, 100., 100 , -100.0, 100.0);
	z_shower = new TH1F("bcal_z_shower","z_shower",450, -50.0, 400);
	E_shower = new TH1F("bcal_E_shower","E_shower", 200, 0.0, 6.0);
	
	E_over_Erec_vs_z = new TH2F("bcal_E_over_Erec_vs_z","E_over_Erec_vs_z", 200, -50.0, 600.0, 200, 0.0, 4.0);
	Ecorr_over_Erec_vs_z = new TH2F("bcal_Ecorr_over_Erec_vs_z","Ecorr_over_Erec_vs_z", 200, -50.0, 600.0, 200, 0.0, 4.0);
	Ereconstructed_vs_Ethrown = new TH2F("Ereconstructed_vs_Ethrown","BCAL total reconstructed E to total thrown E", 200, 0.0, 6.0, 200, 0.0, 6.0);

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
	vector<const DFCALShower*> fcal_showers;
	vector<const DBCALTruthShower*> truthhits;	
	loop->Get(showers, "SIMPLE");
	//loop->Get(fcal_showers);
	loop->Get(truthhits);
	
	double Etot = 0.0;
	double z = 0.0;
	for(unsigned int i=0; i<truthhits.size(); i++){
		Etot += truthhits[i]->E;
		z += truthhits[i]->E*truthhits[i]->z;
	}
	z/=Etot;
	
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
			const DFCALShower *s2 = fcal_showers[j];
			dx = s2->x;
			dy = s2->y;
			dz = FCAL_Z_OFFSET; // shift to coordinate relative to center of target
			R = sqrt(dx*dx + dy*dy + dz*dz);
			double E = s2->E;
			TLorentzVector p2(E*dx/R, E*dy/R, E*dz/R, E);		
			
			TLorentzVector ptot = p1dave+p2;
			bcal_fcal_two_gamma_mass->Fill(ptot.M());

			if(showers.size()==1 && fcal_showers.size()==1){
				bcal_fcal_two_gamma_mass_cut->Fill(ptot.M());
			}
		}
		
	}
	
	// Compare to thrown values
	vector<const DMCThrown*> mcthrowns;
	loop->Get(mcthrowns);
	double Etot_thrown=0.0;
	for(unsigned int i=0; i<mcthrowns.size(); i++){
		Etot_thrown += mcthrowns[i]->E;
		for(unsigned int j=0; j<showers.size(); j++){
			double z = showers[j]->z;
			E_over_Erec_vs_z->Fill(z, mcthrowns[i]->E/showers[j]->Ecorr);

			double Ecorr = showers[j]->Ecorr*(1.106+(z-208.4)*(z-208.4)*6.851E-6);
			Ecorr_over_Erec_vs_z->Fill(z, mcthrowns[i]->E/Ecorr);
		}
	}
	
	Ereconstructed_vs_Ethrown->Fill(Etot_thrown, Etot_reconstructed);

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

