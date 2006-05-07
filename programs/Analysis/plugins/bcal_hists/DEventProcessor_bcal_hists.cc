// $Id$
//
//    File: DEventProcessor_bcal_hists.cc
// Created: Mon Apr  3 11:38:03 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.4.0 powerpc)
//

#include "DEventProcessor_bcal_hists.h"

#include <TLorentzVector.h>

#include <DApplication.h>
#include <DBCALShower.h>
#include <DMCThrown.h>
// Routine used to create our DEventProcessor
extern "C"{
void InitProcessors(DApplication *app){
	app->AddProcessor(new DEventProcessor_bcal_hists());
}
} // "C"


#define BCAL_Z_OFFSET 26.03+123.4+65.0 //convert regina's coordinated to HDGeant's

//------------------
// init
//------------------
derror_t DEventProcessor_bcal_hists::init(void)
{
	// open ROOT file
	ROOTfile = new TFile("bcal_hists.root","RECREATE","Produced by hd_ana");
	cout<<"Opened ROOT file \"bcal_hists.root\""<<endl;
	
	two_gamma_mass = new TH1F("two_gamma_mass","two_gamma_mass",100, 0.0, 0.300);
	two_gamma_mass_corr = new TH1F("two_gamma_mass_corr","two_gamma_mass_corr",100, 0.0, 0.300);
	xy_shower = new TH2F("xy_shower","xy_shower",100, -100.0, 100., 100 , -100.0, 100.0);
	z_shower = new TH1F("z_shower","z_shower",450, -50.0, 400);
	E_shower = new TH1F("E_shower","E_shower", 200, 0.0, 6.0);
	
	E_over_Erec_vs_z = new TH2F("E_over_Erec_vs_z","E_over_Erec_vs_z", 200, -50.0, 600.0, 200, 0.0, 4.0);
	Ecorr_over_Erec_vs_z = new TH2F("Ecorr_over_Erec_vs_z","Ecorr_over_Erec_vs_z", 200, -50.0, 600.0, 200, 0.0, 4.0);

	return NOERROR;
}

//------------------
// brun
//------------------
derror_t DEventProcessor_bcal_hists::brun(DEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
derror_t DEventProcessor_bcal_hists::evnt(DEventLoop *loop, int eventnumber)
{
	vector<const DBCALShower*> showers;
	loop->Get(showers);
	
	// Single shower params
	for(unsigned int i=0; i<showers.size(); i++){
		const DBCALShower *s = showers[i];
		xy_shower->Fill(s->x, s->y);
		z_shower->Fill(s->z);
		E_shower->Fill(s->Ecorr);
	}
	
	// 2-gamma inv. mass
	for(unsigned int i=0; i<showers.size(); i++){
		const DBCALShower *s1 = showers[i];
		double dx = s1->x;
		double dy = s1->y;
		double dz = s1->z+BCAL_Z_OFFSET - 65.0;
		double R = sqrt(dx*dx + dy*dy + dz*dz);
		double E = s1->Ecorr;
		double Edave = s1->Ecorr*(1.106+(dz+65.0-208.4)*(dz+65.0-208.4)*6.851E-6);
		TLorentzVector p1(E*dx/R, E*dy/R, E*dz/R, E);		
		TLorentzVector p1dave(Edave*dx/R, Edave*dy/R, Edave*dz/R, Edave);		
		
		for(unsigned int j=i+1; j<showers.size(); j++){
			const DBCALShower *s2 = showers[j];
			dx = s2->x;
			dy = s2->y;
			dz = s2->z +BCAL_Z_OFFSET- 65.0; // shift to coordinate relative to center of target
			R = sqrt(dx*dx + dy*dy + dz*dz);
			double E = s1->Ecorr;
			double Edave = s2->Ecorr*(1.106+(dz+65.0-208.4)*(dz+65.0-208.4)*6.851E-6);
			TLorentzVector p2(E*dx/R, E*dy/R, E*dz/R, E);		
			TLorentzVector p2dave(Edave*dx/R, Edave*dy/R, Edave*dz/R, Edave);		
			
			TLorentzVector ptot = p1+p2;
			two_gamma_mass->Fill(ptot.M());
			TLorentzVector ptotdave = p1dave+p2dave;
			two_gamma_mass_corr->Fill(ptotdave.M());
		}
	}
	
	// Compare to thrown values
	vector<const DMCThrown*> mcthrowns;
	loop->Get(mcthrowns);
	for(unsigned int i=0; i<mcthrowns.size(); i++){
		for(unsigned int j=0; j<showers.size(); j++){
			double z = showers[j]->z+BCAL_Z_OFFSET;
			E_over_Erec_vs_z->Fill(z, mcthrowns[i]->E/showers[j]->Ecorr);

			double Ecorr = showers[j]->Ecorr*(1.106+(z-208.4)*(z-208.4)*6.851E-6);
			Ecorr_over_Erec_vs_z->Fill(z, mcthrowns[i]->E/Ecorr);
		}
	}

	return NOERROR;
}

//------------------
// erun
//------------------
derror_t DEventProcessor_bcal_hists::erun(void)
{
	// Any final calculations on histograms (like dividing them)
	// should be done here. This may get called more than once.
	return NOERROR;
}

//------------------
// fini
//------------------
derror_t DEventProcessor_bcal_hists::fini(void)
{

	ROOTfile->Write();
	delete ROOTfile;
	cout<<endl<<"Closed ROOT file"<<endl;

	return NOERROR;
}

