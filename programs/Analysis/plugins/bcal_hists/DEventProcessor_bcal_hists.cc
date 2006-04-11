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
// Routine used to create our DEventProcessor
extern "C"{
void InitProcessors(DApplication *app){
	app->AddProcessor(new DEventProcessor_bcal_hists());
}
} // "C"


//------------------
// init
//------------------
derror_t DEventProcessor_bcal_hists::init(void)
{
	// open ROOT file
	ROOTfile = new TFile("bcal_hists.root","RECREATE","Produced by hd_ana");
	cout<<"Opened ROOT file \"bcal_hists.root\""<<endl;
	
	two_gamma_mass = new TH1F("two_gamma_mass","two_gamma_mass",100, 0.0, 0.300);
	xy_shower = new TH2F("xy_shower","xy_shower",100, -100.0, 100., 100 , -100.0, 100.0);
	z_shower = new TH1F("z_shower","z_shower",450, -50.0, 400);

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
	if(showers.size() < 2)return NOERROR;
	
	for(unsigned int i=0; i<showers.size()-1; i++){
		const DBCALShower *s = showers[i];
		xy_shower->Fill(s->x, s->y);
		z_shower->Fill(s->z);
	}
	
	for(unsigned int i=0; i<showers.size()-1; i++){
		const DBCALShower *s1 = showers[i];
		double dx = s1->x;
		double dy = s1->y;
		double dz = s1->z - 65.0;
		double R = sqrt(dx*dx + dy*dy + dz*dz);
		double E = s1->Ecorr;
		TLorentzVector p1(E*dx/R, E*dy/R, E*dz/R, E);
		for(unsigned int j=i+1; j<showers.size(); j++){
			const DBCALShower *s2 = showers[j];
			dx = s2->x;
			dy = s2->y;
			dz = s2->z - 65.0;
			R = sqrt(dx*dx + dy*dy + dz*dz);
			E = s2->Ecorr;
			TLorentzVector p2(E*dx/R, E*dy/R, E*dz/R, E);
			
			TLorentzVector ptot = p1+p2;
			two_gamma_mass->Fill(ptot.M());
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

