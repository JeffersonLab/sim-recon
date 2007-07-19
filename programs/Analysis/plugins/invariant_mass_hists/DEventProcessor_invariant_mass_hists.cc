// $Id: DEventProcessor_invariant_mass_hists.cc 1816 2006-06-06 14:38:18Z davidl $
//
//    File: DEventProcessor_invariant_mass_hists.cc
// Created: Thur Jan 11 15:42:21 EDT 2005
// Creator: davidl 
//

#include <iostream>
#include <cmath>
using namespace std;

#include <TThread.h>
#include <TLorentzVector.h>

#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>

#include "DEventProcessor_invariant_mass_hists.h"
#include "TRACKING/DTrack.h"
#include "TRACKING/DMCThrown.h"


// Routine used to create our JEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_invariant_mass_hists());
}
}

//------------------
// DEventProcessor_invariant_mass_hists
//------------------
DEventProcessor_invariant_mass_hists::DEventProcessor_invariant_mass_hists()
{	
}

//------------------
// ~DEventProcessor_invariant_mass_hists
//------------------
DEventProcessor_invariant_mass_hists::~DEventProcessor_invariant_mass_hists()
{
}

//------------------
// init
//------------------
jerror_t DEventProcessor_invariant_mass_hists::init(void)
{
	// open ROOT file (if needed)
	//if(ROOTfile != NULL) ROOTfile->cd();

	// Create THROWN directory
	TDirectory *dir = new TDirectory("INV_MASS","INV_MASS");
	dir->cd();

	// Create histograms
	mass_1part = new TH1F("mass_1part","Invariant mass of 1 particle GeV/c^{2}",4100, -0.1, 4.0);
	mass_2part = new TH1F("mass_2part","Invariant mass of 2 particles GeV/c^{2}",4100, -0.1, 4.0);
	mass_3part = new TH1F("mass_3part","Invariant mass of 3 particles GeV/c^{2}",4100, -0.1, 4.0);
	mass2_1part = new TH1F("mass2_1part","Invariant mass^{2} of 1 particle GeV/c^{2}",4100, -0.1, 16.0);
	mass2_2part = new TH1F("mass2_2part","Invariant mass^{2} of 2 particles GeV/c^{2}",4100, -0.1, 16.0);
	mass2_3part = new TH1F("mass2_3part","Invariant mass^{2} of 3 particles GeV/c^{2}",4100, -0.1, 16.0);

	thrown_mass_1part = new TH1F("thrown_mass_1part","Invariant mass of 1 thrown particle GeV/c^{2}",4100, -0.1, 4.0);
	thrown_mass_2part = new TH1F("thrown_mass_2part","Invariant mass of 2 thrown particles GeV/c^{2}",4100, -0.1, 4.0);
	thrown_mass_3part = new TH1F("thrown_mass_3part","Invariant mass of 3 thrown particles GeV/c^{2}",4100, -0.1, 4.0);
	thrown_mass2_1part = new TH1F("thrown_mass2_1part","Invariant mass^{2} of 1 thrown particle GeV/c^{2}",4100, -0.1, 16.0);
	thrown_mass2_2part = new TH1F("thrown_mass2_2part","Invariant mass^{2} of 2 thrown particles GeV/c^{2}",4100, -0.1, 16.0);
	thrown_mass2_3part = new TH1F("thrown_mass2_3part","Invariant mass^{2} of 3 thrown particles GeV/c^{2}",4100, -0.1, 16.0);
	
	missing_mass = new TH1F("missing_mass", "Missing mass squared", 200, -1.0, 3.0);
	thrown_missing_mass = new TH1F("thrown_missing_mass", "Missing mass squared from throwns", 200, -1.0, 3.0);
	
	// Go back up to the parent directory
	dir->cd("../");
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_invariant_mass_hists::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DTrack*> tracks;
	vector<const DMCThrown*> mcthrowns;
	loop->Get(tracks);
	loop->Get(mcthrowns);

	double Egamma = 9.0;
	TLorentzVector incident(0.0, 0.0, Egamma, Egamma+0.938);

	// Loop over reconstructed tracks
	vector<TLorentzVector> v;
	for(unsigned int i=0;i<tracks.size();i++){
		TLorentzVector myv;
		if(tracks[i]->chisq>1.0)continue;
		MakeTLorentz(tracks[i], myv);
		v.push_back(myv);
	}

	// Loop over reconstructed tracks
	for(unsigned int i=0;i<v.size();i++){
		TLorentzVector v1 = v[i];
		mass_1part->Fill(v1.M());
		mass2_1part->Fill(v1.M2());
		
		for(unsigned int j=i+1;j<v.size();j++){
			TLorentzVector v2 = v1 + v[j];
			if(tracks[i]->q*tracks[j]->q<0.0){
				// Only fill for opposite charge tracks	
				mass_2part->Fill(v2.M());
				mass2_2part->Fill(v2.M2());
			}

			TLorentzVector mm = incident - v2;
			missing_mass->Fill(mm.M2());
			
			for(unsigned int k=j+1;k<v.size();k++){
				TLorentzVector v3 = v2 + v[k];		
				mass_3part->Fill(v3.M());
				mass2_3part->Fill(v3.M2());
			}
		}
	}

	// Repeat for thrown tracks
	v.clear();
	for(unsigned int i=0;i<mcthrowns.size();i++){
		int type=mcthrowns[i]->type;
		if(type!=8 && type!=9)continue;
		TLorentzVector myv;
		MakeTLorentz(mcthrowns[i], myv);
		v.push_back(myv);
	}

	// Loop over thrown tracks
	for(unsigned int i=0;i<v.size();i++){
		TLorentzVector v1 = v[i];
		thrown_mass_1part->Fill(v1.M());
		thrown_mass2_1part->Fill(v1.M2());

		for(unsigned int j=i+1;j<v.size();j++){
			TLorentzVector v2 = v1 + v[j];		
			if(tracks[i]->q*tracks[j]->q<0.0){
				// Only fill for opposite charge tracks	
				thrown_mass_2part->Fill(v2.M());
				thrown_mass2_2part->Fill(v2.M2());
			}
			
			TLorentzVector mm = incident - v2;
			thrown_missing_mass->Fill(mm.M2());
			
			for(unsigned int k=j+1;k<v.size();k++){
				TLorentzVector v3 = v2 + v[k];		
				thrown_mass_3part->Fill(v3.M());
				thrown_mass2_3part->Fill(v3.M2());
			}
		}
	}

	return NOERROR;
}

//------------------
// MakeTLorentz
//------------------
void DEventProcessor_invariant_mass_hists::MakeTLorentz(const DTrack *track, TLorentzVector &v)
{
	double px = track->p*sin(track->theta)*cos(track->phi);
	double py = track->p*sin(track->theta)*sin(track->phi);
	double pz = track->p*cos(track->theta);
	double E = sqrt(pow((double)0.13957,(double)2.0)+pow((double)track->p, (double)2.0));
	
	v.SetXYZT(px,py,pz,E);
}

//------------------
// MakeTLorentz
//------------------
void DEventProcessor_invariant_mass_hists::MakeTLorentz(const DMCThrown *thrown, TLorentzVector &v)
{
	double p = thrown->p;
	double theta = thrown->theta;
	double phi = thrown->phi;

#if 0	
	// smear values
	p *= (1.0 + SampleGaussian(0.02));
	phi += SampleGaussian(0.1);
	if(phi<0.0)phi+=2.0*M_PI;
	if(phi>=2.0*M_PI)phi-=2.0*M_PI;
	theta += SampleGaussian(0.0325);
	if(theta<0.0)theta=-theta;
	if(theta>=M_PI)theta=2.0*M_PI-theta;
#endif

	double px = p*sin(theta)*cos(phi);
	double py = p*sin(theta)*sin(phi);
	double pz = p*cos(theta);
	double E = sqrt(pow(0.13957,2.0)+pow(p,2.0));
	
	v.SetXYZT(px,py,pz,E);
}

//------------------
// SampleGaussian
//------------------
double DEventProcessor_invariant_mass_hists::SampleGaussian(double sigma)
{
	// We loop to ensure not to return values greater than 3sigma away
	double val;
	do{
		double epsilon = 1.0E-10;
		double r1 = epsilon+((double)random()/(double)RAND_MAX);
		double r2 = (double)random()/(double)RAND_MAX;
		val = sigma*sqrt(-2.0*log(r1))*cos(2.0*M_PI*r2);
	}while(fabs(val/sigma) > 3.0);

	return val;
}


//------------------
// erun
//------------------
jerror_t DEventProcessor_invariant_mass_hists::erun(void)
{

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_invariant_mass_hists::fini(void)
{
	return NOERROR;
}
