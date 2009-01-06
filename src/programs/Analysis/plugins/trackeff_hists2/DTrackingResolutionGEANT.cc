// $Id$
//
//    File: DTrackingResolutionGEANT.cc
// Created: Mon Feb 25 15:06:17 EST 2008
// Creator: davidl (on Darwin fwing-dhcp13.jlab.org 8.11.1 i386)
//

#include <TApplication.h>
#include <TROOT.h>

#include <iostream>
#include <cstdlib>
using namespace std;

#include "DTrackingResolutionGEANT.h"

#define _DBG_ cout<<__FILE__<<":"<<__LINE__<<" "


//---------------------------------
// DTrackingResolutionGEANT    (Constructor)
//---------------------------------
DTrackingResolutionGEANT::DTrackingResolutionGEANT()
{
	
	// Make sure a TApplication exists, otherwise, the following will fail
	if(!gApplication){
		int argc=0;
		new TApplication("myapp", &argc, NULL);
	}

	// Open ROOT file
	file = new TFile("hd_res_charged.root");
	if(!file->IsOpen()){
		cout<<endl;
		cout<<"Couldn't open resolution file \"hd_res_charged.root\"!"<<endl;
		cout<<"Make sure it exists in the current directory and is readable,"<<endl;
		cout<<endl;
		exit(0);
	}
	cout<<"Opened \""<<file->GetName()<<"\""<<endl;

	// Read pt resolution histogram
	file->GetObject("dpt_over_pt_vs_p_vs_theta_2", pt_res_hist);
	if(!pt_res_hist){
		cout<<endl;
		cout<<"Couldn't find resolution histogram \"dpt_over_pt_vs_p_vs_theta_2\""<<endl;
		cout<<"in ROOT file!"<<endl;
		cout<<endl;
		exit(0);
	}

	// Read theta resolution histogram
	file->GetObject("dtheta_vs_p_vs_theta_2", theta_res_hist);
	if(!theta_res_hist){
		cout<<endl;
		cout<<"Couldn't find resolution histogram \"dtheta_vs_p_vs_theta_2\""<<endl;
		cout<<"in ROOT file!"<<endl;
		cout<<endl;
		exit(0);
	}

	// Read phi resolution histogram
	phi_res_hist = (TH2D*)gROOT->FindObject("dphi_vs_p_vs_theta_2");
	if(!phi_res_hist){
		cout<<endl;
		cout<<"Couldn't find resolution histogram \"dphi_vs_p_vs_theta_2\""<<endl;
		cout<<"in ROOT file!"<<endl;
		cout<<endl;
		exit(0);
	}

	// Read in efficiency histogram
	efficiency_hist = (TH2D*)gROOT->FindObject("eff_vs_p_vs_theta");
	if(!efficiency_hist){
		cout<<endl;
		cout<<"Couldn't find efficiency histogram \"eff_vs_p_vs_theta\""<<endl;
		cout<<"in ROOT file!"<<endl;
		cout<<endl;
		exit(0);
	}
}

//---------------------------------
// ~DTrackingResolutionGEANT    (Destructor)
//---------------------------------
DTrackingResolutionGEANT::~DTrackingResolutionGEANT()
{
	delete file;
}

//----------------
// GetResolution
//----------------
void DTrackingResolutionGEANT::GetResolution(int geanttype, const TVector3 &mom, double &pt_res, double &theta_res, double &phi_res)
{
	/// Return the momentum and angular resolutions for a charged
	/// particle based on results from GEANT-based Monte Carlo studies.

	// Find bins for this momentum.
	// Note that we assume the 3 histograms have the same format.
	// Namely, number of bins and range so we only need to calculate
	// the theta and p bins once.
	double p = mom.Mag();
	double theta = mom.Theta()*57.3;
	int pbin = pt_res_hist->GetYaxis()->FindBin(p);
	int thetabin = pt_res_hist->GetXaxis()->FindBin(theta);

//_DBG_<<"pbin="<<pbin<<" thetabin="<<thetabin<<endl;
	if(pbin<1 || pbin>pt_res_hist->GetNbinsY()){pt_res=theta_res=phi_res=0.0; return;}
	if(thetabin<1 || thetabin>pt_res_hist->GetNbinsX()){pt_res=theta_res=phi_res=0.0; return;}

	// Here we should do an interpolation from the surrounding bins.
	// We have fairly small bins though so I can afford to be
	// lazy :)
	pt_res = pt_res_hist->GetBinContent(thetabin, pbin)/100.0;
	theta_res = theta_res_hist->GetBinContent(thetabin, pbin);
	phi_res = phi_res_hist->GetBinContent(thetabin, pbin);
}

//----------------
// GetEfficiency
//----------------
double DTrackingResolutionGEANT::GetEfficiency(int geanttype, const TVector3 &mom)
{
	/// Return the reconstruction efficiency for a charged
	/// particle based on results from GEANT-based Monte Carlo studies.

	// Find bins for this momentum.
	double p = mom.Mag();
	double theta = mom.Theta()*57.3;
	int pbin = efficiency_hist->GetYaxis()->FindBin(p);
	int thetabin = efficiency_hist->GetXaxis()->FindBin(theta);
	
	if(pbin<1 || pbin>efficiency_hist->GetNbinsY())return 0.0;
	if(thetabin<1 || thetabin>efficiency_hist->GetNbinsX())return 0.0;
	
	// Here we should do an interpolation from the surrounding bins.
	// We have fairly small bins though so I can afford to be
	// lazy :)
	return efficiency_hist->GetBinContent(thetabin, pbin);
}



