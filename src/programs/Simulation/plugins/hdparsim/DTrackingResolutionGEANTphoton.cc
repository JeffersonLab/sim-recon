// $Id$
//
//    File: DTrackingResolutionGEANT.cc
// Created: Mon Feb 25 15:06:17 EST 2008
// Creator: davidl (on Darwin fwing-dhcp13.jlab.org 8.11.1 i386)
//

#include <TROOT.h>
#include <TApplication.h>

#include <iostream>
#include <cmath>
using namespace std;

#include <JANA/jerror.h>

#include "DTrackingResolutionGEANTphoton.h"
#include "getwebfile.h"

#define rad2deg (180.0/M_PI)

//---------------------------------
// DTrackingResolutionGEANT    (Constructor)
//---------------------------------
DTrackingResolutionGEANTphoton::DTrackingResolutionGEANTphoton()
{
	//int argc=0;
	//TApplication *app = new TApplication("myapp", &argc, NULL);

	// Get ROOT file from web if it is not already here
	const char *url = "http://www.jlab.org/Hall-D/datatables/hd_res_photon.root";
	getwebfile(url);
	
	TDirectory *savedir = gDirectory;

	//---------------- hd_res_photon ------------------
	// Open ROOT file
	file = new TFile("hd_res_photon.root");
	if(!file->IsOpen()){
		cout<<endl;
		cout<<"Couldn't open resolution file \"hd_res_photon.root\"!"<<endl;
		cout<<"Make sure it exists in the current directory and is readable,"<<endl;
		cout<<endl;
		exit(0);
	}
	//cout<<"Opened \""<<file->GetName()<<"\""<<endl;

	// Read pt resolution histogram
	E_res_hist = (TH2D*)gROOT->FindObject("dE_over_E_vs_p_vs_theta");
	if(!E_res_hist)file->GetObject("dE_over_E_vs_p_vs_theta", E_res_hist);
	if(!E_res_hist){
		cout<<endl;
		cout<<"Couldn't find resolution histogram \"dpt_over_pt_vs_p_vs_theta\""<<endl;
		cout<<"in ROOT file!"<<endl;
		cout<<endl;
		exit(0);
	}

	// Read theta resolution histogram
	theta_res_hist = (TH2D*)gROOT->FindObject("dtheta_vs_p_vs_theta");
	if(!theta_res_hist)file->GetObject("dtheta_vs_p_vs_theta", theta_res_hist);
	if(!theta_res_hist){
		cout<<endl;
		cout<<"Couldn't find resolution histogram \"dtheta_vs_p_vs_theta\""<<endl;
		cout<<"in ROOT file!"<<endl;
		cout<<endl;
		exit(0);
	}

	// Read phi resolution histogram
	phi_res_hist = (TH2D*)gROOT->FindObject("dphi_vs_p_vs_theta");
	if(!phi_res_hist)file->GetObject("dphi_vs_p_vs_theta", phi_res_hist);
	if(!phi_res_hist){
		cout<<endl;
		cout<<"Couldn't find resolution histogram \"dphi_vs_p_vs_theta\""<<endl;
		cout<<"in ROOT file!"<<endl;
		cout<<endl;
		exit(0);
	}

	// Read in efficiency histogram
	efficiency_hist = (TH2D*)gROOT->FindObject("eff_vs_p_vs_theta");
	if(!efficiency_hist)file->GetObject("eff_vs_p_vs_theta", efficiency_hist);
	if(!efficiency_hist){
		cout<<endl;
		cout<<"Couldn't find efficiency histogram \"eff_vs_p_vs_theta\""<<endl;
		cout<<"in ROOT file!"<<endl;
		cout<<endl;
		exit(0);
	}
	
	if(savedir)savedir->cd();
}

//---------------------------------
// ~DTrackingResolutionGEANT    (Destructor)
//---------------------------------
DTrackingResolutionGEANTphoton::~DTrackingResolutionGEANTphoton()
{
	if(file)delete file;
}


//----------------
// GetResolution
//----------------
void DTrackingResolutionGEANTphoton::GetResolution(int geanttype, const TVector3 &mom, double &E_res, double &theta_res, double &phi_res)
{
	/// Return the energy and angular resolutions for a charged
	/// particle based on results from GEANT-based Monte Carlo studies.

	// Find bins for this momentum.
	// Note that we assume the 3 histograms have the same format.
	// Namely, number of bins and range so we only need to calculate
	// the theta and p bins once.
	double p = mom.Mag();
	double theta = mom.Theta()*rad2deg;
	int pbin = E_res_hist->GetYaxis()->FindBin(p);
	int thetabin = E_res_hist->GetXaxis()->FindBin(theta);
	
	if(pbin<1 || pbin>E_res_hist->GetNbinsY()){E_res=theta_res=phi_res=0.0; return;}
	if(thetabin<1 || thetabin>E_res_hist->GetNbinsX()){E_res=theta_res=phi_res=0.0; return;}
	
	// Here we should do an interpolation from the surrounding bins.
	// We have fairly small bins though so I can afford to be
	// lazy :)
	E_res = E_res_hist->GetBinContent(thetabin, pbin); // return as fraction
	theta_res = theta_res_hist->GetBinContent(thetabin, pbin); // return in milliradians
	phi_res = phi_res_hist->GetBinContent(thetabin, pbin); // return in milliradians
}

//----------------
// GetEfficiency
//----------------
double DTrackingResolutionGEANTphoton::GetEfficiency(int geanttype, const TVector3 &mom)
{
	/// Return the reconstruction efficiency for a charged
	/// particle based on results from GEANT-based Monte Carlo studies.

	// Find bins for this momentum.
	double p = mom.Mag();
	double theta = mom.Theta()*rad2deg;
	int pbin = efficiency_hist->GetYaxis()->FindBin(p);
	int thetabin = efficiency_hist->GetXaxis()->FindBin(theta);
	
	if(pbin<1 || pbin>efficiency_hist->GetNbinsY())return 0.0;
	if(thetabin<1 || thetabin>efficiency_hist->GetNbinsX())return 0.0;
	
	// Here we should do an interpolation from the surrounding bins.
	// We have fairly small bins though so I can afford to be
	// lazy :)
	return efficiency_hist->GetBinContent(thetabin, pbin);
}


