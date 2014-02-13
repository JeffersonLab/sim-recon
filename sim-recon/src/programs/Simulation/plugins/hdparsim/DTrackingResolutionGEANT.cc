// $Id$
//
//    File: DTrackingResolutionGEANT.cc
// Created: Mon Feb 25 15:06:17 EST 2008
// Creator: davidl (on Darwin fwing-dhcp13.jlab.org 8.11.1 i386)
//

#include <TROOT.h>
#include <TApplication.h>

#include <iostream>
using namespace std;

#include "DTrackingResolutionGEANT.h"
#include "getwebfile.h"


//---------------------------------
// DTrackingResolutionGEANT    (Constructor)
//---------------------------------
DTrackingResolutionGEANT::DTrackingResolutionGEANT()
{
	//int argc=0;
	//TApplication *app = new TApplication("myapp", &argc, NULL);

	
	TDirectory *savedir = gDirectory;
	
	ReadTableInfo("hd_res_charged_pion.root", pion_info);
	ReadTableInfo("hd_res_charged_proton.root", proton_info);

	if(savedir)savedir->cd();

}

//----------------
// ReadTableInfo
//----------------
void DTrackingResolutionGEANT::ReadTableInfo(const char *fname, TableInfo &ti)
{
	// Get ROOT file from web if it is not already here
	char url[512] = "http://www.jlab.org/Hall-D/datatables/";
	strcat(url, fname);
	getwebfile(url);

	// Open ROOT file
	ti.file = new TFile(fname);
	if(!ti.file->IsOpen()){
		cout<<endl;
		cout<<"Couldn't open resolution file \"hd_res_charged.root\"!"<<endl;
		cout<<"Make sure it exists in the current directory and is readable,"<<endl;
		cout<<endl;
		exit(0);
	}
	//cout<<"Opened \""<<ti.file->GetName()<<"\""<<endl;
	
	// An earlier version used a slightly different naming scheme. We check for the
	// new name first, but if that fails, then look for the old name in case they 
	// are using an older version of the data tables.  Aug. 11, 2009  DL

	// Read pt resolution histogram
	ti.pt_res_hist = (TH2D*)gROOT->FindObject("dpt_over_pt_sigma");
	if(!ti.pt_res_hist)ti.file->GetObject("dpt_over_pt_sigma", ti.pt_res_hist);
	if(!ti.pt_res_hist)ti.pt_res_hist = (TH2D*)gROOT->FindObject("dpt_over_pt_vs_p_vs_theta");
	if(!ti.pt_res_hist)ti.file->GetObject("dpt_over_pt_vs_p_vs_theta", ti.pt_res_hist);
	if(!ti.pt_res_hist){
		cout<<endl;
		cout<<"Couldn't find resolution histogram \"dpt_over_pt_sigma\""<<endl;
		cout<<"in ROOT file!"<<endl;
		cout<<endl;
		exit(0);
	}

	// Read theta resolution histogram
	ti.theta_res_hist = (TH2D*)gROOT->FindObject("dtheta_sigma");
	if(!ti.theta_res_hist)ti.file->GetObject("dtheta_sigma", ti.theta_res_hist);
	if(!ti.theta_res_hist)ti.theta_res_hist = (TH2D*)gROOT->FindObject("dtheta_vs_p_vs_theta");
	if(!ti.theta_res_hist)ti.file->GetObject("dtheta_vs_p_vs_theta", ti.theta_res_hist);
	if(!ti.theta_res_hist){
		cout<<endl;
		cout<<"Couldn't find resolution histogram \"dtheta_sigma\""<<endl;
		cout<<"in ROOT file!"<<endl;
		cout<<endl;
		exit(0);
	}

	// Read phi resolution histogram
	ti.phi_res_hist = (TH2D*)gROOT->FindObject("dphi_sigma");
	if(!ti.phi_res_hist)ti.file->GetObject("dphi_sigma", ti.phi_res_hist);
	if(!ti.phi_res_hist)ti.phi_res_hist = (TH2D*)gROOT->FindObject("dphi_vs_p_vs_theta");
	if(!ti.phi_res_hist)ti.file->GetObject("dphi_vs_p_vs_theta", ti.phi_res_hist);
	if(!ti.phi_res_hist){
		cout<<endl;
		cout<<"Couldn't find resolution histogram \"dphi_sigma\""<<endl;
		cout<<"in ROOT file!"<<endl;
		cout<<endl;
		exit(0);
	}

	// Read in efficiency histogram
	ti.efficiency_hist = (TH2D*)gROOT->FindObject("eff_vs_p_vs_theta");
	if(!ti.efficiency_hist)ti.file->GetObject("eff_vs_p_vs_theta", ti.efficiency_hist);
	if(!ti.efficiency_hist){
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
	if(pion_info.file)delete pion_info.file;
	if(proton_info.file)delete proton_info.file;
}

//----------------
// GetResolution
//----------------
void DTrackingResolutionGEANT::GetResolution(int geanttype, const TVector3 &mom, double &pt_res, double &theta_res, double &phi_res)
{
	switch(geanttype){
		case 14:
			GetResolution(proton_info, geanttype, mom, pt_res, theta_res, phi_res);
			break;
		case 8: // pi+
		case 9: // pi-
		default: // assume everything else is close to pion resolutions
			GetResolution(pion_info, geanttype, mom, pt_res, theta_res, phi_res);
	}
}

//----------------
// GetEfficiency
//----------------
double DTrackingResolutionGEANT::GetEfficiency(int geanttype, const TVector3 &mom)
{
	switch(geanttype){
		case 14:
			return GetEfficiency(proton_info, geanttype, mom);
			break;
		case 8: // pi+
		case 9: // pi-
		default: // assume everything else is close to pion resolutions
			return GetEfficiency(pion_info, geanttype, mom);
	}
}

//----------------
// GetResolution
//----------------
void DTrackingResolutionGEANT::GetResolution(TableInfo &ti, int geanttype, const TVector3 &mom, double &pt_res, double &theta_res, double &phi_res)
{
	/// Return the momentum and angular resolutions for a charged
	/// particle based on results from GEANT-based Monte Carlo studies.

	// Find bins for this momentum.
	// Note that we assume the 3 histograms have the same format.
	// Namely, number of bins and range so we only need to calculate
	// the theta and p bins once.
	double p = mom.Mag();
	double theta = mom.Theta()*57.3;
	int pbin = ti.pt_res_hist->GetYaxis()->FindBin(p);
	int thetabin = ti.pt_res_hist->GetXaxis()->FindBin(theta);
	
	// For tracks with momentum out of the range of our table, use the
	// resolutions for the largest momentum we have
	if(pbin>ti.pt_res_hist->GetNbinsY())pbin=ti.pt_res_hist->GetNbinsY();
	
	if(pbin<1 || pbin>ti.pt_res_hist->GetNbinsY()){pt_res=theta_res=phi_res=0.0; return;}
	if(thetabin<1 || thetabin>ti.pt_res_hist->GetNbinsX()){pt_res=theta_res=phi_res=0.0; return;}
	
	// Here we should do an interpolation from the surrounding bins.
	// We have fairly small bins though so I can afford to be
	// lazy :)
	pt_res = ti.pt_res_hist->GetBinContent(thetabin, pbin); // return as fraction
	theta_res = ti.theta_res_hist->GetBinContent(thetabin, pbin); // return in milliradians
	phi_res = ti.phi_res_hist->GetBinContent(thetabin, pbin); // return in milliradians
}

//----------------
// GetEfficiency
//----------------
double DTrackingResolutionGEANT::GetEfficiency(TableInfo &ti, int geanttype, const TVector3 &mom)
{
	/// Return the reconstruction efficiency for a charged
	/// particle based on results from GEANT-based Monte Carlo studies.

	// Find bins for this momentum.
	double p = mom.Mag();
	double theta = mom.Theta()*57.3;
	int pbin = ti.efficiency_hist->GetYaxis()->FindBin(p);
	int thetabin = ti.efficiency_hist->GetXaxis()->FindBin(theta);
	
	// For tracks with momentum out of the range of our table, use the
	// resolutions for the largest momentum we have
	if(pbin>ti.efficiency_hist->GetNbinsY())pbin=ti.efficiency_hist->GetNbinsY();
	
	if(pbin<1 || pbin>ti.efficiency_hist->GetNbinsY())return 0.0;
	if(thetabin<1 || thetabin>ti.efficiency_hist->GetNbinsX())return 0.0;
	
	// Here we should do an interpolation from the surrounding bins.
	// We have fairly small bins though so I can afford to be
	// lazy :)
	return ti.efficiency_hist->GetBinContent(thetabin, pbin);
}



