
#include <iostream>
#include <cmath>
using namespace std;

#include <TFile.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <TH2.h>

bool SmearCharged(TVector3 &mom, int geanttype);
bool EfficiencyCharged(const TVector3 &mom, int geanttype);
void InitializeResolutionTable(void);
void InitializeEfficiencyTable(void);

// Globals
bool initialized_resolution_table=false;
bool initialized_efficiency_table=false;
TRandom3 rnd;
TH2D* res_charged;
TH2D* eff_charged;

//----------------
// SmearCharged
//----------------
bool SmearCharged(TVector3 &mom, int geanttype)
{
	/// Smear the momentum vector of a charged particle
	/// based on reolutions from GEANT-based Monte Carlo
	/// studies.
	///
	/// The value of geanttype should specify the particle
	/// type using the GEANT particle ids (a few are given
	/// below). At this point, the geanttype is ignored
	/// and is merely here as a place holder.

	// Make sure resolution table is initialized
	if(!initialized_resolution_table)InitializeResolutionTable();
	
	// Base efficiency on input values
	bool reconstructed = EfficiencyCharged(mom, geanttype);

	// Find bin for this momentum
	double p = mom.Mag();
	double theta = mom.Theta()*57.3;
	int pbin = res_charged->GetYaxis()->FindBin(p);
	int thetabin = res_charged->GetXaxis()->FindBin(theta);
	
	// Here we should do an interpolation from the surrounding bins.
	// We have fairly small bins though so I can afford to be
	// lazy :)
	double sigma_dp_over_p = res_charged->GetBinContent(thetabin, pbin)/100.0;

	// Calculate new momentum
	double p_new = p*(1.0 + rnd.Gaus(0.0, sigma_dp_over_p));
	
	// Overwrite input vector with new values
	mom.SetMag(p_new);

	return reconstructed;
}

//----------------
// EfficiencyCharged
//----------------
bool EfficiencyCharged(const TVector3 &mom, int geanttype)
{
	/// Find the efficiency of a charged particle of the 
	/// given momentum and type and return whether or not
	/// the particle would be detected and reconstructed.
	/// The output of this has the efficiency folded in
	/// so when this is called with a momentum vector
	/// which we are 90% efficient at reconstructing, this
	/// will return "true" 90% of the times and "false"
	/// 10% of the times it is called.
	///
	/// The value of geanttype should specify the particle
	/// type using the GEANT particle ids (a few are given
	/// below). At this point, the geanttype is ignored
	/// and is merely here as a place holder.

	// Make sure efficiency table is initialized
	if(!initialized_efficiency_table)InitializeEfficiencyTable();

	// Find bin for this momentum
	double p = mom.Mag();
	double theta = mom.Theta()*57.3;
	int pbin = eff_charged->GetYaxis()->FindBin(p);
	int thetabin = eff_charged->GetXaxis()->FindBin(theta);
	
	// Here we should do an interpolation from the surrounding bins.
	// We have fairly small bins though so I can afford to be
	// lazy :)
	double eff = eff_charged->GetBinContent(thetabin, pbin)/100.0;

	return rnd.Rndm()<=eff;
}

//----------------
// InitializeResolutionTable
//----------------
void InitializeResolutionTable(void)
{
	/// Read in the resolution table for ROOT file
	
	// Open ROOT file
	TFile *f = new TFile("hd_res_charged.root");
	if(!f->IsOpen()){
		cout<<endl;
		cout<<"Couldn't open resolution file \"hd_res_charged.root\"!"<<endl;
		cout<<"Make sure it exists in the current directory and is readable,"<<endl;
		cout<<endl;
		exit(0);
	}
	
	// Read in resolution histogram
	res_charged = (TH2D*)gROOT->FindObject("dpt_over_pt_vs_p_vs_theta_2");
	if(!res_charged){
		cout<<endl;
		cout<<"Couldn't find resolution histogram \"dpt_over_pt_vs_p_vs_theta_2\""<<endl;
		cout<<"in ROOT file!"<<endl;
		cout<<endl;
		exit(0);
	}
	
	initialized_resolution_table=true;
}

//----------------
// InitializeEfficiencyTable
//----------------
void InitializeEfficiencyTable(void)
{
	/// Read in the efficiency table for ROOT file
	
	// Open ROOT file
	TFile *f = new TFile("hd_res_charged.root");
	if(!f->IsOpen()){
		cout<<endl;
		cout<<"Couldn't open resolution file \"hd_res_charged.root\"!"<<endl;
		cout<<"Make sure it exists in the current directory and is readable,"<<endl;
		cout<<endl;
		exit(0);
	}
	
	// Read in resolution histogram
	eff_charged = (TH2D*)gROOT->FindObject("eff_vs_p_vs_theta");
	if(!res_charged){
		cout<<endl;
		cout<<"Couldn't find efficiency histogram \"eff_vs_p_vs_theta\""<<endl;
		cout<<"in ROOT file!"<<endl;
		cout<<endl;
		exit(0);
	}
	
	initialized_efficiency_table=true;
}

