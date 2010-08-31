// $Id: DEventProcessor_radlen_hists.cc 1816 2006-06-06 14:38:18Z davidl $
//
//    File: DEventProcessor_radlen_hists.cc
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include <iostream>
using namespace std;

#include <TThread.h>

#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>
using namespace jana;

extern JApplication *japp;

#include <DANA/DApplication.h>

#include "DEventProcessor_radlen_hists.h"
#include "TRACKING/DMCTrajectoryPoint.h"
#include "TRACKING/DMCThrown.h"

// The executable should define the ROOTfile global variable. It will
// be automatically linked when dlopen is called.
//extern TFile *ROOTfile;

// Routine used to create our JEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_radlen_hists());
}
}

//------------------
// DEventProcessor_radlen_hists
//------------------
DEventProcessor_radlen_hists::DEventProcessor_radlen_hists()
{
	pthread_mutex_init(&mutex, NULL);
}

//------------------
// ~DEventProcessor_radlen_hists
//------------------
DEventProcessor_radlen_hists::~DEventProcessor_radlen_hists()
{
}

//------------------
// init
//------------------
jerror_t DEventProcessor_radlen_hists::init(void)
{
	// open ROOT file (if needed)
	//if(ROOTfile != NULL) ROOTfile->cd();

	pthread_mutex_lock(&mutex);

	// Create THROWN directory
	TDirectory *dir = new TDirectoryFile("RADLEN","RADLEN");
	dir->cd();

	// Create histograms
	dE_vs_r	= new TH2F("dE_vs_r","dE (keV) vs r", 300, 0.0, 65.0,1000, 0.0, 1000.0);
	dE_vs_z	= new TH2F("dE_vs_z","dE (keV) vs z", 650, 0.0, 650.0, 1000, 0.0, 1000.0);

	int Ntheta_bins = 480;
	double theta_min = 0.0;
	double theta_max = 120.0;
	theta_nevents = new TH1F("theta_nevents","Number of events per #theta bin in int_radlen_vs_z_vs_theta", Ntheta_bins, theta_min, theta_max);
	nXo_vs_z_vs_theta = new TH2F("nXo_vs_z_vs_theta","Radiation lengths vs. z and #theta", 650, 0.0, 650.0, Ntheta_bins, theta_min, theta_max);
	nXo_vs_z_vs_theta->SetXTitle("z-position along beamline (cm)");
	nXo_vs_z_vs_theta->SetYTitle("Thrown #theta angle (deg)");
	inXo_vs_z_vs_theta = new TH2F("inXo_vs_z_vs_theta","Integrated radiation lengths vs. z and #theta", 650, 0.0, 650.0, Ntheta_bins, theta_min, theta_max);
	inXo_vs_z_vs_theta->SetXTitle("z-position along beamline (cm)");
	inXo_vs_z_vs_theta->SetYTitle("Thrown #theta angle (deg)");

	nXo_vs_r_vs_theta = new TH2F("nXo_vs_r_vs_theta","Radiation lengths vs. r and #theta", 180, 0.0, 90.0, Ntheta_bins, theta_min, theta_max);
	nXo_vs_r_vs_theta->SetXTitle("Radial distance from beamline (cm)");
	nXo_vs_r_vs_theta->SetYTitle("Thrown #theta angle (deg)");
	inXo_vs_r_vs_theta = new TH2F("inXo_vs_r_vs_theta","Integrated radiation lengths vs. r and #theta", 180, 0.0, 90.0, Ntheta_bins, theta_min, theta_max);
	inXo_vs_r_vs_theta->SetXTitle("Radial distance from beamline (cm)");
	inXo_vs_r_vs_theta->SetYTitle("Thrown #theta angle (deg)");

	nXo_vs_z = new TH1F("nXo_vs_z","Radiation lengths vs. z", 650, 0.0, 650.0);
	inXo_vs_z = new TH1F("inXo_vs_z","Integrated radiation lengths vs. z", 650, 0.0, 650.0);
	nXo_vs_r = new TH1F("nXo_vs_r","Radiation lengths vs. r", 180, 0.0, 90.0);
	inXo_vs_r = new TH1F("inXo_vs_r","Integrated radiation lengths vs. r", 180, 0.0, 90.0);

	nXo_vs_z->SetStats(0);
	nXo_vs_z->SetFillStyle(3000);
	nXo_vs_z->SetFillColor(kMagenta);
	inXo_vs_z->SetStats(0);
	inXo_vs_z->SetFillStyle(3000);
	inXo_vs_z->SetFillColor(kMagenta);
	nXo_vs_r->SetStats(0);
	nXo_vs_r->SetFillStyle(3000);
	nXo_vs_r->SetFillColor(kMagenta);
	inXo_vs_r->SetStats(0);
	inXo_vs_r->SetFillStyle(3000);
	inXo_vs_r->SetFillColor(kMagenta);
	
	// Tree
	rstep_ptr = &rstep;
	tradstep = new TTree("radstep","Radlen steps");
	tradstep->Branch("R","radstep",&rstep_ptr);

	DApplication *dapp = dynamic_cast<DApplication*>(japp);
	bfield = dapp->GetBfield();

	// Go back up to the parent directory
	dir->cd("../");

	pthread_mutex_unlock(&mutex);
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_radlen_hists::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DMCTrajectoryPoint*> trajpoints;
	vector<const DMCThrown*> throwns;
	loop->Get(trajpoints);
	loop->Get(throwns);
	
	// Get theta of thrown particle
	double theta = 0.0;
	if(throwns.size()!=1){
		_DBG_<<"Event doesn't have exactly 1 thrown (has"<<throwns.size()<<"). Skipping..."<<endl;
		return NOERROR;
	}
	theta = throwns[0]->momentum().Theta()*57.3;
	theta_nevents->Fill(theta);
	
	// Loop over trajectory points
	pthread_mutex_lock(&mutex);
	rstep.stot = 0.0;
	rstep.ix_over_Xo = 0.0;
	rstep.iB_cross_p_dl = 0.0;
	rstep.iB_dl = 0.0;
	DVector3 tmp = throwns[0]->momentum();
	rstep.pthrown.SetXYZ(tmp.X(), tmp.Y(), tmp.Z());
	for(unsigned int i=0;i<trajpoints.size();i++){
		const DMCTrajectoryPoint *traj = trajpoints[i];
		
		double dE = traj->dE*1.0E6; // convert to keV
		double r = sqrt(pow((double)traj->x,2.0) + pow((double)traj->y,2.0));
		dE_vs_r->Fill( r, dE);
		dE_vs_z->Fill( traj->z, dE);
		
		double dnXo = traj->step/traj->radlen;
		nXo_vs_r_vs_theta->Fill(r, theta, dnXo);
		nXo_vs_z_vs_theta->Fill(traj->z, theta, dnXo);
		nXo_vs_r->Fill(r, dnXo);
		nXo_vs_z->Fill(traj->z, dnXo);
		
		double Bx, By, Bz;
		bfield->GetField(traj->x, traj->y, traj->z, Bx, By, Bz);
		rstep.B.SetXYZ(Bx, By, Bz);
		TVector3 mom(traj->px, traj->py, traj->pz);
		TVector3 B_cross_p = rstep.B.Cross(mom);
		rstep.iB_cross_p_dl += B_cross_p.Mag() * traj->step/mom.Mag();
		rstep.iB_dl += rstep.B.Mag()*traj->step;

		rstep.pos.SetXYZ(traj->x, traj->y, traj->z);
		rstep.s = traj->step;
		rstep.stot += (double)traj->step;
		rstep.Xo = traj->radlen;
		rstep.ix_over_Xo += dnXo;
		
		tradstep->Fill();
	}
	pthread_mutex_unlock(&mutex);

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_radlen_hists::erun(void)
{
	pthread_mutex_lock(&mutex);

	// Scale the nXo_vs_r_vs_theta and nXo_vs_z_vs_theta histos by the
	// number of events for each theta bin. At the same time, integrate
	// the number of radiation lengths and fill the
	// inXo_vs_r_vs_theta and inXo_vs_z_vs_theta histos.
	TH2F *h = nXo_vs_r_vs_theta;
	for(int i=1; i<=h->GetNbinsY(); i++){
		double N = theta_nevents->GetBinContent(i);
		double nXo_vs_r_sum=0.0;
		for(int j=1; j<=h->GetNbinsX(); j++){
			double v = nXo_vs_r_vs_theta->GetBinContent(j,i);
			double nXo = N==0.0 ? 0.0:v/N;
			nXo_vs_r_sum += nXo;
			nXo_vs_r_vs_theta->SetBinContent(j,i, nXo);
			inXo_vs_r_vs_theta->SetBinContent(j,i, nXo_vs_r_sum);
		}
	}

	h = nXo_vs_z_vs_theta;
	for(int i=1; i<=h->GetNbinsY(); i++){
		double N = theta_nevents->GetBinContent(i);
		double nXo_vs_z=0.0;
		for(int j=1; j<=h->GetNbinsX(); j++){
			double v = nXo_vs_z_vs_theta->GetBinContent(j,i);
			double nXo = N==0.0 ? 0.0:v/N;
			nXo_vs_z += nXo;
			nXo_vs_z_vs_theta->SetBinContent(j,i, nXo);
			inXo_vs_z_vs_theta->SetBinContent(j,i, nXo_vs_z);
		}
	}
	
	double N = theta_nevents->Integral();
_DBG_<<"N="<<N<<endl;
	nXo_vs_r->Scale(1.0/N);
	nXo_vs_z->Scale(1.0/N);
	inXo_vs_r->SetBinContent(1,nXo_vs_r->GetBinContent(1));
	for(int j=2; j<=nXo_vs_r->GetNbinsX(); j++){
		inXo_vs_r->SetBinContent(j, nXo_vs_r->GetBinContent(j)+inXo_vs_r->GetBinContent(j-1));
	}
	inXo_vs_z->SetBinContent(1,nXo_vs_z->GetBinContent(1));
	for(int j=2; j<=nXo_vs_z->GetNbinsX(); j++){
		inXo_vs_z->SetBinContent(j, nXo_vs_z->GetBinContent(j)+inXo_vs_z->GetBinContent(j-1));
	}

	pthread_mutex_unlock(&mutex);
	
	return NOERROR;
}


//------------------
// GapIntegration
//------------------
void DEventProcessor_radlen_hists::GapIntegration(TH1F *hin, TH1F *hout)
{
	//======== THIS IS NOT CURRENTLY USED ==================

	// Generate histogram of integrals of radiation lengths
	// Loop over bins of the histogram looking for gaps where (essentially) zero
	// radiation lengths exist. Integrate the radiation lengths between the
	// gaps.
	int bin =0;
	int Nbins = hin->GetNbinsX();
	double min_nXo = 1.0E-7;
	while(bin<Nbins){
		while(hin->GetBinContent(bin)<min_nXo)if(++bin>=Nbins)break;
		if(bin>=Nbins)break;

		double nXo;
		double nXotot = 0.0;
		double pos = 0.0;
		double left = hin->GetBinCenter(bin);
		double right = left;
		while((nXo=hin->GetBinContent(bin))>=min_nXo){
			nXotot += nXo;
			right = hin->GetBinCenter(bin);
			pos += nXo*hin->GetBinCenter(bin);
			if(++bin>=Nbins)break;
		}

		pos/=nXotot; // center position
		//if(nXotot>1.0E-4)cout<<"r(weighted)="<<pos<<"  r(center)="<<(left+right)/2.0<<"  num_Xo="<<nXotot<<endl;
		
		// Find the bins in hout corrsponding to 
		// left and right positions of this material
		int bin_left = hout->FindBin(left);
		int bin_right = hout->FindBin(right);
		for(int i=bin_left; i<=bin_right; i++)hout->SetBinContent(i, nXotot);
	}
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_radlen_hists::fini(void)
{
	return NOERROR;
}
