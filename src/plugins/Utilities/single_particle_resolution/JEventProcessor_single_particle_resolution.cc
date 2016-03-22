// $Id$
//
//    File: JEventProcessor_single_particle_resolution.cc
// Created: Mon Feb  8 15:12:19 EST 2016
// Creator: dalton (on Linux gluon02.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//

//#include "JEventProcessor_single_particle_resolution.h"
//using namespace jana;


#include <stdint.h>
#include <vector>

#include "JEventProcessor_single_particle_resolution.h"
#include <JANA/JApplication.h>

using namespace std;
using namespace jana;

#include "TRACKING/DMCThrown.h"
#include "TRACKING/DTrackTimeBased.h"
#include "PID/DChargedTrack.h"
#include "PID/DChargedTrackHypothesis.h"

#include <TDirectory.h>
#include <TH2.h>
#include <TH1.h>
#include <TProfile2D.h>
#include <TStyle.h>

// root hist pointers
static TH1I *single_track_num_events;
static TH1I *num_tracks;

static TH2I *mom_recon;
static TH1I *mom_recon_frac;
static TProfile2D *mom_recon_frac_2D;
static TH2F *mom_recon_res_2D;

static TH2I *ang_recon;
static TH1I *ang_recon_diff;
static TProfile2D *ang_recon_diff_2D;
static TH2F *ang_recon_res_2D;

const float PI=3.14159265358;

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_single_particle_resolution());
}
} // "C"


//------------------
// JEventProcessor_single_particle_resolution (Constructor)
//------------------
JEventProcessor_single_particle_resolution::JEventProcessor_single_particle_resolution()
{
	VERBOSE = 0;

	if(gPARMS){
		gPARMS->SetDefaultParameter("single_particle_resolution:VERBOSE", VERBOSE, "single_particle_resolution verbosity level");
	}
}

//------------------
// ~JEventProcessor_single_particle_resolution (Destructor)
//------------------
JEventProcessor_single_particle_resolution::~JEventProcessor_single_particle_resolution()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_single_particle_resolution::init(void)
{
	// This is called once at program startup. If you are creating
	// and filling historgrams in this plugin, you should lock the
	// ROOT mutex like this:
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();
	//
	// First thread to get here makes all histograms. If one pointer is
	// already not NULL, assume all histograms are defined and return now


	// lock all root operations
	japp->RootWriteLock();

	if (single_track_num_events != NULL){
		japp->RootUnLock();
		return NOERROR;
	}
	
	// create root folder for bcal and cd to it, store main dir
	TDirectory *main = gDirectory;
	gDirectory->mkdir("singletrack")->cd();
	gStyle->SetOptStat(1110);

	// set style
	gStyle->SetTitleOffset(1, "Y");
  	gStyle->SetTitleSize(0.05,"xyz");
	gStyle->SetTitleSize(0.08,"h");
	gStyle->SetLabelSize(0.05,"xyz");
	gStyle->SetTitleX(0);
	gStyle->SetTitleAlign(13);
	gStyle->SetNdivisions(505,"xy");

	int nbinsx = 12;
	int nbinsy = 12;

	single_track_num_events = new TH1I("single_track_num_events","Number of single track events",1, 0.5, 1.5);
	num_tracks = new TH1I("num_tracks","Number of tracks",11, -0.5, 10.5);

	mom_recon = new TH2I("mom_recon","Mom recon;Thrown mom;Recon mom",120,0,12,150,0,15);
	mom_recon_frac = new TH1I("mom_recon_frac","Mom recon. frac. diff.;(P_{T} - P_{R})/ P_{T}",200,-1,1);
	mom_recon_frac_2D = new TProfile2D("mom_recon_frac_2D","Mom recon. frac. diff.;P_{T}  (GeV);#theta_{T}   (deg);(P_{T} - P_{R})/ P_{T}",
									   nbinsx,0,12,nbinsy,0,6);

	ang_recon = new TH2I("ang_recon","Ang recon;#theta_{T}  (degrees);#theta_{R}  (degrees)",120,0,6,160,0,8);
	ang_recon_diff = new TH1I("ang_recon_diff","Ang recon. diff.;#theta_{T} - #theta_{R}  (degrees)",200,-2,2);
	ang_recon_diff_2D = new TProfile2D("ang_recon_diff_2D","Ang recon. diff.;P_{T}  (GeV);#theta_{T}   (deg);#theta_{T} - #theta_{R}  (deg)",
									   nbinsx,0,12,nbinsy,0,6);

	
	// int nbinsx =  mom_recon_frac_2D->GetNbinsX();
	// int nbinsy =  mom_recon_frac_2D->GetNbinsY();
	// float lowedgex = mom_recon_frac_2D->GetXaxis()->GetBinLowEdge(1);
	// float lowedgey = mom_recon_frac_2D->GetYaxis()->GetBinLowEdge(1);
	// float hiedgex = mom_recon_frac_2D->GetXaxis()->GetBinUpEdge(nbinsx);
	// float hiedgey = mom_recon_frac_2D->GetYaxis()->GetBinUpEdge(nbinsy);
	mom_recon_res_2D = new TH2F("mom_recon_res_2D","Mom recon. resolution;P_{T}  (GeV);#theta_{T}   (deg);#sigma_{[P_{T} - P_{R})/ P_{T}]}",
								nbinsx,0,12,nbinsy,0,6);
	//nbinsx,lowedgex,hiedgex,nbinsy,lowedgey,hiedgey);
	ang_recon_res_2D = new TH2F("ang_recon_res_2D","Ang recon. resolution;P_{T}  (GeV);#theta_{T}   (deg);#sigma_{[#theta_{T} - #theta_{R}]}",
								nbinsx,0,12,nbinsy,0,6);
	//nbinsx,lowedgex,hiedgex,nbinsy,lowedgey,hiedgey);

	// back to main dir
	main->cd();
	
	// unlock
	japp->RootUnLock();
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_single_particle_resolution::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_single_particle_resolution::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// loop->Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.
	// Here's an example:
	//
	vector<const DChargedTrack*> mytrack;
	vector<const DMCThrown*> mythrown;
	loop->Get(mytrack);
	loop->Get(mythrown);
	
	japp->RootWriteLock();
	
	num_tracks->Fill(mytrack.size());
		
	if ( mytrack.size() == 1 ) { // require exacly 1 detected track
		
		single_track_num_events->Fill(1);

		if ( mythrown.size() == 1 ) {
			
			//const DChargedTrackHypothesis *myhypo = mytrack[0]->Get_Hypothesis(PiPlus);
			const DChargedTrackHypothesis *myhypo = mytrack[0]->Get_BestTrackingFOM();
			float mom_thro = mythrown[0]->pmag();
			float mom_reco = myhypo->pmag();
			float mom_reco_frac = (mom_thro-mom_reco)/mom_thro;

			float ang_thro = 180/PI*mythrown[0]->momentum().Theta();
			float ang_reco = 180/PI*myhypo->momentum().Theta();
			float ang_reco_diff = (ang_thro-ang_reco);

			if (VERBOSE>0 && (mom_reco_frac<-1||mom_reco_frac>1))
				printf("single_particle_resolution   event %5i  thro %7.3f   reco %7.3f   frac %7.3f\n",eventnumber,mom_thro,mom_reco,mom_reco_frac);
			mom_recon->Fill(mom_thro,mom_reco);
			mom_recon_frac->Fill(mom_reco_frac);
			mom_recon_frac_2D->Fill(mom_thro,ang_thro,mom_reco_frac);

			ang_recon->Fill(ang_thro,ang_reco);
			ang_recon_diff->Fill(ang_reco_diff);
			ang_recon_diff_2D->Fill(mom_thro,ang_thro,ang_reco_diff);

		}

	}

	japp->RootUnLock();


	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_single_particle_resolution::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_single_particle_resolution::fini(void)
{
	// Called before program exit after event processing is finished.
	

	int nbinsx =  mom_recon_frac_2D->GetNbinsX();
	int nbinsy =  mom_recon_frac_2D->GetNbinsY();
   
	for (int i=1; i<=nbinsx; i++) {
		for (int j=1; j<=nbinsy; j++) {
			float momres = mom_recon_frac_2D->GetBinError(i,j);
			mom_recon_res_2D->SetBinContent(i,j,momres);
			float angres = ang_recon_diff_2D->GetBinError(i,j);
			ang_recon_res_2D->SetBinContent(i,j,angres);
			//printf("%3i %3i %f %f\n",i,j,momres,angres);
		}
	}
	
	mom_recon_res_2D->SetMaximum(1);
	mom_recon_res_2D->SetMinimum(mom_recon_res_2D->GetMinimum(0.00001));
	mom_recon_res_2D->SetMaximum(1);
	mom_recon_res_2D->SetMinimum(mom_recon_res_2D->GetMinimum(0.00001));
	

	return NOERROR;
}

