// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.cc
//

#include <iostream>
using namespace std;

#include <TThread.h>

#include "MyProcessor.h"

#include "TRACKING/DTrack.h"
#include "TRACKING/DTrackHit.h"
#include "TRACKING/DTrackCandidate.h"
#include "TRACKING/DTrackEfficiency.h"
#include "TRACKING/DMCThrown.h"
#include "GlueX.h"

TFile *ROOTfile = NULL;

TH1F *FDC_z, *FDC_r, *CDC_z, *CDC_r;

//------------------------------------------------------------------
// init   -Open output file here (e.g. a ROOT file)
//------------------------------------------------------------------
jerror_t MyProcessor::init(void)
{
	// open ROOT file
	ROOTfile = new TFile("mctrk_ana.root","RECREATE","Produced by hd_ana");
	cout<<"Opened ROOT file \"mctrk_ana.root\""<<endl;
	
	FDC_z = new TH1F("FDC_z","FDC z-hits", 6510, 0.0, 650.0);
	CDC_z = new TH1F("CDC_z","CDC z-hits", 6510, 0.0, 650.0);
	FDC_r = new TH1F("FDC_r","FDC r-hits", 1100,0.0, 100.0);
	CDC_r = new TH1F("CDC_r","CDC r-hits", 1100, 0.0, 100.0);
	
	//R_vs_theta = new TH2F("R_vs_theta","R_vs_theta", 180, 0.0, M_PI, 1500, 0.0, 1500.0);
	//R_over_sintheta_vs_theta = new TH2F("R_over_sintheta_vs_theta","R_over_sintheta_vs_theta", 180, 0.0, M_PI, 1500, 0.0, 1500.0);
	
	// Create Tree
	fit_parms = new TTree("fitp","Helical Fit parameters");
	fit_parms->Branch("F",val,"p/F:px:py:pz:pcan:pcanx:pcany:pcanz:Ro:phi:theta:p_thrn:px_thrn:py_thrn:pz_thrn:phi_thrn:theta_thrn");

	return NOERROR;
}

//------------------------------------------------------------------
// evnt   -Fill histograms here
//------------------------------------------------------------------
jerror_t MyProcessor::evnt(JEventLoop *loop, int eventnumber)
{
	// Histograms are created and filled in DEventProcessor_TrackHists
	// Automatically since it was added to the app in mctrk_ana.cc

	// Histograms to determine angles from geometry
	vector<const DTrackHit*> trackhits;
	loop->Get(trackhits, "MC");
	TThread::Lock();
	for(unsigned int i=0; i<trackhits.size(); i++){
		const DTrackHit *hit = trackhits[i];
		if(hit->system==SYS_CDC){
			CDC_z->Fill(hit->z);
			CDC_r->Fill(hit->r);
		}
		if(hit->system==SYS_FDC){
			FDC_z->Fill(hit->z);
			FDC_r->Fill(hit->r);
		}
	}
	TThread::UnLock();
	
	vector<const DTrack*> tracks;
	vector<const DTrackEfficiency*> trackeffs;
	vector<const DMCThrown*> mcthrowns;
	vector<const DTrackCandidate*> trackcandidates;
	JFactory<DTrack> *fac_track = loop->Get(tracks);
	loop->Get(trackeffs);
	loop->Get(mcthrowns);
	JFactory<DTrackCandidate> *fac_tc = loop->Get(trackcandidates);
	
	
	// Loop over DMCThrown and DTrackEfficiency objects
	for(unsigned int i=0; i<mcthrowns.size(); i++){
		const DMCThrown *thrown = mcthrowns[i];
		const DTrackEfficiency *trackeff = trackeffs[i];
		const DTrack *track = fac_track->GetByIDT(trackeff->trackid);
		if(!track)continue;
		const DTrackCandidate *tc = fac_tc->GetByIDT(track->candidateid);
		if(!tc)continue;

		LockState();
		val[ 0] = track->p;
		val[ 1] = track->p*sin(track->theta)*cos(track->phi);
		val[ 2] = track->p*sin(track->theta)*sin(track->phi);
		val[ 3] = track->p*cos(track->theta);
		val[ 4] = tc->p;
		val[ 5] = tc->p_trans*cos(tc->phi);
		val[ 6] = tc->p_trans*sin(tc->phi);
		val[ 7] = tc->p*cos(tc->theta);
		val[ 8] = sqrt(tc->x0*tc->x0 + tc->y0*tc->y0);
		val[ 9] = tc->phi;
		val[10] = tc->theta;
		val[11] = thrown->p;
		val[12] = thrown->p*sin(thrown->theta)*cos(thrown->phi);
		val[13] = thrown->p*sin(thrown->theta)*sin(thrown->phi);
		val[14] = thrown->p*cos(thrown->theta);
		val[15] = thrown->phi;
		val[16] = thrown->theta;
		
		fit_parms->Fill();
		UnlockState();
	}
	
#if 0
	for(unsigned int i=0; i<trackcandidates.size(); i++){
		const DTrackCandidate *tc = trackcandidates[i];
		
		double R = sqrt(tc->x0*tc->x0 + tc->y0*tc->y0);
		R_vs_theta->Fill(tc->theta,R);
		R_over_sintheta_vs_theta->Fill(tc->theta,R/sin(tc->theta));
	}
#endif	
	return NOERROR;
}

//------------------------------------------------------------------
// fini   -Close output file here
//------------------------------------------------------------------
jerror_t MyProcessor::fini(void)
{
	ROOTfile->Write();
	delete ROOTfile;
	cout<<endl<<"Closed ROOT file"<<endl;

	return NOERROR;
}

