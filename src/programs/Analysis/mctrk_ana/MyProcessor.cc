// Author: David Lawrence  June 25, 2004
//
//
// MyProcessor.cc
//

#include <iostream>
#include <cmath>
using namespace std;

#include <TThread.h>
#include <TVector3.h>

#include "MyProcessor.h"

#include "TRACKING/DTrack.h"
#include "TRACKING/DTrackHit.h"
#include "TRACKING/DTrackCandidate.h"
//#include "TRACKING/DTrackEfficiency.h"
#include "TRACKING/DMCThrown.h"
#include "FDC/DFDCPseudo.h"
#include "GlueX.h"

TFile *ROOTfile = NULL;

TH1F *FDC_z, *FDC_r, *CDC_z, *CDC_r;
TH3F *FDC_pseudo;

extern const char* OUTPUTFILE;

//------------------------------------------------------------------
// init   -Open output file here (e.g. a ROOT file)
//------------------------------------------------------------------
jerror_t MyProcessor::init(void)
{
	// open ROOT file
	ROOTfile = new TFile(OUTPUTFILE,"RECREATE","Produced by hd_ana");
	cout<<"Opened ROOT file \""<<OUTPUTFILE<<"\""<<endl;
	
	FDC_z = new TH1F("FDC_z","FDC z-hits", 6510, 0.0, 650.0);
	CDC_z = new TH1F("CDC_z","CDC z-hits", 6510, 0.0, 650.0);
	FDC_r = new TH1F("FDC_r","FDC r-hits", 1100,0.0, 100.0);
	CDC_r = new TH1F("CDC_r","CDC r-hits", 1100, 0.0, 100.0);
	
	FDC_pseudo = new TH3F("FDC_pseudo", "",100,-60.0,60.0, 100, -60., 60.,100,200.0, 410.0);
	
	//R_vs_theta = new TH2F("R_vs_theta","R_vs_theta", 180, 0.0, M_PI, 1500, 0.0, 1500.0);
	//R_over_sintheta_vs_theta = new TH2F("R_over_sintheta_vs_theta","R_over_sintheta_vs_theta", 180, 0.0, M_PI, 1500, 0.0, 1500.0);
	
	// Create Tree
	fit_parms = new TTree("fitp","Helical Fit parameters");
	fit_parms->Branch("F",val,"p/F:px:py:pz:x:y:z:pcan:pcanx:pcany:pcanz:canz:Ro:x0:y0:phi:theta:p_thrn:px_thrn:py_thrn:pz_thrn:phi_thrn:theta_thrn:fittable:nhits:chisq");

	gPARMS->GetParameter("TRKFIT:CANDIDATE_TAG", CANDIDATE_TAG);

	return NOERROR;
}

//------------------------------------------------------------------
// evnt   -Fill histograms here
//------------------------------------------------------------------
jerror_t MyProcessor::evnt(JEventLoop *loop, int eventnumber)
{
#if 0
	// Histograms are created and filled in DEventProcessor_TrackHists
	// Automatically since it was added to the app in mctrk_ana.cc

	// Histograms to determine angles from geometry
	vector<const DTrackHit*> trackhits;
	loop->Get(trackhits, "MC");
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
	
	vector<const DTrack*> tracks;
	vector<const DTrackEfficiency*> trackeffs;
	vector<const DMCThrown*> mcthrowns;
	vector<const DTrackCandidate*> trackcandidates;
	JFactory<DTrack> *fac_track = loop->Get(tracks);
	loop->Get(trackeffs);
	loop->Get(mcthrowns);
	JFactory<DTrackCandidate> *fac_tc = loop->Get(trackcandidates, CANDIDATE_TAG.c_str());
	
	
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
		val[ 4] = track->x;
		val[ 5] = track->y;
		val[ 6] = track->z;
		val[ 7] = tc->p;
		val[ 8] = tc->p_trans*cos(tc->phi);
		val[ 9] = tc->p_trans*sin(tc->phi);
		val[10] = tc->p*cos(tc->theta);
		val[11] = tc->z_vertex;
		val[12] = sqrt(tc->x0*tc->x0 + tc->y0*tc->y0);
		val[13] = tc->x0;
		val[14] = tc->y0;
		val[15] = track->phi;
		val[16] = track->theta;
		val[17] = thrown->p;
		val[18] = thrown->p*sin(thrown->theta)*cos(thrown->phi);
		val[19] = thrown->p*sin(thrown->theta)*sin(thrown->phi);
		val[20] = thrown->p*cos(thrown->theta);
		val[21] = thrown->phi;
		val[22] = thrown->theta;
		val[23] = (float)trackeff->fittable;
		val[24] = (float)trackeff->Nhits_thrown;
		val[25] = track->chisq;
		
		fit_parms->Fill();
		UnlockState();
	}
	
	vector<const DFDCPseudo*> fdcs;
	loop->Get(fdcs);
	for(unsigned int i=0; i<fdcs.size(); i++){
		const DFDCPseudo* hit = fdcs[i];
		
		if(!hit->wire)continue;
		TVector3 pos = hit->wire->origin + hit->s*hit->wire->udir;
		FDC_pseudo->Fill(pos.x(), pos.y(), pos.z());
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

