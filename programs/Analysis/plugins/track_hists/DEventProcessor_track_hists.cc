// $Id$
//
//    File: DEventProcessor_track_hists.cc
// Created: Sat Jun  3 13:58:03 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.4.0 powerpc)
//

#include "DEventProcessor_track_hists.h"

#include <TLorentzVector.h>

#include <DApplication.h>
#include <DTrack.h>
#include <DTrackCandidate.h>

static TFile **tfilePtr = NULL;

// Routine used to create our DEventProcessor
extern "C"{
void InitProcessors(DApplication *app){
	app->AddProcessor(new DEventProcessor_track_hists());
}

void SetTFilePtrAddress(TFile **h){
	tfilePtr = h;
}
} // "C"


#define FCAL_Z_OFFSET 640.0-65.0 // I don't know what this value is ???
//#define FCAL_Z_OFFSET 170.0 // I don't know what this value is ???
#define PI_ZERO_MASS 0.13497

//------------------
// init
//------------------
derror_t DEventProcessor_track_hists::init(void)
{
	// open ROOT file (if needed)
	ROOTfile = NULL;
	if(tfilePtr == NULL)tfilePtr = &ROOTfile;
	if(*tfilePtr == NULL){
		*tfilePtr = ROOTfile = new TFile("track_hists.root","RECREATE","Produced by track_hists.so");
		cout<<"Opened ROOT file \"track_hists.root\""<<endl;
	}else{
		(*tfilePtr)->cd();
	}
	
	track_p = new TH1F("track_p","track_p",1000, 0.0, 9.0);
	trackCandidate_p = new TH1F("trackCandidate_p","trackCandidate_p",1000, 0.0, 9.0);
	track_theta = new TH1F("track_theta","track_theta",200, 0.0, M_PI);
	track_phi = new TH1F("track_phi","track_phi",200, 0.0, 2.0*M_PI);

	return NOERROR;
}

//------------------
// brun
//------------------
derror_t DEventProcessor_track_hists::brun(DEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
derror_t DEventProcessor_track_hists::evnt(DEventLoop *loop, int eventnumber)
{
	vector<const DTrackCandidate*> trackcandidates;
	vector<const DTrack*> tracks;
	loop->Get(trackcandidates);
	loop->Get(tracks);
	
	LockState();
	
	// Single shower params
	for(unsigned int i=0; i<trackcandidates.size(); i++){
		const DTrackCandidate *s = trackcandidates[i];
		trackCandidate_p->Fill(s->p);
	}
	
	// 2-gamma inv. mass
	for(unsigned int i=0; i<tracks.size(); i++){
		const DTrack *s = tracks[i];
		track_p->Fill(s->p);
		track_theta->Fill(s->theta);
		track_phi->Fill(s->phi);
	}

	UnlockState();	

	return NOERROR;
}

//------------------
// erun
//------------------
derror_t DEventProcessor_track_hists::erun(void)
{
	// Any final calculations on histograms (like dividing them)
	// should be done here. This may get called more than once.
	return NOERROR;
}

//------------------
// fini
//------------------
derror_t DEventProcessor_track_hists::fini(void)
{

	if(ROOTfile){
		ROOTfile->Write();
		delete ROOTfile;
		cout<<endl<<"Closed ROOT file"<<endl;
	}

	return NOERROR;
}

