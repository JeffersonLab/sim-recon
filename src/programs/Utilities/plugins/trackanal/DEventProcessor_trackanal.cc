// $Id$
//
//    File: DEventProcessor_trackanal.cc
// Created: Wed Jul 14 17:12:56 EDT 2010
// Creator: zihlmann (on Linux chaos.jlab.org 2.6.18-164.2.1.el5 i686)
//

// the local stuff
#include "DEventProcessor_trackanal.h"

// the general stuff
#include <iostream>
#include <iomanip>
using namespace std;

// the cern root stuff
#include <TThread.h>
#include <TDirectoryFile.h>
#include <TLorentzVector.h>

#include <TROOT.h>
#include <TFile.h>

#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>
using namespace jana;

// the DANA stuff
#include <TRACKING/DTrackCandidate.h>
#include <TRACKING/DTrackCandidate_factory.h>
#include <TRACKING/DTrackCandidate_factory_CDC.h>
#include <TRACKING/DTrackCandidate_factory_FDC.h>
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DMCTrackHit.h>
#include <TRACKING/DTrackTimeBased.h>
#include <TRACKING/DTrackTimeBased_factory.h>
#include <FDC/DFDCHit.h>
#include <CDC/DCDCHit.h>
#include <FDC/DFDCPseudo.h>
#include <CDC/DCDCTrackHit.h>
#include <PID/DKinematicData.h>

extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_trackanal());
}
} // "C"

#define DEBUG 0

//------------------
// DEventProcessor_trackanal (Constructor)
//------------------
DEventProcessor_trackanal::DEventProcessor_trackanal()
{

}

//------------------
// ~DEventProcessor_trackanal (Destructor)
//------------------
DEventProcessor_trackanal::~DEventProcessor_trackanal()
{

}

//------------------
// init
//------------------
jerror_t DEventProcessor_trackanal::init(void)
{
  // Create histograms here
  // create Tree to hold TOF data
    
  ROOTFile = new TFile("trackanal.root", "recreate");
  //ROOTFile->cd();
  
  TrackTree = new TTree("TrackTree", "Track raw data");
  
  TrackTree->Branch("EventNum",    &EventNum,    "EventNum/I"); 
  TrackTree->Branch("MaxT",    &MaxT,    "MaxT/I"); 
  TrackTree->Branch("MaxC",    &MaxC,    "MaxC/I"); 
  TrackTree->Branch("MaxF",    &MaxF,    "MaxF/I"); 

  TrackTree->Branch("NTrThrown",   &NTrThrown,   "NTrThrown/I"); 
  TrackTree->Branch("ThrownPType", ThrownPType,  "ThrownPType[NTrThrown]/I"); 
  TrackTree->Branch("ThrownPp",    ThrownPp,     "ThrownPp[NTrThrown]/F"); 
  TrackTree->Branch("ThrownQ",     ThrownQ,      "ThrownQ[NTrThrown]/F"); 

  TrackTree->Branch("NTrCand",     &NTrCand,     "NTrCand/I"); 
  TrackTree->Branch("TrCandQ",      TrCandQ,     "TrCandQ[NTrCand]/F"); 
  TrackTree->Branch("TrCandP",      TrCandP,     "TrCandP[NTrCand]/F"); 
  TrackTree->Branch("TrCandN",      TrCandN,     "TrCandN[NTrCand]/F"); 
  TrackTree->Branch("TrCandM",      TrCandM,     "TrCandM[NTrCand]/F"); 

  // NTRCand: number of track candidates
  // TrCandQ: charge of the track candidates
  // TrCandP: particle momentum of the track candidates
  // TrCandN: number of detrees of freedom of track candidates
  // TrCandM: number of hits in track candidate with found match to TruthPoints

  TrackTree->Branch("NTrFit",      &NTrFit,      "NTrFit/I");
  TrackTree->Branch("trlistcand",  trlistcand,   "trlistcand[NTrFit]/I");
  TrackTree->Branch("trlistPtype", trlistPtype,  "trlistPtype[NTrFit]/I");
  TrackTree->Branch("trlistFOM",   trlistFOM,    "trlistFOM[NTrFit]/F");
  TrackTree->Branch("trlistPp",    trlistPp,     "trlistPp[NTrFit]/F");
  TrackTree->Branch("trlistchisq", trlistchisq,  "trlistchisq[NTrFit]/F");


  TrackTree->Branch("NTrCandHits", &NTrCandHits, "NTrCandHits/I"); 
  TrackTree->Branch("nh",          nh,           "nh[NTrCandHits]/F");
  TrackTree->Branch("ptypes",      ptypes,       "ptypes[NTrCandHits]/F");

  MaxT = MaxTrThrown;
  MaxC = MaxTrCand;
  MaxF = MaxTrFit;

  Int_t MaxArray = MaxF*MaxF;

  // initialize varialbes just for the fun of it
  for (Int_t i=0;i<MaxTrThrown;i++){
    ThrownPType[i] = 0;
    ThrownPp[i]    = 0.0;
    ThrownQ[i]     = 999.;
    TrCandQ[i]     = 999.;
    TrCandP[i]     = 0.0;
    TrCandN[i]     = 0.0;
  }
  for (Int_t i=0;i<MaxTrFit;i++){
    trlistcand[i]  = 0;
    trlistPtype[i] = 0;
    trlistFOM[i]   = 0.0;
    trlistPp[i]    = 0.0;
  }
  for (Int_t i=0;i<MaxArray;i++){
    nh[i] = 0;
    ptypes[i] = 0;
  }

  pthread_mutex_init(&mutex, NULL);

  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_trackanal::brun(JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_trackanal::evnt(JEventLoop *loop, int eventnumber)
{


  // Fill histograms here
  vector<const DMCThrown*> trThrown;
  loop->Get(trThrown);

  vector<const DTrackCandidate*> trCand;
  loop->Get(trCand);

  vector<const DTrackTimeBased*> trFit;
  loop->Get(trFit);
  
  int ntr = trFit.size();

  pthread_mutex_lock(&mutex);

  if (ntr>MaxTrFit)
    ntr = MaxTrFit;

  NTrFit = ntr;

  for (int i=0;i<NTrFit;i++){
    int cid = trFit[i]->candidateid;

    trlistcand[i] = cid; 
    trlistFOM[i] = trFit[i]->FOM;
    trlistchisq[i] = trFit[i]->chisq;
    trlistPp[i] = trFit[i]->pmag();
    trlistPtype[i] = 0;
    float mass = (float)trFit[i]->mass();
    if (mass>.9){
      trlistPtype[i] = 14;
    } else{
      if (trFit[i]->charge() < 0){
	trlistPtype[i] = 9;
      } else{
	trlistPtype[i] = 8;
      }
    }
    
  }
    
  EventNum = eventnumber;
  NTrThrown = trThrown.size();

  for (int i=0;i<NTrThrown;i++){
    ThrownPType[i] = trThrown[i]->type;
    ThrownPp[i] = trThrown[i]->pmag();
    ThrownQ[i] = trThrown[i]->charge();
  }

  ntr  = trCand.size();

  if (ntr>MaxC){
    NTrCand = MaxC;
  } else {
    NTrCand = ntr;
  }
  
  NTrCandHits = NTrCand*MaxTrFit;

  // loop over all track candidates and find how many wire chamber hits contribute
  // and from which particle track they stem.

  for (int j=0;j<MaxF;j++){
    for (int i=0;i<MaxF;i++){
      ptypes[j*MaxF+i] = 0.0 ;
      nh[j*MaxF+i] = 0.0;
    }
  }
  
  int hitcounter;
  for (int k=0;k<NTrCand;k++) {

    if (k>=MaxC){
      continue;
    }

    TrCandQ[k] = trCand[k]->charge();
    TrCandP[k] = trCand[k]->pmag();
    TrCandN[k] = 0.;
    TrCandM[k] = 0.;

    vector<const DFDCPseudo*> fdcPHits;
    trCand[k]->Get(fdcPHits);
    int npfdc = fdcPHits.size();
    TrCandN[k] += (Float_t)npfdc;

    hitcounter = 0;
    for (int j=0;j<npfdc;j++){

      // for each FDC Pseudo hit get associated object DMCTrackHit
      vector<const DMCTrackHit*> mcTrackHits;
      fdcPHits[j]->Get(mcTrackHits);
      Int_t asize = mcTrackHits.size();

      for (int i=0;i<asize;i++) {
	int tr = mcTrackHits[i]->track; // track number if >NTrThrown scondary

	if (tr > NTrThrown){
	  tr = 0; // any track other than thrown is set as zero
	}

	// this means nh contains only hits from initially thrown particles.
	nh[k*MaxF + tr] += 1.0;
	hitcounter++;
	int pt = mcTrackHits[i]->ptype;
	if (pt>14){
	  pt = 0;
	}
	ptypes[k*MaxF + pt] += 1.0;
      }
      TrCandM[k] += (Float_t) asize;
    }

    
    vector<const DCDCTrackHit*> cdcTrackHits;
    trCand[k]->Get(cdcTrackHits);
    int nhitscdc = cdcTrackHits.size();
    TrCandN[k] += (Float_t)nhitscdc;

    // loop over cdc hits used for this track candiate
    for (int j=0;j<nhitscdc;j++){

      // for each hit get the associated hit object DMCTrackHit 
      vector<const DMCTrackHit*> cdcHits;
      cdcTrackHits[j]->Get(cdcHits);

      Int_t asize = cdcHits.size();

      for (int n=0;n<asize;n++) {
	int tr = cdcHits[n]->track;

	if (tr > NTrThrown){
	  tr = 0;
	}
	nh[k*MaxF + tr] += 1.;
	hitcounter++;
	int pt = cdcHits[n]->ptype;
	if (pt>14){
	  pt = 0;
	}
	ptypes[k*MaxF + pt] += 1.;
      }

      TrCandM[k] += (Float_t)asize;
    }

  }

  if (DEBUG){
    cout<<EventNum << "  " << NTrThrown << "  " << NTrCand << "  " << NTrFit 
	<< "  " << NTrCandHits << endl;
      }

  TrackTree->Fill();
  
  pthread_mutex_unlock(&mutex);

  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_trackanal::erun(void)
{
  // Any final calculations on histograms (like dividing them)
  // should be done here. This may get called more than once.
  return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_trackanal::fini(void)
{
	// Called at very end. This will be called only once
  pthread_mutex_lock(&mutex);

  //ROOTFile->cd();
  ROOTFile->Write();
  ROOTFile->Close();
  //delete ROOTFile;
  cout<<endl<<"Close ROOT file"<<endl;

  pthread_mutex_unlock(&mutex);

  return NOERROR;
}

