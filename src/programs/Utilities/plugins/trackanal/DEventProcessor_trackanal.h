// $Id$
//
//    File: DEventProcessor_tofanal.h
// Created: Wed Jul 14 17:12:56 EDT 2010
// Creator: zihlmann (on Linux chaos.jlab.org 2.6.18-164.2.1.el5 i686)
//

#ifndef _DEventProcessor_trackanal_
#define _DEventProcessor_trackanal_

using namespace std;

#include <TTree.h>
#include <TFile.h>

#include <JANA/JEventProcessor.h>
#define MaxTrThrown 25
#define MaxTrFit 30

class DEventProcessor_trackanal:public jana::JEventProcessor{
 public:
  DEventProcessor_trackanal();
  ~DEventProcessor_trackanal();
  const char* className(void){return "DEventProcessor_trackanal";}
  
  TTree *TrackTree;
  TFile *ROOTFile;
  
  Int_t EventNum;
  Int_t NTrThrown; // total number of thrown tracks
  Int_t MaxT;  // Maximum number of thrown tracks
  Int_t MaxF;  // Maximum number of fits
  Int_t MaxC;  // Maximum number of track candidates considered

  Int_t ThrownPType[MaxTrThrown]; // particle type of thrown tracks
  Float_t ThrownPp[MaxTrThrown];  // particle momentum of thrown tracks
  Float_t ThrownQ[MaxTrThrown];   // electric charge of thrown particle
  Int_t NTrCand;
  Float_t TrCandP[MaxTrThrown];   // momentum of track candidate
  Float_t TrCandQ[MaxTrThrown];   // charge of track candidate
  Float_t TrCandN[MaxTrThrown];   // number of hits of track candidate
  Int_t NTrCandHits;
  Int_t NTrFit;
  Int_t   trlistPtype[MaxTrFit];  // particle type of track candidate with best FOM
  Float_t trlistPp[MaxTrFit];   // particle momentum of track candidate with best FOM
  Float_t trlistFOM[MaxTrFit];  // figure of merrit
  Float_t trlistchisq[MaxTrFit];  // chisq
  Int_t   trlistcand[MaxTrFit];   // track candidate number
  Float_t nh[MaxTrFit*MaxTrFit] ; // number of hits for each track candidate
  Float_t ptypes[MaxTrFit*MaxTrFit];    // for each track candidate the chamber hits for each particle type
    
 private:
  jerror_t init(void);	
  jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);
  jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);
  jerror_t erun(void);
  jerror_t fini(void);
  
  
  pthread_mutex_t mutex;
  
};

#endif // _DEventProcessor_trackanal_

