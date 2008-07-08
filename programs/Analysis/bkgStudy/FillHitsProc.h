// $Id$
//
//    File: FillHitsProc.h
// Created: Wed Aug 10 13:17:20 EST 2005
// Creator: shepherd (on Darwin 129-79-159-16.dhcp-bl.indiana.edu 7.9.0 powerpc)
//

#ifndef _FillHitsProc_
#define _FillHitsProc_

#include "JANA/JEventProcessor.h"
#include "TFile.h"
#include "TTree.h"

using namespace jana;

class FillHitsProc:public JEventProcessor{
 public:

  enum { kMaxHits = 2000 };
  
  FillHitsProc(){};
  ~FillHitsProc(){};
  const char* className(void){return "FillHitsProc";}

 private:
  
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.

    TFile* m_fcalFile;
    TTree* m_fcalHitTree;

    TFile* m_bcalFile;
    TTree* m_bcalHitTree;

    TFile* m_cdcFile;
    TTree* m_cdcHitTree;

    TFile* m_fdcFile;
    TTree* m_fdcHitTree;

    TFile* m_tofFile;
    TTree* m_tofHitTree;

    TFile* m_startFile;
    TTree* m_startHitTree;

    TFile* m_ckovFile;
    TTree* m_ckovHitTree;

    int m_nHits;
    float m_x[kMaxHits];
    float m_y[kMaxHits];
    float m_r[kMaxHits];
    float m_E[kMaxHits];
    float m_t[kMaxHits];
    float m_t_north[kMaxHits];
    float m_t_south[kMaxHits];
    float m_dE_north[kMaxHits];
    float m_dE_south[kMaxHits];
    float m_q[kMaxHits];
    float m_module[kMaxHits];
    float m_sector[kMaxHits];
    float m_ring[kMaxHits];
    float m_straw[kMaxHits];
    float m_plane[kMaxHits];
    float m_layer[kMaxHits];
    float m_row[kMaxHits];
    float m_element[kMaxHits];
    float m_bar[kMaxHits];

};

#endif // _FillHitsProc_

