// $Id$
//
//    File: FillHitsProc.cc
// Created: Wed Aug 10 13:17:20 EST 2005
// Creator: shepherd (on Darwin 129-79-159-16.dhcp-bl.indiana.edu 7.9.0 powerpc)
//

#include "FillHitsProc.h"

#include "FCAL/DFCALHit.h"
#include "BCAL/DHDDMBCALHit.h"
#include "CDC/DCDCHit.h"
#include "FDC/DFDCHit.h"
#include "TOF/DHDDMTOFHit.h"
#include "START_COUNTER/DSCHit.h"

//------------------
// init
//------------------
jerror_t FillHitsProc::init(void)
{
  // Create histograms here
	

  m_fcalFile = new TFile( "fcal_hits.root", "RECREATE" );
  m_fcalFile->cd();
  m_fcalHitTree = new TTree( "fcalHits", "FCAL Hits" );
  m_fcalHitTree->Branch( "nHits", &m_nHits, "nHits/I" );
  m_fcalHitTree->Branch( "x", m_x, "x[nHits]/F" );
  m_fcalHitTree->Branch( "y", m_y, "y[nHits]/F" );
  m_fcalHitTree->Branch( "E", m_E, "E[nHits]/F" );
  m_fcalHitTree->Branch( "t", m_t, "t[nHits]/F" );

  m_bcalFile = new TFile( "bcal_hits.root", "RECREATE" );
  m_bcalFile->cd();
  m_bcalHitTree = new TTree( "bcalHits", "BCAL Hits" );
  m_bcalHitTree->Branch( "nHits", &m_nHits, "nHits/I" );
  m_bcalHitTree->Branch( "module", m_module, "module[nHits]/F" );
  m_bcalHitTree->Branch( "layer", m_layer, "layer[nHits]/F" );
  m_bcalHitTree->Branch( "sector", m_sector, "sector[nHits]/F" );
  m_bcalHitTree->Branch( "E", m_E, "E[nHits]/F" );
  m_bcalHitTree->Branch( "t", m_t, "t[nHits]/F" );

  m_cdcFile = new TFile( "cdc_hits.root", "RECREATE" );
  m_cdcFile->cd();
  m_cdcHitTree = new TTree( "cdcHits", "CDC Hits" );
  m_cdcHitTree->Branch( "nHits", &m_nHits, "nHits/I" );
  m_cdcHitTree->Branch( "ring", m_ring, "ring[nHits]/F" );
  m_cdcHitTree->Branch( "straw", m_straw, "straw[nHits]/F" );
  m_cdcHitTree->Branch( "E", m_E, "E[nHits]/F" );
  m_cdcHitTree->Branch( "t", m_t, "t[nHits]/F" );

  m_fdcFile = new TFile( "fdc_hits.root", "RECREATE" );
  m_fdcFile->cd();
  m_fdcHitTree = new TTree( "fdcHits", "FDC Hits" );
  m_fdcHitTree->Branch( "nHits", &m_nHits, "nHits/I" );
  m_fdcHitTree->Branch( "plane", m_plane, "plane[nHits]/F" );
  m_fdcHitTree->Branch( "layer", m_layer, "layer[nHits]/F" );
  m_fdcHitTree->Branch( "element", m_element, "element[nHits]/F" );
  m_fdcHitTree->Branch( "r", m_r, "r[nHits]/F" );
  m_fdcHitTree->Branch( "q", m_q, "q[nHits]/F" );
  m_fdcHitTree->Branch( "t", m_t, "t[nHits]/F" );

  m_tofFile = new TFile( "tof_hits.root", "RECREATE" );
  m_tofFile->cd();
  m_tofHitTree = new TTree( "tofHits", "TOF Hits" );
  m_tofHitTree->Branch( "nHits", &m_nHits, "nHits/I" );
  m_tofHitTree->Branch( "plane", m_plane, "plane[nHits]/F" );
  m_tofHitTree->Branch( "bar", m_bar, "bar[nHits]/F" );
  m_tofHitTree->Branch( "dE_north", m_dE_north, "dE_north[nHits]/F" );
  m_tofHitTree->Branch( "t_north", m_t_north, "t_north[nHits]/F" );
  m_tofHitTree->Branch( "dE_south", m_dE_south, "dE_south[nHits]/F" );
  m_tofHitTree->Branch( "t_south", m_t_south, "t_south[nHits]/F" );

  m_startFile = new TFile( "start_hits.root", "RECREATE" );
  m_startFile->cd();
  m_startHitTree = new TTree( "startHits", "Start Hits" );
  m_startHitTree->Branch( "nHits", &m_nHits, "nHits/I" );
  m_startHitTree->Branch( "E", m_E, "E[nHits]/F" );
  m_startHitTree->Branch( "t", m_t, "t[nHits]/F" );
  m_startHitTree->Branch( "sector", m_sector, "sector[nHits]/F" );

  m_ckovFile = new TFile( "ckov_hits.root", "RECREATE" );
  m_ckovFile->cd();
  m_ckovHitTree = new TTree( "ckovHits", "Cherenkov Hits" );
  
  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t FillHitsProc::brun(JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t FillHitsProc::evnt(JEventLoop *loop, int eventnumber)
{
  
  // FCAL Hits:
  vector<const DFCALHit*> fcalHitVect;
  loop->Get( fcalHitVect );
	
  m_nHits = fcalHitVect.size();
  assert( m_nHits < kMaxHits );

  for( int i = 0; i < m_nHits; ++i ){
	  
    m_x[i] = fcalHitVect[i]->x;
    m_y[i] = fcalHitVect[i]->y;
    m_E[i] = fcalHitVect[i]->E;
    m_t[i] = fcalHitVect[i]->t;
  }
  
  m_fcalHitTree->Fill();

  // BCAL Hits:
  vector<const DHDDMBCALHit*> bcalHitVect;
  loop->Get( bcalHitVect );
	
  m_nHits = bcalHitVect.size();
  assert( m_nHits < kMaxHits );

  for( int i = 0; i < m_nHits; ++i ){
	  
    m_module[i] = bcalHitVect[i]->module;
    m_layer[i] = bcalHitVect[i]->layer;
    m_sector[i] = bcalHitVect[i]->sector;
    m_E[i] = bcalHitVect[i]->E;
    m_t[i] = bcalHitVect[i]->t;
  }
  
  m_bcalHitTree->Fill();

  // CDC Hits:
  vector<const DCDCHit*> cdcHitVect;
  loop->Get( cdcHitVect );
	
  m_nHits = cdcHitVect.size();
  assert( m_nHits < kMaxHits );

  for( int i = 0; i < m_nHits; ++i ){
	  
    m_ring[i] = cdcHitVect[i]->ring;
    m_straw[i] = cdcHitVect[i]->straw;
    m_E[i] = cdcHitVect[i]->dE;
    m_t[i] = cdcHitVect[i]->t;
  }
  
  m_cdcHitTree->Fill();

  // FDC Hits:
  vector<const DFDCHit*> fdcHitVect;
  loop->Get( fdcHitVect );
	
  m_nHits = fdcHitVect.size();
  assert( m_nHits < kMaxHits );

  for( int i = 0; i < m_nHits; ++i ){
	  
    m_plane[i] = fdcHitVect[i]->gPlane;
    m_layer[i] = fdcHitVect[i]->gLayer;
    m_element[i] = fdcHitVect[i]->element;
    m_r[i] = fdcHitVect[i]->r;
    m_q[i] = fdcHitVect[i]->q;
    m_t[i] = fdcHitVect[i]->t;
  }
  
  m_fdcHitTree->Fill();

  // TOF Hits:
  vector<const DHDDMTOFHit*> tofHitVect;
  loop->Get( tofHitVect );
	
  m_nHits = tofHitVect.size();
  assert( m_nHits < kMaxHits );

  for( int i = 0; i < m_nHits; ++i ){
	  
    m_plane[i] = tofHitVect[i]->plane;
    m_bar[i] = tofHitVect[i]->bar;
    m_dE_north[i] = tofHitVect[i]->dE_north;
    m_t_north[i] = tofHitVect[i]->t_north;
    m_dE_south[i] = tofHitVect[i]->dE_south;
    m_t_south[i] = tofHitVect[i]->t_south;
  }
  
  m_tofHitTree->Fill();

  // Start Hits:
  vector<const DSCHit*> scHitVect;
  loop->Get( scHitVect );

  m_nHits = scHitVect.size();
  assert( m_nHits < kMaxHits );

  for( int i = 0; i < m_nHits; ++i ){

    m_sector[i] = scHitVect[i]->sector;
    m_E[i] = scHitVect[i]->dE;
    m_t[i] = scHitVect[i]->t;
  }

  m_startHitTree->Fill();
      
  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t FillHitsProc::erun(void)
{
  // Any final calculations on histograms (like dividing them)
  // should be done here. This may get called more than once.
  return NOERROR;
}

//------------------
// fini
//------------------
jerror_t FillHitsProc::fini(void)
{
  // If another DEventProcessor is in the list ahead of this one, then
  // it will have finished before this is called. e.g. closed the
  // ROOT file!

  m_fcalFile->cd();
  m_fcalHitTree->Write();
  m_fcalFile->Close();

  m_bcalFile->cd();
  m_bcalHitTree->Write();
  m_bcalFile->Close();

  m_cdcFile->cd();
  m_cdcHitTree->Write();
  m_cdcFile->Close();

  m_fdcFile->cd();
  m_fdcHitTree->Write();
  m_fdcFile->Close();

  m_tofFile->cd();
  m_tofHitTree->Write();
  m_tofFile->Close();

  m_startFile->cd();
  m_startHitTree->Write();
  m_startFile->Close();
	
  return NOERROR;
}

