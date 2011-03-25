#ifndef _DBCALCluster_factory_
#define _DBCALCluster_factory_

/*
 *  DBCALCluster_factory.h
 *
 *  Created by Matthew Shepherd on 3/12/11.
 *
 */

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>

using namespace jana;

#include "BCAL/DBCALHit.h"
#include "BCAL/DBCALCluster.h"

#include "TTree.h"
#include "TFile.h"

//#define BCAL_CLUSTER_DIAGNOSTIC

class DBCALCluster_factory : public JFactory< DBCALCluster > {
  
public:
  
  DBCALCluster_factory();
  ~DBCALCluster_factory(){}
  
private:

  jerror_t evnt(JEventLoop *loop, int eventnumber);	
  
  void clearPoints();
  
  // these routines combine points and clusters together

  vector< DBCALCluster > clusterize( vector< const DBCALPoint* > points ) const;
  void merge( vector< DBCALCluster >& clusters ) const;
  
  // these are the routines used for testing whether things should be
  // combined -- right now very basic, but can be fine tuned in the future

  bool overlap( const DBCALCluster& highEClust,
                const DBCALCluster& lowEClust ) const;
  
  bool overlap( const DBCALCluster& clust,
                const DBCALPoint* point ) const;
  
  bool overlap( const DBCALCluster& clust, 
                const DBCALHit* hit ) const; 
  
  float m_mergeSig;
  
  // we may consider a separate factory to provide the BCAL points at
  // a future stage; for now have this factory own and maintain them
  
  vector< DBCALPoint* > m_bcalPoints;
    
#ifdef BCAL_CLUSTER_DIAGNOSTIC
  
#define MAX_POINT 1000
#define MAX_CLUST 50 
  
  jerror_t init();
  jerror_t fini();
  
  TFile* m_rootFile;
  TTree* m_twoEndPtTr;
  TTree* m_firstClustTr;
  TTree* m_ovrlpTr;

  mutable int m_n2EPt;
  mutable float m_rhoPt[MAX_POINT];
  mutable float m_phiPt[MAX_POINT];
  mutable float m_thePt[MAX_POINT];
  mutable float m_rhoSPt[MAX_POINT];
  mutable float m_phiSPt[MAX_POINT];
  mutable float m_theSPt[MAX_POINT];
  mutable float m_ePt[MAX_POINT];
  mutable float m_tPt[MAX_POINT];
  mutable float m_t0Pt[MAX_POINT];
  
  mutable int m_nCl;
  mutable int m_nPts[MAX_CLUST];
  mutable float m_rhoCl[MAX_CLUST];
  mutable float m_phiCl[MAX_CLUST];
  mutable float m_theCl[MAX_CLUST];
  mutable float m_rhoSCl[MAX_CLUST];
  mutable float m_phiSCl[MAX_CLUST];
  mutable float m_theSCl[MAX_CLUST];
  mutable float m_eCl[MAX_CLUST];
  mutable float m_tCl[MAX_CLUST];
  
  mutable float m_dPhi;
  mutable float m_dThe;
  mutable float m_sigPhi;
  mutable float m_sigThe;
  mutable float m_eClus;
  mutable float m_rhoClus;
  mutable float m_theClus;
  mutable float m_phiClus;
  mutable int m_nClClus;
  
#endif
  
};

#endif 

