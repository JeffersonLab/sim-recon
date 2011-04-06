/*
 *  DBCALCluster_factory.cc
 *
 *  Created by Matthew Shepherd on 3/12/11.
 *
 */

#include <iostream>

using namespace std;

#include "BCAL/DBCALGeometry.h"
#include "BCAL/DBCALHit.h"

#include "BCAL/DBCALCluster_factory.h"

#include "units.h"

bool PointSort( const DBCALPoint* p1, const DBCALPoint* p2 ){
  
  return ( p1->E() > p2->E() );
}

bool ClusterSort( const DBCALCluster& c1, const DBCALCluster& c2 ){
  
  return ( c1.E() > c2.E() );
}

DBCALCluster_factory::DBCALCluster_factory() : 
m_mergeSig( 5 ), 
m_moliereRadius( 3.7*k_cm ) {
  
}

#ifdef BCAL_CLUSTER_DIAGNOSTIC

jerror_t
DBCALCluster_factory::init(void){

  m_rootFile = new TFile( "bcal_clust_diag.root", "RECREATE" );
  
  m_twoEndPtTr = new TTree( "twoEndPtTr", "BCAL two-ended points" );
  m_twoEndPtTr->Branch( "n2EPt", &m_n2EPt, "n2EPt/I" );
  m_twoEndPtTr->Branch( "rhoPt", m_rhoPt, "rhoPt[n2EPt]/F" );
  m_twoEndPtTr->Branch( "phiPt", m_phiPt, "phiPt[n2EPt]/F" );
  m_twoEndPtTr->Branch( "thePt", m_thePt, "thePt[n2EPt]/F" );
  m_twoEndPtTr->Branch( "rhoSPt", m_rhoSPt, "rhoSPt[n2EPt]/F" );
  m_twoEndPtTr->Branch( "phiSPt", m_phiSPt, "phiSPt[n2EPt]/F" );
  m_twoEndPtTr->Branch( "theSPt", m_theSPt, "theSPt[n2EPt]/F" );
  m_twoEndPtTr->Branch( "ePt", m_ePt, "ePt[n2EPt]/F" );
  m_twoEndPtTr->Branch( "tPt", m_tPt, "tPt[n2EPt]/F" );
  m_twoEndPtTr->Branch( "t0Pt", m_t0Pt, "t0Pt[n2EPt]/F" );
  
  m_firstClustTr = new TTree( "firstClustTr", "First-pass clusters" );
  m_firstClustTr->Branch( "nCl", &m_nCl, "nCl/I" );
  m_firstClustTr->Branch( "nPts", m_nPts, "nPts[nCl]/I" );
  m_firstClustTr->Branch( "rhoCl", m_rhoCl, "rhoCl[nCl]/F" );
  m_firstClustTr->Branch( "phiCl", m_phiCl, "phiCl[nCl]/F" );
  m_firstClustTr->Branch( "theCl", m_theCl, "theCl[nCl]/F" );
  m_firstClustTr->Branch( "rhoSCl", m_rhoSCl, "rhoSCl[nCl]/F" );
  m_firstClustTr->Branch( "phiSCl", m_phiSCl, "phiSCl[nCl]/F" );
  m_firstClustTr->Branch( "theSCl", m_theSCl, "theSCl[nCl]/F" );
  m_firstClustTr->Branch( "eCl", m_eCl, "eCl[nCl]/F" );
  m_firstClustTr->Branch( "tCl", m_tCl, "tCl[nCl]/F" );
  
  m_ovrlpTr = new TTree( "ovrlpTr", "Point Cluster Overlap" );
  m_ovrlpTr->Branch( "dPhi", &m_dPhi, "dPhi/F" );
  m_ovrlpTr->Branch( "dThe", &m_dThe, "dThe/F" );
  m_ovrlpTr->Branch( "sep", &m_sep, "sep/F" );
  m_ovrlpTr->Branch( "sigPhi", &m_sigPhi, "sigPhi/F" );
  m_ovrlpTr->Branch( "sigThe", &m_sigThe, "sigThe/F" );
  m_ovrlpTr->Branch( "eClus", &m_eClus, "eClus/F" );
  m_ovrlpTr->Branch( "rhoClus", &m_rhoClus, "rhoClus/F" );
  m_ovrlpTr->Branch( "phiClus", &m_phiClus, "phiClus/F" );
  m_ovrlpTr->Branch( "theClus", &m_theClus, "theClus/F" );
  m_ovrlpTr->Branch( "nClClus", &m_nClClus, "nClClus/I" );
  
  return NOERROR;

}

jerror_t
DBCALCluster_factory::fini( void ){
  
  m_rootFile->cd();
  
  m_twoEndPtTr->Write();
  m_firstClustTr->Write();
  m_ovrlpTr->Write();
  m_rootFile->Write();
  m_rootFile->Close();
  
  return NOERROR;
}

#endif

jerror_t
DBCALCluster_factory::evnt( JEventLoop *loop, int eventnumber ){

  clearPoints();
  
  vector< const DBCALHit* > hits;
  loop->Get( hits );
  if( hits.size() <= 0 ) return NOERROR;
  
  // first arrange the list of hits so they are grouped by cell
  map< int, vector< const DBCALHit* > > cellHitMap;
  for( vector< const DBCALHit* >::const_iterator hitPtr = hits.begin();
      hitPtr != hits.end();
      ++hitPtr ){
    
    const DBCALHit& hit = (**hitPtr);
    
    // mcsmear will produce hits with energy zero if the hits do not
    // exceed the threshold -- we want to suppress these hits so we know
    // exactly how many "real" hits we have in a cell
    
    if( hit.E < 0.1*k_MeV ) continue;
    
    int id = DBCALGeometry::cellId( hit.module, hit.layer, hit.sector );
    
    if( cellHitMap.find( id ) == cellHitMap.end() ){
      
      cellHitMap[id] = vector< const DBCALHit* >();
    }
    
    cellHitMap[id].push_back( *hitPtr );
  }
  
  // now go through this list and group hits into BCAL points
  // this combines information from both ends 
  
  vector< const DBCALPoint* > twoEndPoint;
  
  for( map< int, vector< const DBCALHit* > >::iterator mapItr = cellHitMap.begin();
       mapItr != cellHitMap.end();
       ++mapItr ){
    
    if( mapItr->second.size() == 2 && 
        ( mapItr->second[0]->end != mapItr->second[1]->end ) ){
      
      // start with the good stuff -- one hit on each end of a cell
      
      m_bcalPoints.push_back( new DBCALPoint( *(mapItr->second[0]), 
                                              *(mapItr->second[1]) ) );

      twoEndPoint.push_back( m_bcalPoints.back() );      
    }
  }

#ifdef BCAL_CLUSTER_DIAGNOSTIC
  
  m_rootFile->cd();
  
  m_n2EPt = twoEndPoint.size();
  assert( m_n2EPt <= MAX_POINT );
  
  for( int i = 0; i < m_n2EPt; ++i ){
   
    m_rhoPt[i] = twoEndPoint[i]->rho();
    m_phiPt[i] = twoEndPoint[i]->phi();
    m_thePt[i] = twoEndPoint[i]->theta();
    m_rhoSPt[i] = twoEndPoint[i]->sigRho();
    m_phiSPt[i] = twoEndPoint[i]->sigPhi();
    m_theSPt[i] = twoEndPoint[i]->sigTheta();
    m_ePt[i] = twoEndPoint[i]->E();
    m_tPt[i] = twoEndPoint[i]->t();
    m_t0Pt[i] = twoEndPoint[i]->tInnerRadius();
  }
  
  m_twoEndPtTr->Fill();
  
#endif // BCAL_CLUSTER_DIAGNOSTIC
  
  // now try to clusterize the points
  
  vector< DBCALCluster > clusters = clusterize( twoEndPoint );
  
#ifdef BCAL_CLUSTER_DIAGNOSTIC
  
  m_nCl = clusters.size();
  assert( m_nCl <= MAX_CLUST );
  
  for( int i = 0; i < m_nCl; ++i ){
    
    m_nPts[i] = clusters[i].nCells();
    m_rhoCl[i] = clusters[i].rho();
    m_phiCl[i] = clusters[i].phi();
    m_theCl[i] = clusters[i].theta();
    m_rhoSCl[i] = clusters[i].sigRho();
    m_phiSCl[i] = clusters[i].sigPhi();
    m_theSCl[i] = clusters[i].sigTheta();
    m_eCl[i] = clusters[i].E();
    m_tCl[i] = clusters[i].t();
  }
  
  m_firstClustTr->Fill();
  
#endif // BCAL_CLUSTER_DIAGNOSTIC
  
  /*
   
    MRS (6-Apr-11):  needs work!  most single end hits are actually SiPM
    dark noise -- this causes shower to accumulate tons of extra energy.
    This is likely due to a poorly written overlap routine.  The routine
    for example, should take into acount radial separation of the hit from
    the cluster center.
   
   
  // now we should try to add on single hits...
  for( map< int, vector< const DBCALHit* > >::iterator mapItr = cellHitMap.begin();
      mapItr != cellHitMap.end();
      ++mapItr ){
    
    if( mapItr->second.size() == 1 ){
      
      // only one hit in the cell -- loop over clusters and see if there
      // is a cluster in the vicinity of this cell
      
      const DBCALHit* hit = mapItr->second[0];
      
      for( vector< DBCALCluster >::iterator clust = clusters.begin();
          clust != clusters.end();
          ++clust ){
                
        if( overlap( *clust, hit ) ){

          int cellId = 
             DBCALGeometry::cellId( hit->module, hit->layer, hit->sector );
          float r = DBCALGeometry::r( cellId );
          
          // given the location of the cluster, we need the best guess
          // for z with respect to target at this radius
          float z = r / tan( clust->theta() );
          
          // make a point there
          DBCALPoint* pt = new DBCALPoint( *hit, z );
          
          // add it to the cluster and keep track of it for later cleanup
          //
          // note that this point doesn't really provide an independent z
          // measurement for the cluster, but the cluster will treat it as
          // if it does -- probably not a problem since these are low energy
          // and we are really after lost energy not enhanced position
          // resolution
          
          clust->addPoint( pt );
          m_bcalPoints.push_back( pt );
        }
      }
    }
  }
*/
  
  // load our vector of clusters into the factory member data
  for( vector< DBCALCluster >::iterator clust = clusters.begin();
      clust != clusters.end();
      ++clust ){
   
    // put in an energy threshold for clusters
    if( clust->E() > 30*k_MeV ) {
    
      _data.push_back( new DBCALCluster( *clust ) );
    }
  }

  return NOERROR;
}

vector< DBCALCluster >
DBCALCluster_factory::clusterize( vector< const DBCALPoint* > points ) const {

  // first sort the points by energy
  sort( points.begin(), points.end(), PointSort );

  vector< DBCALCluster > clusters( 0 );
  
  // ahh.. more hard coded numbers that should probably
  // come from a database or something
  float seedThresh = 1*k_GeV;
  float minSeed = 20*k_MeV;
  
  while( seedThresh > minSeed ) {
    
    bool usedPoint = false;
    
    for( vector< const DBCALPoint* >::iterator pt = points.begin();
        pt != points.end();
        ++pt ){
    
      // first see if point should be added to an existing
      // cluster
      
      for( vector< DBCALCluster >::iterator clust = clusters.begin();
           clust != clusters.end();
          ++clust ){
        
        if( overlap( *clust, *pt ) ){
                    
          clust->addPoint( *pt );
          points.erase( pt );
          usedPoint = true;
        }
        
        // once we erase a point the iterator is no longer useful
        // and we start the loop over
        if( usedPoint ) break;
      }
      
      if( usedPoint ) break;
        
      // if the point doesn't overlap with a cluster
      // see if it can become a new seed
      
      if( (**pt).E() > seedThresh ){
              
        clusters.push_back( DBCALCluster( *pt ) );
        points.erase( pt );
        usedPoint = true;
      }
      else break;
      
      if( usedPoint ) break;
    }
    
    merge( clusters );
    
    // lower the threshold to look for new seeds if none of 
    // the existing points were used as new clusters or assigned
    // to existing clusters
    if( !usedPoint ) seedThresh /= 2.;
  }
  
  return clusters;
}

void
DBCALCluster_factory::merge( vector< DBCALCluster >& clusters ) const {
  
  if( clusters.size() <= 1 ) return;
  
  sort( clusters.begin(), clusters.end(), ClusterSort );
  
  bool stillMerging = true;
  
  while( stillMerging ){
  
    stillMerging = false;
    for( vector< DBCALCluster >::iterator hClust = clusters.begin();
        hClust != clusters.end() - 1;
        ++hClust ){
    
      for( vector< DBCALCluster >::iterator lClust = hClust + 1;
          lClust != clusters.end();
          ++lClust ){
      
        if( overlap( *hClust, *lClust ) ){
                  
          hClust->mergeClust( *lClust );
          clusters.erase( lClust );
          
          // now iterators are invalid and we need to bail out of loops
          stillMerging = true;
          break;
        }
      }
      if( stillMerging ) break;
    }
  }
}

bool
DBCALCluster_factory::overlap( const DBCALCluster& highEClust,
                               const DBCALCluster& lowEClust ) const {
  
  float sigTheta = fabs( highEClust.theta() - lowEClust.theta() ) / 
    sqrt( highEClust.sigTheta() * highEClust.sigTheta() +
          lowEClust.sigTheta()  * lowEClust.sigTheta() );

  // difference in phi is tricky due to overlap at 0/2pi
  // order based on phi and then take the minimum of the difference
  // and the difference with 2pi added to the smallest
  
  float deltaPhi = highEClust.phi() - lowEClust.phi();
  float deltaPhiAlt = ( highEClust.phi() > lowEClust.phi() ? 
                       highEClust.phi() - lowEClust.phi() - 2*PI :
                       lowEClust.phi() - highEClust.phi() - 2*PI );

  deltaPhi = min( fabs( deltaPhi ), fabs( deltaPhiAlt ) );
  
  float sigPhi = deltaPhi / 
  sqrt( highEClust.sigPhi() * highEClust.sigPhi() +
        lowEClust.sigPhi()  * lowEClust.sigPhi() );
  
  // right now this is a rectangular cut, but it could be circular
  // or elliptical
  
  return( ( sigTheta < m_mergeSig ) && ( sigPhi < m_mergeSig ) );  
}

bool
DBCALCluster_factory::overlap( const DBCALCluster& clust,
                               const DBCALPoint* point ) const {
  
  float deltaTheta = fabs( clust.theta() - point->theta() );
  float sigTheta = deltaTheta / sqrt( clust.sigTheta() * clust.sigTheta() +
                                      point->sigTheta()  * point->sigTheta() );
 
  // difference in phi is tricky due to overlap at 0/2pi
  // order based on phi and then take the minimum of the difference
  // and the difference with 2pi added to the smallest
  
  float deltaPhi = clust.phi() - point->phi();
  float deltaPhiAlt = ( clust.phi() > point->phi() ? 
                        clust.phi() - point->phi() - 2*PI :
                        point->phi() - clust.phi() - 2*PI );

  deltaPhi = min( fabs( deltaPhi ), fabs( deltaPhiAlt ) );
  
  float sigPhi = deltaPhi / 
  sqrt( clust.sigPhi() * clust.sigPhi() +
       point->sigPhi()  * point->sigPhi() );
  
  float rho = ( clust.rho() + point->rho() ) / 2;
  float theta = ( clust.theta() + point->theta() ) / 2;
  
  float sep = sqrt( ( rho * deltaTheta ) * ( rho * deltaTheta ) +
      ( rho * sin( theta ) * deltaPhi ) * ( rho * sin( theta ) * deltaPhi ) );
  
#ifdef BCAL_CLUSTER_DIAGNOSTIC
  
  m_dPhi = deltaPhi;
  m_dThe = fabs( clust.theta() - point->theta() );
  m_sep  = sep;
  m_sigPhi = sigPhi;
  m_sigThe = sigTheta;
  m_eClus = clust.E();
  m_rhoClus = clust.rho();
  m_phiClus = clust.phi();
  m_theClus = clust.theta();
  m_nClClus = clust.nCells();
  
  m_ovrlpTr->Fill();
  
#endif // BCAL_CLUSTER_DIAGNOSTIC
  
  
  if( point->E() / clust.E() < 0.1 ){
    
    return ( sep < 5 * m_moliereRadius );
  }
  else{
    
    return ( sep < 2 * m_moliereRadius );
  }

  // right now this is a rectangular cut, but it could be circular
  // or elliptical

//  return( ( sigTheta < m_mergeSig ) && ( sigPhi < m_mergeSig ) );
}

bool
DBCALCluster_factory::overlap( const DBCALCluster& clust,
                               const DBCALHit* hit ) const {
             
  int cellId = DBCALGeometry::cellId( hit->module, hit->layer, hit->sector );

  float cellPhi = DBCALGeometry::phi( cellId );
  float cellSigPhi = DBCALGeometry::phiSize( cellId );
             
  // annoying +- 2pi business to try to find the best delta phi
               
  float deltaPhi = clust.phi() - cellPhi;
  float deltaPhiAlt = ( clust.phi() > cellPhi ? 
                        clust.phi() - cellPhi - 2*PI :
                        cellPhi - clust.phi() - 2*PI );        
  deltaPhi = min( fabs( deltaPhi ), fabs( deltaPhiAlt ) );
               
  float sigPhi = deltaPhi / 
       sqrt( clust.sigPhi() * clust.sigPhi() + cellSigPhi  * cellSigPhi );
             
  return( sigPhi < m_mergeSig );
}
           
           
void
DBCALCluster_factory::clearPoints() {
 
  for( vector< DBCALPoint* >::iterator pt = m_bcalPoints.begin();
      pt != m_bcalPoints.end();
      ++pt ){
    
    delete (*pt);
  }
  
  m_bcalPoints.clear();
}

