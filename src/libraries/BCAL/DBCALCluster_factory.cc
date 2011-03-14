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
    
    int id = DBCALGeometry::cellId( hit.module, hit.layer, hit.sector );
    
    if( cellHitMap.find( id ) == cellHitMap.end() ){
      
      cellHitMap[id] = vector< const DBCALHit* >();
    }
    
    cellHitMap[id].push_back( *hitPtr );
  }
  
  // now go through this list and group hits into BCAL points
  // this combines information from both ends -- while doing so
  // remove the "consumed" cells from the map
  
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

  /*
  // debugging print of points:
  cout << "BCAL Points: \tE\tphi\tdphi\ttheta\tdtheta" << endl;
  for( vector< const DBCALPoint* >::iterator pt = twoEndPoint.begin();
       pt != twoEndPoint.end();
       ++pt ){
    
    cout << "          \t" << (**pt).E() << "\t" << (**pt).phi() << "\t" << (**pt).sigPhi()
         << "\t" << (**pt).theta() << "\t" << (**pt).sigTheta() << endl;
  }
  */
  
  // now try to clusterize the points
  
  vector< DBCALCluster > clusters = clusterize( twoEndPoint );
  
  // now we should try to add on single hits...
  //
  // MRS note:  14-Mar-11 -- don't seem to see any single ended hits
  // maybe problem in generation of DBCALHit?
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
  
  return( ( sigTheta < m_mergeSig ) && ( sigPhi < m_mergeSig ) );  
}

bool
DBCALCluster_factory::overlap( const DBCALCluster& clust,
                               const DBCALPoint* point ) const {
    
  float sigTheta = fabs( clust.theta() - point->theta() ) / 
  sqrt( clust.sigTheta() * clust.sigTheta() +
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
  
  return( ( sigTheta < m_mergeSig ) && ( sigPhi < m_mergeSig ) );
}

bool
DBCALCluster_factory::overlap( const DBCALCluster& clust,
                               const DBCALHit* hit ) const {
             
  int cellId = DBCALGeometry::cellId( hit->module, hit->layer, hit->sector );

  float cellPhi = DBCALGeometry::phi( cellId );
  float cellSigPhi = DBCALGeometry::phiSize( cellId ) / sqrt( 12 );
             
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

