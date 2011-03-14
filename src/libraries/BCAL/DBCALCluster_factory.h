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

class DBCALCluster_factory : public JFactory< DBCALCluster > {
  
public:
  
  // right now mergeSig is the only parameter that is used to control
  // merging of clusters and points
  
  DBCALCluster_factory() : m_mergeSig( 15 ){}
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
  
};

#endif 

