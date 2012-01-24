//most of the code was orginally written by Matt Shepherd in
//DBCALCluster_factory.cc

#include <vector>

using namespace std;

#include "BCAL/DBCALPoint_factory.h"
#include "BCAL/DBCALHit.h"

#include "units.h"

jerror_t DBCALPoint_factory::evnt(JEventLoop *loop, int eventnumber) {

  vector< const DBCALHit* > hits;
  loop->Get( hits );
  if ( hits.size() <= 0 ) return NOERROR;

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

  for( map< int, vector< const DBCALHit* > >::iterator mapItr = cellHitMap.begin();
       mapItr != cellHitMap.end();
       ++mapItr ){
    
    if( mapItr->second.size() == 2 && 
        ( mapItr->second[0]->end != mapItr->second[1]->end ) ){
      
      // start with the good stuff -- one hit on each end of a cell

      DBCALPoint *point = new DBCALPoint( *(mapItr->second[0]),
                                          *(mapItr->second[1]) );

      point->AddAssociatedObject(mapItr->second[0]);
      point->AddAssociatedObject(mapItr->second[1]);

      _data.push_back(point);
    }
  }

  //Possibly we should also construct points from single-ended hits here.
  //The code for this is currently (commented out) in
  //DBCALCluster_factory.cc

  return NOERROR;

}
