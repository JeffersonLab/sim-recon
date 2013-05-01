//most of the code was orginally written by Matt Shepherd in
//DBCALCluster_factory.cc

#include <vector>
using namespace std;

#include <JANA/JApplication.h>
using namespace jana;

#include "BCAL/DBCALPoint_factory_NEWSMEAR.h"
#include "BCAL/DBCALHit.h"

#include "units.h"

//----------------
// evnt
//----------------
jerror_t DBCALPoint_factory_NEWSMEAR::evnt(JEventLoop *loop, int eventnumber) {
  vector<const DBCALUnifiedHit*> hits;
  loop->Get(hits);
  if (hits.size() <= 0) return NOERROR;

  // first arrange the list of hits so they are grouped by cell
  map< int, cellHits > cellHitMap;
  for( vector< const DBCALUnifiedHit* >::const_iterator hitPtr = hits.begin();
       hitPtr != hits.end();
       ++hitPtr ){
    
    const DBCALUnifiedHit& hit = (**hitPtr);
    
    // mcsmear will produce hits with energy zero if the hits do not
    // exceed the threshold -- we want to suppress these hits so we know
    // exactly how many "real" hits we have in a cell
    
    if( hit.E < 0.1*k_MeV ) continue;
    
    int id = DBCALGeometry::cellId( hit.module, hit.layer, hit.sector );
    if( cellHitMap.find( id ) == cellHitMap.end() ){
      cellHitMap[id] = cellHits();
    }

    // Add hit to appropriate list for this cell
    if(hit.end == DBCALGeometry::kUpstream){
      cellHitMap[id].uphits.push_back( *hitPtr );
    }else{
      cellHitMap[id].dnhits.push_back( *hitPtr );
    }
  }

  // now go through this list and group hits into BCAL points
  // this combines information from both ends
  for( map< int, cellHits >::const_iterator mapItr = cellHitMap.begin();
       mapItr != cellHitMap.end();
       ++mapItr ){
    
    const vector<const DBCALUnifiedHit*> &uphits = mapItr->second.uphits;
    const vector<const DBCALUnifiedHit*> &dnhits = mapItr->second.dnhits;

    //require double-ended hits
    if(uphits.size()==0 || dnhits.size()==0) continue;

    // Each SiPM sum can have multiple hits, some caused purely by
    // dark hits. A more sophisticated algorithm may be needed here
    // to decipher the multi-hit events. For now, we just take the
    // most energetic hit from each end. (Single ended hits are
    // ignored.

    const DBCALUnifiedHit *uphit=uphits[0];
    const DBCALUnifiedHit *dnhit=dnhits[0];

    for(unsigned int i=1; i<uphits.size(); i++){
      if(uphits[i]->E > uphit->E) uphit = uphits[i];
    }

    for(unsigned int i=1; i<dnhits.size(); i++){
      if(dnhits[i]->E > dnhit->E) dnhit = dnhits[i];
    }

    // first check that the hits don't have absurd timing information

    float fibLen = DBCALGeometry::BCALFIBERLENGTH;
    float cEff = DBCALGeometry::C_EFFECTIVE;

    // get the position with respect to the center of the module -- positive
    // z in the downstream direction
    double zLocal = 0.5 * cEff * ( uphit->t - dnhit->t );

    // if the timing information indicates that the z position is more than 80 cm outside the BCAL, likely the hit is contamined by noise or entirely noise, skip this cell
    double tol = 80*k_cm;

    if (zLocal > (0.5*fibLen + tol) || zLocal < (-0.5*fibLen - tol)) continue;

    DBCALPoint *point = new DBCALPoint(*uphit,*dnhit);

    point->AddAssociatedObject(uphit);
    point->AddAssociatedObject(dnhit);

    _data.push_back(point);
  }

  //Possibly we should also construct points from single-ended hits here.
  //The code for this is currently (commented out) in
  //DBCALCluster_factory.cc

  return NOERROR;
}
