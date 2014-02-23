//most of the code was orginally written by Matt Shepherd in
//DBCALCluster_factory.cc

#include <unistd.h>

#include <vector>
using namespace std;

#include <JANA/JApplication.h>
using namespace jana;

#include "BCAL/DBCALPoint_factory_OLDSMEAR.h"
#include "BCAL/DBCALHit.h"

#include "units.h"

static pthread_mutex_t deprecated_warning_mutex = PTHREAD_MUTEX_INITIALIZER;
static bool warning_issued = false;

//----------------
// init
//----------------
jerror_t DBCALPoint_factory_OLDSMEAR::init(void)
{
	pthread_mutex_lock(&deprecated_warning_mutex);
	if(!warning_issued){
		stringstream mess;
		mess << "=========================================================================="<<endl;
		mess << endl;
		mess << "W         W      A      RRRRR     NN    N   IIIIIII   NN    N     GGGGGG  "<<endl;
		mess << "W         W     A A     R    R    N N   N      I      N N   N    G        "<<endl;
		mess << "W         W    A   A    R    R    N  N  N      I      N  N  N   G         "<<endl;
		mess << "W    W    W   AAAAAAA   RRRRR     N   N N      I      N   N N   G    GGGGG"<<endl;
		mess << " W  W W  W    A     A   R    R    N    NN      I      N    NN    G     G  "<<endl;
		mess << "  W    W      A     A   R     R   N     N   IIIIIII   N     N     GGGGG   "<<endl;
		mess << "=========================================================================="<<endl;
		mess << endl;
		mess << "Warning! You are using a deprecated version of the BCAL simulation!"<<endl;
		mess << "The reconstruction code has not yet been updated to handle the newer"<<endl;
		mess << "scheme so it is not yet the default. "<<endl;
		mess << endl;
		mess << "HOWEVER! The results you get for the BCAL may not accurately reflect"<<endl;
		mess << "what the detector is capable of. INTERPRET YOUR RESULTS WITH CAUTION!!"<<endl;
		mess << endl;
		mess << "More importantly, bug the calorimetry group to FIX THIS ISSUE!!!"<<endl;
		mess << "=========================================================================="<<endl;
		cerr << mess.str();

		bool save_monitor_heartbeat = japp->monitor_heartbeat;
		japp->monitor_heartbeat = false;
		japp->SetShowTicker(0);
		for(int i=3; i>=0; i--){
			cout<<".... resumimg in "<<i<<" seconds\r";
			cout.flush();
			sleep(1);
		}
		japp->SetShowTicker(1); // may not be correct! (have no way to check ticker setting beforehand!)
		japp->monitor_heartbeat = save_monitor_heartbeat;
		warning_issued = true;
	}
	pthread_mutex_unlock(&deprecated_warning_mutex);

	return NOERROR;
}

//----------------
// evnt
//----------------
jerror_t DBCALPoint_factory_OLDSMEAR::evnt(JEventLoop *loop, int eventnumber) {
  vector< const DBCALHit* > hits;
  loop->Get( hits );
  if ( hits.size() <= 0 ) return NOERROR;

  // first arrange the list of hits so they are grouped by cell
  map< int, cellHits > cellHitMap;
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

  for( map< int, cellHits >::iterator mapItr = cellHitMap.begin();
       mapItr != cellHitMap.end();
       ++mapItr ){
    
		// Each SiPM sum can have multiple hits, some caused purely by
		// dark hits. A more sophisticated algorithm may be needed here
		// to decipher the multi-hit events. For now, we just take the
		// most energetic hit from each end. (Single ended hits are
		// ignored.
		vector<const DBCALHit *> &uphits = mapItr->second.uphits;
		vector<const DBCALHit *> &dnhits = mapItr->second.dnhits;
		if(uphits.size()==0 || dnhits.size()==0) continue;

		const DBCALHit *uphit=uphits[0];
		const DBCALHit *dnhit=dnhits[0];

		for(unsigned int i=1; i<uphits.size(); i++){
			if(uphits[i]->E > uphit->E) uphit = uphits[i];
		}

		for(unsigned int i=1; i<dnhits.size(); i++){
			if(dnhits[i]->E > dnhit->E) dnhit = dnhits[i];
		}
      
      // start with the good stuff -- one hit on each end of a cell

      // first check that the hits don't have absurd timing information

      float fibLen = DBCALGeometry::BCALFIBERLENGTH;
      float cEff = DBCALGeometry::C_EFFECTIVE;

      double tUp = uphit->t;
      double tDown = dnhit->t;
  
      // get the position with respect to the center of the module -- positive
      // z in the downstream direction
  
      double zLocal = 0.5 * cEff * ( tUp - tDown );

      // if the timing information indicates that the z position is more than 50 cm outside the BCAL, likely the hit is contamined by noise or entirely noise, skip this cell
      double tol = 50*k_cm;
      if ( zLocal > (0.5*fibLen + tol) || zLocal < (-0.5*fibLen - tol) ) continue;
  
      DBCALPoint *point = new DBCALPoint( *uphit,
                                          *dnhit );

      point->AddAssociatedObject(uphit);
      point->AddAssociatedObject(dnhit);

      _data.push_back(point);
  }

  //Possibly we should also construct points from single-ended hits here.
  //The code for this is currently (commented out) in
  //DBCALCluster_factory.cc

  return NOERROR;

}
