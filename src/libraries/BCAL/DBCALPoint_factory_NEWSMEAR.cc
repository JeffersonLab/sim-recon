//most of the code was orginally written by Matt Shepherd in
//DBCALCluster_factory.cc

#include <unistd.h>

#include <vector>
using namespace std;

#include <JANA/JApplication.h>
using namespace jana;

#include "BCAL/DBCALPoint_factory_NEWSMEAR.h"
#include "BCAL/DBCALHit.h"

#include "units.h"

//----------------
// init
//----------------
jerror_t DBCALPoint_factory_NEWSMEAR::init(void)
{

	bcal_points_tree = new TTree("bcal_points_tree","");

	bcal_points_tree->Branch("E",&E_tree,"E/F");
	bcal_points_tree->Branch("t_tdc",&t_tdc_tree,"t_tdc/F");
	bcal_points_tree->Branch("t_adc",&t_adc_tree,"t_adc/F");
	bcal_points_tree->Branch("layer",&layer_tree,"layer/I");
	bcal_points_tree->Branch("end",&end_tree,"end/O");

	return NOERROR;
}

jerror_t DBCALPoint_factory_NEWSMEAR::brun(jana::JEventLoop *eventLoop, int runnumber) {

  //get timewalk corrections from CCDB
  /*
  JCalibration *jcalib = eventLoop->GetJCalibration();
  //these tables hold: module layer sector end c0 c1 c2 c3
  vector<vector<float> > adc_timewalk_table;
  vector<vector<float> > tdc_timewalk_table;
  jcalib->Get("BCAL/timewalk_adc",adc_timewalk_table);
  jcalib->Get("BCAL/timewalk_tdc",tdc_timewalk_table);

  for (vector<vector<float> >::const_iterator iter = adc_timewalk_table.begin();
       iter != adc_timewalk_table.end();
       ++iter) {
    if (iter->size() != 8) {
      cout << "DBCALPoint_factory_NEWSMEAR: Wrong number of values in timewalk_adc table (should be 8)" << endl;
      continue;
    }
    //be really careful about float->int converstions
    int module = (int)((*iter)[0]+0.5);
    int layer = (int)((*iter)[1]+0.5);
    int sector = (int)((*iter)[2]+0.5);
    int endi = (int)((*iter)[3]+0.5);
    DBCALGeometry::End end = (endi==0) ? DBCALGeometry::kUpstream : DBCALGeometry::kDownstream;
    float c0 = (*iter)[4];
    float c1 = (*iter)[5];
    float c2 = (*iter)[6];
    float c3 = (*iter)[7];
    int cellId = DBCALGeometry::cellId(module, layer, sector);
    readout_channel channel(cellId,end);
    adc_timewalk_map[channel] = timewalk_coefficients(c0,c1,c2,c3);
  }

  //check that we have entries in the map for all the expected channels
  for (int module=1; module<=DBCALGeometry::NBCALMODS; module++) {
    //shouldn't be hardcoded
    for (int sector=1; sector<=4; sector++) {
      for (int layer=1; layer<=(DBCALGeometry::NBCALLAYSIN + DBCALGeometry::NBCALLAYSOUT); layer++) {
        int id = DBCALGeometry::cellId(module, layer, sector);
        if (adc_timewalk_map.count(readout_channel(id,DBCALGeometry::kUpstream)) != 1) {
          cout << "DBCALPoint_factory_NEWSMEAR: Channel missing in timewalk_adc_table: "
               << endl << " module " << module << " layer " << layer << " sector " << sector << " upstream" << endl;
        }
        if (adc_timewalk_map.count(readout_channel(id,DBCALGeometry::kDownstream)) != 1) {
          cout << "DBCALPoint_factory_NEWSMEAR: Channel missing in timewalk_adc_table: "
               << endl << " module " << module << " layer " << layer << " sector " << sector << " downstream" << endl;
        }
      }
    }
  }

  for (vector<vector<float> >::const_iterator iter = tdc_timewalk_table.begin();
       iter != tdc_timewalk_table.end();
       ++iter) {
    if (iter->size() != 8) {
      cout << "DBCALPoint_factory_NEWSMEAR: Wrong number of values in timewalk_tdc table (should be 8)" << endl;
      continue;
    }
    //be really careful about float->int conversions
    int module = (int)((*iter)[0]+0.5);
    int layer = (int)((*iter)[1]+0.5);
    int sector = (int)((*iter)[2]+0.5);
    int endi = (int)((*iter)[3]+0.5);
    DBCALGeometry::End end = (endi==0) ? DBCALGeometry::kUpstream : DBCALGeometry::kDownstream;
    float c0 = (*iter)[4];
    float c1 = (*iter)[5];
    float c2 = (*iter)[6];
    float c3 = (*iter)[7];
    int cellId = DBCALGeometry::cellId(module, layer, sector);
    readout_channel channel(cellId,end);
    tdc_timewalk_map[channel] = timewalk_coefficients(c0,c1,c2,c3);
  }

  for (int module=1; module<=DBCALGeometry::NBCALMODS; module++) {
    //shouldn't be hardcoded
    for (int sector=1; sector<=4; sector++) {
      for (int layer=1; layer<=DBCALGeometry::NBCALLAYSIN; layer++) {
        int id = DBCALGeometry::cellId(module, layer, sector);
        if (tdc_timewalk_map.count(readout_channel(id,DBCALGeometry::kUpstream)) != 1) {
          cout << "DBCALPoint_factory_NEWSMEAR: Channel missing in timewalk_tdc_table: "
               << endl << " module " << module << " layer " << layer << " sector " << sector << " upstream" << endl;
        }
        if (tdc_timewalk_map.count(readout_channel(id,DBCALGeometry::kDownstream)) != 1) {
          cout << "DBCALPoint_factory_NEWSMEAR: Channel missing in timewalk_tdc_table: "
               << endl << " module " << module << " layer " << layer << " sector " << sector << " downstream" << endl;
        }
      }
    }
  }*/

  //The code commented out above reads in constants from the CCDB.
  //This is the correct thing to do, but until the timewalk correction
  //is more stable, just put the constants in this file.

  //first index labels the layer, the second labels the coefficient (c0,c1,...)
  //only three layers of TDCs!
  const double tdc_timewalk_array[3][4] = { {15.8251, 28.31, 0.394769, 177.54},
                                            {15.3826, 30.9571, 0.374335, 172.578},
                                            {16.1668, 103.024, 0.599359, 138.334} };

  //In reality, we shouldn't have to do a timewalk correction to the ADC times.
  //But currently the ADC time simulated by mcsmear is a threshold crossing
  //time, so this is necessary.
  const double adc_timewalk_array[4][4] = { {17.6245, 279.018, 0.877655, 32.6738},
                                            {17.4629, 398.125, 0.89339, 12.3497},
                                            {17.5743, 1346.56, 1.08083, -22.4011},
                                            {17.6123, 2225.99, 1.15505, -40.5162} };

  for (int module=1; module<=48; module++) {
    for (int layer=1; layer<=3; layer++) {
      for (int sector=1; sector<=4; sector++) {
        int id = DBCALGeometry::cellId(module, layer, sector);
        readout_channel chan_up(id,DBCALGeometry::kUpstream);
        readout_channel chan_dn(id,DBCALGeometry::kDownstream);

        timewalk_coefficients coeffs(tdc_timewalk_array[layer-1][0],
                                     tdc_timewalk_array[layer-1][1],
                                     tdc_timewalk_array[layer-1][2],
                                     tdc_timewalk_array[layer-1][3]);

        tdc_timewalk_map[chan_up] = coeffs;
        tdc_timewalk_map[chan_dn] = coeffs;
      }
    }
  }

  for (int module=1; module<=48; module++) {
    for (int layer=1; layer<=4; layer++) {
      for (int sector=1; sector<=4; sector++) {
        int id = DBCALGeometry::cellId(module, layer, sector);
        readout_channel chan_up(id,DBCALGeometry::kUpstream);
        readout_channel chan_dn(id,DBCALGeometry::kDownstream);

        timewalk_coefficients coeffs(adc_timewalk_array[layer-1][0],
                                     adc_timewalk_array[layer-1][1],
                                     adc_timewalk_array[layer-1][2],
                                     adc_timewalk_array[layer-1][3]);

        adc_timewalk_map[chan_up] = coeffs;
        adc_timewalk_map[chan_dn] = coeffs;
      }
    }
  }

  return NOERROR;
}

//----------------
// evnt
//----------------
jerror_t DBCALPoint_factory_NEWSMEAR::evnt(JEventLoop *loop, int eventnumber) {
  vector< const DBCALHit* > hits;
  vector<const DBCALTDCHit*> tdc_hits;
  loop->Get( hits );
  loop->Get(tdc_hits);
  if (hits.size() + tdc_hits.size() <= 0) return NOERROR;

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

  //add TDC hits to the same structure
  for(vector<const DBCALTDCHit*>::const_iterator hitPtr = tdc_hits.begin();
      hitPtr != tdc_hits.end();
      ++hitPtr ){
    
    const DBCALTDCHit& hit = (**hitPtr);
    
    int id = DBCALGeometry::cellId (hit.module, hit.layer, hit.sector );

    if( cellHitMap.find(id) == cellHitMap.end() ){
      cellHitMap[id] = cellHits();
    }

    // Add hit to appropriate list for this cell
    if(hit.end == DBCALGeometry::kUpstream){
      cellHitMap[id].tdc_uphits.push_back( *hitPtr );
    }else{
      cellHitMap[id].tdc_dnhits.push_back( *hitPtr );
    }
  }

  // now go through this list and group hits into BCAL points
  // this combines information from both ends
  for( map< int, cellHits >::iterator mapItr = cellHitMap.begin();
       mapItr != cellHitMap.end();
       ++mapItr ){
    
    int cellId = mapItr->first;
    int module = DBCALGeometry::module(cellId);
    int layer = DBCALGeometry::layer(cellId);
    int sector = DBCALGeometry::sector(cellId);

    vector<const DBCALHit *> &uphits = mapItr->second.uphits;
    vector<const DBCALHit *> &dnhits = mapItr->second.dnhits;
    vector<const DBCALTDCHit*> &tdc_uphits = mapItr->second.tdc_uphits;
    vector<const DBCALTDCHit*> &tdc_dnhits = mapItr->second.tdc_dnhits;

    if (uphits.size()==0 && tdc_uphits.size()!=0) {
      cout << "DBCALPoint_factory_NEWSMEAR (event " << eventnumber << "): upstream TDC hits without upstream ADC hits" << endl;
    }
    if (dnhits.size()==0 && tdc_dnhits.size()!=0) {
      cout << "DBCALPoint_factory_NEWSMEAR (event " << eventnumber << "): downstream TDC hits without downstream ADC hits" << endl;
    }

    //if we have insufficient ADC hits, there is nothing to do with the TDC hits either (require double-ended hits)
    if(uphits.size()==0 || dnhits.size()==0) continue;

    // Each SiPM sum can have multiple hits, some caused purely by
    // dark hits. A more sophisticated algorithm may be needed here
    // to decipher the multi-hit events. For now, we just take the
    // most energetic hit from each end. (Single ended hits are
    // ignored.

    const DBCALHit *uphit=uphits[0];
    const DBCALHit *dnhit=dnhits[0];

    for(unsigned int i=1; i<uphits.size(); i++){
      if(uphits[i]->E > uphit->E) uphit = uphits[i];
    }

    for(unsigned int i=1; i<dnhits.size(); i++){
      if(dnhits[i]->E > dnhit->E) dnhit = dnhits[i];
    }

    //keep track of if we have TDC hits that need to be added as associated objects
    bool hasTDCuphit = tdc_uphits.size();
    bool hasTDCdnhit = tdc_dnhits.size();

    //find the TDC hits closest (in time) to the ADC hits
    //this is probably not the most optimal thing to do, as no timewalk corrections have been done yet
    const DBCALTDCHit *tdc_uphit=NULL;
    if (hasTDCuphit) {
      tdc_uphit = tdc_uphits[0];
      for (unsigned int i=1; i<tdc_uphits.size(); i++) {
        if ( fabs(tdc_uphits[i]->t - uphit->t) <
             fabs(tdc_uphit->t - uphit->t) ) {
          tdc_uphit = tdc_uphits[i];
        }
      }
    }

    const DBCALTDCHit *tdc_dnhit=NULL;
    if (hasTDCdnhit) {
      tdc_dnhit = tdc_dnhits[0];
      for (unsigned int i=1; i<tdc_dnhits.size(); i++) {
        if ( fabs(tdc_dnhits[i]->t - dnhit->t) <
             fabs(tdc_dnhit->t - dnhit->t) ) {
          tdc_dnhit = tdc_dnhits[i];
        }
      }
    }

    E_tree = uphit->E;
    t_adc_tree = uphit->t;
    t_tdc_tree = 0;
    layer_tree = uphit->layer;
    end_tree = 0;
    if (hasTDCuphit) t_tdc_tree = tdc_uphit->t;
    bcal_points_tree->Fill();

    E_tree = dnhit->E;
    t_adc_tree = dnhit->t;
    t_tdc_tree = 0;
    layer_tree = dnhit->layer;
    end_tree = 1;
    if (hasTDCdnhit) t_tdc_tree = tdc_dnhit->t;
    bcal_points_tree->Fill();

    float EUp, EDown, tUp, tDown, tUp_ADC, tDown_ADC; //these are values that will be passed to the DBCALPoint constructor

    EUp = uphit->E;
    readout_channel up_channel(cellId,DBCALGeometry::kUpstream);
    timewalk_coefficients up_coeff = adc_timewalk_map[up_channel];
    tUp_ADC = uphit->t;
    tUp_ADC -= up_coeff.c0 + up_coeff.c1/pow(EUp-up_coeff.c3,up_coeff.c2);
    if (hasTDCuphit) {
      timewalk_coefficients tdc_up_coeff = tdc_timewalk_map[up_channel];
      tUp = tdc_uphit->t;
      tUp -= tdc_up_coeff.c0 + tdc_up_coeff.c1/pow(EUp-tdc_up_coeff.c3,tdc_up_coeff.c2);
    } else {
      tUp = tUp_ADC;
    }

    EDown = dnhit->E;
    readout_channel dn_channel(cellId,DBCALGeometry::kDownstream);
    timewalk_coefficients dn_coeff = adc_timewalk_map[dn_channel];
    tDown_ADC = dnhit->t;
    tDown_ADC -= dn_coeff.c0 + dn_coeff.c1/pow(EDown-dn_coeff.c3,dn_coeff.c2);
    if (hasTDCdnhit) {
      timewalk_coefficients tdc_dn_coeff = tdc_timewalk_map[dn_channel];
      tDown = tdc_dnhit->t;
      tDown -= tdc_dn_coeff.c0 + tdc_dn_coeff.c1/pow(EDown-tdc_dn_coeff.c3,tdc_dn_coeff.c2);
    } else {
      tDown = tDown_ADC;
    }

    const float fADC_counts_per_GeV = 15800.0;
    EUp /= fADC_counts_per_GeV;
    EDown /= fADC_counts_per_GeV;

    // start with the good stuff -- one hit on each end of a cell

    // first check that the hits don't have absurd timing information

    float fibLen = DBCALGeometry::BCALFIBERLENGTH;
    float cEff = DBCALGeometry::C_EFFECTIVE;

    // get the position with respect to the center of the module -- positive
    // z in the downstream direction
    double zLocal = 0.5 * cEff * ( tUp - tDown );

    // if the timing information indicates that the z position is more than 80 cm outside the BCAL, likely the hit is contamined by noise or entirely noise, skip this cell
    double tol = 80*k_cm;

    if (zLocal > (0.5*fibLen + tol) || zLocal < (-0.5*fibLen - tol)) continue;
    DBCALPoint *point = new DBCALPoint(module, layer, sector, EUp, EDown, tUp, tDown, tUp_ADC, tDown_ADC);

    point->AddAssociatedObject(uphit);
    point->AddAssociatedObject(dnhit);

    if (hasTDCuphit) point->AddAssociatedObject(tdc_uphit);
    if (hasTDCdnhit) point->AddAssociatedObject(tdc_dnhit);

    _data.push_back(point);
  }

  //Possibly we should also construct points from single-ended hits here.
  //The code for this is currently (commented out) in
  //DBCALCluster_factory.cc

  return NOERROR;
}
