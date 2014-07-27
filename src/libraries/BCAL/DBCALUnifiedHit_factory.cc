#include <vector>
using namespace std;

#include <JANA/JApplication.h>
using namespace jana;

#include "BCAL/DBCALUnifiedHit_factory.h"

#include "units.h"

//----------------
// init
//----------------
jerror_t DBCALUnifiedHit_factory::init(void)
{

  if (enable_debug_output) {
    bcal_points_tree = new TTree("bcal_points_tree","");
    bcal_points_tree->Branch("E",&E_tree,"E/F");
    bcal_points_tree->Branch("t_tdc",&t_tdc_tree,"t_tdc/F");
    bcal_points_tree->Branch("t_adc",&t_adc_tree,"t_adc/F");
    bcal_points_tree->Branch("t_tdc_corrected",&t_tdc_corrected_tree,"t_tdc_corrected/F");
    bcal_points_tree->Branch("t_adc_corrected",&t_adc_corrected_tree,"t_adc_corrected/F");
    bcal_points_tree->Branch("layer",&layer_tree,"layer/I");
    bcal_points_tree->Branch("end",&end_tree,"end/O");
  }

  return NOERROR;
}

jerror_t DBCALUnifiedHit_factory::brun(jana::JEventLoop *eventLoop, int runnumber) {

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
      cout << "DBCALUnifiedHit_factory: Wrong number of values in timewalk_adc table (should be 8)" << endl;
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
          cout << "DBCALUnifiedHit_factory: Channel missing in timewalk_adc_table: "
               << endl << " module " << module << " layer " << layer << " sector " << sector << " upstream" << endl;
        }
        if (adc_timewalk_map.count(readout_channel(id,DBCALGeometry::kDownstream)) != 1) {
          cout << "DBCALUnifiedHit_factory: Channel missing in timewalk_adc_table: "
               << endl << " module " << module << " layer " << layer << " sector " << sector << " downstream" << endl;
        }
      }
    }
  }

  for (vector<vector<float> >::const_iterator iter = tdc_timewalk_table.begin();
       iter != tdc_timewalk_table.end();
       ++iter) {
    if (iter->size() != 8) {
      cout << "DBCALUnifiedHit_factory: Wrong number of values in timewalk_tdc table (should be 8)" << endl;
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
          cout << "DBCALUnifiedHit_factory: Channel missing in timewalk_tdc_table: "
               << endl << " module " << module << " layer " << layer << " sector " << sector << " upstream" << endl;
        }
        if (tdc_timewalk_map.count(readout_channel(id,DBCALGeometry::kDownstream)) != 1) {
          cout << "DBCALUnifiedHit_factory: Channel missing in timewalk_tdc_table: "
               << endl << " module " << module << " layer " << layer << " sector " << sector << " downstream" << endl;
        }
      }
    }
  }*/

  //The code commented out above reads in constants from the CCDB.
  //This is the correct thing to do, but until the timewalk correction
  //is more stable, just put the constants in this file.

  //for now, since we are only dealing with simulated data, assume all
  //sectors and modules are identical and timewalk coefficients only depend
  //on layer

  //first index labels the layer, the second labels the coefficient (c0,c1,...)
  //only three layers of TDCs!
  const double tdc_timewalk_array[3][4] = {{16.0081, 0.365565, 0.54665, 0.00640507},
                                           {15.2826, 0.835660, 0.374335, 0.0111243},
                                           {15.6251, 0.530827, 0.475919, 0.0102508} };

  //In reality, we shouldn't have to do a timewalk correction to the ADC times.
  //But currently the ADC time simulated by mcsmear is a threshold crossing
  //time, so this is necessary.
  const double adc_timewalk_array[4][4] = { {17.6245, 0.0585642, 0.877655, 0.00210614},
                                            {17.4629, 0.0717924, 0.89339, 0.000796056},
                                            {17.5743, 0.0397904, 1.08083, -0.00144397},
                                            {17.6123, 0.0321396, 1.15505, -0.00261166} };

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
jerror_t DBCALUnifiedHit_factory::evnt(JEventLoop *loop, int eventnumber) {
  vector<const DBCALHit*> hits;
  vector<const DBCALTDCHit*> tdc_hits;
  loop->Get(hits);
  loop->Get(tdc_hits);
  if (hits.size() + tdc_hits.size() <= 0) return NOERROR;

  // first arrange the list of hits so they are grouped by cell
  map<readout_channel, cellHits> cellHitMap;
  for( vector<const DBCALHit*>::const_iterator hitPtr = hits.begin();
       hitPtr != hits.end();
       ++hitPtr ){
    
    const DBCALHit& hit = (**hitPtr);
    
    int id = DBCALGeometry::cellId( hit.module, hit.layer, hit.sector );
    readout_channel chan(id, hit.end);

    //this will create cellHitMap[chan] if it doesn't already exist
    cellHitMap[chan].hits.push_back(*hitPtr);
  }

  //add TDC hits to the same structure
  for(vector<const DBCALTDCHit*>::const_iterator hitPtr = tdc_hits.begin();
      hitPtr != tdc_hits.end();
      ++hitPtr ){
    
    const DBCALTDCHit& hit = (**hitPtr);
    
    int id = DBCALGeometry::cellId (hit.module, hit.layer, hit.sector );
    readout_channel chan(id, hit.end);

    if( cellHitMap.find(chan) == cellHitMap.end() ){
      cellHitMap[chan] = cellHits();
    }

    cellHitMap[chan].tdc_hits.push_back(*hitPtr);
  }

  // now go through this list and incorporate TDC hits where available
  for(map<readout_channel,cellHits>::const_iterator mapItr = cellHitMap.begin();
       mapItr != cellHitMap.end();
       ++mapItr) {
    
    readout_channel chan = mapItr->first;
    int cellId = chan.cellId;
    int module = DBCALGeometry::module(cellId);
    int layer = DBCALGeometry::layer(cellId);
    int sector = DBCALGeometry::sector(cellId);

    const vector<const DBCALHit*> &hits = mapItr->second.hits;
    const vector<const DBCALTDCHit*> &tdc_hits = mapItr->second.tdc_hits;

    //if we have no ADC hits, there is nothing to do with the TDC hits either
    if (hits.size()==0) {
      cout << "DBCALUnifiedHit_factory (event " << eventnumber << "): TDC hits without ADC hits" << endl;
      continue;
    }

    // Each SiPM sum can have multiple hits, some caused purely by
    // dark hits. There is a question of how to handle this properly.
    // For now, we ignore TDC info if there are >1 TDC hits per channel.
    // If there are multiple ADC hits, but only one TDC hit, assume that
    // the TDC hit belongs to the highest energy ADC hit. There is probably
    // a better way to do this.

    //Find the index of the highest energy ADC hit.
    unsigned int highEIndex = 0;
    for(unsigned int i=1; i<hits.size(); i++){
      if(hits[i]->E > hits[highEIndex]->E) highEIndex = i;
    }

    //ignore TDC hits unless we have exactly one
    bool hasOneTDCHit = (tdc_hits.size()==1);
    const DBCALTDCHit* tdc_hit=NULL;
    if (hasOneTDCHit) tdc_hit = tdc_hits[0];

    for (unsigned int i=0; i<hits.size(); i++) {
      const DBCALHit* hit=hits[i];

      bool useTDChit = (i==highEIndex) && (hasOneTDCHit);

      if (enable_debug_output) {
        E_tree = hit->E;
        t_adc_tree = hit->t;
        t_tdc_tree = 0;
        layer_tree = hit->layer;
        end_tree = hit->end;
        if (useTDChit) t_tdc_tree = tdc_hit->t;
      }

      float E, t, t_ADC; //these are values that will be assigned to the DBCALUnifiedHit

      E = hit->E;
      timewalk_coefficients coeff = adc_timewalk_map[chan];
      t_ADC = hit->t;

      //In reality, we shouldn't have to do a timewalk correction to the ADC times.
      //But currently the ADC time simulated by mcsmear is a threshold crossing
      //time, so this is necessary.
      //if E < coeff.c3, the correction will be bogus, just skip it (this shouldn't happen)
      if (E > coeff.c3) {
        t_ADC -= coeff.c0 + coeff.c1/pow(E-coeff.c3,coeff.c2);
      }
      if (useTDChit) {
        timewalk_coefficients tdc_coeff = tdc_timewalk_map[chan];
        t = tdc_hit->t;
        //If E < tdc_coeff.c3, this means the energy is low enough we would not
        //expect a TDC hit, so just ignore the TDC hit. Otherwise, use the
        //timewalk-corrected time as the hit time.
        if (E > tdc_coeff.c3) {
          t -= tdc_coeff.c0 + tdc_coeff.c1/pow(E-tdc_coeff.c3,tdc_coeff.c2);
        } else {
          t = t_ADC;
          useTDChit = false;
        }
      } else {
        t = t_ADC;
      }

      if (enable_debug_output) {
        t_adc_corrected_tree = t_ADC;
        t_tdc_corrected_tree = useTDChit ? t : 0;
        bcal_points_tree->Fill();
      }

      DBCALUnifiedHit *uhit = new DBCALUnifiedHit;
      uhit->E = E;
      uhit->t = t;
      uhit->t_ADC = t_ADC;
      uhit->has_TDC_hit = useTDChit;
      uhit->cellId = cellId;
      uhit->module = module;
      uhit->layer = layer;
      uhit->sector = sector;
      uhit->end = chan.end;

      uhit->AddAssociatedObject(hit);
      if (useTDChit) uhit->AddAssociatedObject(tdc_hit);

      _data.push_back(uhit);
    }
  }

  return NOERROR;
}
