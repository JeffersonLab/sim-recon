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

    USE_TDC = false;
    if (gPARMS){
        gPARMS->SetDefaultParameter("BCAL:USE_TDC", USE_TDC, "Set to 1 to use TDC times");
    }

    static pthread_mutex_t print_mutex = PTHREAD_MUTEX_INITIALIZER;
    static bool printMessage = true;
    pthread_mutex_lock(&print_mutex);
    if (printMessage){
        if (USE_TDC){
            jout << "DBCALUnifiedHit_factory: Using TDC times when available." << endl;
        }
        else{
            jout << "DBCALUnifiedHit_factory: Using ADC times only." << endl;
        }
    }
    printMessage = false;
    pthread_mutex_unlock(&print_mutex);

    return NOERROR;
}

jerror_t DBCALUnifiedHit_factory::brun(jana::JEventLoop *eventLoop, int32_t runnumber) {

	// load BCAL geometry
  	vector<const DBCALGeometry *> BCALGeomVec;
  	eventLoop->Get(BCALGeomVec);
  	if(BCALGeomVec.size() == 0)
		throw JException("Could not load DBCALGeometry object!");
	dBCALGeom = BCALGeomVec[0];

    //get timewalk corrections from CCDB
    JCalibration *jcalib = eventLoop->GetJCalibration();
    //these tables hold: module layer sector end c0 c1 c2 c3 threshold
    vector<vector<float> > tdc_timewalk_table;
    
    // jcalib->Get("BCAL/timewalk_tdc",tdc_timewalk_table);
    if(jcalib->Get("BCAL/timewalk_tdc_c4",tdc_timewalk_table)) {
        jerr << "Error loading BCAL/timewalk_tdc !" << endl;
    } else {
        for (vector<vector<float> >::const_iterator iter = tdc_timewalk_table.begin();
             iter != tdc_timewalk_table.end();
             ++iter) {
            //if (iter->size() != 8) {
            //  cout << "DBCALUnifiedHit_factory: Wrong number of values in timewalk_tdc table (should be 8)" << endl;
            if (iter->size() != 9) {
                cout << "DBCALUnifiedHit_factory: Wrong number of values in timewalk_tdc table (should be 9)" << endl;
                continue;
            }
            //be really careful about float->int conversions
            int module = (int)((*iter)[0]+0.5);
            int layer  = (int)((*iter)[1]+0.5);
            int sector = (int)((*iter)[2]+0.5);
            int endi   = (int)((*iter)[3]+0.5);
            DBCALGeometry::End end = (endi==0) ? DBCALGeometry::kUpstream : DBCALGeometry::kDownstream;
            float c0 = (*iter)[4];
            float c1 = (*iter)[5];
            float c2 = (*iter)[6];
            float c3 = (*iter)[7];
            float thresh = (*iter)[8];
            int cellId = dBCALGeom->cellId(module, layer, sector);
            readout_channel channel(cellId,end);
            tdc_timewalk_map_c4[channel] = timewalk_coefficients_c4(c0,c1,c2,c3,thresh);
            //tdc_timewalk_map[channel] = timewalk_coefficients(c0,c1,c2,a_thresh);
        }
        
        for (int module=1; module<=dBCALGeom->GetBCAL_Nmodules(); module++) {
            //shouldn't be hardcoded
            for (int sector=1; sector<=4; sector++) {
                for (int layer=1; layer<=dBCALGeom->GetBCAL_NInnerLayers(); layer++) {
                    int id = dBCALGeom->cellId(module, layer, sector);
                    //if (tdc_timewalk_map.count(readout_channel(id,dBCALGeom->kUpstream)) != 1) {
                    if (tdc_timewalk_map_c4.count(readout_channel(id,dBCALGeom->kUpstream)) != 1) {
                        cout << "DBCALUnifiedHit_factory: Channel missing in timewalk_tdc_table: "
                             << endl << " module " << module << " layer " << layer << " sector " << sector << " upstream" << endl;
                    }
                    //if (tdc_timewalk_map.count(readout_channel(id,dBCALGeom->kDownstream)) != 1) {
                    if (tdc_timewalk_map_c4.count(readout_channel(id,dBCALGeom->kDownstream)) != 1) {
                        cout << "DBCALUnifiedHit_factory: Channel missing in timewalk_tdc_table: "
                             << endl << " module " << module << " layer " << layer << " sector " << sector << " downstream" << endl;
                    }
                }
            }
        }
    }
    
    return NOERROR;
}

//----------------
// evnt
//----------------
jerror_t DBCALUnifiedHit_factory::evnt(JEventLoop *loop, uint64_t eventnumber) {
    vector<const DBCALHit*> hits;
    vector<const DBCALTDCHit*> tdc_hits;
    loop->Get(hits);
    loop->Get(tdc_hits);
    if (hits.size() + tdc_hits.size() <= 0) return NOERROR;
	if (VERBOSE>=3)  jout << eventnumber << " " << hits.size() << " ADC hits, " << tdc_hits.size() << " TDC hits" << endl;

    // first arrange the list of hits so they are grouped by cell
    map<readout_channel, cellHits> cellHitMap;
    for( vector<const DBCALHit*>::const_iterator hitPtr = hits.begin();
            hitPtr != hits.end();
            ++hitPtr ){

        const DBCALHit& hit = (**hitPtr);

        int id = dBCALGeom->cellId( hit.module, hit.layer, hit.sector );
        readout_channel chan(id, hit.end);

        //this will create cellHitMap[chan] if it doesn't already exist
        cellHitMap[chan].hits.push_back(*hitPtr);
    }

    //add TDC hits to the same structure
    for(vector<const DBCALTDCHit*>::const_iterator hitPtr = tdc_hits.begin();
            hitPtr != tdc_hits.end();
            ++hitPtr ){

        const DBCALTDCHit& hit = (**hitPtr);

        int id = dBCALGeom->cellId (hit.module, hit.layer, hit.sector );
        readout_channel chan(id, hit.end);

        if( cellHitMap.find(chan) == cellHitMap.end() ){
            cellHitMap[chan] = cellHits();
        }

        cellHitMap[chan].tdc_hits.push_back(*hitPtr);
    }

    // now go through this list cell by cell and creat a UnifiedHit if appropriate
    for(map<readout_channel,cellHits>::const_iterator mapItr = cellHitMap.begin();
            mapItr != cellHitMap.end();
            ++mapItr) {

        readout_channel chan = mapItr->first;
        int cellId = chan.cellId;
        int module = dBCALGeom->module(cellId);
        int layer = dBCALGeom->layer(cellId);
        int sector = dBCALGeom->sector(cellId);

        const vector<const DBCALHit*> &hits = mapItr->second.hits;
        const vector<const DBCALTDCHit*> &tdc_hits = mapItr->second.tdc_hits;

        //if we have no ADC hits in the cell, there is nothing to do with the TDC hits either
		if (hits.size()==0) {
			if (VERBOSE>=1) {
				if (tdc_hits.size()!=0) {
					static uint64_t Nwarnings = 0;
					if(++Nwarnings <= 10) cout << "DBCALUnifiedHit_factory (event " << eventnumber << "): TDC hits without ADC hits" << endl;
					if(  Nwarnings == 10) cout << "DBCALUnifiedHit_factory: LAST WARNING (others will be suppressed)" <<endl;
				}
			}
			continue;
        }

        // Find the index of the earliest ADC hit.
        unsigned int firstIndex = 0;
        for(unsigned int i=1; i<hits.size(); i++){
            if (hits[i]->t < hits[firstIndex]->t) firstIndex = i;
			if (VERBOSE>=4 && hits.size()>=2) {
				printf("DBCALUnifiedHit_factory  event %5llu (%2i,%i,%i,%i) peak %4i, ADC %i %6.1f\n",
					   (unsigned long long int)eventnumber, module, layer, sector, chan.end, hits[i]->pulse_peak, i, hits[i]->t);
			}
        }

		const DBCALHit* hit=hits[firstIndex];
		float pulse_peak, E, t, t_ADC, t_TDC=0; //these are values that will be assigned to the DBCALUnifiedHit
		pulse_peak = hit->pulse_peak; 
		E = hit->E;
		t_ADC = hit->t;

		// Loop through the TDC hits, apply timewalk correction and find the TDC time closest to the ADC time
		int goodTDCindex=-1;
		//timewalk_coefficients tdc_coeff = tdc_timewalk_map[chan];
		timewalk_coefficients_c4 tdc_coeff;
        bool good_timewalk_params = false;
        if( tdc_timewalk_map_c4.find(chan) != tdc_timewalk_map_c4.end() ) {
            // really we should probably print out some more errors here, if we can't find the timewalk correction factor
            // but since we complain enough above, it is probably fine...
            tdc_coeff = tdc_timewalk_map_c4[chan];
            good_timewalk_params = true;
        }
		float t_diff=100000;
		bool haveTDChit = false;
		for (unsigned int i=0; i<tdc_hits.size(); i++) {
			haveTDChit = true;
			const DBCALTDCHit* tdc_hit=tdc_hits[i];
			float tdc_hit_t = tdc_hit->t;

            if(good_timewalk_params) {
                //tdc_hit_t -= tdc_coeff.c0 + tdc_coeff.c1/pow(pulse_peak/tdc_coeff.a_thresh, tdc_coeff.c2);
                float shifted_peak = pulse_peak+tdc_coeff.c2; // only apply formula if time is above TDC zero
                if (shifted_peak>0) tdc_hit_t -= tdc_coeff.c0 + tdc_coeff.c1*pow(shifted_peak,tdc_coeff.c3);
                if (VERBOSE>=4) printf("tamewalk %f -> %f: (%f,%f,%f,%f)\n",tdc_hit->t,tdc_hit_t,tdc_coeff.c0,tdc_coeff.c1,tdc_coeff.c2,tdc_coeff.c3);
            }

			float tdc_adc_diff = tdc_hit_t - t_ADC;
			if (fabs(tdc_adc_diff) < fabs(t_diff)) {
				goodTDCindex=i;
				t_diff=tdc_adc_diff;
				t_TDC=tdc_hit_t;
			}
			if (VERBOSE>=4 || (VERBOSE>=3&&tdc_hits.size()>1)) {
				printf("DBCALUnifiedHit_factory  event %5llu (%2i,%i,%i,%i) peak %4.0f, ADC 0 %6.1f, TDC %i %6.1f  diff %6.1f    best %2i %6.1f\n",
					   (unsigned long long int)eventnumber, module, layer, sector, chan.end, pulse_peak, t_ADC, i, tdc_hit_t, tdc_adc_diff, goodTDCindex, t_diff);
			}
		}
		if (VERBOSE>=4 && !haveTDChit) {
			// printf("DBCALUnifiedHit_factory  event %5llu (%2i,%i,%i,%i) peak %4.0f, ADC 0 %6.1f, TDC           diff           best %2i %6.1f\n",
			// 	   (unsigned long long int)eventnumber, module, layer, sector, chan.end, pulse_peak, hit->t, goodTDCindex, t_diff);
			printf("DBCALUnifiedHit_factory  event %5llu (%2i,%i,%i,%i) peak %4.0f, ADC 0 %6.1f\n",
				   (unsigned long long int)eventnumber, module, layer, sector, chan.end, pulse_peak, hit->t);
		}
		// Decide whether to use TDC time
		if (pulse_peak>tdc_coeff.thresh && USE_TDC==1 && goodTDCindex>=0 && fabs(t_diff)<2 && good_timewalk_params) {
			t = t_TDC;
		} else {
			t = t_ADC;
		}

		// Create the DBCALUnifiedHit
		DBCALUnifiedHit *uhit = new DBCALUnifiedHit;
		uhit->E = E;
		uhit->t = t;
		uhit->t_ADC = t_ADC;
		uhit->t_TDC = t_TDC;
		uhit->has_TDC_hit = haveTDChit;
		uhit->cellId = cellId;
		uhit->module = module;
		uhit->layer = layer;
		uhit->sector = sector;
		uhit->end = chan.end;

		uhit->AddAssociatedObject(hit);
		if (goodTDCindex>=0) uhit->AddAssociatedObject(tdc_hits[goodTDCindex]);

		_data.push_back(uhit);
    } // end loop over cells

    return NOERROR;
}
