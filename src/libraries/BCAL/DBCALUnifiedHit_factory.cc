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

    if (USE_TDC){
        jout << "DBCALUnifiedHit_factory: Using TDC times when available." << endl;
    }
    else{
        jout << "DBCALUnifiedHit_factory: Using ADC times only." << endl;
    }

    return NOERROR;
}

jerror_t DBCALUnifiedHit_factory::brun(jana::JEventLoop *eventLoop, int32_t runnumber) {

    if (USE_TDC){
        //get timewalk corrections from CCDB
        JCalibration *jcalib = eventLoop->GetJCalibration();
        //these tables hold: module layer sector end c0 c1 c2 c3
        vector<vector<float> > tdc_timewalk_table;
        jcalib->Get("BCAL/timewalk_tdc",tdc_timewalk_table);

        for (vector<vector<float> >::const_iterator iter = tdc_timewalk_table.begin();
                iter != tdc_timewalk_table.end();
                ++iter) {
            if (iter->size() != 8) {
                cout << "DBCALUnifiedHit_factory: Wrong number of values in timewalk_tdc table (should be 8)" << endl;
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
            float a_thresh = (*iter)[7];
            int cellId = DBCALGeometry::cellId(module, layer, sector);
            readout_channel channel(cellId,end);
            tdc_timewalk_map[channel] = timewalk_coefficients(c0,c1,c2,a_thresh);
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

    // now go through this list cell by cell and creat a UnifiedHit if appropriate
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

        //if we have no ADC hits in the cell, there is nothing to do with the TDC hits either
        if (hits.size()==0) {
            static uint64_t Nwarnings = 0;
            if(++Nwarnings <= 10) cout << "DBCALUnifiedHit_factory (event " << eventnumber << "): TDC hits without ADC hits" << endl;
            if(  Nwarnings == 10) cout << "DBCALUnifiedHit_factory: LAST WARNING (others will be suppressed)" <<endl;
            continue;
        }

        // At the moment we only allow 1 ADC hit in the firmware.
        // Need to revisit handling of multiple ADC events for a hypothetical future.

        // Find the index of the highest energy ADC hit.
        unsigned int highEIndex = 0;
        for(unsigned int i=1; i<hits.size(); i++){
            if(hits[i]->E > hits[highEIndex]->E) highEIndex = i;
        }

        // We trust the ADC hit times.  Choose the TDC time that is closest to the ADC time
        bool haveTDChit = false;
        if (tdc_hits.size() > 0) haveTDChit = true;

        // For each ADC hit 
        for (unsigned int i=0; i<hits.size(); i++) {
            const DBCALHit* hit=hits[i];

            float pulse_peak, E, t, t_ADC, t_TDC=0; //these are values that will be assigned to the DBCALUnifiedHit

            pulse_peak = hit->pulse_peak; 
            E = hit->E;
            t_ADC = hit->t;

            int goodTDCindex=0;
            if (haveTDChit) {
                // Loop through the TDC hits and find the closest to the ADC hit
                float t_diff=100000;
                for (unsigned int i=0; i<tdc_hits.size(); i++) {
                    const DBCALTDCHit* tdc_hit=tdc_hits[i];
                    float tdc_adc_diff = tdc_hit->t - hit->t;
                    if (fabs(tdc_adc_diff) < fabs(t_diff)) {
                        goodTDCindex=i;
                        t_diff=tdc_adc_diff;
                    }
                    if (VERBOSE>5) {
                        printf("DBCALUnifiedHit_factory  event %5llu (%2i,%i,%i,%i) TDC %i %6.1f  ADC %6.1f diff %6.1f    best %2i %6.1f\n",
                                (unsigned long long int)eventnumber, module, layer, sector, chan.end, i, tdc_hit->t, hit->t, tdc_adc_diff, goodTDCindex, t_diff);
                    }
                }
                t_TDC = tdc_hits[goodTDCindex]->t;
                if (USE_TDC){
                    // Apply the timewalk correction
                    timewalk_coefficients tdc_coeff = tdc_timewalk_map[chan];
                    t_TDC -= tdc_coeff.c0 + tdc_coeff.c1/pow(pulse_peak/tdc_coeff.a_thresh, tdc_coeff.c2);
                }
            }

            // Decide which time to use for further analysis
            t = t_ADC;
            if ( USE_TDC && haveTDChit){
                t = t_TDC;
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
            if (haveTDChit) uhit->AddAssociatedObject(tdc_hits[goodTDCindex]);

            _data.push_back(uhit);
        } // end loop over ADC hist
    } // end loop over cells

    return NOERROR;
}
