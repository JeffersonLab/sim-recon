#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

#include <TRIGGER/DL1Trigger.h>
#include <TPOL/DTPOLSectorDigiHit.h>
#include <TPOL/DTPOLRingDigiHit.h>
#include <DAQ/Df250PulsePedestal.h>
#include <DAQ/Df250Config.h>
#include "DTPOLHit_factory.h"

using namespace std;
using namespace jana;

// static consts need initialization
const double DTPOLHit_factory::SECTOR_DIVISION = 360. / DTPOLHit_factory::NSECTORS;
const double DTPOLHit_factory::INNER_RADIUS = 22. / 2;       // From "ACTIVE INNER DIAMETER" in catalog
const double DTPOLHit_factory::OUTER_RADIUS = 70. / 2;       // From "ACTIVE OUTER DIAMETER" in catalog
const double DTPOLHit_factory::RING_DIVISION   = (OUTER_RADIUS - INNER_RADIUS ) / DTPOLHit_factory::NRINGS;

// Order hits by sector/ring. If same sector/ring, use ADC time to order hits.
bool DTPOLSectorHit_fadc_cmp(const DTPOLSectorDigiHit *a,const DTPOLSectorDigiHit *b){
    if (a->sector==b->sector) return (a->pulse_time<b->pulse_time);
    return (a->sector<b->sector);
}
bool DTPOLRingHit_fadc_cmp(const DTPOLRingDigiHit *a,const DTPOLRingDigiHit *b){
    if (a->ring==b->ring) return (a->pulse_time<b->pulse_time);
    return (a->ring<b->ring);
}

//------------------
// init
//------------------
jerror_t DTPOLHit_factory::init(void)
{
    ADC_THRESHOLD = 50.0;
    gPARMS->SetDefaultParameter("TPOL:ADC_THRESHOLD", ADC_THRESHOLD,
    "ADC pulse-height threshold");

    /// set the base conversion scales
    a_scale    = 0.0001;
    t_scale    = 0.0625; // 62.5 ps/count

    return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTPOLHit_factory::brun(jana::JEventLoop *eventLoop, int32_t runnumber)
{
    /// Read in calibration constants
    //jout << "In DTPOLHit_factory, loading constants..." << endl;
    // load scale factors
    /*map<string,double> scale_factors;
    // a_scale (SC_ADC_SCALE)
    if (eventLoop->GetCalib("/TPOL/digi_scales", scale_factors))
    jout << "Error loading /TPOL/digi_scales !" << endl;
    if (scale_factors.find("TPOL_ADC_ASCALE") != scale_factors.end())
    a_scale = scale_factors["TPOL_ADC_ASCALE"];
    else
    jerr << "Unable to get TPOL_ADC_ASCALE from /TPOL/digi_scales !"
    << endl;
    // t_scale (SC_ADC_SCALE)
    if (scale_factors.find("TPOL_ADC_TSCALE") != scale_factors.end())
    t_scale = scale_factors["TPOL_ADC_TSCALE"];
    else
    jerr << "Unable to get TPOL_ADC_TSCALE from /TPOL/digi_scales !"
    << endl;

    // load base time offset
    map<string,double> base_time_offset;
    if (eventLoop->GetCalib("/TPOL/base_time_offset",base_time_offset))
    jout << "Error loading /TPOL/base_time_offset !" << endl;
    if (base_time_offset.find("TPOL_BASE_TIME_OFFSET") != base_time_offset.end())
    t_base = base_time_offset["TPOL_BASE_TIME_OFFSET"];
    else
    jerr << "Unable to get TPOL_BASE_TIME_OFFSET from /TPOL/base_time_offset !" << endl;

    // load constant tables
    // a_gains (gains)
    if (eventLoop->GetCalib("/TPOL/gains", a_gains))
    jout << "Error loading /TPOL/gains !" << endl;
    // a_pedestals (pedestals)
    if (eventLoop->GetCalib("/TPOL/pedestals", a_pedestals))
    jout << "Error loading /TPOL/pedestals !" << endl;
    // adc_time_offsets (adc_timing_offsets)
    if (eventLoop->GetCalib("/TPOL/adc_timing_offsets", adc_time_offsets))
    jout << "Error loading /TPOL/adc_timing_offsets !" << endl;*/
    return NOERROR;
}



//------------------
// evnt
//------------------
jerror_t DTPOLHit_factory::evnt(JEventLoop *loop, uint64_t eventnumber)
{
    /// Generate DTPOLHit object for each DTPOLSectorDigiHit
    /// and DTPOLRingDigiHit object.
    /// This is where the first set of calibration constants
    /// is applied to convert from digitzed units into natural
    /// units.
    ///
    /// Note that this code does NOT get called for simulated
    /// data in HDDM format. The HDDM event source will copy
    /// the precalibrated values directly into the _data vector.
    //
    // Get fADC250 hits
    /*vector<const DTPOLSectorDigiHit*> sectordigihits;
    loop->Get(sectordigihits);
    sort(sectordigihits.begin(),sectordigihits.end(),DTPOLSectorHit_fadc_cmp);
    vector<const DTPOLRingDigiHit*> ringdigihits;
    loop->Get(ringdigihits);
    sort(ringdigihits.begin(),ringdigihits.end(),DTPOLRingHit_fadc_cmp);
    char str[256];
    // Loop over SECTOR hits
    for (unsigned int i = 0; i < sectordigihits.size(); i++){
        //
        const DTPOLSectorDigiHit *sectordigihit = sectordigihits[i];
        const Df250PulsePedestal* PPobj = NULL;
        sectordigihit->GetSingle(PPobj);
        double pulse_peak = 0.0; double pedestal = 0.0;
        if (PPobj != NULL){
            if (PPobj->pedestal == 0 || PPobj->pulse_peak == 0) continue;
            pedestal = PPobj->pedestal;
            pulse_peak = PPobj->pulse_peak - PPobj->pedestal;
        }
        if (pulse_peak < ADC_THRESHOLD) continue;
        // Make sure sector is in valid range
        if( sectordigihit->sector <= 0 || sectordigihit->sector > DTPOLHit_factory::NSECTORS){
            sprintf(str, "DTPOLSectorDigiHit sector out of range! sector=%d (should be 1-%d)",
            sectordigihit->sector, DTPOLHit_factory::NSECTORS);
            throw JException(str);
        }
        // Initialize pedestal to one found in CCDB, but override it
        // with one found in event if is available
        //double pedestal = a_pedestals[sectordigihit->sector-1];
        const Df250PulseIntegral *pulse_integral = NULL;
        sectordigihit->GetSingle(pulse_integral);
        double pedestal_sum = 0.0;
        if (pulse_integral != NULL) {
            double single_sample_ped = (double)pulse_integral->pedestal;
            double nsamples_integral = (double)pulse_integral->nsamples_integral;
            double nsamples_pedestal = (double)pulse_integral->nsamples_pedestal;
            pedestal_sum = single_sample_ped * nsamples_integral/nsamples_pedestal;
        }
        // Apply calibration constants here
        double A = (double)sectordigihit->pulse_integral;
        double T = (double)sectordigihit->pulse_time;
        if (T == 0.0 || pedestal == 0.0) continue;
        double dA = A - pedestal_sum;
        DTPOLHit *hit = new DTPOLHit;
        hit->sector = sectordigihit->sector;
        hit->ring = 0;
        hit->pulse_peak = pulse_peak;
        hit->integral = dA;
        hit->dE = dA; // This will be scaled to energy units later
        //hit->t = t_scale * T - adc_time_offsets[hit->sector-1] + t_base;
        hit->t = t_scale*T;
        hit->AddAssociatedObject(sectordigihit);
        _data.push_back(hit);
    }*/
    /*
    // Apply calibration constants to convert pulse integrals to energy
    // units
    for (unsigned int i=0;i<_data.size();i++)
    {
    _data[i]->dE*=a_scale * a_gains[_data[i]->sector-1];
}
*/
    // get trigger type
    const DL1Trigger *trig_words = NULL;
    uint32_t trig_mask, fp_trig_mask;
    try {
        loop->GetSingle(trig_words);
    } catch(...) {};
    if (trig_words) {
        trig_mask = trig_words->trig_mask;
        fp_trig_mask = trig_words->fp_trig_mask;
    }
    else {
        trig_mask = 0;
        fp_trig_mask = 0;
    }
    int trig_bits = fp_trig_mask > 0 ? 10 + fp_trig_mask:trig_mask;
    // skim PS triggers
    if (trig_bits!=8) {
        return NOERROR;
    }
    // get raw samples and make TPOL hits
    vector<const Df250WindowRawData*> windowraws;
    loop->Get(windowraws);
    for(unsigned int i=0; i< windowraws.size(); i++){
        const Df250WindowRawData *windowraw = windowraws[i];
        if (windowraw->rocid!=84) continue; // choose rocPS2
        if (!(windowraw->slot==13||windowraw->slot==14)) continue; // azimuthal sectors: 13,14; rings: 15,16
        int slot = windowraw->slot;
        int channel = windowraw->channel;
        int sector = GetSector(slot,channel);
        // Get a vector of the samples for this channel
        const vector<uint16_t> &samplesvector = windowraw->samples;
        unsigned int nsamples=samplesvector.size();
        // loop over the samples to calculate integral, min, max
        double w_samp1 = 0.0; double w_min = 0.0; double w_max = 0.0; double w_integral = 0.0;
        for (uint16_t c_samp=0; c_samp<nsamples; c_samp++) {
            if (c_samp==0) {  // use first sample for initialization
                w_integral = samplesvector[0];
                w_min = samplesvector[0];
                w_max = samplesvector[0];
                w_samp1 = samplesvector[0];
            } else {
                w_integral += samplesvector[c_samp];
                if (w_min > samplesvector[c_samp]) w_min = samplesvector[c_samp];
                if (w_max < samplesvector[c_samp]) w_max = samplesvector[c_samp];
            }
        }
        double pulse_height = w_max - w_min;
        if (w_samp1 > 133.0) continue; // require first sample above readout threshold
        if (pulse_height < ADC_THRESHOLD) continue;
        DTPOLHit *hit = new DTPOLHit;
        hit->sector = sector;
        hit->phi = GetPhi(sector);
        hit->ring = 0;
        hit->theta = 0;
        hit->pulse_peak = pulse_height;
        hit->integral = w_integral - nsamples*w_samp1;
        hit->dE = pulse_height;
        hit->t = t_scale*GetPulseTime(samplesvector,w_min,w_max,ADC_THRESHOLD);
        _data.push_back(hit);
    }
    return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DTPOLHit_factory::erun(void)
{
    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTPOLHit_factory::fini(void)
{
    return NOERROR;
}

int DTPOLHit_factory::GetSector(int slot,int channel)
{
    int sector = 0;
    if (slot == 13) sector = 25 - channel;
    if (slot == 14) {
        if (channel <= 8) sector = 9 - channel;
        else sector = NSECTORS + 9 - channel;
    }
    // fix cable swap
    if (sector == 9) sector = 6;
    else if (sector == 6) sector = 9;
    if (sector == 0) cerr << "sector did not change from initial value (0)." << endl;
    return sector;
}
double DTPOLHit_factory::GetPhi(int sector)
{
    double phi = -10.0;
    if(sector <= 8) phi = (sector + 23)*SECTOR_DIVISION + 0.5*SECTOR_DIVISION;
    if(sector >= 9) phi = (sector - 9)*SECTOR_DIVISION + 0.5*SECTOR_DIVISION;
    return phi;
}
double DTPOLHit_factory::GetPulseTime(const vector<uint16_t> waveform,double w_min,double w_max,double minpeakheight)
{
    // find the time to cross half peak height
    int lastbelowsamp=0; double peakheight = w_max-w_min;
    double threshold = w_min + peakheight/2.0;
    double  firstaboveheight=0, lastbelowheight=0;
    double w_time=0;
    if (peakheight > minpeakheight) {
        for (uint16_t c_samp=0; c_samp<waveform.size(); c_samp++) {
            if (waveform[c_samp]>threshold) {
                firstaboveheight = waveform[c_samp];
                lastbelowsamp = c_samp-1;
                lastbelowheight = waveform[c_samp-1];
                break;
            }
        }
        w_time =  lastbelowsamp + (threshold-lastbelowheight)/(firstaboveheight-lastbelowheight);
    }
    return 64.0*w_time;
}
//------------------------------------
// GetConstant
//   Allow a few different interfaces
//------------------------------------
const double DTPOLHit_factory::GetConstant(const vector<double> &the_table,const int in_sector) const{
    char str[256];

    if ( (in_sector < 0) || (in_sector >= DTPOLHit_factory::NSECTORS)){
        sprintf(str, "Bad sector # requested in DTPOLHit_factory::GetConstant()!"
        " requested=%d , should be %ud", in_sector, DTPOLHit_factory::NSECTORS);
        cerr << str << endl;
        throw JException(str);
    }

    return the_table[in_sector];
}

const double DTPOLHit_factory::GetConstant(const vector<double> &the_table,const DTPOLHit *in_hit) const {
    char str[256];

    if ( (in_hit->sector < 0) || (in_hit->sector >= DTPOLHit_factory::NSECTORS)){
        sprintf(str, "Bad sector # requested in DTPOLHit_factory::GetConstant()!"
        " requested=%d , should be %ud",
        in_hit->sector, DTPOLHit_factory::NSECTORS);
        cerr << str << endl;
        throw JException(str);
    }

    return the_table[in_hit->sector];
}

