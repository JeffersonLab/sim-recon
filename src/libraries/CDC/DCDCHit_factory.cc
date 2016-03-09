// $Id$
//
//    File: DCDCHit_factory.cc
// Created: Tue Aug  6 11:29:56 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include <CDC/DCDCDigiHit.h>
#include "DCDCHit_factory.h"
#include "DCDCWire.h"
#include <DAQ/Df125PulseIntegral.h>
#include <DAQ/Df125Config.h>
#include <DAQ/Df125CDCPulse.h>

using namespace jana;

static double DIGI_THRESHOLD = -1000000.0;
//#define ENABLE_UPSAMPLING

//------------------
// init
//------------------
jerror_t DCDCHit_factory::init(void)
{
    gPARMS->SetDefaultParameter("CDC:DIGI_THRESHOLD",DIGI_THRESHOLD,
            "Do not convert CDC digitized hits into DCDCHit objects"
            " that would have q less than this");

    // default values
    Nrings = 0;
    a_scale = 0.;
    t_scale = 0.;
    t_base = 0.;

    // Set default number of number of detector channels
    maxChannels = 3522;

    /// set the base conversion scales
    a_scale = 4.0E3/1.0E2; 
    t_scale = 8.0/10.0;    // 8 ns/count and integer time is in 1/10th of sample
    t_base  = 0.;       // ns

    return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DCDCHit_factory::brun(jana::JEventLoop *eventLoop, int32_t runnumber)
{
    // Only print messages for one thread whenever run number change
    static pthread_mutex_t print_mutex = PTHREAD_MUTEX_INITIALIZER;
    static set<int> runs_announced;
    pthread_mutex_lock(&print_mutex);
    bool print_messages = false;
    if(runs_announced.find(runnumber) == runs_announced.end()){
        print_messages = true;
        runs_announced.insert(runnumber);
    }
    pthread_mutex_unlock(&print_mutex);

    // calculate the number of straws in each ring
    CalcNstraws(eventLoop, runnumber, Nstraws);
    Nrings = Nstraws.size();

    /// Read in calibration constants
    vector<double> raw_gains;
    vector<double> raw_pedestals;
    vector<double> raw_time_offsets;

    if(print_messages) jout << "In DCDCHit_factory, loading constants..." << std::endl;

    // load scale factors
    map<string,double> scale_factors;
    if (eventLoop->GetCalib("/CDC/digi_scales", scale_factors))
        jout << "Error loading /CDC/digi_scales !" << endl;
    if (scale_factors.find("CDC_ADC_ASCALE") != scale_factors.end())
        a_scale = scale_factors["CDC_ADC_ASCALE"];
    else
        jerr << "Unable to get CDC_ADC_ASCALE from /CDC/digi_scales !" << endl;

#ifdef ENABLE_UPSAMPLING
    //t_scale=1.;
#else
    if (scale_factors.find("CDC_ADC_TSCALE") != scale_factors.end())
        t_scale = scale_factors["CDC_ADC_TSCALE"];
    else
        jerr << "Unable to get CDC_ADC_TSCALE from /CDC/digi_scales !" << endl;
#endif

    // load base time offset
    map<string,double> base_time_offset;
    if (eventLoop->GetCalib("/CDC/base_time_offset",base_time_offset))
        jout << "Error loading /CDC/base_time_offset !" << endl;
    if (base_time_offset.find("CDC_BASE_TIME_OFFSET") != base_time_offset.end())
        t_base = base_time_offset["CDC_BASE_TIME_OFFSET"];
    else
        jerr << "Unable to get CDC_BASE_TIME_OFFSET from /CDC/base_time_offset !" << endl;

    // load constant tables
    if (eventLoop->GetCalib("/CDC/wire_gains", raw_gains))
        jout << "Error loading /CDC/wire_gains !" << endl;
    if (eventLoop->GetCalib("/CDC/pedestals", raw_pedestals))
        jout << "Error loading /CDC/pedestals !" << endl;
    if (eventLoop->GetCalib("/CDC/timing_offsets", raw_time_offsets))
        jout << "Error loading /CDC/timing_offsets !" << endl;

    // fill the tables
    FillCalibTable(gains, raw_gains, Nstraws);
    FillCalibTable(pedestals, raw_pedestals, Nstraws);
    FillCalibTable(time_offsets, raw_time_offsets, Nstraws);

    // Verify that the right number of rings was read for each set of constants
    char str[256];
    if (gains.size() != Nrings) {
        sprintf(str, "Bad # of rings for CDC gain from CCDB! CCDB=%zu , should be %d", gains.size(), Nrings);
        std::cerr << str << std::endl;
        throw JException(str);
    }
    if (pedestals.size() != Nrings) {
        sprintf(str, "Bad # of rings for CDC pedestal from CCDB! CCDB=%zu , should be %d", pedestals.size(), Nrings);
        std::cerr << str << std::endl;
        throw JException(str);
    }
    if (time_offsets.size() != Nrings) {
        sprintf(str, "Bad # of rings for CDC time_offset from CCDB!"
                " CCDB=%zu , should be %d", time_offsets.size(), Nrings);
        std::cerr << str << std::endl;
        throw JException(str);
    }

    // Verify the right number of straws was read for each ring for each set of constants
    for (unsigned int i=0; i < Nrings; i++) {
        if (gains[i].size() != Nstraws[i]) {
            sprintf(str, "Bad # of straws for CDC gain from CCDB!"
                    " CCDB=%zu , should be %d for ring %d",
                    gains[i].size(), Nstraws[i], i+1);
            std::cerr << str << std::endl;
            throw JException(str);
        }
        if (pedestals[i].size() != Nstraws[i]) {
            sprintf(str, "Bad # of straws for CDC pedestal from CCDB!"
                    " CCDB=%zu , should be %d for ring %d",
                    pedestals[i].size(), Nstraws[i], i+1);
            std::cerr << str << std::endl;
            throw JException(str);
        }
        if (time_offsets[i].size() != Nstraws[i]) {
            sprintf(str, "Bad # of straws for CDC time_offset from CCDB!"
                    " CCDB=%zu , should be %d for ring %d",
                    time_offsets[i].size(), Nstraws[i], i+1);
            std::cerr << str << std::endl;
            throw JException(str);
        }
    }

    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DCDCHit_factory::evnt(JEventLoop *loop, uint64_t eventnumber)
{
    /// Generate DCDCHit object for each DCDCDigiHit object.
    /// This is where the first set of calibration constants
    /// is applied to convert from digitzed units into natural
    /// units.
    ///
    /// Note that this code does NOT get called for simulated
    /// data in HDDM format. The HDDM event source will copy
    /// the precalibrated values directly into the _data vector.

    /// In order to use the new Flash125 data types and maintain compatibility with the old code, what is below is a bit of a mess

    vector<const DCDCDigiHit*> digihits;
    loop->Get(digihits);
    char str[256];
    for (unsigned int i=0; i < digihits.size(); i++) {
        const DCDCDigiHit *digihit = digihits[i];

        if ( (digihit->QF & 0x1) != 0 ) continue; // Cut bad timing quality factor hits... (should check effect on efficiency)

        const int &ring  = digihit->ring;
        const int &straw = digihit->straw;

        // Make sure ring and straw are in valid range
        if ( (ring < 1) || (ring > (int)Nrings)) {
            sprintf(str, "DCDCDigiHit ring out of range!"
                    " ring=%d (should be 1-%d)", ring, Nrings);
            throw JException(str);
        }
        if ( (straw < 1) || (straw > (int)Nstraws[ring-1])) {
            sprintf(str, "DCDCDigiHit straw out of range!"
                    " straw=%d for ring=%d (should be 1-%d)",
                    straw, ring, Nstraws[ring-1]);
            throw JException(str);
        }

        // Get pedestal. Preference is given to pedestal measured
        // for event. Otherwise, use statistical one from CCDB
        double pedestal = pedestals[ring-1][straw-1];

        // Grab the pedestal from the digihit since this should be consistent between the old and new formats
        uint32_t raw_ped           = digihit->pedestal;
        uint32_t nsamples_integral = digihit->nsamples_integral;

        // There are a few values from the new data type that are critical for the interpretation of the data
        uint16_t IBIT = 0; // 2^{IBIT} Scale factor for integral
        //        uint16_t ABIT = 0; // 2^{ABIT} Scale factor for amplitude
        uint16_t PBIT = 0; // 2^{PBIT} Scale factor for pedestal
        uint16_t NW   = 0;

        // This is the place to make quality cuts on the data. 
        // Try to get the new data type, if that fails, try to get the old...
        const Df125CDCPulse *CDCPulseObj = NULL;
        digihit->GetSingle(CDCPulseObj);
        if (CDCPulseObj != NULL){
            const Df125Config *config = NULL;
            CDCPulseObj->GetSingle(config);

            // Set some constants to defaults until they appear correctly in the config words in the future
            if(config){
                IBIT = config->IBIT == 0xffff ? 4 : config->IBIT;
                //            ABIT = config->ABIT == 0xffff ? 3 : config->ABIT;
                PBIT = config->PBIT == 0xffff ? 0 : config->PBIT;
                NW   = config->NW   == 0xffff ? 180 : config->NW;
            }else{
                static int Nwarnings = 0;
                if(Nwarnings<10){
                    _DBG_ << "NO Df125Config object associated with Df125FDCPulse object!" << endl;
                    Nwarnings++;
                    if(Nwarnings==10) _DBG_ << " --- LAST WARNING!! ---" << endl;
                }
            }
            if(NW==0) NW=180; // some data was taken (<=run 4700) where NW was written as 0 to file

            // The integration window in the CDC should always extend past the end of the window
            // Only true after about run 4100
            nsamples_integral = (NW - (digihit->pulse_time / 10));  

        }
        else{ // Use the old format
            // This code will at some point become deprecated in the future...
            // This applies to the firmware for data taken until the fall of 2015.
            // Mode 8: Raw data and processed data (except pulse integral).
            // Mode 7: Processed data only.

            // This error state is only present in mode 8
            if (digihit->pulse_time==0.) continue;

            // There is a slight difference between Mode 7 and 8 data
            // The following condition signals an error state in the flash algorithm
            // Do not make hits out of these
            const Df125PulsePedestal* PPobj = NULL;
            digihit->GetSingle(PPobj);
            if (PPobj != NULL){
                if (PPobj->pedestal == 0 || PPobj->pulse_peak == 0) continue;
                if (PPobj->pulse_number == 1) continue; // Unintentionally had 2 pulses found in fall 2014 data (0-1 counting issue)
            }

            const Df125PulseIntegral* PIobj = NULL;
            digihit->GetSingle(PIobj);
            if (PPobj == NULL || PIobj == NULL) continue; // We don't want hits where ANY of the associated information is missing
        }

        // Complete the pedestal subtracion here since we should know the correct number of samples.
        uint32_t scaled_ped = raw_ped << PBIT;
        pedestal = double(scaled_ped * nsamples_integral);

        // Apply calibration constants here
        double integral = double(digihit->pulse_integral << IBIT);
        double t_raw = double(digihit->pulse_time);

        double q = a_scale * gains[ring-1][straw-1] * (integral - pedestal);
        double t = t_scale * t_raw - time_offsets[ring-1][straw-1] + t_base;

        if (q < DIGI_THRESHOLD) 
            continue;

        DCDCHit *hit = new DCDCHit;
        hit->ring  = ring;
        hit->straw = straw;

        // Values for d, itrack, ptype only apply to MC data
        // note that ring/straw counting starts at 1
        hit->q = q;
        hit->t = t;
        hit->d = 0.0;
        hit->itrack = -1;
        hit->ptype = 0;

        hit->AddAssociatedObject(digihit);

        _data.push_back(hit);
    }

    return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DCDCHit_factory::erun(void)
{
    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DCDCHit_factory::fini(void)
{
    return NOERROR;
}

//------------------
// CalcNstraws
//------------------
void DCDCHit_factory::CalcNstraws(jana::JEventLoop *eventLoop, int32_t runnumber, vector<unsigned int> &Nstraws)
{
    DGeometry *dgeom;
    vector<vector<DCDCWire *> >cdcwires;

    // Get pointer to DGeometry object
    DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
    dgeom  = dapp->GetDGeometry(runnumber);

    // Get the CDC wire table from the XML
    dgeom->GetCDCWires(cdcwires);

    // Fill array with the number of straws for each layer
    // Also keep track of the total number of straws, i.e., the total number of detector channels
    maxChannels = 0;
    Nstraws.clear();
    for (unsigned int i=0; i<cdcwires.size(); i++) {
        Nstraws.push_back( cdcwires[i].size() );
        maxChannels += cdcwires[i].size();
    }

    // clear up all of the wire information
    for (unsigned int i=0; i<cdcwires.size(); i++) {
        for (unsigned int j=0; j<cdcwires[i].size(); j++) {
            delete cdcwires[i][j];
        }
    }    
    cdcwires.clear();
}


//------------------
// FillCalibTable
//------------------
void DCDCHit_factory::FillCalibTable(vector< vector<double> > &table, vector<double> &raw_table, 
        vector<unsigned int> &Nstraws)
{
    int ring = 0;
    int straw = 0;

    // reset table before filling it
    table.clear();
    table.resize( Nstraws.size() );

    for (unsigned int channel=0; channel<raw_table.size(); channel++,straw++) {
        // make sure that we don't try to load info for channels that don't exist
        if (channel == maxChannels) break;

        // if we've hit the end of the ring, move on to the next
        if (straw == (int)Nstraws[ring]) {
            ring++;
            straw = 0;
        }

        table[ring].push_back( raw_table[channel] );
    }

}


//------------------------------------
// GetConstant
//   Allow a few different interfaces
//------------------------------------
const double DCDCHit_factory::GetConstant(const cdc_digi_constants_t &the_table,
        const int in_ring, const int in_straw) const {

    char str[256];

    if ( (in_ring <= 0) || (static_cast<unsigned int>(in_ring) > Nrings)) {
        sprintf(str, "Bad ring # requested in DCDCHit_factory::GetConstant()!"
                " requested=%d , should be %ud", in_ring, Nrings);
        std::cerr << str << std::endl;
        throw JException(str);
    }
    if ( (in_straw <= 0) || 
            (static_cast<unsigned int>(in_straw) > Nstraws[in_ring]))
    {
        sprintf(str, "Bad straw # requested in DCDCHit_factory::GetConstant()!"
                " requested=%d , should be %ud", in_straw, Nstraws[in_ring]);
        std::cerr << str << std::endl;
        throw JException(str);
    }

    return the_table[in_ring-1][in_straw-1];
}

const double DCDCHit_factory::GetConstant(const cdc_digi_constants_t &the_table,
        const DCDCDigiHit *in_digihit) const {

    char str[256];

    if ( (in_digihit->ring <= 0) || 
            (static_cast<unsigned int>(in_digihit->ring) > Nrings))
    {
        sprintf(str, "Bad ring # requested in DCDCHit_factory::GetConstant()!"
                " requested=%d , should be %ud", in_digihit->ring, Nrings);
        std::cerr << str << std::endl;
        throw JException(str);
    }
    if ( (in_digihit->straw <= 0) || 
            (static_cast<unsigned int>(in_digihit->straw) > Nstraws[in_digihit->ring]))
    {
        sprintf(str, "Bad straw # requested in DCDCHit_factory::GetConstant()!"
                " requested=%d , should be %ud",
                in_digihit->straw, Nstraws[in_digihit->ring]);
        std::cerr << str << std::endl;
        throw JException(str);
    }

    return the_table[in_digihit->ring-1][in_digihit->straw-1];
}

const double DCDCHit_factory::GetConstant(const cdc_digi_constants_t &the_table,
        const DCDCHit *in_hit) const {

    char str[256];

    if ( (in_hit->ring <= 0) || (static_cast<unsigned int>(in_hit->ring) > Nrings)) {
        sprintf(str, "Bad ring # requested in DCDCHit_factory::GetConstant()!"
                " requested=%d , should be %ud", in_hit->ring, Nrings);
        std::cerr << str << std::endl;
        throw JException(str);
    }
    if ( (in_hit->straw <= 0) || 
            (static_cast<unsigned int>(in_hit->straw) > Nstraws[in_hit->ring])) {
        sprintf(str, "Bad straw # requested in DCDCHit_factory::GetConstant()!"
                " requested=%d , should be %ud",
                in_hit->straw, Nstraws[in_hit->ring]);
        std::cerr << str << std::endl;
        throw JException(str);
    }

    return the_table[in_hit->ring-1][in_hit->straw-1];
}
/*
   const double DCDCHit_factory::GetConstant(const cdc_digi_constants_t &the_table,
   const DTranslationTable *ttab,
   const int in_rocid, const int in_slot, const int in_channel) const {

   char str[256];

   DTranslationTable::csc_t daq_index = { in_rocid, in_slot, in_channel };
   DTranslationTable::DChannelInfo channel_info = ttab->GetDetectorIndex(daq_index);

   if ( (channel_info.cdc.ring <= 0) || (channel_info.cdc.ring > Nrings)) {
   sprintf(str, "Bad ring # requested in DCDCHit_factory::GetConstant()!"
   " requested=%d , should be %ud", channel_info.cdc.ring, Nrings);
   std::cerr << str << std::endl;
   throw JException(str);
   }
   if ( (channel_info.cdc.straw <= 0) || 
   (channel_info.cdc.straw > Nstraws[channel_info.cdc.ring]))
   {
   sprintf(str, "Bad straw # requested in DCDCHit_factory::GetConstant()!"
   " requested=%d , should be %ud",
   channel_info.cdc.ring, Nstraws[channel_info.cdc.straw]);
   std::cerr << str << std::endl;
   throw JException(str);
   }

   return the_table[channel_info.cdc.ring-1][channel_info.cdc.straw-1];
   }
   */
