// $Id$
//
//    File: JEventProcessor_HLDetectorTiming.cc
// Created: Mon Jan 12 14:37:56 EST 2015
// Creator: mstaib (on Linux egbert 2.6.32-431.20.3.el6.x86_64 x86_64)
//

#include "JEventProcessor_HLDetectorTiming.h"
using namespace jana;

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>

#include "PID/DDetectorMatches.h"
#include "HistogramTools.h"

extern "C"{
void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_HLDetectorTiming());
}
} // "C"


//------------------
// JEventProcessor_HLDetectorTiming (Constructor)
//------------------
JEventProcessor_HLDetectorTiming::JEventProcessor_HLDetectorTiming()
{

}

//------------------
// ~JEventProcessor_HLDetectorTiming (Destructor)
//------------------
JEventProcessor_HLDetectorTiming::~JEventProcessor_HLDetectorTiming()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_HLDetectorTiming::init(void)
{
    BEAM_CURRENT = 50; // Assume that there is beam until first EPICs event. Set from EPICS evio data, can override on command line

    fBeamEventCounter = 0;
    REQUIRE_BEAM = 1;
    BEAM_EVENTS_TO_KEEP = 1000000000; // Set enormously high
    DO_ROUGH_TIMING = 0;
    DO_CDC_TIMING = 0;
    DO_TDC_ADC_ALIGN = 0;
    DO_TRACK_BASED = 0;
    DO_VERIFY = 1;
    DO_OPTIONAL = 0;

    if(gPARMS){
        gPARMS->SetDefaultParameter("HLDETECTORTIMING:DO_ROUGH_TIMING", DO_ROUGH_TIMING, "Set to > 0 to do rough timing of all detectors");
        gPARMS->SetDefaultParameter("HLDETECTORTIMING:DO_CDC_TIMING", DO_CDC_TIMING, "Set to > 0 to do CDC Per channel Alignment");
        gPARMS->SetDefaultParameter("HLDETECTORTIMING:DO_TDC_ADC_ALIGN", DO_TDC_ADC_ALIGN, "Set to > 0 to do TDC/ADC alignment of SC,TOF,TAGM,TAGH");
        gPARMS->SetDefaultParameter("HLDETECTORTIMING:DO_TRACK_BASED", DO_TRACK_BASED, "Set to > 0 to do Track Based timing corrections");
        gPARMS->SetDefaultParameter("HLDETECTORTIMING:DO_VERIFY", DO_VERIFY, "Set to > 0 to verify timing with current constants");
        gPARMS->SetDefaultParameter("HLDETECTORTIMING:REQUIRE_BEAM", REQUIRE_BEAM, "Set to 0 to skip beam current check");
        gPARMS->SetDefaultParameter("HLDETECTORTIMING:BEAM_EVENTS_TO_KEEP", BEAM_EVENTS_TO_KEEP, "Set to the number of beam on events to use");
        gPARMS->SetDefaultParameter("HLDETECTORTIMING:DO_OPTIONAL", BEAM_EVENTS_TO_KEEP, "Set to >0 to enable optional histograms ");
    }
    
    // Would like the code with no arguments to simply verify the current status of the calibration
    if (DO_ROUGH_TIMING > 0 || DO_CDC_TIMING > 0 || DO_TDC_ADC_ALIGN > 0 || DO_TRACK_BASED > 0) DO_VERIFY = 0;

    return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_HLDetectorTiming::brun(JEventLoop *eventLoop, int runnumber)
{
    // This is called whenever the run number changes
    // Get the particleID object for each run
    // Thanks to Justin for the How-To ;)

    vector<const DParticleID *> dParticleID_algos;
    eventLoop->Get(dParticleID_algos);
    if(dParticleID_algos.size()<1){
        _DBG_<<"Unable to get a DParticleID object! NO PID will be done!"<<endl;
        return RESOURCE_UNAVAILABLE;
    }
    dParticleID = dParticleID_algos[0];

    // load base time offsets
    map<string,double> base_time_offset;
    //CDC
    if (eventLoop->GetCalib("/CDC/base_time_offset",base_time_offset))
        jout << "Error loading /CDC/base_time_offset !" << endl;
    if (base_time_offset.find("CDC_BASE_TIME_OFFSET") != base_time_offset.end())
        cdc_t_base = base_time_offset["CDC_BASE_TIME_OFFSET"];
    else
        jerr << "Unable to get CDC_BASE_TIME_OFFSET from /CDC/base_time_offset !" << endl;

    //FCAL
    if (eventLoop->GetCalib("/FCAL/base_time_offset",base_time_offset))
        jout << "Error loading /FCAL/base_time_offset !" << endl;
    if (base_time_offset.find("FCAL_BASE_TIME_OFFSET") != base_time_offset.end())
        fcal_t_base = base_time_offset["FCAL_BASE_TIME_OFFSET"];
    else
        jerr << "Unable to get FCAL_BASE_TIME_OFFSET from /FCAL/base_time_offset !" << endl;
    // BCAL
    if (eventLoop->GetCalib("/BCAL/base_time_offset",base_time_offset))
        jout << "Error loading /BCAL/base_time_offset !" << endl;
    if (base_time_offset.find("BCAL_BASE_TIME_OFFSET") != base_time_offset.end()) {
        bcal_t_base = base_time_offset["BCAL_BASE_TIME_OFFSET"];
    }
    else
        jerr << "Unable to get BCAL_BASE_TIME_OFFSET from /BCAL/base_time_offset !" << endl;
    //TOF
    if (eventLoop->GetCalib("/TOF/base_time_offset",base_time_offset))
        jout << "Error loading /TOF/base_time_offset !" << endl;
    if (base_time_offset.find("TOF_BASE_TIME_OFFSET") != base_time_offset.end())
        tof_t_base_fadc = base_time_offset["TOF_BASE_TIME_OFFSET"];
    else
        jerr << "Unable to get TOF_BASE_TIME_OFFSET from /TOF/base_time_offset !" << endl;  

    if (base_time_offset.find("TOF_TDC_BASE_TIME_OFFSET") != base_time_offset.end())
        tof_t_base_tdc = base_time_offset["TOF_TDC_BASE_TIME_OFFSET"];
    else
        jerr << "Unable to get TOF_TDC_BASE_TIME_OFFSET from /TOF/base_time_offset !" << endl;

    // The SC, TOF, TAGM and TAGH we will need the per channel offsets to adjust
    // SC
    if (eventLoop->GetCalib("/START_COUNTER/tdc_timing_offsets", sc_tdc_time_offsets))
        jout << "Error loading /START_COUNTER/tdc_timing_offsets !" << endl;  

    // TOF
    if(eventLoop->GetCalib("TOF/timing_offsets", tof_tdc_time_offsets))
        jout << "Error loading /TOF/timing_offsets !" << endl; 

    // TAGM
    std::vector< std::map<std::string, double> > table;
    if (eventLoop->GetCalib("/PHOTON_BEAM/microscope/fadc_time_offsets", table))
        jout << "Error loading /PHOTON_BEAM/microscope/fadc_time_offsets from ccdb!" << std::endl;
    for (unsigned int i=0; i < table.size(); ++i) {
        //int row = (table[i])["row"];
        int col = (table[i])["column"];
        tagm_fadc_time_offsets[col] = (table[i])["offset"]; // Only aligning per column for now...
    }    

    if (eventLoop->GetCalib("/PHOTON_BEAM/microscope/tdc_time_offsets", table))
        jout << "Error loading /PHOTON_BEAM/microscope/tdc_time_offsets from ccdb!" << std::endl;
    for (unsigned int i=0; i < table.size(); ++i) {
        //int row = (table[i])["row"];
        int col = (table[i])["column"];
        tagm_tdc_time_offsets[col] = (table[i])["offset"]; // Only aligning per column for now...
    }
    // TAGH
    if (eventLoop->GetCalib("/PHOTON_BEAM/hodoscope/fadc_time_offsets", table))
        jout << "Error loading /PHOTON_BEAM/hodoscope/fadc_time_offsets from ccdb!" << std::endl;
    for (unsigned int i=0; i < table.size(); ++i) {
        int counter = (table[i])["id"];
        tagh_fadc_time_offsets[counter] = (table[i])["offset"];
    }

    if (eventLoop->GetCalib("/PHOTON_BEAM/hodoscope/tdc_time_offsets", table))
        jout << "Error loading /PHOTON_BEAM/hodoscope/tdc_time_offsets from ccdb!" << std::endl;
    for (unsigned int i=0; i < table.size(); ++i) {
        int counter = (table[i])["id"];
        tagh_tdc_time_offsets[counter] = (table[i])["offset"];
    }

    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_HLDetectorTiming::evnt(JEventLoop *loop, int eventnumber)
{

    // Get the EPICs events and update beam current. Skip event if current too low (<10 nA).
    if (REQUIRE_BEAM){
        vector<const DEPICSvalue *> epicsValues;
        loop->Get(epicsValues);
        for(unsigned int j = 0; j < epicsValues.size(); j++){
            const DEPICSvalue *thisValue = epicsValues[j];
            if (strcmp((thisValue->name).c_str(), "IBCAD00CRCUR6") == 0){
                BEAM_CURRENT = thisValue->fval;
                Fill1DHistogram("HLDetectorTiming", "", "Beam Current",
                        BEAM_CURRENT,
                        "Beam Current; Beam Current [nA]; Entries",
                        100, 0, 200);
            }
            //cout << "EPICS Name " <<  (thisValue->name).c_str() << " Value " << thisValue->fval << endl;
        }
        if (BEAM_CURRENT < 10.0) {
            Fill1DHistogram("HLDetectorTiming", "" , "Beam Events",
                    0, "Beam On Events (0 = no beam, 1 = beam > 10nA)",
                    2, -0.5, 1.5);
            return NOERROR; // Skip events where we can't verify the beam current
        }
        Fill1DHistogram("HLDetectorTiming", "" , "Beam Events",
                1, "Beam On Events (0 = no beam, 1 = beam > 10nA)",
                2, -0.5, 1.5);
        fBeamEventCounter++;
        if (fBeamEventCounter >= BEAM_EVENTS_TO_KEEP) {
            japp->Quit();
            return NOERROR;
        }
    }

    vector<const DCDCHit *> cdcHitVector;
    loop->Get(cdcHitVector);

    vector<const DFDCHit *> fdcHitVector;
    loop->Get(fdcHitVector);

    vector<const DSCHit *> scHitVector;
    loop->Get(scHitVector);

    vector<const DBCALHit *> bcalHitVector;
    loop->Get(bcalHitVector);

    vector<const DTOFHit *> tofHitVector;
    loop->Get(tofHitVector);

    vector<const DFCALHit *> fcalHitVector;
    loop->Get(fcalHitVector);

    vector<const DTAGMHit *> tagmHitVector;
    loop->Get(tagmHitVector);

    vector<const DTAGHHit *> taghHitVector;
    loop->Get(taghHitVector);

    unsigned int i = 0;

    int nBins = 2000;
    float xMin = -500, xMax = 1500;
    for (i = 0; i < cdcHitVector.size(); i++){
        Fill1DHistogram ("HLDetectorTiming", "CDC", "CDCHit time", cdcHitVector[i]->t, 
                "CDCHit time", nBins, xMin, xMax);
        if(DO_VERIFY || DO_CDC_TIMING){
            int nStraws = 3522;
            Fill2DHistogram("HLDetectorTiming", "CDC", "CDCHit time per Straw Raw", 
                    cdcHitVector[i]->t, GetCCDBIndexCDC(cdcHitVector[i]),
                    "Hit time for each CDC wire; t [ns]; CCDB Index",
                    750, -500, 1000, nStraws, 0.5, nStraws + 0.5);
            if (DO_OPTIONAL){
                Fill2DHistogram("HLDetectorTiming", "CDC", "CDC Energy Vs. Drift Time",
                        cdcHitVector[i]->t, cdcHitVector[i]->q , 
                        "CDC Energy Vs. Drift Time; Drift Time [ns]; Integral [arb.]",
                        800, 0, 800, 500, -10000, 7000000);
            }
        }
    }

    for (i = 0; i < fdcHitVector.size(); i++){
        if(fdcHitVector[i]->type != 1 ) continue; // ADC hits only
        Fill1DHistogram ("HLDetectorTiming", "FDC", "FDCHit time", fdcHitVector[i]->t,
                "FDCHit time", nBins, xMin, xMax);
    }

    for (i = 0; i < scHitVector.size(); i++){
        //if(!scHitVector[i]->has_fADC || !scHitVector[i]->has_TDC) continue;
        Fill1DHistogram ("HLDetectorTiming", "SC", "SCHit time", scHitVector[i]->t,
                "SCHit time", nBins, xMin, xMax);
    }
    for (i = 0; i < bcalHitVector.size(); i++){
        Fill1DHistogram ("HLDetectorTiming", "BCAL", "BCALHit time", bcalHitVector[i]->t,
                "BCALHit time", nBins, xMin, xMax);
    }
    for (i = 0; i < tofHitVector.size(); i++){
        Fill1DHistogram ("HLDetectorTiming", "TOF", "TOFHit time", tofHitVector[i]->t,
                "TOFHit time", nBins, xMin, xMax);
    }
    for (i = 0; i < fcalHitVector.size(); i++){
        Fill1DHistogram ("HLDetectorTiming", "FCAL", "FCALHit time", fcalHitVector[i]->t,
                "FCALHit time", nBins, xMin, xMax);
    }
    for (i = 0; i < tagmHitVector.size(); i++){
        Fill1DHistogram ("HLDetectorTiming", "TAGM", "TAGMHit time", tagmHitVector[i]->t,
                "TAGMHit time", nBins, xMin, xMax);
    }
    for (i = 0; i < taghHitVector.size(); i++){
        Fill1DHistogram ("HLDetectorTiming", "TAGH", "TAGHHit time", taghHitVector[i]->t,
                "TAGHHit time", nBins, xMin, xMax);
    }

    // The detectors with both TDCs and ADCs need these two to be aligned
    // These detectors are the SC,TAGM,TAGH,TOF

    // Break these histograms up into hits coming from the TDC and hits coming from the ADC
    for (i = 0; i < scHitVector.size(); i++){
        //if(!scHitVector[i]->has_fADC || !scHitVector[i]->has_TDC) continue;
        const DSCHit *thisSCHit = scHitVector[i];
        if (thisSCHit->has_fADC && !thisSCHit->has_TDC){
            Fill1DHistogram ("HLDetectorTiming", "SC", "SCHit ADC time", scHitVector[i]->t,
                    "SCHit ADC only time", nBins, xMin, xMax);
        }
        else if (!thisSCHit->has_fADC && thisSCHit->has_TDC){
            Fill1DHistogram ("HLDetectorTiming", "SC", "SCHit TDC time", scHitVector[i]->t,
                    "SCHit TDC only time", nBins, xMin, xMax);
        }
        else{
            Fill1DHistogram ("HLDetectorTiming", "SC", "SCHit Matched time", scHitVector[i]->t,
                    "SCHit Matched ADC/TDC time", nBins, xMin, xMax);
            Fill1DHistogram ("HLDetectorTiming", "SC", "SCHit ADC time", scHitVector[i]->t_fADC,
                    "SCHit ADC only time", nBins, xMin, xMax);
            Fill1DHistogram ("HLDetectorTiming", "SC", "SCHit TDC time", scHitVector[i]->t_TDC,
                    "SCHit TDC only time", nBins, xMin, xMax);

            int nSCCounters = 30;
            Fill2DHistogram("HLDetectorTiming", "SC", "SCHit TDC_ADC Difference",
                    scHitVector[i]->sector, scHitVector[i]->t_TDC - scHitVector[i]->t_fADC,
                    "SC #Deltat TDC-ADC; Sector ;t_{TDC} - t_{ADC} [ns]", nSCCounters, 0.5, nSCCounters + 0.5, 100, -10, 10);
            if (DO_TDC_ADC_ALIGN){
                char name [200];
                char title[500];
                sprintf(name, "Sector %.2i", scHitVector[i]->sector);
                sprintf(title, "SC Sector %i #Deltat Vs dE; dE [Arb.]; t_{TDC} - t_{ADC} [ns]", scHitVector[i]->sector);
                Fill2DHistogram("HLDetectorTiming", "SC_TimeWalk", name,
                        scHitVector[i]->dE, scHitVector[i]->t_TDC - scHitVector[i]->t_fADC,
                        title,
                        50, 0, 30000, 50, -10, 10, false);
            }
        }

    }
    for (i = 0; i < tagmHitVector.size(); i++){
        const DTAGMHit *thisTAGMHit = tagmHitVector[i];
        if(thisTAGMHit->has_fADC && !thisTAGMHit->has_TDC){
            Fill1DHistogram ("HLDetectorTiming", "TAGM", "TAGMHit ADC time", tagmHitVector[i]->t,
                    "TAGMHit ADC only time", nBins, xMin, xMax);
        }
        else if (!thisTAGMHit->has_fADC && thisTAGMHit->has_TDC){
            Fill1DHistogram ("HLDetectorTiming", "TAGM", "TAGMHit TDC time", tagmHitVector[i]->t,
                    "TAGMHit TDC only time", nBins, xMin, xMax);
        }
        else{
            Fill1DHistogram ("HLDetectorTiming", "TAGM", "TAGMHit Matched time", tagmHitVector[i]->t,
                    "TAGMHit Matched ADC/TDC time", nBins, xMin, xMax);
            Fill1DHistogram ("HLDetectorTiming", "TAGM", "TAGMHit ADC time", tagmHitVector[i]->time_fadc,
                    "TAGMHit ADC only time", nBins, xMin, xMax);
            Fill1DHistogram ("HLDetectorTiming", "TAGM", "TAGMHit TDC time", tagmHitVector[i]->time_tdc,
                    "TAGMHit TDC only time", nBins, xMin, xMax);

            int nTAGMCounters = 102;
            // We want to look at the timewalk within these ADC/TDC detectors
            Fill2DHistogram("HLDetectorTiming", "TAGM", "TAGMHit TDC_ADC Difference",
                    tagmHitVector[i]->column, tagmHitVector[i]->time_tdc - tagmHitVector[i]->time_fadc,
                    "TAGM #Deltat TDC-ADC; Column ;t_{TDC} - t_{ADC} [ns]", nTAGMCounters, 0.5, nTAGMCounters + 0.5, 100, -10, 10);
            if (DO_TDC_ADC_ALIGN){
                // TString::Form() is not thread safe! (I learned this the hard way)
                // Do it the old fashioned way
                char name [200];
                char title[500];
                sprintf(name, "Column %.3i", tagmHitVector[i]->column);
                sprintf(title, "TAGM Column %i #Deltat Vs Integral; Integral [ADC Units]; t_{TDC} - t_{ADC} [ns]", tagmHitVector[i]->column);
                Fill2DHistogram("HLDetectorTiming", "TAGM_TimeWalk", name,
                        tagmHitVector[i]->integral, tagmHitVector[i]->time_tdc - tagmHitVector[i]->time_fadc,
                        title,
                        50, 0, 5000, 50, -10, 10, false);
            }
        }

    }
    for (i = 0; i < taghHitVector.size(); i++){
        const DTAGHHit *thisTAGHHit = taghHitVector[i];
        if(thisTAGHHit->has_fADC && !thisTAGHHit->has_TDC){
            Fill1DHistogram ("HLDetectorTiming", "TAGH", "TAGHHit ADC time", taghHitVector[i]->t,
                    "TAGHHit ADC only time", nBins, xMin, xMax);
        }
        else if (!thisTAGHHit->has_fADC && thisTAGHHit->has_TDC){
            Fill1DHistogram ("HLDetectorTiming", "TAGH", "TAGHHit TDC time", taghHitVector[i]->t,
                    "TAGHHit TDC only time", nBins, xMin, xMax);
        }
        else{
            Fill1DHistogram ("HLDetectorTiming", "TAGH", "TAGHHit Matched time", taghHitVector[i]->t,
                    "TAGHHit Matched ADC/TDC time", nBins, xMin, xMax);
            Fill1DHistogram ("HLDetectorTiming", "TAGH", "TAGHHit ADC time", taghHitVector[i]->time_fadc,
                    "TAGHHit ADC only time", nBins, xMin, xMax);
            Fill1DHistogram ("HLDetectorTiming", "TAGH", "TAGHHit TDC time", taghHitVector[i]->time_tdc,
                    "TAGHHit TDC only time", nBins, xMin, xMax);

            int nTAGHCounters = 274;
            // We want to look at the timewalk within these ADC/TDC detectors
            Fill2DHistogram("HLDetectorTiming", "TAGH", "TAGHHit TDC_ADC Difference",
                    taghHitVector[i]->counter_id, taghHitVector[i]->time_tdc - taghHitVector[i]->time_fadc,
                    "TAGH #Deltat TDC-ADC; Counter ID ;t_{TDC} - t_{ADC} [ns]", nTAGHCounters, 0.5, nTAGHCounters + 0.5, 100, -10, 10);

            if (DO_TDC_ADC_ALIGN){
                char name [200];
                char title[500];
                sprintf(name, "Counter ID %.3i", taghHitVector[i]->counter_id);
                sprintf(title, "TAGH Counter ID %i #Deltat Vs Integral; Integral [ADC Units]; t_{TDC} - t_{ADC} [ns]", taghHitVector[i]->counter_id);
                Fill2DHistogram("HLDetectorTiming", "TAGH_TimeWalk", name,
                        taghHitVector[i]->integral, taghHitVector[i]->time_tdc - taghHitVector[i]->time_fadc,
                        title,
                        50, 0, 5000, 50, -10, 10, false);
            }
        }
    }
    for (i = 0; i < tofHitVector.size(); i++){
        const DTOFHit *thisTOFHit = tofHitVector[i];
        if(thisTOFHit->has_fADC && !thisTOFHit->has_TDC){
            Fill1DHistogram ("HLDetectorTiming", "TOF", "TOFHit ADC time", tofHitVector[i]->t,
                    "TOFHit ADC only time", nBins, xMin, xMax);
        }
        else if (!thisTOFHit->has_fADC && thisTOFHit->has_TDC){
            Fill1DHistogram ("HLDetectorTiming", "TOF", "TOFHit TDC time", tofHitVector[i]->t,
                    "TOFHit TDC only time", nBins, xMin, xMax);
        }
        else{
            Fill1DHistogram ("HLDetectorTiming", "TOF", "TOFHit Matched time", tofHitVector[i]->t,
                    "TOFHit Matched ADC/TDC time", nBins, xMin, xMax);
            Fill1DHistogram ("HLDetectorTiming", "TOF", "TOFHit ADC time", tofHitVector[i]->t_fADC,
                    "TOFHit ADC only time", nBins, xMin, xMax);
            Fill1DHistogram ("HLDetectorTiming", "TOF", "TOFHit TDC time", tofHitVector[i]->t_TDC,
                    "TOFHit TDC only time", nBins, xMin, xMax);

            int nTOFCounters = 176;
            Fill2DHistogram("HLDetectorTiming", "TOF", "TOFHit TDC_ADC Difference",
                    GetCCDBIndexTOF(tofHitVector[i]), tofHitVector[i]->t_TDC - tofHitVector[i]->t_fADC,
                    "TOF #Deltat TDC-ADC; CDCB Index ;t_{TDC} - t_{ADC} [ns]", nTOFCounters, 0.5, nTOFCounters + 0.5, 100, -10, 10, true);

            if (DO_TDC_ADC_ALIGN){
                // TString->Form(...) is not thread safe, trying sprintf
                char name [200];
                char title[500];
                sprintf(name, "Counter ID %.3i", GetCCDBIndexTOF(tofHitVector[i]));
                sprintf(title, "TOF CCDB Counter ID %i #Deltat Vs dE; dE [Arb.]; t_{TDC} - t_{ADC} [ns]", GetCCDBIndexTOF(tofHitVector[i]));
                Fill2DHistogram("HLDetectorTiming", "TOF_TimeWalk", name,
                        tofHitVector[i]->dE, tofHitVector[i]->t_TDC - tofHitVector[i]->t_fADC,
                        title,
                        50, 0, 30000, 50, -10, 10, false);
            }
        }
    }

    // Next the relative times between detectors using tracking
    // By the time we get to this point, our first guess at the timing should be fairly good. 
    // Certainly good enough to take a pass at the time based tracking
    // This will be the final alignment step for now

    if (!DO_TRACK_BASED && !DO_VERIFY && !DO_CDC_TIMING) return NOERROR; // Before this stage we aren't really ready yet, so just return

    vector<const DTrackTimeBased *> timeBasedTrackVector;
    loop->Get(timeBasedTrackVector);

    const DDetectorMatches* locDetectorMatches = NULL;
    loop->GetSingle(locDetectorMatches);

    for (i = 0; i < timeBasedTrackVector.size(); i++){

        const DTrackTimeBased * thisTrack = timeBasedTrackVector[i];

        // Some quality cuts for the tracks we will use
        // Keep this minimal for now and investigate later
        //float trackingFOMCut = 0.01;
        float trackingFOMCut =0.0027;
        //int trackingNDFCut = 5;
        int trackingNDFCut = 7;

        if(thisTrack->FOM < trackingFOMCut) continue;
        if(thisTrack->Ndof < trackingNDFCut) continue;


        // Do CDC Fine timing with these tracks (Requires a lot of statistics)
        if (DO_VERIFY || DO_CDC_TIMING){
            const vector<DTrackFitter::pull_t> *pulls = &thisTrack->pulls;

            // Since there is no smoothing of the track, the first value in the pulls is tied to the 
            // tracking result. Dont use this one since the residual is always zero.

            for (unsigned int iPull = 1; iPull < pulls->size(); iPull++){
                if ((*pulls)[iPull].cdc_hit == NULL) continue;
                int ring = (*pulls)[iPull].cdc_hit->wire->ring;
                int straw = (*pulls)[iPull].cdc_hit->wire->straw;

                double tdrift = (*pulls)[iPull].tdrift;
                double residual = (*pulls)[iPull].resi;
                int nStraws = 3522;
                int CCDBIndex = GetCCDBIndexCDC(ring, straw);
                Fill2DHistogram("HLDetectorTiming", "CDC", "CDC Drift time per Straw",
                        tdrift, CCDBIndex,
                        "Drift time for each CDC wire (flight time/B-field corr.); t [ns]; CCDB Index",
                        750, -500, 1000, nStraws, 0.5, nStraws + 0.5);
                Fill2DHistogram("HLDetectorTiming", "CDC", "CDC Residual Vs Drift Time",
                        tdrift, residual,
                        "CDC Residual Vs Drift time; Drift time [ns]; Residual [cm]",
                        1020,-20, 1000, 200, -0.1, 0.1);
            }
        }


        // For calibration of the timing, we might as well just look at the negatively charged tracks 
        // Should reject proton candidates
        // Will move this cut to later, statistics are lacking for Tagger SC alignment.
        // if(thisTrack->charge() > 0) continue; 

        // At this point the timing alignment is good enough to make some use of the PID libraries
        // Copied from Justin's physics_detcom plugin

        //////////////////////////////////////////
        // get best matches to SC/TOF/FCAL/BCAL //
        //////////////////////////////////////////
        DSCHitMatchParams locSCHitMatchParams;
        DTOFHitMatchParams locTOFHitMatchParams;
        DFCALShowerMatchParams locFCALShowerMatchParams;
        DBCALShowerMatchParams locBCALShowerMatchParams;

        bool foundSC = dParticleID->Get_BestSCMatchParams(thisTrack, locDetectorMatches, locSCHitMatchParams);
        bool foundTOF = dParticleID->Get_BestTOFMatchParams(thisTrack, locDetectorMatches, locTOFHitMatchParams);
        bool foundFCAL = dParticleID->Get_BestFCALMatchParams(thisTrack, locDetectorMatches, locFCALShowerMatchParams);
        bool foundBCAL = dParticleID->Get_BestBCALMatchParams(thisTrack, locDetectorMatches, locBCALShowerMatchParams);

        // We will only use tracks matched to the start counter for our calibration since this will be our reference for t0
        if (!foundSC) continue;

        // the idea will be to fix the SC time and reference the other PID detectors off of this

        // We want to plot the delta t at the target between the SC hit and the tagger hits
        float nBinsE = 160, EMin = 2.0, EMax = 10;
        float nBinsdT = 100, dTMin = -50.0, dTMax = 50.0;

        // These "flightTime" corrected time are essentially that detector's estimate of the target time
        float flightTimeCorrectedSCTime = locSCHitMatchParams.dHitTime - locSCHitMatchParams.dFlightTime; 
        // There appears to be some problem with the CDC times when the SC is used as t0
        // This histogram should help diagnose the issue
        vector < const DCDCTrackHit *> cdcTrackHitVector;
        thisTrack->Get(cdcTrackHitVector);
        if (cdcTrackHitVector.size() != 0){
            float earliestTime = 10000; // Initialize high
            for (unsigned int iCDC = 0; iCDC < cdcTrackHitVector.size(); iCDC++){
                if (cdcTrackHitVector[iCDC]->tdrift < earliestTime) earliestTime = cdcTrackHitVector[iCDC]->tdrift;
            }

            Fill1DHistogram("HLDetectorTiming", "TRACKING", "Earliest CDC Time Minus Matched SC Time",
                    earliestTime - locSCHitMatchParams.dHitTime,
                    "Earliest CDC Time Minus Matched SC Time; t_{CDC} - t_{SC} [ns];",
                    200, -50, 150);
        }
        // Loop over TAGM hits
        for (unsigned int j = 0 ; j < tagmHitVector.size(); j++){
            int nTAGMColumns = 102;
            // We want to look at the timewalk within these ADC/TDC detectors
            Fill2DHistogram("HLDetectorTiming", "TRACKING", "TAGM - SC Target Time",
                    tagmHitVector[j]->column, tagmHitVector[j]->t - flightTimeCorrectedSCTime,
                    "#Deltat TAGM-SC; Column ;t_{TAGM} - t_{SC @ target} [ns]", nTAGMColumns, 0.5, nTAGMColumns + 0.5, 100, -50, 50);
            Fill2DHistogram("HLDetectorTiming", "TRACKING", "Tagger - SC Target Time",
                    tagmHitVector[j]->t - flightTimeCorrectedSCTime, tagmHitVector[j]->E,
                    "Tagger - SC Target Time; #Deltat_{Tagger - SC} [ns]; Energy [GeV]",
                    nBinsdT, dTMin, dTMax, nBinsE, EMin, EMax);   
            Fill1DHistogram("HLDetectorTiming", "TRACKING", "Tagger - SC 1D Target Time",
                    tagmHitVector[j]->t - flightTimeCorrectedSCTime,
                    "Tagger - SC Time at Target; #Deltat_{Tagger - SC} [ns]; Entries",
                    40, -20, 20);
        }
        // Loop over TAGH hits
        for (unsigned int j = 0 ; j < taghHitVector.size(); j++){
            int nTAGHCounters = 274;
            // We want to look at the timewalk within these ADC/TDC detectors
            Fill2DHistogram("HLDetectorTiming", "TRACKING", "TAGH - SC Target Time",
                    taghHitVector[j]->counter_id, taghHitVector[j]->t - flightTimeCorrectedSCTime,
                    "#Deltat TAGH-SC; Counter ID ;t_{TAGH} - t_{SC @ target} [ns]", nTAGHCounters, 0.5, nTAGHCounters + 0.5, 100, -50, 50);

            Fill2DHistogram("HLDetectorTiming", "TRACKING", "Tagger - SC Target Time",
                    taghHitVector[j]->t - flightTimeCorrectedSCTime, taghHitVector[j]->E,
                    "Tagger - SC Target Time; #Deltat_{Tagger - SC} [ns]; Energy [GeV]",
                    nBinsdT, dTMin, dTMax, nBinsE, EMin, EMax);

            Fill1DHistogram("HLDetectorTiming", "TRACKING", "Tagger - SC 1D Target Time",
                    taghHitVector[j]->t - flightTimeCorrectedSCTime,
                    "Tagger - SC Time at Target; #Deltat_{Tagger - SC} [ns]; Entries",
                    40, -20, 20);
        }

        if(thisTrack->charge() > 0) continue; // Cut on the charge to kill protons

        if (foundTOF){
            // Now check the TOF matching. Do this on a full detector level.
            float flightTimeCorrectedTOFTime = locTOFHitMatchParams.dHitTime - locTOFHitMatchParams.dFlightTime;
            Fill1DHistogram("HLDetectorTiming", "TRACKING", "TOF - SC Target Time",
                    flightTimeCorrectedTOFTime - flightTimeCorrectedSCTime,
                    "t_{TOF} - t_{SC} at Target; t_{TOF} - t_{SC} at Target [ns]; Entries",
                    200, -100, 100);
        }
        if (foundBCAL){
            float flightTimeCorrectedBCALTime = locBCALShowerMatchParams.dBCALShower->t - locBCALShowerMatchParams.dFlightTime;
            Fill1DHistogram("HLDetectorTiming", "TRACKING", "BCAL - SC Target Time",
                    flightTimeCorrectedBCALTime - flightTimeCorrectedSCTime,
                    "t_{BCAL} - t_{SC} at Target; t_{BCAL} - t_{SC} [ns]; Entries",
                    600, -300, 300);
            // Add histogram suggested by Mark Dalton
            Fill2DHistogram("HLDetectorTiming", "TRACKING", "BCAL - SC Target Time Vs Correction",
                    locBCALShowerMatchParams.dFlightTime, flightTimeCorrectedBCALTime - flightTimeCorrectedSCTime,
                    "t_{BCAL} - t_{SC} at Target; Flight time [ns]; t_{BCAL} - t_{SC} [ns]",
                    100, 0, 20, 50, -10, 10);
        }
        if (foundFCAL){
            float flightTimeCorrectedFCALTime = locFCALShowerMatchParams.dFCALShower->getTime() - locFCALShowerMatchParams.dFlightTime;
            Fill1DHistogram("HLDetectorTiming", "TRACKING", "FCAL - SC Target Time",
                    flightTimeCorrectedFCALTime - flightTimeCorrectedSCTime,
                    "t_{FCAL} - t_{SC} at Target; t_{FCAL} - t_{SC} [ns]; Entries",
                    600, -300, 300);
        }


    } // End of loop over time based tracks

    return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_HLDetectorTiming::erun(void)
{
    // This is called whenever the run number changes, before it is
    // changed to give you a chance to clean up before processing
    // events from the next run number.

    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_HLDetectorTiming::fini(void)
{
    // Called before program exit after event processing is finished.
    //Here is where we do the fits to the data to see if we have a reasonable alignment
    SortDirectories(); //Defined in HistogramTools.h

    if (DO_ROUGH_TIMING) DoRoughTiming();
    if (DO_TDC_ADC_ALIGN) DoTDCADCAlign();
    if (DO_TRACK_BASED) DoTrackBased();
    return NOERROR;
}
// Do the roughest level of timing for this run
void JEventProcessor_HLDetectorTiming::DoRoughTiming(void){

    //Fit all plots with expected funtional form, output files for CCDB input

    float TAGH_ADC_Offset = 0, TAGH_TDC_Offset = 0; 
    TH1I *thisHist = Get1DHistogram("HLDetectorTiming", "TAGH", "TAGHHit TDC time");
    if(thisHist != NULL){
        Double_t maximum = thisHist->GetBinCenter(thisHist->GetMaximumBin());
        TFitResultPtr fr = thisHist->Fit("gaus", "S", "", maximum - 50, maximum + 50);
        float mean = fr->Parameter(1);
        TAGH_TDC_Offset=mean;
    }

    thisHist = Get1DHistogram("HLDetectorTiming", "TAGH", "TAGHHit ADC time");
    if(thisHist != NULL){
        // In the tagger hodoscope the last few bins with data can be problematic. Ignore these last few in the fit
        Int_t lastBin = thisHist->FindLastBinAbove( 1 , 1); // Find last bin with content above 1 in the histogram
        for (int i = 0; i <= 12; i++){
            if ((lastBin - i) > 0) thisHist->SetBinContent((lastBin - i), 0);
        }
        Double_t maximum = thisHist->GetBinCenter(thisHist->GetMaximumBin());
        TFitResultPtr fr = thisHist->Fit("gaus", "S", "", maximum - 50, maximum + 50);
        float mean = fr->Parameter(1);
        TAGH_ADC_Offset=mean;

    }

    // Send offset to output file
    TAGH_ADC_Offset *= -1; TAGH_TDC_Offset *= -1; // Additive constants
    ofstream outFile;
    outFile.open("tagh_base_time.txt");
    outFile << TAGH_ADC_Offset << " " << TAGH_TDC_Offset << endl;
    outFile.close();

    float TAGM_ADC_Offset = 0, TAGM_TDC_Offset = 0;
    thisHist = Get1DHistogram("HLDetectorTiming", "TAGM", "TAGMHit ADC time");
    if(thisHist != NULL){
        // In the tagger microscope the last few bins with data can be problomatic. Ignore these last few in the fit
        Int_t lastBin = thisHist->FindLastBinAbove( 1 , 1); // Find last bin with content above 1 in the histogram
        for (int i = 0; i <= 12; i++){
            if ((lastBin - i) > 0) thisHist->SetBinContent((lastBin - i), 0);
        }
        Double_t maximum = thisHist->GetBinCenter(thisHist->GetMaximumBin());
        TFitResultPtr fr = thisHist->Fit("gaus", "S", "", maximum - 50, maximum + 50);
        float mean = fr->Parameter(1);
        TAGM_ADC_Offset=mean;
    }

    thisHist = Get1DHistogram("HLDetectorTiming", "TAGM", "TAGMHit TDC time");
    if(thisHist != NULL){
        Double_t maximum = thisHist->GetBinCenter(thisHist->GetMaximumBin());
        TFitResultPtr fr = thisHist->Fit("gaus", "S", "", maximum - 50, maximum + 50);
        float mean = fr->Parameter(1);
        TAGM_TDC_Offset=mean;
    }

    TAGM_ADC_Offset *= -1; TAGM_TDC_Offset *= -1;
    outFile.open("tagm_base_time.txt");
    outFile << TAGM_ADC_Offset << " " << TAGM_TDC_Offset << endl;
    outFile.close();

    float SC_ADC_Offset = 0, SC_TDC_Offset = 0;
    thisHist = Get1DHistogram("HLDetectorTiming", "SC", "SCHit ADC time");
    if(thisHist != NULL){
        //Gaussian
        Double_t maximum = thisHist->GetBinCenter(thisHist->GetMaximumBin());
        TFitResultPtr fr = thisHist->Fit("gaus", "S", "", maximum - 50, maximum + 50);
        float mean = fr->Parameter(1);
        SC_ADC_Offset = mean;
    }

    thisHist = Get1DHistogram("HLDetectorTiming", "SC", "SCHit TDC time");
    if(thisHist != NULL){
        //Gaussian
        Double_t maximum = thisHist->GetBinCenter(thisHist->GetMaximumBin());
        TFitResultPtr fr = thisHist->Fit("gaus", "S", "", maximum - 50, maximum + 50);
        float mean = fr->Parameter(1);
        SC_TDC_Offset = mean;
    }

    SC_ADC_Offset *= -1; SC_TDC_Offset*= -1;
    outFile.open("sc_base_time.txt");
    outFile << SC_ADC_Offset << " " << SC_TDC_Offset << endl;
    outFile.close();
    /*
       float BCAL_ADC_Offset = 0, BCAL_TDC_Offset = 0;
       thisHist = Get1DHistogram("HLDetectorTiming", "BCAL", "BCALHit time");
       if(thisHist != NULL){
    //Fit a gaussian, choose start time mean - 3sigma    
    Double_t maximum = thisHist->GetBinCenter(thisHist->GetMaximumBin());
    TFitResultPtr fr = thisHist->Fit("gaus", "S", "", maximum - 50, maximum + 50);
    float mean = fr->Parameter(1);
    float sigma = fr->Parameter(2);
    BCAL_ADC_Offset = mean - 3*sigma;
    }

    BCAL_ADC_Offset *= -1;
    outFile.open("bcal_base_time.txt");
    outFile << BCAL_ADC_Offset << " " << BCAL_TDC_Offset << endl;
    outFile.close();
    */
    float CDC_ADC_Offset = 0.0;
    thisHist = Get1DHistogram("HLDetectorTiming", "CDC", "CDCHit time");
    if(thisHist != NULL){
        //Fit a gaussian to the left of the main peak
        Double_t maximum = thisHist->GetBinCenter(thisHist->GetMaximumBin());
        TF1 *f = new TF1("f", "gaus");
        f->SetParameters(100, maximum, 20);
        f->FixParameter(1 , maximum);
        TFitResultPtr fr = thisHist->Fit(f, "S", "", maximum - 40, maximum + 10); // Cant fix value at end of range
        float mean = fr->Parameter(1);
        float sigma = fr->Parameter(2);
        CDC_ADC_Offset = mean - 3*sigma;
        delete f;
    }

    CDC_ADC_Offset *= -1;
    outFile.open("cdc_base_time.txt");
    outFile << CDC_ADC_Offset << endl;
    outFile.close();

    float FDC_ADC_Offset = 0.0, FDC_TDC_Offset = 0.0;
    thisHist = Get1DHistogram("HLDetectorTiming", "FDC", "FDCHit time");
    if(thisHist != NULL){
        //Fit a gaussian to the left of the main peak
        Double_t maximum = thisHist->GetBinCenter(thisHist->GetMaximumBin());
        TF1 *f = new TF1("f", "gaus");
        f->SetParameters(100, maximum, 20);
        f->FixParameter(1 , maximum);
        TFitResultPtr fr = thisHist->Fit(f, "S", "", maximum - 40, maximum + 10); // Cant fix value at end of range
        float mean = fr->Parameter(1);
        float sigma = fr->Parameter(2);
        FDC_ADC_Offset = mean - 3*sigma;
        delete f;
    }

    FDC_ADC_Offset *= -1;
    outFile.open("fdc_base_time.txt");
    outFile << FDC_ADC_Offset << " " << FDC_TDC_Offset << endl;
    outFile.close();
    /*
       float FCAL_ADC_Offset = 0.0;
       thisHist = Get1DHistogram("HLDetectorTiming", "FCAL", "FCALHit time");
       if(thisHist != NULL){
       Double_t maximum = thisHist->GetBinCenter(thisHist->GetMaximumBin());
       TFitResultPtr fr = thisHist->Fit("gaus", "S", "", maximum - 20, maximum + 20);
       float mean = fr->Parameter(1);
       float sigma = fr->Parameter(2);
       FCAL_ADC_Offset = mean - 3*sigma;
       }

       FCAL_ADC_Offset *= -1;
       outFile.open("fcal_base_time.txt");
       outFile << FCAL_ADC_Offset << endl;
       outFile.close();
       */
    float TOF_ADC_Offset = 0.0, TOF_TDC_Offset = 0.0;
    thisHist = Get1DHistogram("HLDetectorTiming", "TOF", "TOFHit ADC time");
    if(thisHist != NULL){
        //Fit a gaussian
        Double_t maximum = thisHist->GetBinCenter(thisHist->GetMaximumBin());
        TFitResultPtr fr = thisHist->Fit("gaus", "S", "", maximum - 20, maximum + 20);
        float mean = fr->Parameter(1);
        //float sigma = fr->Parameter(2);
        TOF_ADC_Offset = mean - 20; //Guess at 20ns
    }

    thisHist = Get1DHistogram("HLDetectorTiming", "TOF", "TOFHit TDC time");
    if(thisHist != NULL){
        Double_t maximum = thisHist->GetBinCenter(thisHist->GetMaximumBin());
        TFitResultPtr fr = thisHist->Fit("gaus", "S", "", maximum - 20, maximum + 20);
        float mean = fr->Parameter(1);
        //float sigma = fr->Parameter(2);
        TOF_TDC_Offset = mean - 20; //Guess at 20ns
    }

    TOF_ADC_Offset *= -1; TOF_TDC_Offset *= -1;
    outFile.open("tof_base_time.txt");
    outFile << TOF_ADC_Offset << " " << TOF_TDC_Offset << endl;
    outFile.close();

}

void JEventProcessor_HLDetectorTiming::DoTDCADCAlign(void){

    //Do a finer alignement of the TDC and ADC's in the detectors that have both
    int minHits = 7;

    // Use incremental shifts in case you want to run the alignment twice
    ofstream outFile;
    TH2I *thisHist = Get2DHistogram("HLDetectorTiming", "SC", "SCHit TDC_ADC Difference");
    if(thisHist != NULL){
        TObjArray histArray;
        thisHist->FitSlicesY(0, 0, -1, minHits, "QNR", &histArray); // require at least minHits entries to do the fit
        TH1D * meanHist = (TH1D *) histArray[1];
        int nbins = meanHist->GetNbinsX();
        outFile.open("sc_tdc_timing_offsets.txt");
        for (int i = 1; i <= nbins; i++){
            double mean = meanHist->GetBinContent(i);
            // Don't allow any large outliers if the fit fails
            if (fabs(mean) > 10.0) mean = 0.0;
            outFile << mean + sc_tdc_time_offsets[i - 1] << endl; // Vector indexed from zero
        }
        outFile.close();
    }

    thisHist = Get2DHistogram("HLDetectorTiming", "TOF", "TOFHit TDC_ADC Difference");
    if(thisHist != NULL){
        // The FitSlicesY routine isn't working great here. Will handle this like the TAGM/TAGH.
        /*
           TObjArray histArray;
           thisHist->FitSlicesY(0, 0, -1, minHits, "QNR", &histArray); 
           TH1D * meanHist = (TH1D *) histArray[1];
           int nbins = meanHist->GetNbinsX();
           outFile.open("tof_tdc_timing_offsets.txt");
           for (int i = 1; i <= nbins; i++){
           double mean = meanHist->GetBinContent(i);
        // Don't allow any large outliers if the fit fails
        if (fabs(mean) > 10.0) mean = 0.0;
        outFile << mean + tof_tdc_time_offsets[i - 1]<< endl; // Vector indexed from zero
        }
        outFile.close();
        */
        int nBinsX = thisHist->GetNbinsX();
        int nBinsY = thisHist->GetNbinsY();
        TH1D * selected_TOF_TDCADCOffset = new TH1D("selected_TOF_TDCADCOffset", "Selected TOF TDC/ADC Offset; CCDB Index; Offset [ns]", nBinsX, 0.5, nBinsX + 0.5);

        outFile.open("tof_tdc_timing_offsets.txt");
        outFile.close(); // Clear Output File

        for (int i = 1 ; i <= nBinsX; i++){
            TH1D *projY = thisHist->ProjectionY("temp", i, i);
            // Scan over the histogram
            double maxEntries = 0;
            double maxMean = 0;
            for (int j = 1 ; j <= projY->GetNbinsX();j++){
                int minBin = j;
                int maxBin = (j + (projY->GetNbinsX() / 15) <= projY->GetNbinsX()) ? (j + projY->GetNbinsX() / 15) : projY->GetNbinsX();
                double sum = 0, nEntries = 0;
                for (int bin = minBin; bin <= maxBin; bin++){
                    sum += projY->GetBinContent(bin) * projY->GetBinCenter(bin);
                    nEntries += projY->GetBinContent(bin);
                    if (bin == maxBin){
                        if (nEntries > maxEntries) {
                            maxMean = sum / nEntries;
                            maxEntries = nEntries;
                        }
                    }
                }
            }
            outFile.open("tof_tdc_timing_offsets.txt", ios::out | ios::app);
            outFile << maxMean + tof_tdc_time_offsets[i - 1]<< endl;
            outFile.close();
            // Some histograms for tracking the magnitude of the shifts
            selected_TOF_TDCADCOffset->SetBinContent(i, maxMean);
            Fill2DHistogram("HLDetectorTiming", "TOF", "Selected TOF TDC/ADC Offset",
                    i, maxMean,
                    "Selected Mean position for TOF TDC/ADC offset",
                    nBinsX, 0.5, nBinsX + 0.5, nBinsY, projY->GetBinCenter(0), projY->GetBinCenter(projY->GetNbinsX()));
        }
    }

    thisHist = Get2DHistogram("HLDetectorTiming", "TAGM", "TAGMHit TDC_ADC Difference");
    if(thisHist != NULL){
        TObjArray histArray;
        thisHist->FitSlicesY(0, 0, -1, minHits, "QNR", &histArray); 
        TH1D * meanHist = (TH1D *) histArray[1];
        int nbins = meanHist->GetNbinsX();
        outFile.open("tagm_tdc_timing_offsets.txt");
        for (int i = 1; i <= nbins; i++){
            double mean = meanHist->GetBinContent(i);
            outFile << "0 " << i << " " << mean << endl;
            if (i == 7 || i == 25 || i == 79 || i == 97){
                for(int j = 1; j <= 5; j++){
                    outFile << j << " " << i << " " << mean + tagm_tdc_time_offsets[i] << endl;
                }
            }
        }
        outFile.close();
    }

    thisHist = Get2DHistogram("HLDetectorTiming", "TAGH", "TAGHHit TDC_ADC Difference");
    if(thisHist != NULL){
        TObjArray histArray;
        thisHist->FitSlicesY(0, 0, -1, minHits, "QNR", &histArray); 
        TH1D * meanHist = (TH1D *) histArray[1];
        int nbins = meanHist->GetNbinsX();
        outFile.open("tagh_tdc_timing_offsets.txt");
        for (int i = 1; i <= nbins; i++){
            double mean = meanHist->GetBinContent(i);
            outFile << i << " " << mean + tagh_tdc_time_offsets[i] << endl;
        }
        outFile.close();
    }

    return;
}

void JEventProcessor_HLDetectorTiming::DoTrackBased(void){

    // Do our final step in the timing alignment with tracking
    ofstream outFile;
    TH2I *thisHist = Get2DHistogram("HLDetectorTiming", "TRACKING", "TAGM - SC Target Time");
    if (thisHist != NULL){
        //Statistics on these histograms are really quite low we will have to rebin and do some interpolation
        outFile.open("tagm_tdc_timing_offsets.txt", ios::out | ios::trunc);
        outFile.close(); // clear file
        outFile.open("tagm_adc_timing_offsets.txt", ios::out | ios::trunc);
        outFile.close(); // clear file
        int nBinsX = thisHist->GetNbinsX();
        int nBinsY = thisHist->GetNbinsY();
        TH1D * selectedTAGMOffset = new TH1D("selectedTAGMOffset", "Selected TAGM Offset; Column; Offset [ns]", nBinsX, 0.5, nBinsX + 0.5);
        for (int i = 1 ; i <= nBinsX; i++){ 
            TH1D *projY = thisHist->ProjectionY("temp", i, i);
            // Scan over the histogram
            double maxEntries = 0;
            double maxMean = 0;
            for (int j = 1 ; j <= projY->GetNbinsX();j++){
                int minBin = j;
                int maxBin = (j + (projY->GetNbinsX() / 15) <= projY->GetNbinsX()) ? (j + projY->GetNbinsX() / 15) : projY->GetNbinsX();
                double sum = 0, nEntries = 0;
                for (int bin = minBin; bin <= maxBin; bin++){
                    sum += projY->GetBinContent(bin) * projY->GetBinCenter(bin);
                    nEntries += projY->GetBinContent(bin);
                    if (bin == maxBin){
                        if (nEntries > maxEntries) {
                            maxMean = sum / nEntries;
                            maxEntries = nEntries;
                        }
                    } 
                }
            }
            selectedTAGMOffset->SetBinContent(i, maxMean);
            Fill2DHistogram("HLDetectorTiming", "TRACKING", "Selected TAGM Offset",
                    i, maxMean,
                    "Selected Mean position for TAGM offset",
                    nBinsX, 0.5, nBinsX + 0.5, nBinsY, projY->GetBinCenter(0), projY->GetBinCenter(projY->GetNbinsX()));
        }

        TFitResultPtr fr1 = selectedTAGMOffset->Fit("pol2", "SQ", "", 0.5, nBinsX + 0.5);

        for (int i = 1 ; i <= nBinsX; i++){
            double x0 = fr1->Parameter(0);
            double x1 = fr1->Parameter(1);
            double x2 = fr1->Parameter(2);
            double fitResult = x0 + i*x1 + i*i*x2;

            double outlierCut = 7;
            double valueToUse = selectedTAGMOffset->GetBinContent(i);
            if (fabs(selectedTAGMOffset->GetBinContent(i) - fitResult) > outlierCut && valueToUse != 0.0){
                valueToUse = fitResult;
            }

            selectedTAGMOffset->SetBinContent(i, valueToUse);

            outFile.open("tagm_tdc_timing_offsets.txt", ios::out | ios::app);
            outFile << "0 " << i << " " << valueToUse + tagm_tdc_time_offsets[i] << endl;
            if (i == 7 || i == 25 || i == 79 || i == 97){
                for(int j = 1; j <= 5; j++){
                    outFile << j << " " << i << " " << valueToUse + tagm_tdc_time_offsets[i] << endl;
                }
            }
            outFile.close();
            // Apply the same shift to the adc offsets
            outFile.open("tagm_adc_timing_offsets.txt", ios::out | ios::app);
            outFile << "0 " << i << " " << valueToUse + tagm_fadc_time_offsets[i] << endl;
            if (i == 7 || i == 25 || i == 79 || i == 97){
                for(int j = 1; j <= 5; j++){
                    outFile << j << " " << i << " " << valueToUse + tagm_fadc_time_offsets[i] << endl;
                }
            }
            outFile.close();
        }
    }

    thisHist = Get2DHistogram("HLDetectorTiming", "TRACKING", "TAGH - SC Target Time");
    if(thisHist != NULL){
        outFile.open("tagh_tdc_timing_offsets.txt", ios::out | ios::trunc);
        outFile.close(); // clear file
        outFile.open("tagh_adc_timing_offsets.txt", ios::out | ios::trunc);
        outFile.close(); // clear file

        int nBinsX = thisHist->GetNbinsX();
        int nBinsY = thisHist->GetNbinsY();
        TH1D * selectedTAGHOffset = new TH1D("selectedTAGHOffset", "Selected TAGH Offset; ID; Offset [ns]", nBinsX, 0.5, nBinsX + 0.5);
        for (int i = 1 ; i <= nBinsX; i++){
            TH1D *projY = thisHist->ProjectionY("temp", i, i);
            // Scan over the histogram
            double maxEntries = 0;
            double maxMean = 0;
            for (int j = 1 ; j <= projY->GetNbinsX();j++){
                int minBin = j;
                int maxBin = (j + (projY->GetNbinsX() / 15) <= projY->GetNbinsX()) ? (j + projY->GetNbinsX() / 15) : projY->GetNbinsX();
                double sum = 0; 
                double nEntries = 0;
                for (int bin = minBin; bin <= maxBin; bin++){
                    sum += projY->GetBinContent(bin) * projY->GetBinCenter(bin);
                    nEntries += projY->GetBinContent(bin);
                    if (bin == maxBin){
                        if (nEntries > maxEntries){
                            maxMean = sum / nEntries;
                            maxEntries = nEntries;
                        }
                    }
                }
            }

            Fill2DHistogram("HLDetectorTiming", "TRACKING", "Selected TAGH Offset",
                    i, maxMean,
                    "Selected Mean position for offset",
                    nBinsX, 0.5, nBinsX + 0.5, nBinsY, projY->GetBinCenter(0), projY->GetBinCenter(projY->GetNbinsX()));
            selectedTAGHOffset->SetBinContent(i, maxMean);
            /*
               outFile.open("tagh_tdc_timing_offsets.txt", ios::out | ios::app);
               outFile << i << " " << maxMean + tagh_tdc_time_offsets[i] << endl;
               outFile.close();
               outFile.open("tagh_adc_timing_offsets.txt", ios::out | ios::app);
               outFile << i << " " << maxMean + tagh_fadc_time_offsets[i] << endl;
               outFile.close();
               */
        }

        // Fit 1D histogram. If value is far from the fit use the fitted value
        // Two behaviors above and below microscope
        TFitResultPtr fr1 = selectedTAGHOffset->Fit("pol2", "SQ", "", 0.5, 131.5);
        TFitResultPtr fr2 = selectedTAGHOffset->Fit("pol2", "SQ", "", 182.5, 274.5);        

        for (int i = 1 ; i <= nBinsX; i++){
            double fitResult = 0.0;
            if (i < 150){
                double x0 = fr1->Parameter(0);
                double x1 = fr1->Parameter(1);
                double x2 = fr1->Parameter(2);
                fitResult = x0 + i*x1 + i*i*x2;
            }
            else{
                double x0 = fr2->Parameter(0);
                double x1 = fr2->Parameter(1);
                double x2 = fr2->Parameter(2);
                fitResult = x0 + i*x1 + i*i*x2;
            }

            double outlierCut = 7;
            double valueToUse = selectedTAGHOffset->GetBinContent(i);
            if (fabs(selectedTAGHOffset->GetBinContent(i) - fitResult) > outlierCut && valueToUse != 0.0){
                valueToUse = fitResult;
            }

            selectedTAGHOffset->SetBinContent(i, valueToUse);

            outFile.open("tagh_tdc_timing_offsets.txt", ios::out | ios::app);
            outFile << i << " " << valueToUse + tagh_tdc_time_offsets[i] << endl;
            outFile.close();
            outFile.open("tagh_adc_timing_offsets.txt", ios::out | ios::app);
            outFile << i << " " << valueToUse + tagh_fadc_time_offsets[i] << endl;
            outFile.close();
        }
    }

    TH1I *this1DHist = Get1DHistogram("HLDetectorTiming", "TRACKING", "TOF - SC Target Time");
    if(this1DHist != NULL){
        //Gaussian
        Double_t maximum = this1DHist->GetBinCenter(this1DHist->GetMaximumBin());
        TFitResultPtr fr = this1DHist->Fit("gaus", "S", "", maximum - 10, maximum + 10);
        float mean = fr->Parameter(1);
        outFile.open("tof_base_time.txt");
        outFile << tof_t_base_fadc - mean << " " << tof_t_base_tdc - mean << endl;
        outFile.close();
    }

    this1DHist = Get1DHistogram("HLDetectorTiming", "TRACKING", "BCAL - SC Target Time");
    if(this1DHist != NULL){
        //Gaussian
        Double_t maximum = this1DHist->GetBinCenter(this1DHist->GetMaximumBin());
        TFitResultPtr fr = this1DHist->Fit("gaus", "S", "", maximum - 10, maximum + 10);
        float mean = fr->Parameter(1);
        outFile.open("bcal_base_time.txt");
        outFile << bcal_t_base - mean << " 0.0" << endl; // TDC info not used
        outFile.close();
    }

    this1DHist = Get1DHistogram("HLDetectorTiming", "TRACKING", "FCAL - SC Target Time");
    if(this1DHist != NULL){
        //Gaussian
        Double_t maximum = this1DHist->GetBinCenter(this1DHist->GetMaximumBin());
        TFitResultPtr fr = this1DHist->Fit("gaus", "S", "", maximum - 10, maximum + 10);
        float mean = fr->Parameter(1);
        outFile.open("fcal_base_time.txt");
        outFile << fcal_t_base - mean << endl; 
        outFile.close();
    }

    this1DHist = Get1DHistogram("HLDetectorTiming", "TRACKING", "Earliest CDC Time Minus Matched SC Time");
    if(this1DHist != NULL){
        //Gaussian
        Double_t maximum = this1DHist->GetBinCenter(this1DHist->GetMaximumBin());
        TFitResultPtr fr = this1DHist->Fit("gaus", "S", "", maximum - 15, maximum + 10);
        float mean = fr->Parameter(1);
        outFile.open("cdc_base_time.txt");
        outFile << cdc_t_base - mean << endl;
        outFile.close();
    }

}

int JEventProcessor_HLDetectorTiming::GetCCDBIndexTOF(const DTOFHit *thisHit){
    // Returns the CCDB index of a particular hit
    // This 
    int plane = thisHit->plane;
    int bar = thisHit->bar;
    int end = thisHit->end;
    // 44 bars per plane
    int CCDBIndex = plane * 88 + end * 44 + bar; 
    return CCDBIndex;
}

int JEventProcessor_HLDetectorTiming::GetCCDBIndexBCAL(const DBCALHit *thisHit){
    return 0;
}

int JEventProcessor_HLDetectorTiming::GetCCDBIndexCDC(const DCDCHit *thisHit){

    int ring = thisHit->ring;
    int straw = thisHit->straw;

    int CCDBIndex = GetCCDBIndexCDC(ring, straw);
    return CCDBIndex;
}

int JEventProcessor_HLDetectorTiming::GetCCDBIndexCDC(int ring, int straw){

    //int Nstraws[28] = {42, 42, 54, 54, 66, 66, 80, 80, 93, 93, 106, 106, 123, 123, 135, 135, 146, 146, 158, 158, 170, 170, 182, 182, 197, 197, 209, 209};
    int StartIndex[28] = {0, 42, 84, 138, 192, 258, 324, 404, 484, 577, 670, 776, 882, 1005, 1128, 1263, 1398, 1544, 1690, 1848, 2006, 2176, 2346, 2528, 2710, 2907, 3104, 3313};

    int CCDBIndex = StartIndex[ring - 1] + straw;
    return CCDBIndex;
}

