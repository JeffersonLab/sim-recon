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
#include "PID/DEventRFBunch.h"
#include "TTAB/DTTabUtilities.h"
#include "TTAB/DTranslationTable.h"
#include "FCAL/DFCALGeometry.h"
#include "BCAL/DBCALGeometry.h"
#include "HistogramTools.h"

extern "C"{
void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_HLDetectorTiming());
    app->AddFactoryGenerator(new DFactoryGenerator_p2pi()); //register the factory generator
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
    REQUIRE_BEAM = 0;
    BEAM_EVENTS_TO_KEEP = 1000000000; // Set enormously high
    DO_ROUGH_TIMING = 0;
    DO_CDC_TIMING = 0;
    DO_TDC_ADC_ALIGN = 0;
    DO_TRACK_BASED = 0;
    DO_VERIFY = 1;
    DO_OPTIONAL = 0;
    DO_REACTION = 0;
    USE_RF_BUNCH = 1;

    if(gPARMS){
        gPARMS->SetDefaultParameter("HLDETECTORTIMING:DO_ROUGH_TIMING", DO_ROUGH_TIMING, "Set to > 0 to do rough timing of all detectors");
        gPARMS->SetDefaultParameter("HLDETECTORTIMING:DO_CDC_TIMING", DO_CDC_TIMING, "Set to > 0 to do CDC Per channel Alignment");
        gPARMS->SetDefaultParameter("HLDETECTORTIMING:DO_TDC_ADC_ALIGN", DO_TDC_ADC_ALIGN, "Set to > 0 to do TDC/ADC alignment of SC,TOF,TAGM,TAGH");
        gPARMS->SetDefaultParameter("HLDETECTORTIMING:DO_TRACK_BASED", DO_TRACK_BASED, "Set to > 0 to do Track Based timing corrections");
        gPARMS->SetDefaultParameter("HLDETECTORTIMING:DO_VERIFY", DO_VERIFY, "Set to > 0 to verify timing with current constants");
        gPARMS->SetDefaultParameter("HLDETECTORTIMING:REQUIRE_BEAM", REQUIRE_BEAM, "Set to 0 to skip beam current check");
        gPARMS->SetDefaultParameter("HLDETECTORTIMING:BEAM_EVENTS_TO_KEEP", BEAM_EVENTS_TO_KEEP, "Set to the number of beam on events to use");
        gPARMS->SetDefaultParameter("HLDETECTORTIMING:DO_OPTIONAL", DO_OPTIONAL, "Set to >0 to enable optional histograms ");
        gPARMS->SetDefaultParameter("HLDETECTORTIMING:DO_REACTION", DO_REACTION, "Set to >0 to run DReaction");
        gPARMS->SetDefaultParameter("HLDETECTORTIMING:USE_RF_BUNCH", USE_RF_BUNCH, "Set to 0 to disable use of 2 vote RF Bunch");
    }

    // Would like the code with no arguments to simply verify the current status of the calibration
    if (DO_ROUGH_TIMING > 0 || DO_CDC_TIMING > 0 || DO_TDC_ADC_ALIGN > 0 || DO_TRACK_BASED > 0) DO_VERIFY = 0;

    // Increase range for initial search
    if(DO_TDC_ADC_ALIGN){
        NBINS_TDIFF = 1400; MIN_TDIFF = -150.0; MAX_TDIFF = 550.0;
    }
    else{
        NBINS_TDIFF = 50; MIN_TDIFF = -10.0; MAX_TDIFF = 10.0;
    }

    if (DO_TRACK_BASED){
        NBINS_TAGGER_TIME = 1600; MIN_TAGGER_TIME = 0; MAX_TAGGER_TIME = 400;
        NBINS_MATCHING = 500; MIN_MATCHING_T = -100; MAX_MATCHING_T = 400;
    }
    else{
        NBINS_TAGGER_TIME = 100; MIN_TAGGER_TIME = -50; MAX_TAGGER_TIME = 50;
        NBINS_MATCHING = 100; MIN_MATCHING_T = -25; MAX_MATCHING_T = 25;
    }
    NBINS_RF_COMPARE = 200; MIN_RF_COMPARE = -2.2; MAX_RF_COMPARE = 2.2;
    return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_HLDetectorTiming::brun(JEventLoop *eventLoop, int runnumber)
{
    // This is called whenever the run number changes
    // Get the particleID object for each run
    vector<const DParticleID *> dParticleID_algos;
    eventLoop->Get(dParticleID_algos);
    if(dParticleID_algos.size()<1){
        _DBG_<<"Unable to get a DParticleID object! NO PID will be done!"<<endl;
        return RESOURCE_UNAVAILABLE;
    }
    dParticleID = dParticleID_algos[0];

    // We want to be use some of the tools available in the RFTime factory 
    // Specifivally steping the RF back to a chosen time
    dRFTimeFactory = static_cast<DRFTime_factory*>(eventLoop->GetFactory("DRFTime"));

    //be sure that DRFTime_factory::init() and brun() are called
    vector<const DRFTime*> locRFTimes;
    eventLoop->Get(locRFTimes);

    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_HLDetectorTiming::evnt(JEventLoop *loop, int eventnumber)
{

    // Get the EPICs events and update beam current. Skip event if current too low (<10 nA).
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
    // There is a caveat here when running multithreaded
    // Another thread might be the one to catch the EPICS event
    // and there is no way to reject events that may have come from earilier
    // Expect number of entries in histograms to vary slightly over the same file with many threads
    if (BEAM_CURRENT < 10.0) {
        Fill1DHistogram("HLDetectorTiming", "" , "Beam Events",
                0, "Beam On Events (0 = no beam, 1 = beam > 10nA)",
                2, -0.5, 1.5);
        if (REQUIRE_BEAM){
            return NOERROR; // Skip events where we can't verify the beam current
        }
    }
    else{
        Fill1DHistogram("HLDetectorTiming", "" , "Beam Events",
                1, "Beam On Events (0 = no beam, 1 = beam > 10nA)",
                2, -0.5, 1.5);
        fBeamEventCounter++;
    }

    if (fBeamEventCounter >= BEAM_EVENTS_TO_KEEP) { // Able to specify beam ON events instead of just events
        cout<< "Maximum number of Beam Events reached" << endl;
        japp->Quit();
        return NOERROR;
    }


    vector<const DCDCHit *> cdcHitVector;
    loop->Get(cdcHitVector);

    vector<const DFDCHit *> fdcHitVector;
    loop->Get(fdcHitVector);

    vector<const DSCHit *> scHitVector;
    loop->Get(scHitVector);

    // Grab the unified hits since we are also interested in the TDC times
    //vector<const DBCALHit *> bcalHitVector;
    //loop->Get(bcalHitVector);

    vector<const DBCALUnifiedHit *> bcalUnifiedHitVector;
    loop->Get(bcalUnifiedHitVector);

    vector<const DTOFHit *> tofHitVector;
    loop->Get(tofHitVector);

    vector<const DFCALHit *> fcalHitVector;
    loop->Get(fcalHitVector);

    vector<const DTAGMHit *> tagmHitVector;
    loop->Get(tagmHitVector);

    vector<const DTAGHHit *> taghHitVector;
    loop->Get(taghHitVector);

    vector<const DRFTDCDigiTime *> tdcRFTimeVector;
    loop->Get(tdcRFTimeVector);

    // TTabUtilities object used for RF time conversion
    const DTTabUtilities* locTTabUtilities = NULL;
    loop->GetSingle(locTTabUtilities);

    const DRFTime* thisRFTime = NULL;
    vector <const DRFTime*> RFTimeVector;
    loop->Get(RFTimeVector);
    if (RFTimeVector.size() != 0){
        thisRFTime = RFTimeVector[0];
    }

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
                /*
                   Fill2DHistogram("HLDetectorTiming", "CDC", "CDC Energy Vs. Drift Time",
                   cdcHitVector[i]->t, cdcHitVector[i]->q , 
                   "CDC Energy Vs. Drift Time; Drift Time [ns]; Integral [arb.]",
                   800, 0, 800, 500, -10000, 7000000);
                   */
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
    for (i = 0; i < bcalUnifiedHitVector.size(); i++){
        int the_cell = (bcalUnifiedHitVector[i]->module - 1) * 16 + (bcalUnifiedHitVector[i]->layer - 1) * 4 + bcalUnifiedHitVector[i]->sector;
        // There is one less layer of TDCs so the numbering relects this
        int the_tdc_cell = (bcalUnifiedHitVector[i]->module - 1) * 12 + (bcalUnifiedHitVector[i]->layer - 1) * 4 + bcalUnifiedHitVector[i]->sector;
        // Get the underlying associated objects
        const DBCALHit * thisADCHit;
        const DBCALTDCHit * thisTDCHit;
        bcalUnifiedHitVector[i]->GetSingle(thisADCHit);
        bcalUnifiedHitVector[i]->GetSingle(thisTDCHit);

        if (thisADCHit != NULL){ //This should never be NULL but might as well check
            Fill1DHistogram ("HLDetectorTiming", "BCAL", "BCALHit ADC time", thisADCHit->t,
                    "BCALHit ADC time; t_{ADC} [ns]; Entries", nBins, xMin, xMax);

            if (DO_OPTIONAL){
                if (bcalUnifiedHitVector[i]->end == 0){
                    Fill2DHistogram ("HLDetectorTiming", "BCAL", "BCALHit Upstream Per Channel ADC Hit Time",
                            the_cell, thisADCHit->t,
                            "BCALHit Upstream Per Channel Hit Time; cellID; t_{ADC} [ns] ",
                            768, 0.5, 768.5, 250, -50, 50);
                }
                else{
                    Fill2DHistogram ("HLDetectorTiming", "BCAL", "BCALHit Downstream Per Channel ADC Hit Time",
                            the_cell, thisADCHit->t,
                            "BCALHit Downstream Per Channel Hit Time; cellID; t_{ADC} [ns] ",
                            768, 0.5, 768.5, 250, -50, 50);
                }
            }
        }

        if (thisTDCHit != NULL){
            Fill1DHistogram ("HLDetectorTiming", "BCAL", "BCALHit TDC time", thisTDCHit->t,
                    "BCALHit TDC time; t_{TDC} [ns]; Entries", nBins, xMin, xMax);

            if (DO_OPTIONAL){
                if (bcalUnifiedHitVector[i]->end == 0){
                    Fill2DHistogram ("HLDetectorTiming", "BCAL", "BCALHit Upstream Per Channel TDC Hit Time",
                            the_tdc_cell, thisTDCHit->t,
                            "BCALHit Upstream Per Channel TDC Hit Time; cellID; t_{TDC} [ns] ",
                            576, 0.5, 576.5, 350, -50, 300);
                }
                else{
                    Fill2DHistogram ("HLDetectorTiming", "BCAL", "BCALHit Downstream Per Channel TDC Hit Time",
                            the_tdc_cell, thisTDCHit->t,
                            "BCALHit Downstream Per Channel TDC Hit Time; cellID; t_{TDC} [ns] ",
                            576, 0.5, 576.5, 350, -50, 300);
                }
            }
        }

        if (thisADCHit != NULL && thisTDCHit != NULL){
            if (bcalUnifiedHitVector[i]->end == 0){
                Fill2DHistogram ("HLDetectorTiming", "BCAL", "BCALHit Upstream Per Channel TDC-ADC Hit Time",
                        the_tdc_cell, thisTDCHit->t - thisADCHit->t,
                        "BCALHit Upstream Per Channel TDC-ADC Hit Time; cellID; t_{TDC} - t_{ADC} [ns] ",
                        576, 0.5, 576.5, NBINS_TDIFF, MIN_TDIFF, MAX_TDIFF);
            }
            else{
                Fill2DHistogram ("HLDetectorTiming", "BCAL", "BCALHit Downstream Per Channel TDC-ADC Hit Time",
                        the_tdc_cell, thisTDCHit->t - thisADCHit->t,
                        "BCALHit Downstream Per Channel TDC-ADC Hit Time; cellID; t_{TDC} - t_{ADC} [ns] ",
                        576, 0.5, 576.5, NBINS_TDIFF, MIN_TDIFF, MAX_TDIFF);
            }
        }
    }

    for (i = 0; i < tofHitVector.size(); i++){
        Fill1DHistogram ("HLDetectorTiming", "TOF", "TOFHit time", tofHitVector[i]->t,
                "TOFHit time", nBins, xMin, xMax);
    }
    for (i = 0; i < fcalHitVector.size(); i++){
        Fill1DHistogram ("HLDetectorTiming", "FCAL", "FCALHit time", fcalHitVector[i]->t,
                "FCALHit time", nBins, xMin, xMax);
        if (DO_OPTIONAL){
            // extract the FCAL Geometry
            vector<const DFCALGeometry*> fcalGeomVect;
            loop->Get( fcalGeomVect );
            if (fcalGeomVect.size() < 1){
                cout << "FCAL Geometry not available?" << endl;
                return OBJECT_NOT_AVAILABLE;
            }
            const DFCALGeometry& fcalGeom = *(fcalGeomVect[0]);
            Fill2DHistogram("HLDetectorTiming", "FCAL", "FCALHit Per Channel Time",
                    fcalGeom.channel(fcalHitVector[i]->row, fcalHitVector[i]->column), fcalHitVector[i]->t,
                    "FCAL Per Channel Hit time; channel; t [ns]",
                    fcalGeom.numActiveBlocks(), 0.5, fcalGeom.numActiveBlocks() + 0.5, 250, -50, 50); 
        }
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
        if(thisSCHit->has_fADC && thisRFTime != NULL){
            char name [200];
            char title[500];
            sprintf(name, "Sector %.2i", scHitVector[i]->sector);
            sprintf(title, "SC Sector %i t_{ADC} - t_{RF}; t_{ADC} - t_{RF} [ns]; Entries", scHitVector[i]->sector);
            double locShiftedTime = dRFTimeFactory->Step_TimeToNearInputTime(thisRFTime->dTime, scHitVector[i]->t_fADC);
            Fill1DHistogram("HLDetectorTiming", "SC_ADC_RF_Compare", name,
                    scHitVector[i]->t_fADC - locShiftedTime,
                    title,
                    NBINS_RF_COMPARE, MIN_RF_COMPARE, MAX_RF_COMPARE);
        }

        if(thisSCHit->has_TDC && thisRFTime != NULL){
            char name [200];
            char title[500];
            sprintf(name, "Sector %.2i", scHitVector[i]->sector);
            sprintf(title, "SC Sector %i t_{TDC} - t_{RF}; t_{TDC} - t_{RF} [ns]; Entries", scHitVector[i]->sector);
            double locShiftedTime = dRFTimeFactory->Step_TimeToNearInputTime(thisRFTime->dTime, scHitVector[i]->t_TDC);
            Fill1DHistogram("HLDetectorTiming", "SC_TDC_RF_Compare", name,
                    scHitVector[i]->t_TDC - locShiftedTime,
                    title,
                    NBINS_RF_COMPARE, MIN_RF_COMPARE, MAX_RF_COMPARE);
        }
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
                    "SC #Deltat TDC-ADC; Sector ;t_{TDC} - t_{ADC} [ns]", nSCCounters, 0.5, nSCCounters + 0.5, NBINS_TDIFF, MIN_TDIFF, MAX_TDIFF);
            if (DO_TRACK_BASED){
                char name [200];
                char title[500];
                sprintf(name, "Sector %.2i", scHitVector[i]->sector);
                sprintf(title, "SC Sector %i #Deltat Vs dE; dE [Arb.]; t_{TDC} - t_{ADC} [ns]", scHitVector[i]->sector);
                Fill2DHistogram("HLDetectorTiming", "SC_TimeWalk", name,
                        scHitVector[i]->dE, scHitVector[i]->t_TDC - scHitVector[i]->t_fADC,
                        title,
                        100, 0, 0.0075, 100, -10, 10, false);
            }
            if (DO_OPTIONAL){
                Fill2DHistogram("HLDetectorTiming", "SC", "SC Per Channel TDC Time",
                        scHitVector[i]->sector, scHitVector[i]->t_TDC,
                        "SC t_{TDC} ; Sector; t_{TDC} [ns]",
                        nSCCounters, 0.5, nSCCounters + 0.5, 100, -50, 50);
            }
        }

    }
    for (i = 0; i < tagmHitVector.size(); i++){
        const DTAGMHit *thisTAGMHit = tagmHitVector[i];
        if(thisTAGMHit->has_fADC && thisRFTime != NULL){
            char name [200];
            char title[500];
            sprintf(name, "Column %.3i Row %.1i", tagmHitVector[i]->column, tagmHitVector[i]->row);
            sprintf(title, "TAGM Column %i t_{ADC} - t_{RF}; t_{ADC} - t_{RF} [ns]; Entries", tagmHitVector[i]->column);
            double locShiftedTime = dRFTimeFactory->Step_TimeToNearInputTime(thisRFTime->dTime, tagmHitVector[i]->time_fadc);
            Fill1DHistogram("HLDetectorTiming", "TAGM_ADC_RF_Compare", name,
                    tagmHitVector[i]->time_fadc - locShiftedTime,
                    title,
                    NBINS_RF_COMPARE, MIN_RF_COMPARE, MAX_RF_COMPARE);
        }

        if(thisTAGMHit->has_TDC && thisRFTime != NULL){
            char name [200];
            char title[500];
            sprintf(name, "Column %.3i Row %.1i", tagmHitVector[i]->column, tagmHitVector[i]->row);
            sprintf(title, "TAGM Column %i t_{TDC} - t_{RF}; t_{TDC} - t_{RF} [ns]; Entries", tagmHitVector[i]->column);
            double locShiftedTime = dRFTimeFactory->Step_TimeToNearInputTime(thisRFTime->dTime, tagmHitVector[i]->time_tdc);
            Fill1DHistogram("HLDetectorTiming", "TAGM_TDC_RF_Compare", name,
                    tagmHitVector[i]->time_tdc - locShiftedTime,
                    title,
                    NBINS_RF_COMPARE, MIN_RF_COMPARE, MAX_RF_COMPARE);
        }

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

            int nTAGMCounters = 122; // 102 + 20 including 4 fully read out columns
            Fill2DHistogram("HLDetectorTiming", "TAGM", "TAGMHit TDC_ADC Difference",
                    GetCCDBIndexTAGM(tagmHitVector[i]), tagmHitVector[i]->time_tdc - tagmHitVector[i]->time_fadc,
                    "TAGM #Deltat TDC-ADC; Column ;t_{TDC} - t_{ADC} [ns]", nTAGMCounters, 0.5, nTAGMCounters + 0.5, NBINS_TDIFF, MIN_TDIFF, MAX_TDIFF);
            if (DO_OPTIONAL){
                Fill2DHistogram("HLDetectorTiming", "TAGM", "TAGM Per Channel TDC Time",
                        GetCCDBIndexTAGM(tagmHitVector[i]), tagmHitVector[i]->time_tdc,
                        "TAGM Per Channel TDC time; Column ;t_{TDC} [ns]", nTAGMCounters, 0.5, nTAGMCounters + 0.5, 100, -50, 50);
            }

            if (DO_TRACK_BASED){
                // TString::Form() is not thread safe! (I learned this the hard way)
                // Do it the old fashioned way
                char name [200];
                char title[500];
                sprintf(name, "Column %.3i", tagmHitVector[i]->column);
                sprintf(title, "TAGM Column %i #Deltat Vs Integral; Integral [ADC Units]; t_{TDC} - t_{ADC} [ns]", tagmHitVector[i]->column);
                Fill2DHistogram("HLDetectorTiming", "TAGM_TimeWalk", name,
                        tagmHitVector[i]->integral, tagmHitVector[i]->time_tdc - tagmHitVector[i]->time_fadc,
                        title,
                        100, 0, 2000, 100, -10, 10, false);
            }
        }

    }
    for (i = 0; i < taghHitVector.size(); i++){
        const DTAGHHit *thisTAGHHit = taghHitVector[i];
        if(thisTAGHHit->has_fADC && thisRFTime != NULL){
            char name [200];
            char title[500];
            sprintf(name, "Counter ID %.3i", taghHitVector[i]->counter_id);
            sprintf(title, "TAGH Counter ID %i t_{ADC} - t_{RF}; t_{ADC} - t_{RF} [ns]; Entries", taghHitVector[i]->counter_id);
            double locShiftedTime = dRFTimeFactory->Step_TimeToNearInputTime(thisRFTime->dTime, taghHitVector[i]->time_fadc);
            Fill1DHistogram("HLDetectorTiming", "TAGH_ADC_RF_Compare", name,
                    taghHitVector[i]->time_fadc - locShiftedTime,
                    title,
                    NBINS_RF_COMPARE, MIN_RF_COMPARE, MAX_RF_COMPARE);
        }
        if(thisTAGHHit->has_TDC && thisRFTime != NULL){
            char name [200];
            char title[500];
            sprintf(name, "Counter ID %.3i", taghHitVector[i]->counter_id);
            sprintf(title, "TAGH Counter ID %i t_{TDC} - t_{RF}; t_{TDC} - t_{RF} [ns]; Entries", taghHitVector[i]->counter_id);
            double locShiftedTime = dRFTimeFactory->Step_TimeToNearInputTime(thisRFTime->dTime, taghHitVector[i]->time_tdc);
            Fill1DHistogram("HLDetectorTiming", "TAGH_TDC_RF_Compare", name,
                    taghHitVector[i]->time_tdc - locShiftedTime,
                    title,
                    NBINS_RF_COMPARE, MIN_RF_COMPARE, MAX_RF_COMPARE);
        }

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
                    "TAGH #Deltat TDC-ADC; Counter ID ;t_{TDC} - t_{ADC} [ns]", nTAGHCounters, 0.5, nTAGHCounters + 0.5, NBINS_TDIFF, MIN_TDIFF, MAX_TDIFF);
            if (DO_OPTIONAL){
                Fill2DHistogram("HLDetectorTiming", "TAGH", "TAGM Per Channel TDC Time",
                        taghHitVector[i]->counter_id, taghHitVector[i]->time_tdc,
                        "TAGH Per Channel TDC Time; Counter ID ;t_{TDC} [ns]", nTAGHCounters, 0.5, nTAGHCounters + 0.5, 100, -50, 50);
            }

            if (DO_TRACK_BASED){
                char name [200];
                char title[500];
                sprintf(name, "Counter ID %.3i", taghHitVector[i]->counter_id);
                sprintf(title, "TAGH Counter ID %i #Deltat Vs Integral; Integral [ADC Units]; t_{TDC} - t_{ADC} [ns]", taghHitVector[i]->counter_id);
                Fill2DHistogram("HLDetectorTiming", "TAGH_TimeWalk", name,
                        taghHitVector[i]->integral, taghHitVector[i]->time_tdc - taghHitVector[i]->time_fadc,
                        title,
                        100, 1000, 7000, 100, -10, 10, false);
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
                    "TOF #Deltat TDC-ADC; CDCB Index ;t_{TDC} - t_{ADC} [ns]", nTOFCounters, 0.5, nTOFCounters + 0.5, NBINS_TDIFF, MIN_TDIFF, MAX_TDIFF);
            if (DO_OPTIONAL){
                Fill2DHistogram("HLDetectorTiming", "TOF", "TOFHit Per Channel TDC Time",
                        GetCCDBIndexTOF(tofHitVector[i]), tofHitVector[i]->t_TDC,
                        "TOF Per Channel TDC Time; CDCB Index ;t_{TDC} [ns]", nTOFCounters, 0.5, nTOFCounters + 0.5, 100, -50, 50);
            }

            if (DO_TRACK_BASED){
                // TString->Form(...) is not thread safe, trying sprintf
                char name [200];
                char title[500];
                sprintf(name, "Counter ID %.3i", GetCCDBIndexTOF(tofHitVector[i]));
                sprintf(title, "TOF CCDB Counter ID %i #Deltat Vs dE; dE [Arb.]; t_{TDC} - t_{ADC} [ns]", GetCCDBIndexTOF(tofHitVector[i]));
                Fill2DHistogram("HLDetectorTiming", "TOF_TimeWalk", name,
                        tofHitVector[i]->dE, tofHitVector[i]->t_TDC - tofHitVector[i]->t_fADC,
                        title,
                        150, 0, 0.01, 100, -10, 10, false);
            }
        }
    }

    // Next the relative times between detectors using tracking
    // By the time we get to this point, our first guess at the timing should be fairly good. 
    // Certainly good enough to take a pass at the time based tracking
    // This will be the final alignment step for now

    if (!DO_TRACK_BASED && !DO_VERIFY && !DO_CDC_TIMING) return NOERROR; // Before this stage we aren't really ready yet, so just return

    // We want to plot the delta t at the target between the SC hit and the tagger hits
    // Some limits for these
    float nBinsE = 160, EMin = 2.0, EMax = 10;

    const DEventRFBunch *thisRFBunch = NULL;
    loop->GetSingle(thisRFBunch, "Calibrations"); // SC only hits

    if (USE_RF_BUNCH){
        if (thisRFBunch == NULL || thisRFBunch->dNumParticleVotes < 2) return NOERROR; // Require 2 votes on this RF bunch
        // Loop over TAGM hits
        for (unsigned int j = 0 ; j < tagmHitVector.size(); j++){
            int nTAGMColumns = 122;
            // We want to look at the timewalk within these ADC/TDC detectors
            Fill2DHistogram("HLDetectorTiming", "TRACKING", "TAGM - RFBunch Time",
                    GetCCDBIndexTAGM(tagmHitVector[j]), tagmHitVector[j]->t - thisRFBunch->dTime,
                    "#Deltat TAGM-RFBunch; Column ;t_{TAGM} - t_{SC @ target} [ns]",
                    nTAGMColumns, 0.5, nTAGMColumns + 0.5, NBINS_TAGGER_TIME,MIN_TAGGER_TIME,MAX_TAGGER_TIME);
            Fill2DHistogram("HLDetectorTiming", "TRACKING", "Tagger - RFBunch Time",
                    tagmHitVector[j]->t - thisRFBunch->dTime, tagmHitVector[j]->E,
                    "Tagger - RFBunch Time; #Deltat_{Tagger - SC} [ns]; Energy [GeV]",
                    NBINS_TAGGER_TIME,MIN_TAGGER_TIME,MAX_TAGGER_TIME, nBinsE, EMin, EMax);
            Fill1DHistogram("HLDetectorTiming", "TRACKING", "Tagger - RFBunch 1D Time",
                    tagmHitVector[j]->t - thisRFBunch->dTime,
                    "Tagger - RFBunch Time; #Deltat_{Tagger - RFBunch} [ns]; Entries",
                    160, -20, 20);
        }
        // Loop over TAGH hits
        for (unsigned int j = 0 ; j < taghHitVector.size(); j++){
            int nTAGHCounters = 274;
            Fill2DHistogram("HLDetectorTiming", "TRACKING", "TAGH -  RFBunch Time",
                    taghHitVector[j]->counter_id, taghHitVector[j]->t - thisRFBunch->dTime,
                    "#Deltat TAGH-RFBunch; Counter ID ;t_{TAGH} - t_{RFBunch} [ns]",
                    nTAGHCounters, 0.5, nTAGHCounters + 0.5, NBINS_TAGGER_TIME,MIN_TAGGER_TIME,MAX_TAGGER_TIME);

            Fill2DHistogram("HLDetectorTiming", "TRACKING", "Tagger - RFBunch Time",
                    taghHitVector[j]->t - thisRFBunch->dTime, taghHitVector[j]->E,
                    "Tagger - RFBunch Time; #Deltat_{Tagger - RFBunch} [ns]; Energy [GeV]",
                    NBINS_TAGGER_TIME,MIN_TAGGER_TIME,MAX_TAGGER_TIME, nBinsE, EMin, EMax);

            Fill1DHistogram("HLDetectorTiming", "TRACKING", "Tagger - RFBunch 1D Time",
                    taghHitVector[j]->t - thisRFBunch->dTime,
                    "Tagger - RFBunch Time; #Deltat_{Tagger - RFBunch} [ns]; Entries",
                    160, -20, 20);
        }
    }

    if (DO_REACTION){
        // Trigger the analysis
        vector<const DAnalysisResults*> locAnalysisResultsVector;
        loop->Get(locAnalysisResultsVector);
        // Get the time from the results and fill histograms
        deque<const DParticleCombo*> locPassedParticleCombos;
        locAnalysisResultsVector[0]->Get_PassedParticleCombos(locPassedParticleCombos);

        for (i=0; i < locPassedParticleCombos.size(); i++){
            double taggerTime = locPassedParticleCombos[i]->Get_ParticleComboStep(0)->Get_InitialParticle_Measured()->time();
            //cout << "We have a tagger time of " << taggerTime << endl;
            //Find matching hit by time
            for (unsigned int j = 0 ; j < tagmHitVector.size(); j++){
                if (taggerTime == tagmHitVector[j]->t){
                    int nTAGMColumns = 122;
                    // We want to look at the timewalk within these ADC/TDC detectors
                    Fill2DHistogram("HLDetectorTiming", "TRACKING", "TAGM - RFBunch Time p2pi",
                            GetCCDBIndexTAGM(tagmHitVector[j]), tagmHitVector[j]->t - thisRFBunch->dTime,
                            "#Deltat TAGM-RFBunch; Column ;t_{TAGM} - t_{SC @ target} [ns]",
                            nTAGMColumns, 0.5, nTAGMColumns + 0.5, NBINS_TAGGER_TIME,MIN_TAGGER_TIME,MAX_TAGGER_TIME);
                    Fill2DHistogram("HLDetectorTiming", "TRACKING", "Tagger - RFBunch Time p2pi",
                            tagmHitVector[j]->t - thisRFBunch->dTime, tagmHitVector[j]->E,
                            "Tagger - RFBunch Time; #Deltat_{Tagger - SC} [ns]; Energy [GeV]",
                            NBINS_TAGGER_TIME,MIN_TAGGER_TIME,MAX_TAGGER_TIME, nBinsE, EMin, EMax);
                    Fill1DHistogram("HLDetectorTiming", "TRACKING", "Tagger - RFBunch 1D Time p2pi",
                            tagmHitVector[j]->t - thisRFBunch->dTime,
                            "Tagger - RFBunch Time; #Deltat_{Tagger - RFBunch} [ns]; Entries",
                            160, -20, 20);
                }
            }
            // Loop over TAGH hits
            for (unsigned int j = 0 ; j < taghHitVector.size(); j++){
                if (taggerTime == taghHitVector[j]->t){
                    int nTAGHCounters = 274;
                    Fill2DHistogram("HLDetectorTiming", "TRACKING", "TAGH -  RFBunch Time p2pi",
                            taghHitVector[j]->counter_id, taghHitVector[j]->t - thisRFBunch->dTime,
                            "#Deltat TAGH-RFBunch; Counter ID ;t_{TAGH} - t_{RFBunch} [ns]",
                            nTAGHCounters, 0.5, nTAGHCounters + 0.5, NBINS_TAGGER_TIME,MIN_TAGGER_TIME,MAX_TAGGER_TIME);

                    Fill2DHistogram("HLDetectorTiming", "TRACKING", "Tagger - RFBunch Time p2pi",
                            taghHitVector[j]->t - thisRFBunch->dTime, taghHitVector[j]->E,
                            "Tagger - RFBunch Time; #Deltat_{Tagger - RFBunch} [ns]; Energy [GeV]",
                            NBINS_TAGGER_TIME,MIN_TAGGER_TIME,MAX_TAGGER_TIME, nBinsE, EMin, EMax);

                    Fill1DHistogram("HLDetectorTiming", "TRACKING", "Tagger - RFBunch 1D Time p2pi",
                            taghHitVector[j]->t - thisRFBunch->dTime,
                            "Tagger - RFBunch Time; #Deltat_{Tagger - RFBunch} [ns]; Entries",
                            160, -20, 20);
                }
            }
        }
    }

    // Try using the detector matches
    // Loop over the charged tracks
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
        if (DO_OPTIONAL || DO_CDC_TIMING){
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
            int nTAGMColumns = 122;
            // We want to look at the timewalk within these ADC/TDC detectors
            Fill2DHistogram("HLDetectorTiming", "TRACKING", "TAGM - SC Target Time",
                    GetCCDBIndexTAGM(tagmHitVector[j]), tagmHitVector[j]->t - flightTimeCorrectedSCTime,
                    "#Deltat TAGM-SC; Column ;t_{TAGM} - t_{SC @ target} [ns]", nTAGMColumns, 0.5, nTAGMColumns + 0.5, NBINS_TAGGER_TIME,MIN_TAGGER_TIME,MAX_TAGGER_TIME);
            Fill2DHistogram("HLDetectorTiming", "TRACKING", "Tagger - SC Target Time",
                    tagmHitVector[j]->t - flightTimeCorrectedSCTime, tagmHitVector[j]->E,
                    "Tagger - SC Target Time; #Deltat_{Tagger - SC} [ns]; Energy [GeV]",
                    NBINS_TAGGER_TIME,MIN_TAGGER_TIME,MAX_TAGGER_TIME, nBinsE, EMin, EMax);   
            Fill1DHistogram("HLDetectorTiming", "TRACKING", "Tagger - SC 1D Target Time",
                    tagmHitVector[j]->t - flightTimeCorrectedSCTime,
                    "Tagger - SC Time at Target; #Deltat_{Tagger - SC} [ns]; Entries",
                    160, -20, 20);
        }
        // Loop over TAGH hits
        for (unsigned int j = 0 ; j < taghHitVector.size(); j++){
            int nTAGHCounters = 274;
            Fill2DHistogram("HLDetectorTiming", "TRACKING", "TAGH - SC Target Time",
                    taghHitVector[j]->counter_id, taghHitVector[j]->t - flightTimeCorrectedSCTime,
                    "#Deltat TAGH-SC; Counter ID ;t_{TAGH} - t_{SC @ target} [ns]", nTAGHCounters, 0.5, nTAGHCounters + 0.5, NBINS_TAGGER_TIME,MIN_TAGGER_TIME,MAX_TAGGER_TIME);

            Fill2DHistogram("HLDetectorTiming", "TRACKING", "Tagger - SC Target Time",
                    taghHitVector[j]->t - flightTimeCorrectedSCTime, taghHitVector[j]->E,
                    "Tagger - SC Target Time; #Deltat_{Tagger - SC} [ns]; Energy [GeV]",
                    NBINS_TAGGER_TIME,MIN_TAGGER_TIME,MAX_TAGGER_TIME, nBinsE, EMin, EMax);

            Fill1DHistogram("HLDetectorTiming", "TRACKING", "Tagger - SC 1D Target Time",
                    taghHitVector[j]->t - flightTimeCorrectedSCTime,
                    "Tagger - SC Time at Target; #Deltat_{Tagger - SC} [ns]; Entries",
                    160, -20, 20);
        }

        if(thisTrack->charge() > 0) continue; // Cut on the charge to kill protons

        if (foundTOF){
            // Now check the TOF matching. Do this on a full detector level.
            float flightTimeCorrectedTOFTime = locTOFHitMatchParams.dHitTime - locTOFHitMatchParams.dFlightTime;
            Fill1DHistogram("HLDetectorTiming", "TRACKING", "TOF - SC Target Time",
                    flightTimeCorrectedTOFTime - flightTimeCorrectedSCTime,
                    "t_{TOF} - t_{SC} at Target; t_{TOF} - t_{SC} at Target [ns]; Entries",
                    NBINS_MATCHING, MIN_MATCHING_T, MAX_MATCHING_T);
        }
        if (foundBCAL){
            float flightTimeCorrectedBCALTime = locBCALShowerMatchParams.dBCALShower->t - locBCALShowerMatchParams.dFlightTime;
            Fill1DHistogram("HLDetectorTiming", "TRACKING", "BCAL - SC Target Time",
                    flightTimeCorrectedBCALTime - flightTimeCorrectedSCTime,
                    "t_{BCAL} - t_{SC} at Target; t_{BCAL} - t_{SC} [ns]; Entries",
                    NBINS_MATCHING, MIN_MATCHING_T, MAX_MATCHING_T);
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
                    NBINS_MATCHING, MIN_MATCHING_T, MAX_MATCHING_T);
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

    return NOERROR;
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

int JEventProcessor_HLDetectorTiming::GetCCDBIndexTAGM(const DTAGMHit *thisHit){
    // Since there are a few counters where each row is read out seperately this is a bit of a mess
    int row = thisHit->row;
    int column = thisHit->column;

    int CCDBIndex = column + row;
    if (column > 9) CCDBIndex += 5;
    if (column > 27) CCDBIndex += 5;
    if (column > 81) CCDBIndex += 5;
    if (column > 99) CCDBIndex += 5;

    return CCDBIndex;
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

