// $Id$
//
//    File: JEventProcessor_BCAL_TDC_Timing.cc
// Created: Tue Jul 28 10:55:56 EDT 2015
// Creator: mstaib (on Linux egbert 2.6.32-504.30.3.el6.x86_64 x86_64)
//

/*************************************
  This plugin is designed to calibrate the TDC times for the BCAL.
  These calibrations include the time-walk for each TDC channel, effective velocities and relative offsets between ends of the detector.
  Creation of constants from histograms is controlled by ROOT scripts in the FitScripts directory.
 *************************************/

#include "JEventProcessor_BCAL_TDC_Timing.h"
#include "HistogramTools.h" // To make my life easier
#include "BCAL/DBCALHit.h"
#include "BCAL/DBCALTDCHit.h"
#include "BCAL/DBCALCluster.h"
#include "BCAL/DBCALPoint.h"
#include "BCAL/DBCALUnifiedHit.h"
#include "BCAL/DBCALGeometry.h"
#include "PID/DChargedTrack.h"
#include "TRACKING/DTrackTimeBased.h"
#include "PID/DEventRFBunch.h"
#include "PID/DDetectorMatches.h"
#include "DAQ/Df250PulsePedestal.h" // Needed for pulse peak information
#include "BCAL/DBCALDigiHit.h"

using namespace jana;

#include "DANA/DApplication.h"
// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_BCAL_TDC_Timing());
}
} // "C"


//------------------
// JEventProcessor_BCAL_TDC_Timing (Constructor)
//------------------
JEventProcessor_BCAL_TDC_Timing::JEventProcessor_BCAL_TDC_Timing()
{

}

//------------------
// ~JEventProcessor_BCAL_TDC_Timing (Destructor)
//------------------
JEventProcessor_BCAL_TDC_Timing::~JEventProcessor_BCAL_TDC_Timing()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_BCAL_TDC_Timing::init(void)
{
    // This is called once at program startup on a single thread.

    return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_BCAL_TDC_Timing::brun(JEventLoop *loop, int runnumber)
{
    // This is called whenever the run number changes
    DApplication* app = dynamic_cast<DApplication*>(loop->GetJApplication());
    DGeometry* geom = app->GetDGeometry(runnumber);
    geom->GetTargetZ(Z_TARGET);
    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_BCAL_TDC_Timing::evnt(JEventLoop *loop, int eventnumber)
{
    vector<const DBCALUnifiedHit *> bcalUnifiedHitVector;
    loop->Get(bcalUnifiedHitVector);

    /**********************************************
     _____ _                   __        __    _ _    
    |_   _(_)_ __ ___   ___    \ \      / /_ _| | | __
      | | | | '_ ` _ \ / _ \____\ \ /\ / / _` | | |/ /
      | | | | | | | | |  __/_____\ V  V / (_| | |   < 
      |_| |_|_| |_| |_|\___|      \_/\_/ \__,_|_|_|\_\
     ********************************************/
    // The very first thing to do is correct for the timewalk. If this calibration is screwed up, the calibrations that follow will be wrong...
    // The following plots can be used for the calibration or to check the calibration if it has already been done

    int NBINS_TDIFF = 100;
    double MIN_TDIFF = -10.0, MAX_TDIFF = 10.0;

    for (unsigned int i = 0; i < bcalUnifiedHitVector.size(); i++){
        //int the_cell = (bcalUnifiedHitVector[i]->module - 1) * 16 + (bcalUnifiedHitVector[i]->layer - 1) * 4 + bcalUnifiedHitVector[i]->sector;
        // There is one less layer of TDCs so the numbering relects this
        int the_tdc_cell = (bcalUnifiedHitVector[i]->module - 1) * 12 + (bcalUnifiedHitVector[i]->layer - 1) * 4 + bcalUnifiedHitVector[i]->sector;
        // Get the underlying associated objects
        const DBCALHit * thisADCHit;
        const DBCALTDCHit * thisTDCHit;
        bcalUnifiedHitVector[i]->GetSingle(thisADCHit);
        bcalUnifiedHitVector[i]->GetSingle(thisTDCHit);

        // From the ADC hit we can extract the pulse peak information, just have to go down the chain a little...
        const DBCALDigiHit *thisDigiHit;
        thisADCHit->GetSingle(thisDigiHit);
        const Df250PulsePedestal *pp;
        thisDigiHit->GetSingle(pp);
        int pulse_peak = 0;
        if (pp != NULL){
            pulse_peak = pp->pulse_peak;
            pulse_peak -= pp->pedestal;
        }

        // The raw information from the DBCALHit and DBCALTDCHit is not corrected for timewalk yet, so we can always plot the before and after.
        if (thisADCHit != NULL && thisTDCHit != NULL){
            char name[200];
            sprintf(name, "Module %.2i Layer %.2i Sector %.2i", bcalUnifiedHitVector[i]->module, bcalUnifiedHitVector[i]->layer, bcalUnifiedHitVector[i]->sector);
            if (bcalUnifiedHitVector[i]->end == 0){
                Fill2DHistogram ("BCAL_TDC_Timing", "BCAL_Upstream_Timewalk_NoCorrection_E", name,
                        bcalUnifiedHitVector[i]->E, thisTDCHit->t - thisADCHit->t,
                        "Timewalk; E [GeV]; t_{TDC} - t_{ADC} [ns]",
                        250, 0.0, 0.35, NBINS_TDIFF, MIN_TDIFF , MAX_TDIFF );
                Fill2DHistogram ("BCAL_TDC_Timing", "BCAL_Upstream_Timewalk_NoCorrection_PP", name,
                        pulse_peak, thisTDCHit->t - thisADCHit->t,
                        "Timewalk; Pulse Peak [ADC Counts]; t_{TDC} - t_{ADC} [ns]",
                        500, 0.0, 1500.0, NBINS_TDIFF, MIN_TDIFF , MAX_TDIFF );
                Fill2DHistogram ("BCAL_TDC_Timing", "BCAL", "BCALHit Upstream Per Channel TDC-ADC Hit Time - No Timewalk",
                        the_tdc_cell, thisTDCHit->t - thisADCHit->t,
                        "BCALHit Upstream Per Channel TDC-ADC Hit Time; cellID; t_{TDC} - t_{ADC} [ns] ",
                        576, 0.5, 576.5, NBINS_TDIFF, MIN_TDIFF, MAX_TDIFF);
            }

            else{
                Fill2DHistogram ("BCAL_TDC_Timing", "BCAL_Downstream_Timewalk_NoCorrection_E", name,
                        bcalUnifiedHitVector[i]->E, thisTDCHit->t - thisADCHit->t,
                        "Timewalk; E [GeV]; t_{TDC} - t_{ADC} [ns]",
                        250, 0.0, 0.35, NBINS_TDIFF, MIN_TDIFF, MAX_TDIFF);
                Fill2DHistogram ("BCAL_TDC_Timing", "BCAL_Downstream_Timewalk_NoCorrection_PP", name,
                        pulse_peak, thisTDCHit->t - thisADCHit->t,
                        "Timewalk; Pulse Peak [ADC Counts]; t_{TDC} - t_{ADC} [ns]",
                        500, 0.0, 1500.0, NBINS_TDIFF, MIN_TDIFF, MAX_TDIFF);
                Fill2DHistogram ("BCAL_TDC_Timing", "BCAL", "BCALHit Downstream Per Channel TDC-ADC Hit Time - No Timewalk",
                        the_tdc_cell, thisTDCHit->t - thisADCHit->t,
                        "BCALHit Downstream Per Channel TDC-ADC Hit Time; cellID; t_{TDC} - t_{ADC} [ns] ",
                        576, 0.5, 576.5, NBINS_TDIFF, MIN_TDIFF, MAX_TDIFF);
            }
        }
        // Next look directly at the DBCALUnifiedHit to get the corrected times and plot those seperately
        if (bcalUnifiedHitVector[i]->has_TDC_hit){
            char name[200];
            sprintf(name, "Module %.2i Layer %.2i Sector %.2i", bcalUnifiedHitVector[i]->module, bcalUnifiedHitVector[i]->layer, bcalUnifiedHitVector[i]->sector);
            if (bcalUnifiedHitVector[i]->end == 0){
                Fill2DHistogram ("BCAL_TDC_Timing", "BCAL_Upstream_Timewalk_WithCorrection_E", name,
                        bcalUnifiedHitVector[i]->E, bcalUnifiedHitVector[i]->t_TDC - bcalUnifiedHitVector[i]->t_ADC,
                        "Timewalk; E [GeV]; t_{TDC} - t_{ADC} [ns]",
                        250, 0.0, 0.35, NBINS_TDIFF, MIN_TDIFF , MAX_TDIFF );
                Fill2DHistogram ("BCAL_TDC_Timing", "BCAL_Upstream_Timewalk_WithCorrection_PP", name,
                        pulse_peak, bcalUnifiedHitVector[i]->t_TDC - bcalUnifiedHitVector[i]->t_ADC,
                        "Timewalk; Pulse Peak [ADC Counts]; t_{TDC} - t_{ADC} [ns]",
                        500, 0.0, 1500.0, NBINS_TDIFF, MIN_TDIFF , MAX_TDIFF );
                Fill2DHistogram ("BCAL_TDC_Timing", "BCAL", "BCALHit Upstream Per Channel TDC-ADC Hit Time - With Timewalk",
                        the_tdc_cell, bcalUnifiedHitVector[i]->t_TDC - bcalUnifiedHitVector[i]->t_ADC,
                        "BCALHit Upstream Per Channel TDC-ADC Hit Time - TW Corrected; cellID; t_{TDC} - t_{ADC} [ns] ",
                        576, 0.5, 576.5, NBINS_TDIFF, MIN_TDIFF, MAX_TDIFF);
            }
            else{
                Fill2DHistogram ("BCAL_TDC_Timing", "BCAL_Downstream_Timewalk_WithCorrection_E", name,
                        bcalUnifiedHitVector[i]->E, bcalUnifiedHitVector[i]->t_TDC - bcalUnifiedHitVector[i]->t_ADC,
                        "Timewalk; E [GeV]; t_{TDC} - t_{ADC} [ns]",
                        250, 0.0, 0.35, NBINS_TDIFF, MIN_TDIFF, MAX_TDIFF);
                Fill2DHistogram ("BCAL_TDC_Timing", "BCAL_Downstream_Timewalk_WithCorrection_PP", name,
                        pulse_peak, bcalUnifiedHitVector[i]->t_TDC - bcalUnifiedHitVector[i]->t_ADC,
                        "Timewalk; Pulse Peak [ADC Counts]; t_{TDC} - t_{ADC} [ns]",
                        500, 0.0, 1500.0, NBINS_TDIFF, MIN_TDIFF, MAX_TDIFF);
                Fill2DHistogram ("BCAL_TDC_Timing", "BCAL", "BCALHit Downstream Per Channel TDC-ADC Hit Time - With Timewalk",
                        the_tdc_cell, bcalUnifiedHitVector[i]->t_TDC - bcalUnifiedHitVector[i]->t_ADC,
                        "BCALHit Downstream Per Channel TDC-ADC Hit Time - TW Corrected; cellID; t_{TDC} - t_{ADC} [ns] ",
                        576, 0.5, 576.5, NBINS_TDIFF, MIN_TDIFF, MAX_TDIFF);
            }
        }
    }
    /*************************************************
        _________  _____  _______       _
       /_  __/ _ \/ ___/ /_  __(_)_ _  (_)__  ___ _
        / / / // / /__    / / / /  ' \/ / _ \/ _ `/
       /_/ /____/\___/   /_/ /_/_/_/_/_/_//_/\_, /
                                            /___/
    **************************************************/

    // Now we will use the track matching in order to work out the time offsets and effective velocities
    // We just need to grab the charged particles in order to get their detector matches
    // Note that this method is only really applicable for runs with the solenoid on...

    // We need the RF bunch for the event in order to check the global alignemtn of the timing
    const DEventRFBunch *thisRFBunch = NULL;
    loop->GetSingle(thisRFBunch);

    vector <const DChargedTrack *> chargedTrackVector;
    loop->Get(chargedTrackVector);

    for (unsigned int iTrack = 0; iTrack < chargedTrackVector.size(); iTrack++){
        // Pick out the best charged track hypothesis for this charged track based only on the Tracking FOM
        const DChargedTrackHypothesis* bestHypothesis = chargedTrackVector[iTrack]->Get_BestTrackingFOM();

        // Now from this hypothesis we can get the detector matches to the BCAL
        const DBCALShowerMatchParams* bcalMatch = bestHypothesis->Get_BCALShowerMatchParams();
        const DSCHitMatchParams* scMatch = bestHypothesis->Get_SCHitMatchParams(); // Needed for quality cut later
        if (bcalMatch == NULL) continue; 

        // We also need the reference trajectory, which is buried deep in there
        const DTrackTimeBased *timeBasedTrack = (const DTrackTimeBased *) bcalMatch->dTrack;
        const DReferenceTrajectory *rt = timeBasedTrack->rt;

        // Get the shower from the match
        const DBCALShower *thisShower = bcalMatch->dBCALShower;

        // Get the clusters from the showers
        vector <const DBCALCluster *> clusterVector;
        thisShower->Get(clusterVector);

        // Loop over clusters within the shower
        for (unsigned int iCluster = 0; iCluster < clusterVector.size(); iCluster++){
            // Get the points
            vector <const DBCALPoint*> pointVector = clusterVector[iCluster]->points();
            // Loop over the points within the cluster
            for (unsigned int iPoint = 0; iPoint < pointVector.size(); iPoint++){
                const DBCALPoint *thisPoint = pointVector[iPoint];
                if (thisPoint->E() < 0.05) continue; // The timing is known not to be great for very low energy, so only use our best info 
                DVector3 proj_pos = rt->GetLastDOCAPoint();
                double pathLength, flightTime;
                double rpoint = thisPoint->r();
                if (rt->GetIntersectionWithRadius(rpoint,proj_pos, &pathLength, &flightTime)==NOERROR){
                    // Now proj_pos contains the projected position of the track at this particular point within the BCAL
                    // We can plot the difference of the projected position and the BCAL position as a function of the channel
                    char name[200];
                    sprintf(name , "Module %.2i Layer %.2i Sector %.2i", thisPoint->module(), thisPoint->layer(), thisPoint->sector());
                    // These results are in slightly different coordinate systems. We want one where the center of the BCAL is z=0
                    double localTrackHitZ = proj_pos.z() - DBCALGeometry::GLOBAL_CENTER;
                    double localBCALHitZ = thisPoint->z() - DBCALGeometry::GLOBAL_CENTER + Z_TARGET;
                    Fill2DHistogram ("BCAL_TDC_Offsets", "Z Position", name,
                            localTrackHitZ, localBCALHitZ,
                            "Z_{point} Vs. Z_{Track}; Z_{Track} [cm]; Z_{Point} [cm]",
                            500, -250, 250, 500, -250, 250); 

                    // Now fill some histograms that are useful for aligning the BCAL with the rest of the detector systems
                    if (thisRFBunch->dNumParticleVotes >= 2 && scMatch != NULL){ // Require good RF bunch and this track match the SC
                        // Get the time of the BCAL point
                        double pointTime = thisPoint->t();
                        // We have the flight time to our BCAL point, so we can get the target time
                        double targetCenterTime = pointTime - flightTime - ((timeBasedTrack->position()).Z() - Z_TARGET) / SPEED_OF_LIGHT;

                        // Now we just plot the difference in from the RF Time to get out the correction
                        int the_cell = (thisPoint->module() - 1) * 16 + (thisPoint->layer() - 1) * 4 + thisPoint->sector();
                        Fill2DHistogram("BCAL_Global_Offsets", "Target Time", "Target Time Minus RF Time Vs. Cell Number",
                                the_cell, targetCenterTime - thisRFBunch->dTime,
                                "Target Time Minus RF Time Vs. Cell Number; CCDB Index; t_{Target} - t_{RF} [ns]",
                                768, 0.5, 768.5, 200, -10, 10);
                    }
                }
            }
        }
    }

    return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_BCAL_TDC_Timing::erun(void)
{
    // This is called whenever the run number changes, before it is
    // changed to give you a chance to clean up before processing
    // events from the next run number.
    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_BCAL_TDC_Timing::fini(void)
{
    // Called before program exit after event processing is finished.
    SortDirectories(); // Sort the histogram directories by name
    return NOERROR;
}

