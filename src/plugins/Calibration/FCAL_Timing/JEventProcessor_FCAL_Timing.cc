// $Id$
//
//    File: JEventProcessor_FCAL_Timing.cc
// Created: Mon Jan  4 13:03:20 EST 2016
// Creator: mstaib (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "JEventProcessor_FCAL_Timing.h"
using namespace jana;


// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
#include "HistogramTools.h"
#include "FCAL/DFCALShower.h"
#include "FCAL/DFCALCluster.h"
#include "PID/DDetectorMatches.h"
#include "PID/DChargedTrack.h"
#include "TRACKING/DTrackTimeBased.h"
#include "PID/DEventRFBunch.h"
#include "DANA/DApplication.h"

extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_FCAL_Timing());
}
} // "C"


//------------------
// JEventProcessor_FCAL_Timing (Constructor)
//------------------
JEventProcessor_FCAL_Timing::JEventProcessor_FCAL_Timing()
{

}

//------------------
// ~JEventProcessor_FCAL_Timing (Destructor)
//------------------
JEventProcessor_FCAL_Timing::~JEventProcessor_FCAL_Timing()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_FCAL_Timing::init(void)
{
	// This is called once at program startup. If you are creating
	// and filling historgrams in this plugin, you should lock the
	// ROOT mutex like this:
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();
	//

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_FCAL_Timing::brun(JEventLoop *eventLoop, int32_t runnumber)
{
    // This is called whenever the run number changes
    DApplication* app = dynamic_cast<DApplication*>(eventLoop->GetJApplication());
    DGeometry* geom = app->GetDGeometry(runnumber);
    geom->GetTargetZ(Z_TARGET);
    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_FCAL_Timing::evnt(JEventLoop *loop, uint64_t eventnumber)
{
    // Now we will use the track matching in order to work out the time offsets and effective velocities
    // We just need to grab the charged particles in order to get their detector matches
    // Note that this method is only really applicable for runs with the solenoid on...

    // We need the RF bunch for the event in order to check the global alignemtn of the timing
    const DEventRFBunch *thisRFBunch = NULL;
    loop->GetSingle(thisRFBunch);
    if (thisRFBunch->dNumParticleVotes < 2 ) return NOERROR; // Require the RF bunch to be fairly well determined

    vector <const DChargedTrack *> chargedTrackVector;
    loop->Get(chargedTrackVector);

    for (unsigned int iTrack = 0; iTrack < chargedTrackVector.size(); iTrack++){
        // Pick out the best charged track hypothesis for this charged track based only on the Tracking FOM
        const DChargedTrackHypothesis* bestHypothesis = chargedTrackVector[iTrack]->Get_BestTrackingFOM();

        // Now from this hypothesis we can get the detector matches to the BCAL
        const DFCALShowerMatchParams* fcalMatch = bestHypothesis->Get_FCALShowerMatchParams();
        const DSCHitMatchParams* scMatch = bestHypothesis->Get_SCHitMatchParams(); 
        if (fcalMatch == NULL || scMatch == NULL) continue; // Skip this track if it doesn't match SC and FCAL

        // We also need the reference trajectory, which is buried deep in there
        const DTrackTimeBased *timeBasedTrack = (const DTrackTimeBased *) fcalMatch->dTrack;

        // Get the shower from the match
        const DFCALShower *thisShower = fcalMatch->dFCALShower;
        double flightTime = fcalMatch->dFlightTime;

        // Get the clusters from the showers
        vector <const DFCALCluster *> clusterVector;
        thisShower->Get(clusterVector);
        
        // Loop over clusters within the shower
        for (unsigned int iCluster = 0; iCluster < clusterVector.size(); iCluster++){
            // Get the hits
            const vector<DFCALCluster::DFCALClusterHit_t> hitVector = clusterVector[iCluster]->GetHits();

            // Get the time correction applied to the shower so that we can apply it to the hits also...
            // Essentially accounts for light propogation to the readout
            double tCorr = thisShower->getTime() - clusterVector[iCluster]->getTime();

            //Loop over hits
            //for (unsigned int iHit = 0; iHit < hitVector.size(); iHit++){
            for (unsigned int iHit = 0; iHit < 1; iHit++){ // Test only using the highest energy hit...
                float hitTime = hitVector[iHit].t;
                hitTime += tCorr; // Apply the t corection used for the cluster/shower conversion (tCorr will be negative)

                int ChannelNumber = hitVector[iHit].ch;

                double targetCenterTime = hitTime - flightTime - ((timeBasedTrack->position()).Z() - Z_TARGET) / SPEED_OF_LIGHT;
                double tDiff = targetCenterTime - thisRFBunch->dTime;

                // Fill the only important histogram
                Fill2DHistogram("FCAL_Timing", "", "Target Time minus RF Time Vs. Channel Number",
                        ChannelNumber, tDiff,
                        "Target Time minus RF Time Vs. Channel Number;Channel Number; Target Time minus RF Time",
                        2800,-0.5,2799.5,250,-20, 20); // Channels are numbered from zero...AGAINST CONVENTION
            }
        }
    }

    return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_FCAL_Timing::erun(void)
{
    // This is called whenever the run number changes, before it is
    // changed to give you a chance to clean up before processing
    // events from the next run number.
    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_FCAL_Timing::fini(void)
{
    // Called before program exit after event processing is finished.
    return NOERROR;
}

