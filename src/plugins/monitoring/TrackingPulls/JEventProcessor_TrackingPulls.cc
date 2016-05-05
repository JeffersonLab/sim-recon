// $Id$
//
//    File: JEventProcessor_TrackingPulls.cc
// Created: Tue Apr 26 09:29:02 EDT 2016
// Creator: mstaib (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "JEventProcessor_TrackingPulls.h"
#include "HistogramTools.h"

using namespace jana;


// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_TrackingPulls());
}
} // "C"


//------------------
// JEventProcessor_TrackingPulls (Constructor)
//------------------
JEventProcessor_TrackingPulls::JEventProcessor_TrackingPulls()
{

}

//------------------
// ~JEventProcessor_TrackingPulls (Destructor)
//------------------
JEventProcessor_TrackingPulls::~JEventProcessor_TrackingPulls()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_TrackingPulls::init(void)
{
    EXCLUDERING=0;
    EXCLUDEPLANE=0;
    if (gPARMS){
        gPARMS->SetDefaultParameter("KALMAN:RING_TO_SKIP", EXCLUDERING);
        gPARMS->SetDefaultParameter("KALMAN:PLANE_TO_SKIP", EXCLUDEPLANE);
    }
    if(EXCLUDERING == 0 ){
        jout << "TrackingPulls::Did not set KALMAN:RING_TO_SKIP on the command line -- Using Biased CDC fits" << endl;
    }
    if(EXCLUDEPLANE == 0 ){
        jout << "TrackingPulls::Did not set KALMAN:PLANE_TO_SKIP on the command line -- Using Biased FDC fits" << endl;
    }
    return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_TrackingPulls::brun(JEventLoop *eventLoop, int32_t runnumber)
{
    // This is called whenever the run number changes
    DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
    JCalibration *jcalib = dapp->GetJCalibration(runnumber);
    // This is called whenever the run number changes
    // Get the straw sag parameters from the database
    unsigned int numstraws[28]={42,42,54,54,66,66,80,80,93,93,106,106,123,123,
        135,135,146,146,158,158,170,170,182,182,197,197,
        209,209};
    max_sag.clear();
    sag_phi_offset.clear();
    vector< map<string, double> > tvals;
    unsigned int straw_count=0,ring_count=0;
    if (jcalib->Get("CDC/sag_parameters", tvals)==false){
        vector<double>temp,temp2;
        for(unsigned int i=0; i<tvals.size(); i++){
            map<string, double> &row = tvals[i];

            temp.push_back(row["offset"]);
            temp2.push_back(row["phi"]);

            straw_count++;
            if (straw_count==numstraws[ring_count]){
                max_sag.push_back(temp);
                sag_phi_offset.push_back(temp2);
                temp.clear();
                temp2.clear();
                straw_count=0;
                ring_count++;
            }
        }
    }
    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_TrackingPulls::evnt(JEventLoop *loop, uint64_t eventnumber)
{
    // For now we eill grab tracks using the higest FOM
    vector <const DChargedTrack *> chargedTrackVector;
    loop->Get(chargedTrackVector);

    for (unsigned int iTrack = 0; iTrack < chargedTrackVector.size(); iTrack++){

        const DChargedTrackHypothesis* bestHypothesis = chargedTrackVector[iTrack]->Get_BestTrackingFOM();

        // Require Single track events
        //if (trackCandidateVector.size() != 1) return NOERROR;
        //const DTrackCandidate* thisTrackCandidate = trackCandidateVector[0];
        // Cut very loosely on the track quality
        const DTrackTimeBased *thisTimeBasedTrack;
        bestHypothesis->GetSingle(thisTimeBasedTrack);
        if (thisTimeBasedTrack->FOM < 1E-20) return NOERROR;
        vector<DTrackFitter::pull_t> pulls = thisTimeBasedTrack->pulls;
        // Loop over the pulls to get the appropriate information for our ring
        for (unsigned int i = 0; i < pulls.size(); i++){
            if (pulls[i].cdc_hit != NULL) FillCDCPulls(pulls[i], thisTimeBasedTrack);
            else if (pulls[i].fdc_hit != NULL) FillFDCPulls(pulls[i]);
        }
    }
    return NOERROR;
}

void JEventProcessor_TrackingPulls::FillCDCPulls( DTrackFitter::pull_t thisPull , const DTrackTimeBased *thisTimeBasedTrack){
    // Funtion to fill histograms for CDC
    const DCDCTrackHit* thisCDCHit = thisPull.cdc_hit;

    int ring = thisCDCHit->wire->ring;
    //int straw = thisCDCHit->wire->straw;
    // Allow for unbiased fits
    if ( EXCLUDERING != 0 && ring != EXCLUDERING) return;

    double residual = thisPull.resi;
    double error = thisPull.err;
    double time = thisPull.tdrift;
    double docaphi = thisPull.docaphi;
    if (docaphi > TMath::Pi()) docaphi -= 2 * TMath::Pi();
    double docaz = thisPull.z;
    double dz = docaz - 92.0;
    //if (docaz < 70.0 || docaz > 110.0) continue; // Only focus on the center of the chamber
    //if (docaz < 140.0) continue; // Only focus on downstream end of chamber
    double distance = thisPull.d; // This is the distance from the lookup table
    double predictedDistance = distance - residual;

    if( EXCLUDERING != 0){
        Fill1DHistogram("TrackingPulls", "CDC", "Residual", residual, "Residual; Residual [cm]; Entries", 100, -0.05, 0.05);
        Fill1DHistogram("TrackingPulls", "CDC", "Pull", residual/error, "Pull; Pull; Entries", 100, -3., 3.);
    }

    char folder[100];
    sprintf(folder, "CDC Ring %.2i", ring);

    Fill1DHistogram("TrackingPulls", folder, "Residual", residual, "Residual; Residual [cm]; Entries", 100, -0.05, 0.05);
    Fill1DHistogram("TrackingPulls", folder, "Pull", residual/error, "Pull; Pull; Entries", 100, -3., 3.);
    Fill2DHistogram("TrackingPulls", folder, "Residual Vs. Momentum", 
            thisTimeBasedTrack->pmag(), residual,
            "Residual Vs. Momentum; Momentum [GeV/c]; Residual [cm]",
            50, 0.0, 12.0, 100, -0.05, 0.05);
    Fill2DHistogram("TrackingPulls", folder, "Residual Vs. Theta",
            thisTimeBasedTrack->momentum().Theta()*TMath::RadToDeg(), residual,
            "Residual Vs. Theta; Theta [deg]; Residual [cm]",
            60, 0.0, 180.0, 100, -0.05, 0.05);
    Fill2DHistogram("TrackingPulls", folder, "Residual Vs. Z",
            dz, residual,
            "Residual Vs. Z; Z (Measured from CDC center) [cm]; Residual [cm]",
            100, -75.0, 75.0, 100, -0.05, 0.05);
    Fill2DHistogram("TrackingPulls", folder, "Residual Vs. Tracking FOM",
            thisTimeBasedTrack->FOM, residual,
            "Residual Vs. Tracking FOM; Tracking FOM; Residual [cm]",
            100, 0.0, 1.0, 100, -0.05, 0.05);
    Fill2DHistogram("TrackingPulls", folder, "Residual Vs. Drift Time",
            time, residual,
            "Residual Vs. Drift Time; Drift Time [ns]; Residual [cm]",
            500, -10, 1500 , 100, -0.05, 0.05);
    Fill2DHistogram("TrackingPulls", folder, "Pull Vs. Drift Time",
            time, residual/error,
            "Pull Vs. Drift Time; Drift Time [ns]; Pull",
            500, -10, 1500 , 100, -3., 3.);
    Fill1DHistogram("TrackingPulls", folder, "Drift Time", time, "Drift Time; Drift Time [ns]; Entries", 500, -10, 1500);
    Fill1DHistogram("TrackingPulls", folder, "Drift Distance", distance, "Drift Distance; Drift Distance [cm]; Entries", 50, 0.0, 1.2);
    Fill1DHistogram("TrackingPulls", folder, "Predicted Drift Distance", predictedDistance, "Predicted Drift Distance; Drift Distance [cm]; Entries", 50, 0.0, 1.2);
    return;
}

void JEventProcessor_TrackingPulls::FillFDCPulls( DTrackFitter::pull_t thisPull ){
    double residual = thisPull.resi;
    double error = thisPull.err;

    Fill1DHistogram("TrackingPulls", "FDC", "Residual", residual, "Residual; Residual [cm]; Entries", 100, -0.05, 0.05);
    Fill1DHistogram("TrackingPulls", "FDC", "Pull", residual/error, "Pull; Pull; Entries", 100, -3., 3.);

    return;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_TrackingPulls::erun(void)
{
    // This is called whenever the run number changes, before it is
    // changed to give you a chance to clean up before processing
    // events from the next run number.
    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_TrackingPulls::fini(void)
{
    // Called before program exit after event processing is finished.
    return NOERROR;
}

