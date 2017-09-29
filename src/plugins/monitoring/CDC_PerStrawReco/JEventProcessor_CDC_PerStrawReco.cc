// $Id$
//
//    File: JEventProcessor_CDC_PerStrawReco.cc
// Created: Mon Jul  6 13:00:51 EDT 2015
// Creator: mstaib (on Linux egbert 2.6.32-504.16.2.el6.x86_64 x86_64)
//

#include "JEventProcessor_CDC_PerStrawReco.h"
#include "PID/DChargedTrack.h"
#include "TRACKING/DTrackTimeBased.h"
#include "HistogramTools.h"

using namespace jana;


// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_CDC_PerStrawReco());
}
} // "C"


//------------------
// JEventProcessor_CDC_PerStrawReco (Constructor)
//------------------
JEventProcessor_CDC_PerStrawReco::JEventProcessor_CDC_PerStrawReco()
{

}

//------------------
// ~JEventProcessor_CDC_PerStrawReco (Destructor)
//------------------
JEventProcessor_CDC_PerStrawReco::~JEventProcessor_CDC_PerStrawReco()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_CDC_PerStrawReco::init(void)
{
    EXCLUDERING=0;
    if (gPARMS){
        gPARMS->SetDefaultParameter("CDCCOSMIC:EXCLUDERING", EXCLUDERING, "Ring Excluded from the fit");
        gPARMS->SetDefaultParameter("KALMAN:RING_TO_SKIP", EXCLUDERING);
    }
    if(EXCLUDERING == 0 ){
        jout << "Did not set CDCCOSMIC:EXCLUDERING on the command line -- Using Biased fits" << endl;
    }
    return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_CDC_PerStrawReco::brun(JEventLoop *eventLoop, int32_t runnumber)
{
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
jerror_t JEventProcessor_CDC_PerStrawReco::evnt(JEventLoop *loop, uint64_t eventnumber)
{
    // Getting the charged tracks will allow us to use the field on data
    vector <const DChargedTrack *> chargedTrackVector;
    loop->Get(chargedTrackVector);

    for (unsigned int iTrack = 0; iTrack < chargedTrackVector.size(); iTrack++){

        const DChargedTrackHypothesis* bestHypothesis = chargedTrackVector[iTrack]->Get_BestTrackingFOM();

        // Require Single track events
        //if (trackCandidateVector.size() != 1) return NOERROR;
        //const DTrackCandidate* thisTrackCandidate = trackCandidateVector[0];
        // Cut very loosely on the track quality
        auto thisTimeBasedTrack = bestHypothesis->Get_TrackTimeBased();
        if (thisTimeBasedTrack->FOM < 1E-20) continue;
        if (!thisTimeBasedTrack->IsSmoothed) continue;
        vector<DTrackFitter::pull_t> pulls = thisTimeBasedTrack->pulls;
        // Loop over the pulls to get the appropriate information for our ring
        for (unsigned int i = 0; i < pulls.size(); i++){
            DTrackFitter::pull_t thisPull = pulls[i];
            double residual = thisPull.resi;
            //double error = thisPull.err;
            double time = thisPull.tcorr;
            double docaphi = thisPull.docaphi;
            if (docaphi > TMath::Pi()) docaphi -= 2 * TMath::Pi();
            double docaz = thisPull.z;
            double dz = docaz - 92.0;
            bool isMiddle = false, isDownstream = false;
            if (docaz > 70.0 && docaz < 110.0) isMiddle = true; 
            if (docaz > 140.0) isDownstream = true;
            if (!isMiddle && !isDownstream) continue;
            double predictedDistance = thisPull.d; // This is the DOCA from the track
            double distance = residual + predictedDistance; // This is the distance from the T-D lookup
            const DCDCTrackHit* thisCDCHit = thisPull.cdc_hit;

            if (thisCDCHit == NULL) continue;

            int ring = thisCDCHit->wire->ring;
            int straw = thisCDCHit->wire->straw;
            // Allow for unbiased fits
            if ( EXCLUDERING != 0 && ring != EXCLUDERING) continue;
            // Now we have just the unbiased information for the ring we have chosen
            // Now just make a bunch of histograms to display all of the information
            char folder[100];
            sprintf(folder, "Ring %.2i", ring);

            const char * regionString = "";
            const char * regionStringPerStraw = "";
            if (isMiddle) {
               regionString = "CDCReco_Middle";
               regionStringPerStraw = "CDCPerStrawReco_Middle";
            }
            else if (isDownstream){
               regionString = "CDCReco_Downstream";
               regionStringPerStraw = "CDCPerStrawReco_Downstream";
            }


            Fill1DHistogram(regionString, folder, "Residuals", residual, "Residuals; Residual [cm]; Entries", 200, -0.05, 0.05);
            Fill2DHistogram(regionString, folder, "Residual Vs. Momentum", 
                  thisTimeBasedTrack->pmag(), residual,
                  "Residual Vs. Momentum; Momentum [GeV/c]; Residual [cm]",
                  50, 0.0, 12.0, 100, -0.05, 0.05);
            Fill2DHistogram(regionString, folder, "Residual Vs. Theta",
                  thisTimeBasedTrack->momentum().Theta()*TMath::RadToDeg(), residual,
                  "Residual Vs. Theta; Theta [deg]; Residual [cm]",
                  60, 0.0, 180.0, 100, -0.05, 0.05);
            Fill2DHistogram(regionString, folder, "Residual Vs. Z",
                  dz, residual,
                  "Residual Vs. Z; Z (Measured from CDC center) [cm]; Residual [cm]",
                  100, -75.0, 75.0, 100, -0.05, 0.05);
            Fill2DHistogram(regionString, folder, "Residual Vs. Tracking FOM",
                  thisTimeBasedTrack->FOM, residual,
                  "Residual Vs. Tracking FOM; Tracking FOM; Residual [cm]",
                  100, 0.0, 1.0, 100, -0.05, 0.05);
            Fill1DHistogram(regionString, folder, "Drift Time", time, "Drift Time; Drift Time [ns]; Entries", 500, -10, 1500);
            Fill1DHistogram(regionString, folder, "Drift Distance", distance, "Drift Distance; Drift Distance [cm]; Entries", 250, 0.0, 1.2);
            Fill1DHistogram(regionString, folder, "Predicted Drift Distance", predictedDistance, "Predicted Drift Distance; Drift Distance [cm]; Entries", 250, 0.0, 1.2);

            char strawname[100];
            char strawtitle[256];
            sprintf(strawname,"Straw %.3i Drift time Vs phi_DOCA", straw);
            sprintf(strawtitle,"Ring %.2i Straw %.3i Drift time Vs phi_DOCA;#phi_{DOCA};Drift Time [ns]", ring, straw);
            Fill2DHistogram(regionStringPerStraw,folder,strawname, docaphi, time,
                  strawtitle, 8, -3.14, 3.14,  500, -10, 1500);
            sprintf(strawname,"Straw %.3i Predicted Drift Distance Vs phi_DOCA", straw);
            sprintf(strawtitle,"Ring %.2i Straw %.3i Predicted Drift Distance Vs phi_DOCA; #phi_{DOCA};Predicted Distance [cm]", ring, straw);
            Fill2DHistogram(regionStringPerStraw,folder,strawname, docaphi, predictedDistance,
                  strawtitle, 16, -3.14, 3.14,  400, 0.0, 1.2);
            char residualname[100];
            char residualtitle[256];
            sprintf(residualname,"Straw %.3i Residual", straw);
            sprintf(residualtitle,"Ring %.2i Straw %.3i Residual;Residual [cm]", ring, straw);
            Fill1DHistogram(regionStringPerStraw,folder,residualname, residual, residualtitle, 200, -0.05, 0.05);
            sprintf(residualname,"Straw %.3i Residual Vs. Z", straw);
            sprintf(residualtitle,"Ring %.2i Straw %.3i Residual;Z [cm]; Residual [cm]", ring, straw);
            Fill2DHistogram(regionStringPerStraw,folder,residualname, 
                  dz,residual, 
                  residualtitle, 
                  30, -75.0,75.0,200, -0.05, 0.05);

            //Time to distance relation in bins
            // Calcuate delta
            double delta = max_sag[ring - 1][straw - 1]*(1.-dz*dz/5625.)
               *cos(docaphi + sag_phi_offset[ring - 1][straw - 1]);
            sprintf(strawname,"Straw %.3i residual Vs delta", straw);
            sprintf(strawtitle,"Ring %.2i Straw %.3i Residual Vs #delta; #delta [cm]; Residual [cm]",ring,  straw);
            double binwidth = 0.005;
            if ( 2 * max_sag[ring - 1][straw - 1] > binwidth){
               Fill2DHistogram(regionStringPerStraw,folder,strawname, delta, residual,
                     strawtitle, Int_t(2 * max_sag[ring - 1][straw - 1] / binwidth), -1 * max_sag[ring - 1][straw - 1], max_sag[ring - 1][straw - 1], 100, -0.05, 0.05);
            }

            char binname[150];
            char bintitle[150];
            sprintf(binname,"Straw %.3i Predicted Drift Distance Vs. Drift Time", straw);
            sprintf(bintitle,"Ring %.2i Straw %.3i Predicted Drift Distance Vs. Drift Time", ring, straw);
            Fill2DHistogram(regionStringPerStraw,folder,binname, time, predictedDistance,
                  bintitle, 250, -50, 200, 250, 0.0, 0.4);

            sprintf(binname,"Straw %.3i Predicted Drift Distance Vs. delta", straw);
            sprintf(bintitle,"Ring %.2i Straw %.3i Predicted Drift Distance Vs. #delta;#delta [cm]; Predicted Drift Distance - Nominal Radius [cm]", ring, straw);
            Fill2DHistogram(regionStringPerStraw,folder,binname, delta, predictedDistance - 0.78,
                  bintitle, 20, -0.25, 0.25, 250, -0.25, 0.25);

            if (delta > 0){ // Long side of straw
               sprintf(binname,"Straw %.3i Predicted Drift Distance Vs. Drift Time Positive Delta", straw);
               sprintf(bintitle,"Ring %.2i Straw %.3i Predicted Drift Distance Vs. Drift Time (Positive Delta)", ring, straw);
               Fill2DHistogram(regionStringPerStraw,folder,binname, time, predictedDistance,
                     bintitle, 250, -10, 1500, 50, 0.0, 1.2);
               sprintf(binname,"Straw %.3i Residual Vs. Drift Time Positive Delta", straw);
               sprintf(bintitle,"Ring %.2i Straw %.3i Residual Vs. Drift Time (Positive Delta)", ring, straw);
               Fill2DHistogram(regionStringPerStraw,folder,binname, time, residual,
                     bintitle, 100, -10, 1500, 100, -0.05, 0.05);
               sprintf(binname,"Straw %.3i Residual Positive Delta", straw);
               sprintf(bintitle,"Ring %.2i Straw %.3i Residual (Positive Delta); Residual [cm]; Entries", ring, straw);
               Fill1DHistogram(regionStringPerStraw,folder,binname, residual, bintitle, 200, -0.05, 0.05);
            }
            else { // Short side of straw
               sprintf(binname,"Straw %.3i Predicted Drift Distance Vs. Drift Time Negative Delta", straw);
               sprintf(bintitle,"Ring %.2i Straw %.3i Predicted Drift Distance Vs. Drift Time (Negative Delta)", ring, straw);
               Fill2DHistogram(regionStringPerStraw,folder,binname, time, predictedDistance,
                     bintitle, 250, -10, 1500, 50, 0.0, 1.2);
               sprintf(binname,"Straw %.3i Residual Vs. Drift Time Negative Delta", straw);
               sprintf(bintitle,"Ring %.2i Straw %.3i Residual Vs. Drift Time (Negative Delta)", ring, straw);
               Fill2DHistogram(regionStringPerStraw,folder,binname, time, residual,
                     bintitle, 100, -10, 1500, 100, -0.05, 0.05);
               sprintf(binname,"Straw %.3i Residual Negative Delta", straw);
               sprintf(bintitle,"Ring %.2i Straw %.3i Residual (Negative Delta); Residual [cm]; Entries", ring, straw);
               Fill1DHistogram(regionStringPerStraw,folder,binname, residual, bintitle, 200, -0.05, 0.05);
            }

            Fill2DHistogram(regionString, folder, "Residual Vs. Drift Time", time, residual,
                  "Residual Vs. Drift Time; Drift Time [ns]; Residual [cm]", 
                  500, -10, 1500, 100, -0.05, 0.05);
            Fill2DHistogram(regionString, folder, "Residual Vs. Drift Distance", distance, residual,
                  "Residual Vs. Drift Distance; Drift Distance [cm]; Residual [cm]",
                  50, 0.0, 1.0, 100, -0.05, 0.05);
            Fill2DHistogram(regionString, folder, "Residual Vs. Predicted Drift Distance", predictedDistance, residual,
                  "Residual Vs. Predicted Drift Distance; Predicted Drift Distance [cm]; Residual [cm]",
                  50, 0.0, 1.0, 100, -0.05, 0.05);

            Fill2DHistogram(regionString, folder, "Predicted Drift Distance Vs. Drift Time", time, predictedDistance,
                  "Predicted Drift Distance Vs. Drift Time; Drift Time [ns]; Predicted Drift Distance [cm]",
                  500, -10, 1500, 100, 0.0, 1.0);

        } 
    }
    return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_CDC_PerStrawReco::erun(void)
{
   // This is called whenever the run number changes, before it is
   // changed to give you a chance to clean up before processing
   // events from the next run number.
   return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_CDC_PerStrawReco::fini(void)
{
   // Called before program exit after event processing is finished.
   //SortDirectories();
   return NOERROR;
}

