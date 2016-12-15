// $Id$
//
//    File: JEventProcessor_TrackingPulls.cc
// Created: Thu Nov  3 14:30:19 EDT 2016
// Creator: mstaib (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "JEventProcessor_TrackingPulls.h"
#include "HistogramTools.h"
#include "TRACKING/DTrackTimeBased.h"
#include "PID/DChargedTrack.h"

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
   // This is called once at program startup. 

   return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_TrackingPulls::brun(JEventLoop *eventLoop, int32_t runnumber)
{
   // This is called whenever the run number changes
   return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_TrackingPulls::evnt(JEventLoop *loop, uint64_t eventnumber)
{
   unsigned int numstraws[28]={42,42,54,54,66,66,80,80,93,93,106,106,123,123,
      135,135,146,146,158,158,170,170,182,182,197,197,
      209,209};
   // Loop over the tracks, get the tracking pulls, and fill some histograms. Easy peasy

   vector<const DChargedTrack *> chargedTrackVector;
   loop->Get(chargedTrackVector);

   for (size_t i = 0; i < chargedTrackVector.size(); i++){
      // TODO: Should be changed to use PID FOM when ready
      const DChargedTrackHypothesis *bestHypothesis = chargedTrackVector[i]->Get_BestTrackingFOM();

      if (bestHypothesis == NULL) continue;

      double trackingFOM = TMath::Prob(bestHypothesis->dChiSq_Track, bestHypothesis->dNDF_Track);

      // Some quality cuts for the tracks we will use
      // Keep this minimal for now and investigate later
      float trackingFOMCut = 0.001;
      unsigned int trackingNDFCut = 5;

      if(trackingFOM < trackingFOMCut) continue;
      if( bestHypothesis->dNDF_Track < trackingNDFCut) continue;

      double phi = bestHypothesis->momentum().Phi()*TMath::RadToDeg();
      double theta = bestHypothesis->momentum().Theta()*TMath::RadToDeg();
      double pmag = bestHypothesis->momentum().Mag();

      if (pmag < 0.5) continue;

      // Fill some track information
      Fill1DHistogram("TrackingPulls", "TrackInfo", "Tracking FOM",
            trackingFOM,
            "Tracking FOM", 200, 0.0, 1.0);
      Fill2DHistogram("TrackingPulls", "TrackInfo", "P Vs. Theta",
            theta,  pmag,
            "P Vs. #theta; #theta [deg.]; |P| [GeV/c]", 70, 0.0, 140.0, 50, 0.0, 10.0);
      Fill2DHistogram("TrackingPulls", "TrackInfo", "Phi Vs. Theta",
            theta,  phi,
            "#phi Vs. #theta; #theta [deg.];  #phi [deg.]", 70, 0.0, 140.0, 180, -180.0, 180.0);
      Fill2DHistogram("TrackingPulls", "TrackInfo", "P Vs. Phi",
            phi,  pmag,
            "P Vs. #phi; #phi [deg.]; |P| [GeV/c]", 180, -180, 180.0, 50, 0.0, 10.0);

      // Get the pulls vector from the track
      const DTrackTimeBased *thisTimeBasedTrack;
      bestHypothesis->GetSingle(thisTimeBasedTrack);

      vector<DTrackFitter::pull_t> pulls = thisTimeBasedTrack->pulls;
      for (size_t iPull = 0; iPull < pulls.size(); iPull++){
         // Here is all of the information currently stored in the pulls from the fit
         // From TRACKING/DTrackFitter.h
         double resi                 = pulls[iPull].resi;   // residual of measurement
         double err                  = pulls[iPull].err;      // estimated error of measurement
         //double s                    = pulls[iPull].s;
         double tdrift               = pulls[iPull].tdrift;      // drift time of this measurement
         //double d                    = pulls[iPull].d;  // doca to wire
         const DCDCTrackHit *cdc_hit = pulls[iPull].cdc_hit;
         const DFDCPseudo *fdc_hit   = pulls[iPull].fdc_hit;
         //double docaphi              = pulls[iPull].docaphi; // phi of doca in CDC straws
         double z                    = pulls[iPull].z;// z position at doca
         //double tcorr                = pulls[iPull].tcorr; // drift time with correction for B
         double resic                = pulls[iPull].resic; // residual for FDC cathode measuremtns
         double errc                 = pulls[iPull].errc;

         Fill1DHistogram("TrackingPulls", "TrackPulls","All Pulls",
               resi/err,
               "Residual/Error", 100, -5.0, 5.0);
         Fill2DHistogram("TrackingPulls", "TrackPulls","All Pulls Vs. P",
               pmag, resi/err,
               ";|P| ;Residual/Error", 100, 0.0, 10.0, 100, -5.0, 5.0);
         Fill2DHistogram("TrackingPulls", "TrackPulls","All Pulls Vs. Phi",
               phi, resi/err,
               ";#phi ;Residual/Error", 180, -180.0, 180.0, 100, -5.0, 5.0);
         Fill2DHistogram("TrackingPulls", "TrackPulls","All Pulls Vs. Theta",
               theta, resi/err,
               ";#theta ;Residual/Error", 140, 0.0, 140.0, 100, -5.0, 5.0);

         // Fill some detector specific info
         // Fill them in order = super-hacked
         static int nextPlane = 1;
         static int nextRing = 1;

         if (fdc_hit != nullptr && fdc_hit->wire->layer <= nextPlane){
            if(fdc_hit->wire->layer == nextPlane) nextPlane++;
            Fill1DHistogram("TrackingPulls", "FDCPulls","All Wire Pulls",
                  resi/err,
                  "Residual/Error", 100, -5.0, 5.0);
            Fill1DHistogram("TrackingPulls", "FDCPulls","All Cathode Pulls",
                  resic/errc,
                  "Residual/Error", 100, -5.0, 5.0);
            Fill1DHistogram("TrackingPulls", "FDCPulls","All Wire Residuals",
                  resi,
                  "Residual", 100, -0.1, 0.1);
            Fill2DHistogram("TrackingPulls", "FDCPulls","All Wire Residuals Vs. Plane",
                  fdc_hit->wire->layer, resi,
                  ";plane ;Residual", 24, 0.5, 24.5, 100, -0.1, 0.1);
            Fill2DHistogram("TrackingPulls", "FDCPulls","All Cathode Residuals Vs. Plane",
                  fdc_hit->wire->layer, resic,
                  ";plane ;Residual", 24, 0.5, 24.5, 100, -0.1, 0.1);
            Fill2DHistogram("TrackingPulls", "FDCPulls","All Wire Pulls Vs. Plane",
                  fdc_hit->wire->layer, resi/err,
                  ";plane ;Residual/Error", 24, 0.5, 24.5, 100, -5.0, 5.0);
            Fill2DHistogram("TrackingPulls", "FDCPulls","All Wire Residuals Vs Drift Time",
                  tdrift, resi,
                  ";Drift Time;Residual", 170, -20.0, 150.0, 100, -0.1, 0.1);
            Fill2DHistogram("TrackingPulls", "FDCPulls","All Wire Pulls Vs Drift Time",
                  tdrift, resi/err,
                  ";Drift Time;Residual/Error", 170, -20.0, 150.0, 100, -5.0, 5.0);
            Fill1DHistogram("TrackingPulls", "FDCPulls","All Cathode Residuals",
                  resic,
                  "Residual", 100, -0.1, 0.1);
            Fill2DHistogram("TrackingPulls", "FDCPulls","All Wire Pulls Vs. P",
                  pmag, resic/errc,
                  ";|P| ;Residual/Error", 100, 0.0, 10.0, 100, -5.0, 5.0);
            Fill2DHistogram("TrackingPulls", "FDCPulls","All Wire Pulls Vs. Phi",
                  phi, resi/err,
                  ";#phi ;Residual/Error", 180, -180.0, 180.0, 100, -5.0, 5.0);
            Fill2DHistogram("TrackingPulls", "FDCPulls","All Wire Pulls Vs. Theta",
                  theta, resi/err,
                  ";#theta ;Residual/Error", 50, 0.0, 25.0, 100, -5.0, 5.0);
            Fill2DHistogram("TrackingPulls", "FDCPulls","All Cathode Pulls Vs. P",
                  pmag, resi/err,
                  ";|P| ;Residual/Error", 100, 0.0, 10.0, 100, -5.0, 5.0);
            Fill2DHistogram("TrackingPulls", "FDCPulls","All Cathode Pulls Vs. Phi",
                  phi, resi/err,
                  ";#phi ;Residual/Error", 180, -180.0, 180.0, 100, -5.0, 5.0);
            Fill2DHistogram("TrackingPulls", "FDCPulls","All Cathode Pulls Vs. Theta",
                  theta, resi/err,
                  ";#theta ;Residual/Error", 50, 0.0, 25.0, 100, -5.0, 5.0);

            // Make the Per-Plane Histograms
            char planeName[256];
            sprintf(planeName,"FDCPulls_Plane%.2i", fdc_hit->wire->layer);

            Fill1DHistogram("TrackingPulls", planeName,"All Wire Pulls",
                  resi/err,
                  "Residual/Error", 100, -5.0, 5.0);
            Fill1DHistogram("TrackingPulls", planeName,"All Cathode Pulls",
                  resic/errc,
                  "Residual/Error", 100, -5.0, 5.0);
            Fill1DHistogram("TrackingPulls", planeName,"All Wire Residuals",
                  resi,
                  "Residual", 100, -0.1, 0.1);
            Fill2DHistogram("TrackingPulls", planeName,"All Wire Residuals Vs Drift Time",
                  tdrift, resi,
                  ";Drift Time;Residual", 170, -20.0, 150.0, 100, -0.1, 0.1);
            Fill2DHistogram("TrackingPulls", planeName,"All Wire Pulls Vs Drift Time",
                  tdrift, resi/err,
                  ";Drift Time;Residual/Error", 170, -20.0, 150.0, 100, -5.0, 5.0);
            Fill1DHistogram("TrackingPulls", planeName,"All Cathode Residuals",
                  resic,
                  "Residual", 100, -0.1, 0.1);
            Fill2DHistogram("TrackingPulls", planeName,"All Wire Pulls Vs. P",
                  pmag, resi/err,
                  ";|P| ;Residual/Error", 100, 0.0, 10.0, 100, -5.0, 5.0);
            Fill2DHistogram("TrackingPulls", planeName,"All Wire Pulls Vs. Phi",
                  phi, resi/err,
                  ";#phi ;Residual/Error", 180, -180.0, 180.0, 100, -5.0, 5.0);
            Fill2DHistogram("TrackingPulls", planeName,"All Wire Pulls Vs. Theta",
                  theta, resi/err,
                  ";#theta ;Residual/Error", 50, 0.0, 25.0, 100, -5.0, 5.0);
            Fill2DHistogram("TrackingPulls", planeName,"All Wire Residuals Vs. P",
                  pmag, resi,
                  ";|P| ;Residual", 100, 0.0, 10.0, 100, -0.1, 0.1);
            Fill2DHistogram("TrackingPulls", planeName,"All Wire Residuals Vs. Phi",
                  phi, resi,
                  ";#phi ;Residual", 180, -180.0, 180.0, 100, -0.1, 0.1);
            Fill2DHistogram("TrackingPulls", planeName,"All Wire Residuals Vs. Theta",
                  theta, resi,
                  ";#theta ;Residual", 50, 0.0, 25.0, 100, -0.1, 0.1);
            Fill2DHistogram("TrackingPulls", planeName,"All Cathode Pulls Vs. P",
                  pmag, resic/errc,
                  ";|P| ;Residual/Error", 100, 0.0, 10.0, 100, -5.0, 5.0);
            Fill2DHistogram("TrackingPulls", planeName,"All Cathode Pulls Vs. Phi",
                  phi, resic/errc,
                  ";#phi ;Residual/Error", 180, -180.0, 180.0, 100, -5.0, 5.0);
            Fill2DHistogram("TrackingPulls", planeName,"All Cathode Pulls Vs. Theta",
                  theta, resic/errc,
                  ";#theta ;Residual/Error", 140, 0.0, 140.0, 100, -5.0, 5.0);
            Fill2DHistogram("TrackingPulls", planeName,"All Cathode Residuals Vs. P",
                  pmag, resic,
                  ";|P| ;Residual", 100, 0.0, 10.0, 100, -0.1, 0.1);
            Fill2DHistogram("TrackingPulls", planeName,"All Cathode Residuals Vs. Phi",
                  phi, resic,
                  ";#phi ;Residual", 180, -180.0, 180.0, 100, -0.1, 0.1);
            Fill2DHistogram("TrackingPulls", planeName,"All Cathode Residuals Vs. Theta",
                  theta, resic,
                  ";#theta ;Residual", 140, 0.0, 140.0, 100, -0.1, 0.1);
            Fill2DHistogram("TrackingPulls", planeName,"Wire Pulls",
                  fdc_hit->wire->wire,resi/err,
                  ";Wire Number ;Residual/Error", 96, 0.5, 96.5, 100, -5.0, 5.0);
            Fill2DHistogram("TrackingPulls", planeName,"Wire Residuals",
                  fdc_hit->wire->wire,resi,
                  ";Wire Number ;Residual", 96, 0.5, 96.5, 100, -0.1, 0.1);
            if(fabs(resi/err) < 5.0){
               Fill2DProfile("TrackingPulls", planeName,"2D Wire Hit Pulls",
                     fdc_hit->xy.X(), fdc_hit->xy.Y(), resi/err,
                     "Mean of Wire Pulls vs. PseudoHit XY",
                     100, -50., 50., 100, -50., 50.);
            }
            if(fabs(resi) < 0.1){
               Fill2DProfile("TrackingPulls", planeName,"2D Wire Hit Residuals",
                     fdc_hit->xy.X(), fdc_hit->xy.Y(), resi,
                     "Mean of Wire Residuals vs. PseudoHit XY",
                     100, -50., 50., 100, -50., 50.);
               Fill2DProfile("TrackingPulls", planeName,"2D Wire Hit Residuals Local",
                     fdc_hit->w, fdc_hit->s, resi,
                     "Mean of Wire Residuals vs. PseudoHit WS;Perpendicular Distance to Wire; Distance Along the Wire",
                     100, -50., 50., 100, -50., 50.);
            }
            if(fabs(resic/errc) < 5.0){
               Fill2DProfile("TrackingPulls", planeName,"2D Cathode Hit Pulls",
                     fdc_hit->xy.X(), fdc_hit->xy.Y(), resic/errc,
                     "Mean of Cathode Pulls vs. PseudoHit XY",
                     100, -50., 50., 100, -50., 50.);
            }
            if(fabs(resic) < 0.1){
               Fill2DProfile("TrackingPulls", planeName,"2D Cathode Hit Residuals",
                     fdc_hit->xy.X(), fdc_hit->xy.Y(), resic,
                     "Mean of Cathode Residuals vs. PseudoHit XY",
                     100, -50., 50., 100, -50., 50.);
               Fill2DProfile("TrackingPulls", planeName,"2D Cathode Hit Residuals Local",
                     fdc_hit->w, fdc_hit->s, resic,
                     "Mean of Cathode Residuals vs. PseudoHit WS;Perpendicular Distance to Wire; Distance Along the Wire",
                     100, -50., 50., 100, -50., 50.);
            }
         }

         // Once we are done with the FDC, move on to the CDC.
         if (cdc_hit != nullptr && cdc_hit->wire->ring <= nextRing && nextPlane == 25){
            if(cdc_hit->wire->ring == nextRing) nextRing++;

            Fill2DHistogram("TrackingPulls", "CDCPulls","All Pulls Vs. Ring",
                  cdc_hit->wire->ring, resi/err,
                  ";ring ;Residual/Error", 28, 0.5, 28.5, 100, -5.0, 5.0);
            Fill2DHistogram("TrackingPulls", "CDCPulls","All Residuals Vs. Ring",
                  cdc_hit->wire->ring, resi,
                  ";ring ;Residual", 28, 0.5, 28.5, 100, -0.1, 0.1);
            Fill2DHistogram("TrackingPulls", "CDCPulls","All Pulls Vs. tdrift",
                  tdrift, resi/err,
                  ";tdrift [ns] ;Residual/Error", 200, 0.0, 1000.0, 100, -5.0, 5.0);
            Fill2DHistogram("TrackingPulls", "CDCPulls","All Residuals Vs. tdrift",
                  tdrift, resi,
                  ";tdrift [ns] ;Residual", 200, 0.0, 1000.0, 100, -0.1, 0.1);
            Fill2DHistogram("TrackingPulls", "CDCPulls","All Pulls Vs. P",
                  pmag, resi/err,
                  ";|P| ;Residual/Error", 50, 0.0, 10.0, 100, -5.0, 5.0);
            Fill2DHistogram("TrackingPulls", "CDCPulls","All Pulls Vs. Phi",
                  phi, resi/err,
                  ";#phi ;Residual/Error", 180, -180.0, 180.0, 100, -5.0, 5.0);
            Fill2DHistogram("TrackingPulls", "CDCPulls","All Pulls Vs. Theta",
                  theta, resi/err,
                  ";#theta ;Residual/Error", 140, 0.0, 140.0, 100, -5.0, 5.0);
            Fill2DHistogram("TrackingPulls", "CDCPulls","All Residuals Vs. P",
                  pmag, resi,
                  ";|P| ;Residual", 50, 0.0, 10.0, 100, -0.1, 0.1);
            Fill2DHistogram("TrackingPulls", "CDCPulls","All Residuals Vs. Phi",
                  phi, resi,
                  ";#phi ;Residual", 180, -180.0, 180.0, 100, -0.1, 0.1);
            Fill2DHistogram("TrackingPulls", "CDCPulls","All Residuals Vs. Theta",
                  theta, resi,
                  ";#theta ;Residual", 140, 0.0, 140.0, 100, -0.1, 0.1);

            // Make the Per-Ring Histograms
            char ringName[256];
            sprintf(ringName,"CDCPulls_Ring%.2i", cdc_hit->wire->ring);

            Fill2DHistogram("TrackingPulls", ringName,"All Pulls Vs. tdrift",
                  tdrift, resi/err,
                  ";tdrift [ns] ;Residual/Error", 200, 0.0, 1000.0, 100, -5.0, 5.0);
            Fill2DHistogram("TrackingPulls", ringName,"All Pulls Vs. z",
                  z, resi/err,
                  ";z [cm] ;Residual/Error", 200, -30., 200., 100, -5.0, 5.0);
            Fill2DHistogram("TrackingPulls", ringName,"All Residuals Vs. tdrift",
                  tdrift, resi,
                  ";tdrift [ns] ;Residual", 200, 0.0, 1000.0, 100, -0.1, 0.1);
            Fill2DHistogram("TrackingPulls", ringName,"All Residuals Vs. z",
                  z, resi,
                  ";z [cm] ;Residual", 200, -30., 200.0, 100, -0.1, 0.1);
            Fill2DHistogram("TrackingPulls", ringName,"All Pulls Vs. P",
                  pmag, resi/err,
                  ";|P| ;Residual/Error", 50, 0.0, 10.0, 100, -5.0, 5.0);
            Fill2DHistogram("TrackingPulls", ringName,"All Pulls Vs. Phi",
                  phi, resi/err,
                  ";#phi ;Residual/Error", 180, -180.0, 180.0, 100, -5.0, 5.0);
            Fill2DHistogram("TrackingPulls", ringName,"All Pulls Vs. Theta",
                  theta, resi/err,
                  ";#theta ;Residual/Error", 140, 0.0, 140.0, 100, -5.0, 5.0);
            Fill2DHistogram("TrackingPulls", ringName,"All Residuals Vs. P",
                  pmag, resi,
                  ";|P| ;Residual", 50, 0.0, 10.0, 100, -0.1, 0.1);
            Fill2DHistogram("TrackingPulls", ringName,"All Residuals Vs. Phi",
                  phi, resi,
                  ";#phi ;Residual", 180, -180.0, 180.0, 100, -0.1, 0.1);
            Fill2DHistogram("TrackingPulls", ringName,"All Residuals Vs. Theta",
                  theta, resi,
                  ";#theta ;Residual", 140, 0.0, 140.0, 100, -0.1, 0.1);

            double nStraw = numstraws[cdc_hit->wire->ring-1];
            double phiIntersect = (cdc_hit->wire->origin + (z-92.0)*cdc_hit->wire->udir).Phi() * TMath::RadToDeg();

            Fill2DHistogram("TrackingPulls", ringName,"Per Straw Pulls",
                  cdc_hit->wire->straw, resi/err,
                  ";Straw Number ;Residual/Error", nStraw, 0.5, nStraw+0.5, 100, -5.0, 5.0);
            Fill2DHistogram("TrackingPulls", ringName,"Per Straw Residuals",
                  cdc_hit->wire->straw, resi,
                  ";Straw Number ;Residual", nStraw, 0.5, nStraw+0.5, 100, -0.1, 0.1);

            if (fabs(resi) < 0.1){
               Fill2DProfile("TrackingPulls", ringName, "Residual Vs Phi-Theta",
                     theta, phi, resi,
                     ";#theta;#phi", 70, 0.0, 140.0, 180, -180.0, 180.0); 
               Fill2DProfile("TrackingPulls", ringName, "Residual Vs Phi-z",
                     z, phi, resi,
                     ";z;#phi", 200, 0.0, 200.0,  180, -180.0, 180.0);
               Fill2DProfile("TrackingPulls", ringName, "Residual Vs PhiIntersect-z",
                     z, phiIntersect, resi,
                     ";z;#phi Intersect", 200, 0.0, 200.0, nStraw, -180.0, 180.0);
               Fill2DProfile("TrackingPulls", ringName, "Residual Vs P-Theta",
                     theta, pmag, resi,
                     ";#theta;|P|", 70, 0.0, 140.0, 50, 0.0, 10.0);
            }

         }
      }
   }

   return NOERROR;
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

