// $Id$
//
//    File: JEventProcessor_BCAL_ADC_4ns.cc
// Created: Fri Jul 21 10:41:38 EDT 2017
// Creator: dalton (on Linux gluon106.jlab.org 2.6.32-642.3.1.el6.x86_64 x86_64)
//
/*************************************
  This plugin is designed to calibrate the TDC times for the BCAL.make the histograms necessary to 
  in order to extract the ADC 4 ns offsets for each channel.
 *************************************/

#include "JEventProcessor_BCAL_ADC_4ns.h"
#include "HistogramTools.h" 
#include "BCAL/DBCALHit.h"
#include "BCAL/DBCALTDCHit.h"
#include "BCAL/DBCALCluster.h"
#include "BCAL/DBCALDigiHit.h"
#include "BCAL/DBCALPoint.h"
#include "BCAL/DBCALUnifiedHit.h"
#include "BCAL/DBCALGeometry.h"
#include "DANA/DStatusBits.h"
#include "PID/DChargedTrack.h"
#include "PID/DEventRFBunch.h"
#include "PID/DDetectorMatches.h"
#include "PID/DNeutralShower.h"
#include "PID/DVertex.h"
#include "TRACKING/DTrackTimeBased.h"
#include "TRIGGER/DL1Trigger.h"

using namespace jana;


// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_BCAL_ADC_4ns());
}
} // "C"


//------------------
// JEventProcessor_BCAL_ADC_4ns (Constructor)
//------------------
JEventProcessor_BCAL_ADC_4ns::JEventProcessor_BCAL_ADC_4ns()
{

}

//------------------
// ~JEventProcessor_BCAL_ADC_4ns (Destructor)
//------------------
JEventProcessor_BCAL_ADC_4ns::~JEventProcessor_BCAL_ADC_4ns()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_BCAL_ADC_4ns::init(void)
{
	// This is called once at program startup. 

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_BCAL_ADC_4ns::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_BCAL_ADC_4ns::evnt(JEventLoop *loop, uint64_t eventnumber)
{
   // First check that this is not a font panel trigger or no trigger
   bool goodtrigger=1;
   const DL1Trigger *trig = NULL;
   try {
       loop->GetSingle(trig);
   } catch (...) {}
   if (trig) {
       if (trig->fp_trig_mask){
           goodtrigger=0;
       }
   }
   if (!goodtrigger) {
       return NOERROR;
   }

   const DEventRFBunch *thisRFBunch = NULL;
   loop->GetSingle(thisRFBunch);

   vector <const DChargedTrack *> chargedTrackVector;
   loop->Get(chargedTrackVector);

   for (unsigned int iTrack = 0; iTrack < chargedTrackVector.size(); iTrack++){
      // get charge and choose pion hypothesis as most likely
      int charge = chargedTrackVector[iTrack]->Get_Charge();
      char q[2];
      const DChargedTrackHypothesis* bestHypothesis;
      if (charge>0) {
          bestHypothesis = chargedTrackVector[iTrack]->Get_Hypothesis(PiPlus);
          sprintf(q,"+");
      } else {
          bestHypothesis = chargedTrackVector[iTrack]->Get_Hypothesis(PiMinus);
          sprintf(q,"-");
      }
      if (bestHypothesis == NULL) continue;
      // Now from this hypothesis we can get the detector matches to the BCAL
      const DBCALShowerMatchParams* bcalMatch = bestHypothesis->Get_BCALShowerMatchParams();
      if (bcalMatch == NULL) continue; 
      const DSCHitMatchParams* scMatch = bestHypothesis->Get_SCHitMatchParams();
      if (scMatch == NULL) continue;
      DVector3 position = bestHypothesis->position();

      // We also need the reference trajectory, which is buried deep in there
      const DTrackTimeBased *timeBasedTrack = nullptr;
      bestHypothesis->GetSingle(timeBasedTrack);
      const DReferenceTrajectory *rt = timeBasedTrack->rt;
      if (timeBasedTrack->FOM < 0.0027) continue; // 3-sigma cut on tracking FOM

      if (timeBasedTrack->Ndof < 10) continue; // CDC: 5 params in fit, 10 dof => [15 hits]; FDC [10 hits]

      // Get the shower from the match
      const DBCALShower *thisShower = bcalMatch->dBCALShower;

      // Fill histograms based on the shower
      char  title[200];
      DVector3 proj_pos = rt->GetLastDOCAPoint();
      double pathLength, flightTime;

      // Get the points from the shower
      vector <const DBCALPoint*> pointVector;
      thisShower->Get(pointVector);

      // Loop over the points within the cluster
      for (unsigned int iPoint = 0; iPoint < pointVector.size(); iPoint++){
         const DBCALPoint *thisPoint = pointVector[iPoint];
         if (thisPoint->E() < 0.05) continue; // The timing is known not to be great for very low energy, so only use our best info 
         double rpoint = thisPoint->r();

         float zminhall = 0;
         float zmaxhall = 450; 
         if (rt->GetIntersectionWithRadius(rpoint,proj_pos, &pathLength, &flightTime)==NOERROR){
            // Now proj_pos contains the projected position of the track at this particular point within the BCAL
            char channame[255];
            sprintf(channame, "M%02iL%iS%i", thisPoint->module(), thisPoint->layer(), thisPoint->sector());
            double trackHitZ = proj_pos.z();

            vector <const DBCALHit*> hitVector;
            thisPoint->Get(hitVector);

                //double Deltat = thisPoint->t_US() - thisPoint->t_DS();
            double Deltat = hitVector[0]->t_raw - hitVector[1]->t_raw;
            if (hitVector[0]->end==1) Deltat = -Deltat;

            sprintf(title, "%s  Z_{Track} vs #Delta t;#Delta t = t_{US}-t_{DS};Z_{Track} [cm]", channame);
            Fill2DHistogram ("BCAL_ADC_Deltat", "ZvsDeltat", channame, Deltat, trackHitZ, title,
                             480, -30, 30, 250, zminhall, zmaxhall); 
         }
      }
   }
   return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_BCAL_ADC_4ns::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_BCAL_ADC_4ns::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}

