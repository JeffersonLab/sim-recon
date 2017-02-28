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
#include "TRIGGER/DL1Trigger.h"
#include "DANA/DStatusBits.h"

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
jerror_t JEventProcessor_BCAL_TDC_Timing::brun(JEventLoop *loop, int32_t runnumber)
{
    // This is called whenever the run number changes
    DApplication* app = dynamic_cast<DApplication*>(loop->GetJApplication());
    DGeometry* geom = app->GetDGeometry(runnumber);
    geom->GetTargetZ(Z_TARGET);

    //////
    // THIS NEEDS TO CHANGE IF THERE IS A NEW TABLE
    //////
    //get timewalk corrections from CCDB
    JCalibration *jcalib = loop->GetJCalibration();
    //these tables hold: module layer sector end c0 c1 c2 c3
    vector<vector<float> > tdc_timewalk_table;
    jcalib->Get("BCAL/timewalk_tdc",tdc_timewalk_table);

    for (vector<vector<float> >::const_iterator iter = tdc_timewalk_table.begin();
          iter != tdc_timewalk_table.end();
          ++iter) {
       if (iter->size() != 8) {
          cout << "DBCALUnifiedHit_factory: Wrong number of values in timewalk_tdc table (should be 8)" << endl;
          continue;
       }
       //be really careful about float->int conversions
       int module = (int)((*iter)[0]+0.5);
       int layer  = (int)((*iter)[1]+0.5);
       int sector = (int)((*iter)[2]+0.5);
       int endi   = (int)((*iter)[3]+0.5);
       DBCALGeometry::End end = (endi==0) ? DBCALGeometry::kUpstream : DBCALGeometry::kDownstream;
       float c0 = (*iter)[4];
       float c1 = (*iter)[5];
       float c2 = (*iter)[6];
       float a_thresh = (*iter)[7];
       int cellId = DBCALGeometry::cellId(module, layer, sector);
       readout_channel channel(cellId,end);
       tdc_timewalk_map[channel] = timewalk_coefficients(c0,c1,c2,a_thresh);
    }

	/// Read in calibration constants
	vector<double> raw_channel_global_offset;
    //if(print_messages) jout << "In BCAL_TDC_Timing, loading constants..." << endl;
	if(loop->GetCalib("/BCAL/channel_global_offset", raw_channel_global_offset))
		jout << "Error loading /BCAL/channel_global_offset !" << endl;

	japp->RootFillLock(this);
	TH1D *CCDB_raw_channel_global_offset = new TH1D("CCDB_raw_channel_global_offset","Offsets at time of running;channel;offset",
													768,0.5,768.5);
	int counter = 1;
	for (vector<double>::iterator iter = raw_channel_global_offset.begin(); iter != raw_channel_global_offset.end(); ++iter) {
		CCDB_raw_channel_global_offset->SetBinContent(counter,*iter);
		counter++;
	}
	japp->RootFillUnLock(this);

    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_BCAL_TDC_Timing::evnt(JEventLoop *loop, uint64_t eventnumber)
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
   } else {
	   // HDDM files are from simulation, so keep them even though they have no trigger
	   bool locIsHDDMEvent = loop->GetJEvent().GetStatusBit(kSTATUS_HDDM);
	   if (!locIsHDDMEvent) goodtrigger=0;
   }

   if (!goodtrigger) {
	   return NOERROR;
   }

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
      int cellId = DBCALGeometry::cellId(bcalUnifiedHitVector[i]->module, bcalUnifiedHitVector[i]->layer, bcalUnifiedHitVector[i]->sector);
      // Get the underlying associated objects
      const DBCALHit * thisADCHit;
      const DBCALTDCHit * thisTDCHit;
      bcalUnifiedHitVector[i]->GetSingle(thisADCHit);
      bcalUnifiedHitVector[i]->GetSingle(thisTDCHit);

      // From the ADC hit we can extract the pulse peak information...
      int pulse_peak = thisADCHit->pulse_peak;

      // The raw information from the DBCALHit and DBCALTDCHit is not corrected for timewalk yet, so we can always plot the before and after.
      if (thisADCHit != NULL && thisTDCHit != NULL){
         char name[200];
         sprintf(name, "Module%.2iLayer%.2iSector%.2i", bcalUnifiedHitVector[i]->module, bcalUnifiedHitVector[i]->layer, bcalUnifiedHitVector[i]->sector);
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
         double correctedTDCTime = bcalUnifiedHitVector[i]->t_TDC;
         if (bcalUnifiedHitVector[i]->t_TDC != bcalUnifiedHitVector[i]->t){ // Time walk correction has not been applied
            // Apply the timewalk correction
            readout_channel channel(cellId,bcalUnifiedHitVector[i]->end);
            timewalk_coefficients tdc_coeff = tdc_timewalk_map[channel];
            // HERE IS THE TIMEWALK CORRECTION
            correctedTDCTime -= tdc_coeff.c0 + tdc_coeff.c1/pow(pulse_peak/tdc_coeff.a_thresh, tdc_coeff.c2);
         }
         char name[200];
         sprintf(name, "Module %.2i Layer %.2i Sector %.2i", bcalUnifiedHitVector[i]->module, bcalUnifiedHitVector[i]->layer, bcalUnifiedHitVector[i]->sector);
         if (bcalUnifiedHitVector[i]->end == 0){
            Fill2DHistogram ("BCAL_TDC_Timing", "BCAL_Upstream_Timewalk_WithCorrection_E", name,
                  bcalUnifiedHitVector[i]->E, correctedTDCTime - bcalUnifiedHitVector[i]->t_ADC,
                  "Timewalk; E [GeV]; t_{TDC} - t_{ADC} [ns]",
                  250, 0.0, 0.35, NBINS_TDIFF, MIN_TDIFF , MAX_TDIFF );
            Fill2DHistogram ("BCAL_TDC_Timing", "BCAL_Upstream_Timewalk_WithCorrection_PP", name,
                  pulse_peak, correctedTDCTime - bcalUnifiedHitVector[i]->t_ADC,
                  "Timewalk; Pulse Peak [ADC Counts]; t_{TDC} - t_{ADC} [ns]",
                  500, 0.0, 1500.0, NBINS_TDIFF, MIN_TDIFF , MAX_TDIFF );
            Fill2DHistogram ("BCAL_TDC_Timing", "BCAL", "BCALHit Upstream Per Channel TDC-ADC Hit Time - With Timewalk",
                  the_tdc_cell, correctedTDCTime - bcalUnifiedHitVector[i]->t_ADC,
                  "BCALHit Upstream Per Channel TDC-ADC Hit Time - TW Corrected; cellID; t_{TDC} - t_{ADC} [ns] ",
                  576, 0.5, 576.5, NBINS_TDIFF, MIN_TDIFF, MAX_TDIFF);
         }
         else{
            Fill2DHistogram ("BCAL_TDC_Timing", "BCAL_Downstream_Timewalk_WithCorrection_E", name,
                  bcalUnifiedHitVector[i]->E, correctedTDCTime - bcalUnifiedHitVector[i]->t_ADC,
                  "Timewalk; E [GeV]; t_{TDC} - t_{ADC} [ns]",
                  250, 0.0, 0.35, NBINS_TDIFF, MIN_TDIFF, MAX_TDIFF);
            Fill2DHistogram ("BCAL_TDC_Timing", "BCAL_Downstream_Timewalk_WithCorrection_PP", name,
                  pulse_peak, correctedTDCTime - bcalUnifiedHitVector[i]->t_ADC,
                  "Timewalk; Pulse Peak [ADC Counts]; t_{TDC} - t_{ADC} [ns]",
                  500, 0.0, 1500.0, NBINS_TDIFF, MIN_TDIFF, MAX_TDIFF);
            Fill2DHistogram ("BCAL_TDC_Timing", "BCAL", "BCALHit Downstream Per Channel TDC-ADC Hit Time - With Timewalk",
                  the_tdc_cell, correctedTDCTime - bcalUnifiedHitVector[i]->t_ADC,
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
      // const DChargedTrackHypothesis* bestHypothesis = chargedTrackVector[iTrack]->Get_BestTrackingFOM();

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
      const DSCHitMatchParams* scMatch = bestHypothesis->Get_SCHitMatchParams(); // Needed for quality cut later
      if (bcalMatch == NULL) continue; 

      // We also need the reference trajectory, which is buried deep in there
      const DTrackTimeBased *timeBasedTrack = nullptr;
      bestHypothesis->GetSingle(timeBasedTrack);
      const DReferenceTrajectory *rt = timeBasedTrack->rt;

	  // Use CDC dEdx to help reject protons
	  double dEdx=1e6*timeBasedTrack->ddEdx_CDC;
	  double P=timeBasedTrack->momentum().Mag();
	  bool dEdx_pion = 0;
	  if (dEdx<2.5) dEdx_pion = 1;

      // Get the shower from the match
      const DBCALShower *thisShower = bcalMatch->dBCALShower;

	  // Fill histograms based on the shower
	  char name[200], title[200];
	  DVector3 proj_pos = rt->GetLastDOCAPoint();
	  double pathLength, flightTime;
	  double r_shower = sqrt(thisShower->x*thisShower->x+thisShower->y*thisShower->y);
	  double t_shower = thisShower->t;
	  double E_shower = thisShower->E;
	  if (rt->GetIntersectionWithRadius(r_shower,proj_pos, &pathLength, &flightTime)==NOERROR){
		  if (thisRFBunch->dNumParticleVotes >= 2 && scMatch != NULL){ // Require good RF bunch and this track match the SC
			  // We have the flight time to our BCAL point, so we can get the target time
			  double targetCenterTime = t_shower - flightTime - ((timeBasedTrack->position()).Z() - Z_TARGET) / SPEED_OF_LIGHT;
			  sprintf(name, "AllShowers_q%s", q);
			  sprintf(title, "Timing resolution; E [GeV]; t_{Target} - t_{RF} [ns]");
			  Fill2DHistogram("BCAL_Global_Offsets", "Showers", name,
							  E_shower, targetCenterTime - thisRFBunch->dTime, title,
							  1000, 0.0, 5.0, 200, -10, 10);
			  sprintf(name, "dEdxVsP_q%s", q);
			  sprintf(title, "CDC dE/dx vs P; P [GeV]; dE/dx [keV/cm]");
			  Fill2DHistogram("BCAL_Global_Offsets", "Showers_PID", name,
							  P, dEdx, title,
							  200, 0.0, 5.0, 200, 0.0, 5.0);
			  sprintf(name, "deltaTVsdEdx_q%s", q);
			  sprintf(title, "PID; dE/dx [keV/cm]; t_{Target} - t_{RF} [ns]");
			  Fill2DHistogram("BCAL_Global_Offsets", "Showers_PID", name,
							  dEdx, targetCenterTime - thisRFBunch->dTime, title,
							  200, 0.0, 5.0, 200, -10, 10);
			  sprintf(name, "EoverPVsP_q%s", q);
			  sprintf(title, "PID; dE/dx [keV/cm]; t_{Target} - t_{RF} [ns]");
			  Fill2DHistogram("BCAL_Global_Offsets", "Showers_PID", name,
							  P, E_shower/P, title,
							  200, 0.0, 5.0, 200, 0, 2);
			  if (dEdx_pion) {
				  sprintf(name, "PionShowers_q%s", q);
				  sprintf(title, "Timing resolution; E [GeV]; t_{Target} - t_{RF} [ns]");
				  Fill2DHistogram("BCAL_Global_Offsets", "Showers", name,
								  E_shower, targetCenterTime - thisRFBunch->dTime, title,
								  1000, 0.0, 5.0, 200, -10, 10);
			  }
		  }
	  }


      // Get the points from the shower
      vector <const DBCALPoint*> pointVector;
      thisShower->Get(pointVector);

      // Loop over the points within the cluster
      for (unsigned int iPoint = 0; iPoint < pointVector.size(); iPoint++){
         const DBCALPoint *thisPoint = pointVector[iPoint];
         //if (thisPoint->E() < 0.05) continue; // The timing is known not to be great for very low energy, so only use our best info 
         double rpoint = thisPoint->r();
         if (rt->GetIntersectionWithRadius(rpoint,proj_pos, &pathLength, &flightTime)==NOERROR){
            // Now proj_pos contains the projected position of the track at this particular point within the BCAL
            // We can plot the difference of the projected position and the BCAL position as a function of the channel
            sprintf(name , "Module%.2iLayer%.2iSector%.2i", thisPoint->module(), thisPoint->layer(), thisPoint->sector());
            // These results are in slightly different coordinate systems. We want one where the center of the BCAL is z=0
            double localTrackHitZ = proj_pos.z() - DBCALGeometry::GetBCAL_center();
            double localBCALHitZ = thisPoint->z() - DBCALGeometry::GetBCAL_center() + Z_TARGET;
			Fill2DHistogram ("BCAL_TDC_Offsets", "Z Position", "AllPoints",
							 localTrackHitZ, localBCALHitZ,
							 "Z_{point} Vs. Z_{Track}; Z_{Track} [cm]; Z_{Point} [cm]",
							 500, -250, 250, 500, -250, 250); 
            Fill2DHistogram ("BCAL_TDC_Offsets", "Z Position", name,
                  localTrackHitZ, localBCALHitZ,
                  "Z_{point} Vs. Z_{Track}; Z_{Track} [cm]; Z_{Point} [cm]",
                  500, -250, 250, 500, -250, 250); 

            // Now fill some histograms that are useful for aligning the BCAL with the rest of the detector systems
            if (thisRFBunch->dNumParticleVotes >= 2 && scMatch != NULL && dEdx_pion){ // Require good RF bunch and this track match the SC
               // Get the time of the BCAL point
               double pointTime = thisPoint->t();
               // We have the flight time to our BCAL point, so we can get the target time
               double targetCenterTime = pointTime - flightTime - ((timeBasedTrack->position()).Z() - Z_TARGET) / SPEED_OF_LIGHT;

               // Now we just plot the difference in from the RF Time to get out the correction
			   if (thisPoint->E() > 0.05) { // The timing is known not to be great for very low energy, so only use our best info 
				   int the_cell = (thisPoint->module() - 1) * 16 + (thisPoint->layer() - 1) * 4 + thisPoint->sector();
				   Fill2DHistogram("BCAL_Global_Offsets", "Target Time", "deltaTVsCell",
								   the_cell, targetCenterTime - thisRFBunch->dTime,
								   "Target Time Minus RF Time Vs. Cell Number; CCDB Index; t_{Target} - t_{RF} [ns]",
								   768, 0.5, 768.5, 200, -10, 10);
				   sprintf(name , "deltaTVsCell_q%s", q);
				   Fill2DHistogram("BCAL_Global_Offsets", "Target Time", name,
								   the_cell, targetCenterTime - thisRFBunch->dTime,
								   "Target Time Minus RF Time Vs. Cell Number; CCDB Index; t_{Target} - t_{RF} [ns]",
								   768, 0.5, 768.5, 200, -10, 10);
			   }
			   // Get the unifiedhits
			   vector <const DBCALUnifiedHit*> unifiedhitVector;
			   thisPoint->Get(unifiedhitVector);

			   const DBCALUnifiedHit *thisUnifiedhit1 = unifiedhitVector[0];
			   const DBCALUnifiedHit *thisUnifiedhit2 = unifiedhitVector[1];
			   float t1 = thisUnifiedhit1->t;
			   float t_ADC1 = thisUnifiedhit1->t_ADC;
			   float t_TDC1 = thisUnifiedhit1->t_TDC;
			   float t2 = thisUnifiedhit2->t;
			   float t_ADC2 = thisUnifiedhit2->t_ADC;
			   float t_TDC2 = thisUnifiedhit2->t_TDC;
			   char type[10];
			   sprintf(type,"Mixed");
			   if (t1 == t_ADC1 && t2 == t_ADC2) sprintf(type,"ADC");
			   if (t1 == t_TDC1 && t2 == t_TDC2) sprintf(type,"TDC");

			   const DBCALHit * thisADCHit1;
			   thisUnifiedhit1->GetSingle(thisADCHit1);
			   const DBCALHit * thisADCHit2;
			   thisUnifiedhit2->GetSingle(thisADCHit2);

			   float pulse_peak_max = max(thisADCHit1->pulse_peak,thisADCHit2->pulse_peak);
			   float pulse_peak_min = min(thisADCHit1->pulse_peak,thisADCHit2->pulse_peak);
			   float E_point = thisPoint->E();

			   sprintf(title, "Timing resolution; E [GeV]; t_{Target} - t_{RF} [ns]");
			   sprintf(name, "AllHits_%s_q%s",type,q);
			   Fill2DHistogram("BCAL_Global_Offsets", "Hits_deltaTVsE", name,
							   E_point, targetCenterTime - thisRFBunch->dTime, title,
							   1000, 0.0, 2.0, 200, -10, 10);

			   sprintf(title, "Timing resolution; peak [counts]; t_{Target} - t_{RF} [ns]");
			   sprintf(name, "AllHits_%s_q%s",type,q);
			   Fill2DHistogram("BCAL_Global_Offsets", "Hits_deltaTVsPPmax", name,
							   pulse_peak_max, targetCenterTime - thisRFBunch->dTime, title,
							   1000, 0.0, 4000, 200, -10, 10);
			   sprintf(name, "AllHits_%s_q%s",type,q);
			   Fill2DHistogram("BCAL_Global_Offsets", "Hits_deltaTVsPPmin", name,
							   pulse_peak_min, targetCenterTime - thisRFBunch->dTime, title,
							   1000, 0.0, 4000, 200, -10, 10);

			   int layer = thisPoint->layer();
			   sprintf(name, "AllPoints_q%s", q);
			   sprintf(title, "Timing resolution, all points; E [GeV]; t_{Target} - t_{RF} [ns]");
			   Fill2DHistogram("BCAL_Global_Offsets", "Points_deltaTVsEnergy", name,
							   E_point, targetCenterTime - thisRFBunch->dTime, title,
							   1000, 0.0, 2.0, 200, -10, 10);
			   sprintf(name, "Layer%i_q%s", layer, q);
			   sprintf(title, "Timing resolution, layer %i; E [GeV]; t_{Target} - t_{RF} [ns]", layer);
			   Fill2DHistogram("BCAL_Global_Offsets", "Points_deltaTVsEnergy", name,
							   E_point, targetCenterTime - thisRFBunch->dTime, title,
							   1000, 0.0, 2.0, 200, -10, 10);
			   sprintf(name, "AllPoints_%s_q%s", type, q);
			   sprintf(title, "Timing resolution, all points; E [GeV]; t_{Target} - t_{RF} [ns]");
			   Fill2DHistogram("BCAL_Global_Offsets", "Points_deltaTVsEnergy", name,
							   E_point, targetCenterTime - thisRFBunch->dTime, title,
							   1000, 0.0, 2.0, 200, -10, 10);
			   sprintf(name, "Layer%i_%s_q%s", layer, type, q);
			   sprintf(title, "Timing resolution, layer %i; E [GeV]; t_{Target} - t_{RF} [ns]", layer);
			   Fill2DHistogram("BCAL_Global_Offsets", "Points_deltaTVsEnergy", name,
							   E_point, targetCenterTime - thisRFBunch->dTime, title,
							   1000, 0.0, 2.0, 200, -10, 10);
			   sprintf(name, "AllPoints_q%s", q);
			   sprintf(title, "Timing resolution, all points; E [GeV]; t_{Target} - t_{RF} [ns]");
			   Fill2DHistogram("BCAL_Global_Offsets", "Points_deltaTVsShowerEnergy", name,
							   E_shower, targetCenterTime - thisRFBunch->dTime, title,
							   1000, 0.0, 2.0, 200, -10, 10);
			   sprintf(name, "Layer%i_q%s", layer, q);
			   sprintf(title, "Timing resolution, layer %i; E [GeV]; t_{Target} - t_{RF} [ns]", layer);
			   Fill2DHistogram("BCAL_Global_Offsets", "Points_deltaTVsShowerEnergy", name,
							   E_shower, targetCenterTime - thisRFBunch->dTime, title,
							   1000, 0.0, 2.0, 200, -10, 10);
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
   //SortDirectories(); // Sort the histogram directories by name
   return NOERROR;
}

