// $Id$
//
//    File: JEventProcessor_ST_Propagation_Time.cc
// Created: Thu Jan 14 12:16:20 EST 2016
// Creator: mkamel (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "JEventProcessor_ST_Propagation_Time.h"
#include "TRIGGER/DTrigger.h"
using namespace jana;

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_ST_Propagation_Time());
}
} // "C"

//------------------
// JEventProcessor_ST_Propagation_Time (Constructor)
//------------------
JEventProcessor_ST_Propagation_Time::JEventProcessor_ST_Propagation_Time()
{

}

//------------------
// ~JEventProcessor_ST_Propagation_Time (Destructor)
//------------------
JEventProcessor_ST_Propagation_Time::~JEventProcessor_ST_Propagation_Time()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_ST_Propagation_Time::init(void)
{
	// This is called once at program startup. If you are creating
	// and filling historgrams in this plugin, you should lock the
	// ROOT mutex like this:
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();
	//
  //****** Define Some Constants*********************************
  int NoBins_time = 200;
  int NoBins_z = 60;
  double time_lower_limit = -10.0;      
  double time_upper_limit =  10.0;
  double z_lower_limit = 0.0;
  double z_upper_limit = 60.0;
  Photonspeed = 29.9792458;
  
  // **************** define histograms *************************

  //Create root folder and cd to it, store main dir
  TDirectory *main = gDirectory;
  gDirectory->mkdir("ST_Propagation_Time")->cd();

  h2_PropTime_z_SS_chan = new TH2I*[NCHANNELS];
  h2_PropTime_z_BS_chan = new TH2I*[NCHANNELS];
  h2_PropTime_z_NS_chan = new TH2I*[NCHANNELS];
  h2_PpropTime_z = new TH2I*[NCHANNELS];  
  h2_PropTimeCorr_z_SS_chan = new TH2I*[NCHANNELS];
  h2_PropTimeCorr_z_BS_chan = new TH2I*[NCHANNELS];
  h2_PropTimeCorr_z_NS_chan = new TH2I*[NCHANNELS];
  h2_CorrectedTime_z = new TH2I*[NCHANNELS];
  for (Int_t i = 0; i < NCHANNELS; i++)
    { 
      h2_PropTime_z_SS_chan[i] = new TH2I(Form("h2_PropTime_z_SS_chan_%i", i+1), "Prop Time vs. Path length along the paddle; Path length along the paddle (cm); Propagation Time (ns)", NoBins_z,z_lower_limit,z_upper_limit, NoBins_time, time_lower_limit, time_upper_limit);
      h2_PropTime_z_BS_chan[i] = new TH2I(Form("h2_PropTime_z_BS_chan_%i", i+1), "Prop Time vs. Path length along the paddle; Path length along the paddle (cm); Propagation Time (ns)", NoBins_z,z_lower_limit,z_upper_limit, NoBins_time, time_lower_limit, time_upper_limit);
      h2_PropTime_z_NS_chan[i] = new TH2I(Form("h2_PropTime_z_NS_chan_%i", i+1), "Prop Time vs. Path length along the paddle; Path length along the paddle (cm); Propagation Time (ns)", NoBins_z,z_lower_limit,z_upper_limit, NoBins_time, time_lower_limit, time_upper_limit);
      h2_PpropTime_z[i] = new TH2I(Form("h2_PpropTime_z_%i", i+1), "Propagation Time vs. Path length along the paddle; Path length along the paddle (cm); Propagation Time (ns)", NoBins_z,z_lower_limit,z_upper_limit, NoBins_time, time_lower_limit, time_upper_limit);


      h2_PropTimeCorr_z_SS_chan[i] = new TH2I(Form("h2_PropTimeCorr_z_SS_chan_%i", i+1), "Prop Time vs. Path length along the paddle; Path length along the paddle (cm); Propagation Time (ns)", NoBins_z,z_lower_limit,z_upper_limit, NoBins_time, time_lower_limit, time_upper_limit);
      h2_PropTimeCorr_z_BS_chan[i] = new TH2I(Form("h2_PropTimeCorr_z_BS_chan_%i", i+1), "Prop Time vs. Path length along the paddle; Path length along the paddle (cm); Propagation Time (ns)", NoBins_z,z_lower_limit,z_upper_limit, NoBins_time, time_lower_limit, time_upper_limit);
      h2_PropTimeCorr_z_NS_chan[i] = new TH2I(Form("h2_PropTimeCorr_z_NS_chan_%i", i+1), "Prop Time vs. Path length along the paddle; Path length along the paddle (cm); Propagation Time (ns)", NoBins_z,z_lower_limit,z_upper_limit, NoBins_time, time_lower_limit, time_upper_limit);
      h2_CorrectedTime_z[i] = new TH2I(Form("h2_CorrectedTime_z_%i", i+1), "Corrected Time vs. Path length along the paddle; Path length along the paddle (cm); Propagation Time (ns)", NoBins_z,z_lower_limit,z_upper_limit, NoBins_time, time_lower_limit, time_upper_limit);

    }
  // cd back to main directory
  main->cd();

  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_ST_Propagation_Time::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes
  // Get the particleID object for each run
  vector<const DParticleID *> dParticleID_algos;
  eventLoop->Get(dParticleID_algos);
  if(dParticleID_algos.size() < 1)
    {
      _DBG_<<"Unable to get a DParticleID object! NO PID will be done!"<<endl;
      return RESOURCE_UNAVAILABLE;
    }
  dParticleID = dParticleID_algos[0];
  // We want to be use some of the tools available in the RFTime factory 
  // Specifically steping the RF back to a chosen time
  dRFTimeFactory = static_cast<DRFTime_factory*>(eventLoop->GetFactory("DRFTime"));
  // Be sure that DRFTime_factory::init() and brun() are called
  vector<const DRFTime*> locRFTimes;
  eventLoop->Get(locRFTimes);
 //RF Period
  vector<double> locRFPeriodVector;
  eventLoop->GetCalib("PHOTON_BEAM/RF/rf_period", locRFPeriodVector);
  dRFBunchPeriod = locRFPeriodVector[0];
  
  // Obtain the target center along z;
  map<string,double> target_params;
  if (eventLoop->GetCalib("/TARGET/target_parms", target_params))
    jout << "Error loading /TARGET/target_parms/ !" << endl;
  if (target_params.find("TARGET_Z_POSITION") != target_params.end())
    z_target_center = target_params["TARGET_Z_POSITION"];
  else
    jerr << "Unable to get TARGET_Z_POSITION from /TARGET/target_parms !" << endl;
  // Obtain the Start Counter geometry
  DApplication* dapp = dynamic_cast<DApplication*>(eventLoop->GetJApplication());
  if(!dapp)
    _DBG_<<"Cannot get DApplication from JEventLoop! (are you using a JApplication based program?)"<<endl; 
  DGeometry* locGeometry = dapp->GetDGeometry(eventLoop->GetJEvent().GetRunNumber());
  sc_angle_corr = 1.;
  if(locGeometry->GetStartCounterGeom(sc_pos, sc_norm)) {
      double theta = sc_norm[0][sc_norm[0].size()-2].Theta(); 
      sc_angle_corr = 1./cos(M_PI_2 - theta);
  }  

  // Propagation Time constant
  if(eventLoop->GetCalib("START_COUNTER/propagation_time_corr", propagation_time_corr))
    jout << "Error loading /START_COUNTER/propagation_time_corr !" << endl;
  // Propagation Time fit Boundaries
  if(eventLoop->GetCalib("START_COUNTER/PTC_Boundary", PTC_Boundary))
    jout << "Error loading /START_COUNTER/PTC_Boundary !" << endl;

  // configure parameters
  trackingFOMCut = 0.0027;  // 3 sigma cut

  return NOERROR;

}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_ST_Propagation_Time::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// loop->Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.
	// Here's an example:
	//
	// vector<const MyDataClass*> mydataclasses;
	// loop->Get(mydataclasses);
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();

  
        // select events with physics events, i.e., not LED and other front panel triggers
        const DTrigger* locTrigger = NULL; 
	loop->GetSingle(locTrigger); 
	if(locTrigger->Get_L1FrontPanelTriggerBits() != 0) 
	  return NOERROR;

  // SC hits
  vector<const DSCHit *> scHitVector;
  loop->Get(scHitVector);
  
  // RF time object
  const DRFTime* thisRFTime = NULL;
  vector <const DRFTime*> RFTimeVector;
  loop->Get(RFTimeVector);
  if (RFTimeVector.size() != 0)
    thisRFTime = RFTimeVector[0];

  // Grab charged tracks
  vector<const DChargedTrack *> chargedTrackVector;
  loop->Get(chargedTrackVector);
  
  // Grab the associated detector matches object
  const DDetectorMatches* locDetectorMatches = NULL;
  loop->GetSingle(locDetectorMatches);
  
  // Grab the associated RF bunch object
  const DEventRFBunch *thisRFBunch = NULL;
  loop->GetSingle(thisRFBunch);

  double locTOFRFShiftedTime = 9.9E+9;
  // Loop over the charged tracks
  for (uint32_t i = 0; i < chargedTrackVector.size(); i++)
    {   
      shared_ptr<const DSCHitMatchParams>  locSCHitMatchParams; 
      shared_ptr<const DTOFHitMatchParams> locTOFHitMatchParams;
        
      // Grab the charged track and declare time based track object
      const DChargedTrack   *thisChargedTrack = chargedTrackVector[i];
      // Grab associated time based track object by selecting charged track with best FOM
      const DTrackTimeBased *timeBasedTrack = thisChargedTrack->Get_BestTrackingFOM()->Get_TrackTimeBased();

      // Implement quality cuts for the time based tracks 
      if(timeBasedTrack->FOM  < trackingFOMCut) continue;
  
      // Define vertex vector and cut on target/scattering chamber geometry
      DVector3 vertex = timeBasedTrack->position();
      double z_v = vertex.z();
      double r_v = vertex.Perp();
      bool z_vertex_cut = fabs(z_target_center - z_v) <= 15.0;
      bool r_vertex_cut = r_v < 0.5;
      // Apply  vertex cut
      if (!z_vertex_cut) continue;
      if (!r_vertex_cut) continue;

      // Grab the TOF hit match params object and cut on tracks matched to the TOF
      // Want to use TOF/RF time for t0 for SC hit time
      bool foundTOF = dParticleID->Get_BestTOFMatchParams(timeBasedTrack, locDetectorMatches, locTOFHitMatchParams);
      if (!foundTOF) continue;
      
      // If there is a matched track to the SC then skip this track (avoid bias in calibration)
      bool foundSCandTOF = dParticleID->Get_BestSCMatchParams(timeBasedTrack, locDetectorMatches, locSCHitMatchParams);
      if (foundSCandTOF) continue;
      
      // Cut on the number of particle votes to find the best RF time
      if (thisRFBunch->dNumParticleVotes < 1) continue;
      // Calculate the TOF estimate of the target time
      double locTOFHitTime                 = locTOFHitMatchParams->dHitTime;    // Corrected for timewalk and propagation time 
      double locTOFTrackFlightTime         = locTOFHitMatchParams->dFlightTime;
      double locFlightTimeCorrectedTOFTime = locTOFHitTime - locTOFTrackFlightTime;
      
      // Calculate the RF estimate of the target time
      double locCenteredRFTime       = thisRFTime->dTime;                                         // RF time at center of target
      double locCenterToVertexRFTime = (timeBasedTrack->z() - z_target_center)*(1.0/Photonspeed);  // Time correction for photon from target center to vertex of track
      double locVertexRFTime         = locCenteredRFTime + locCenterToVertexRFTime;               // RF time progated to vertex time
      // Correlate the proper RF target time with the TOF "target time"
      locTOFRFShiftedTime = dRFTimeFactory->Step_TimeToNearInputTime(locVertexRFTime, locFlightTimeCorrectedTOFTime);
      if (locTOFRFShiftedTime != 9.9E+9)
      	break;	
    }
  if (locTOFRFShiftedTime != 9.9E+9)
    {
      // Loop over the charged tracks
      for (uint32_t i = 0; i < chargedTrackVector.size(); i++)
      {
          shared_ptr<DSCHitMatchParams>  locSCHitMatchParams;
          shared_ptr<const DSCHitMatchParams>  locBestSCHitMatchParams;
          //DTOFHitMatchParams locTOFHitMatchParams;

          // Grab the charged track and declare time based track object
          const DChargedTrack   *thisChargedTrack = chargedTrackVector[i];
          // Grab associated time based track object by selecting charged track with best FOM
	      const DTrackTimeBased *timeBasedTrack = thisChargedTrack->Get_BestTrackingFOM()->Get_TrackTimeBased();

	  
          // Implement quality cuts for the time based tracks 
          if(timeBasedTrack->FOM  < trackingFOMCut) continue;
          
          // Define vertex vector and cut on target/scattering chamber geometry
          DVector3 vertex = timeBasedTrack->position();
          double z_v = vertex.z();
          double r_v = vertex.Perp();
          bool z_vertex_cut = fabs(z_target_center - z_v) <= 15.0;
          bool r_vertex_cut = r_v < 0.5;
          // Apply  vertex cut
          if (!z_vertex_cut) continue;
          if (!r_vertex_cut) continue;
          // Grab the ST hit match params object and cut on tracks matched to the ST
          bool foundSC = dParticleID->Get_BestSCMatchParams(timeBasedTrack, locDetectorMatches, locBestSCHitMatchParams);
          if (!foundSC) continue;
	  
          // Define sector array index
          int sc_index = locBestSCHitMatchParams->dSCHit->sector - 1;
          // Start Counter geometry in hall coordinates 
          double sc_pos_soss = sc_pos[sc_index][0].z();   // Start of straight section
          double sc_pos_eoss = sc_pos[sc_index][1].z();   // End of straight section
          double sc_pos_eobs = sc_pos[sc_index][11].z();  // End of bend section
          double sc_pos_eons = sc_pos[sc_index][12].z();  // End of nose section
          
          vector<shared_ptr<const DSCHitMatchParams>> st_params;
          bool sc_match = locDetectorMatches->Get_SCMatchParams(timeBasedTrack, st_params); 
          // If st_match = true, there is a match between this track and the ST
          if (!sc_match) continue;
          DVector3 IntersectionPoint, IntersectionMomentum;
	  vector<DTrackFitter::Extrapolation_t>extrapolations=timeBasedTrack->extrapolations.at(SYS_START);
          bool sc_match_pid = dParticleID->Cut_MatchDistance(extrapolations, st_params[0]->dSCHit, st_params[0]->dSCHit->t, locSCHitMatchParams, true, &IntersectionPoint, &IntersectionMomentum);

          if(!sc_match_pid) continue;
          // For each paddle calculate the hit time, flight time, intersection point (z), and t0 from TOF
          // For each hit we want to calculate thit - tflight - t0 from TOF
          // Then correlate with hit position along z
	  
          double locSCHitTime                 = st_params[0]->dSCHit->t; //Only corrected for time walk
          double locSCTrackFlightTime         = locSCHitMatchParams->dFlightTime;  
          double locFlightTimeCorrectedSCTime = locSCHitTime - locSCTrackFlightTime;
	  
          double locSCPropTime = locFlightTimeCorrectedSCTime - locTOFRFShiftedTime;
          // Z intersection of charged track and SC 
          double locSCzIntersection = IntersectionPoint.z();
          // Calculate the path along the paddle
	  double path_ss=0.,path_bs=0.,path_ns=0.; 
          double SS_Length = sc_pos_eoss - sc_pos_soss;// same for along z or along the paddle
          
          /////////////////////////////////////////////
          // Fill the histograms before corrections////
          ////////////////////////////////////////////
          // FILL HISTOGRAMS
          // Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
          japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
          
          // Straight Sections
          if (locSCzIntersection > sc_pos_soss && locSCzIntersection <= sc_pos_eoss)
          {
              path_ss = IntersectionPoint.z() - sc_pos_soss;
              h2_PropTime_z_SS_chan[sc_index]->Fill(path_ss,locSCPropTime);
              h2_PpropTime_z[sc_index]->Fill(path_ss,locSCPropTime);
          }
          // Bend Sections
          if(locSCzIntersection > sc_pos_eoss && locSCzIntersection <= sc_pos_eobs)
          {
	    path_bs = SS_Length +  (locSCzIntersection - sc_pos_eoss)*sc_angle_corr;;
              h2_PropTime_z_BS_chan[sc_index]->Fill(path_bs,locSCPropTime);
              h2_PpropTime_z[sc_index]->Fill(path_bs,locSCPropTime);
          }
          // Nose Sections
          if(locSCzIntersection > sc_pos_eobs && locSCzIntersection <= sc_pos_eons) 
          {
	    path_ns = SS_Length +  (locSCzIntersection - sc_pos_eoss)*sc_angle_corr;;
              h2_PropTime_z_NS_chan[sc_index]->Fill(path_ns,locSCPropTime);
              h2_PpropTime_z[sc_index]->Fill(path_ns,locSCPropTime);
          }
          ///////////////////////////////////////////////////////////////////
          /// Fill the propagation time corrected histograms////////////////
          /////////////////////////////////////////////////////////////////
          // Read the constants from CCDB
          double incpt_ss   = propagation_time_corr[sc_index][0];
          double slope_ss   = propagation_time_corr[sc_index][1];
          double incpt_bs   = propagation_time_corr[sc_index][2];
          double slope_bs   = propagation_time_corr[sc_index][3];
          double incpt_ns   = propagation_time_corr[sc_index][4];
          double slope_ns   = propagation_time_corr[sc_index][5];
	  //Read fit boundary from CCDB
	  double Bound1 = PTC_Boundary[0][0];
	  double Bound2 = PTC_Boundary[1][0];
          // Straight Sections
          if (sc_pos_soss < locSCzIntersection  && locSCzIntersection <= sc_pos_eoss)
          {
              double Corr_Time_ss = locSCPropTime - (incpt_ss + (slope_ss *  path_ss));
              h2_PropTimeCorr_z_SS_chan[sc_index]->Fill(path_ss,Corr_Time_ss);
              h2_CorrectedTime_z[sc_index]->Fill(path_ss,Corr_Time_ss);
          }
          // Bend Sections:
	  //In new calibration we are using SS constants up to 44 cm.
          if( sc_pos_eoss < locSCzIntersection  && locSCzIntersection <= sc_pos_eobs)
          {
              double Corr_Time_bs = locSCPropTime - (incpt_ss + (slope_ss *  path_bs));
              h2_PropTimeCorr_z_BS_chan[sc_index]->Fill(path_bs,Corr_Time_bs);
              h2_CorrectedTime_z[sc_index]->Fill(path_bs,Corr_Time_bs);
          }
          // Nose Sections
          if(sc_pos_eobs < locSCzIntersection   && locSCzIntersection <= sc_pos_eons)
          { 

	    if (path_ns <= Bound1)
	      {        
		double Corr_Time_ns = locSCPropTime - (incpt_ss + (slope_ss *  path_ns));
              h2_PropTimeCorr_z_NS_chan[sc_index]->Fill(path_ns,Corr_Time_ns);
              h2_CorrectedTime_z[sc_index]->Fill(path_ns,Corr_Time_ns);
	      }

	    else if ((Bound1 < path_ns)&&(path_ns <= Bound2))
	      {
		double Corr_Time_ns = locSCPropTime - (incpt_bs + (slope_bs *  path_ns));
		h2_PropTimeCorr_z_NS_chan[sc_index]->Fill(path_ns,Corr_Time_ns);
		h2_CorrectedTime_z[sc_index]->Fill(path_ns,Corr_Time_ns);
	      }
	    else
	      {
		double Corr_Time_ns = locSCPropTime - (incpt_ns + (slope_ns *  path_ns));
		h2_PropTimeCorr_z_NS_chan[sc_index]->Fill(path_ns,Corr_Time_ns);
		h2_CorrectedTime_z[sc_index]->Fill(path_ns,Corr_Time_ns);
	      }
          }
          japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK
      } // sc charged tracks
    }// TOF reference time


  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_ST_Propagation_Time::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_ST_Propagation_Time::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}

