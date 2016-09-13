// $Id$
//
//    File: JEventProcessor_ST_Tresolution.cc
// Created: Fri Jan  8 09:07:34 EST 2016
// Creator: mkamel (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "JEventProcessor_ST_Tresolution.h"
using namespace jana;


// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_ST_Tresolution());
}
} // "C"


//------------------
// JEventProcessor_ST_Tresolution (Constructor)
//------------------
JEventProcessor_ST_Tresolution::JEventProcessor_ST_Tresolution()
{

}

//------------------
// ~JEventProcessor_ST_Tresolution (Destructor)
//------------------
JEventProcessor_ST_Tresolution::~JEventProcessor_ST_Tresolution()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_ST_Tresolution::init(void)
{
	// This is called once at program startup. If you are creating
	// and filling historgrams in this plugin, you should lock the
	// ROOT mutex like this:
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();
	//
  // **************** define histograms *************************

  //Create root folder and cd to it, store main dir
  TDirectory *main = gDirectory;
  gDirectory->mkdir("ST_Tresolution")->cd();

  h2_CorrectedTime_z = new TH2I*[NCHANNELS];
  // All my Calculations in 2015 were using the binning below
  int NoBins_time = 100;
  int NoBins_z = 300;
  double time_lower_limit = -5.0;      
  double time_upper_limit = 5.0;
  double z_lower_limit = 0.0;
  double z_upper_limit = 60.0;
  for (Int_t i = 0; i < NCHANNELS; i++)
    { 
      h2_CorrectedTime_z[i] = new TH2I(Form("h2_CorrectedTime_z_%i", i+1), "Corrected Time vs. Z; Z (cm); Propagation Time (ns)", NoBins_z,z_lower_limit,z_upper_limit, NoBins_time, time_lower_limit, time_upper_limit);
    }

  // cd back to main directory
  main->cd();
  
  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_ST_Tresolution::brun(JEventLoop *eventLoop, int32_t runnumber)
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

  // set some parameters
  trackingFOMCut = 0.0027;  // 3 sigma cut

  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_ST_Tresolution::evnt(JEventLoop *loop, uint64_t eventnumber)
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
  double speed_light = 29.9792458;
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

  for (uint32_t i = 0; i < chargedTrackVector.size(); i++)
    {   
      // Grab the charged track and declare time based track object
      const DChargedTrack   *thisChargedTrack = chargedTrackVector[i];
      const DTrackTimeBased *timeBasedTrack;
      // Grab associated time based track object by selecting charged track with best FOM
      thisChargedTrack->Get_BestTrackingFOM()->GetSingle(timeBasedTrack);
      // Implement quality cuts for the time based tracks 
      //trackingFOMCut = 0.0027;  // 3 sigma cut
      //trackingFOMCut = 0.0001;  // 5 sigma cut
      if(timeBasedTrack->FOM  < trackingFOMCut) continue;

      // Grab the ST hit match params object and cut on only tracks matched to the ST
      DSCHitMatchParams locSCHitMatchParams;
      bool foundSC = dParticleID->Get_BestSCMatchParams(timeBasedTrack, locDetectorMatches, locSCHitMatchParams);
      if (!foundSC) continue;
      
      // Define vertex vector and cut on target/scattering chamber geometry
      DVector3 vertex = timeBasedTrack->position();
      double z_v = vertex.z();
      double r_v = vertex.Perp();
      
      bool z_vertex_cut = fabs(z_target_center - z_v) <= 15.0;
      bool r_vertex_cut = r_v < 0.5;
      // Apply  vertex cut
      if (!z_vertex_cut) continue;
      if (!r_vertex_cut) continue;
      vector<DSCHitMatchParams> st_params;
      bool st_match = locDetectorMatches->Get_SCMatchParams(timeBasedTrack, st_params); 
      // If st_match = true, there is a match between this track and the ST
      if (!st_match) continue;

      DVector3 IntersectionPoint;
      DVector3 IntersectionDir;
      bool sc_match_pid = dParticleID->MatchToSC(timeBasedTrack, 
                                                 timeBasedTrack->rt, 
                                                 st_params[0].dSCHit, 
                                                 st_params[0].dSCHit->t, 
                                                 locSCHitMatchParams, 
                                                 true, NULL,
                                                 &IntersectionPoint, &IntersectionDir); 
      if(!sc_match_pid) continue; 
      // Cut on the number of particle votes to find the best RF time
      if (thisRFBunch->dNumParticleVotes < 2) continue;
      // Calculate the TOF estimate of the target time

      // Calculate the RF estimate of the target time
      double locCenteredRFTime       = thisRFTime->dTime;
      // RF time at center of target
      double locCenterToVertexRFTime = (timeBasedTrack->z() - z_target_center)*(1.0/speed_light);  // Time correction for photon from target center to vertex of track
      double locVertexRFTime         = locCenteredRFTime + locCenterToVertexRFTime;
      int sc_index= locSCHitMatchParams.dSCHit->sector - 1;
      // Start Counter geometry in hall coordinates 
      double sc_pos_soss = sc_pos[sc_index][0].z();   // Start of straight section
      double sc_pos_eoss = sc_pos[sc_index][1].z();   // End of straight section
      double sc_pos_eobs = sc_pos[sc_index][11].z();  // End of bend section
      double sc_pos_eons = sc_pos[sc_index][12].z();  // End of nose section
      //Get the ST time walk corrected time
      double st_time = st_params[0].dSCHit->t;
      // Get the Flight time 
      double FlightTime = locSCHitMatchParams.dFlightTime; 
      //St time corrected for the flight time
      double st_corr_FlightTime =  st_time - FlightTime;
      // SC_RFShiftedTime = dRFTimeFactory->Step_TimeToNearInputTime(locVertexRFTime,  st_corr_FlightTime);
      // Z intersection of charged track and SC 
      double locSCzIntersection = IntersectionPoint.z();
      ////////////////////////////////////////////////////////////////////
      /// Fill the sc time histograms corrected for walk and propagation//
      ////////////////////////////////////////////////////////////////////
      // Read the constants from CCDB
      double incpt_ss   = propagation_time_corr[sc_index][0];
      double slope_ss   = propagation_time_corr[sc_index][1];
      double incpt_bs   = propagation_time_corr[sc_index][2];
      double slope_bs   = propagation_time_corr[sc_index][3];
      double incpt_ns   = propagation_time_corr[sc_index][4];
      double slope_ns   = propagation_time_corr[sc_index][5];
      ///////////////////////////////////////
      //Calculate the path along the paddle
      /////////////////////////////////////
      //Define some parameters
      //double Radius = 12.0;
      //double theta  = 18.5 * pi/180.0;
      double SS_Length = sc_pos_eoss - sc_pos_soss;// same for along z or along the paddle
      //double BS_Length = Radius *  theta ; // along the paddle
      //double NS_Length = (sc_pos_eons - sc_pos_eobs)/cos(theta);// along the paddle
      
      // FILL HISTOGRAMS
      // Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
      japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

      // Straight Sections
      if (sc_pos_soss < locSCzIntersection && locSCzIntersection <= sc_pos_eoss)
      {
          double path_ss = locSCzIntersection - sc_pos_soss;
          double Corr_Time_ss = st_corr_FlightTime  - (incpt_ss + (slope_ss *  path_ss));
          double SC_RFShiftedTime = dRFTimeFactory->Step_TimeToNearInputTime(locVertexRFTime,  Corr_Time_ss);
          h2_CorrectedTime_z[sc_index]->Fill(path_ss, Corr_Time_ss -SC_RFShiftedTime);
      }
      // Bend Sections
      if(sc_pos_eoss < locSCzIntersection && locSCzIntersection <= sc_pos_eobs)
      {
          //double path_bs = SS_Length + Radius * asin((locSCzIntersection - sc_pos_eoss)/Radius);
          double path_bs = SS_Length +  (locSCzIntersection - sc_pos_eoss)*sc_angle_corr;
          double Corr_Time_bs =  st_corr_FlightTime  - (incpt_bs + (slope_bs *  path_bs));
          double SC_RFShiftedTime = dRFTimeFactory->Step_TimeToNearInputTime(locVertexRFTime,  Corr_Time_bs);
          h2_CorrectedTime_z[sc_index]->Fill(path_bs,Corr_Time_bs - SC_RFShiftedTime);
      }
      // Nose Sections
      if(sc_pos_eobs < locSCzIntersection && locSCzIntersection <= sc_pos_eons)
      { 
          //double path_ns = SS_Length + BS_Length +((locSCzIntersection - sc_pos_eobs)/cos(theta));
          double path_ns = SS_Length +  (locSCzIntersection - sc_pos_eoss)*sc_angle_corr;
          double Corr_Time_ns =  st_corr_FlightTime  - (incpt_ns + (slope_ns *  path_ns));
          double SC_RFShiftedTime = dRFTimeFactory->Step_TimeToNearInputTime(locVertexRFTime,  Corr_Time_ns);
          h2_CorrectedTime_z[sc_index]->Fill(path_ns,Corr_Time_ns - SC_RFShiftedTime);
      }
      japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

    } // sc charged tracks


  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_ST_Tresolution::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_ST_Tresolution::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}

