// $Id$
//
//    File: JEventProcessor_ST_ZEff.cc
// Created: Wed Aug 26 17:18:47 EDT 2015
// Creator: mkamel (on Linux ifarm1102 2.6.32-431.el6.x86_64 x86_64)
//
#include "JEventProcessor_ST_ZEff.h"
#include "TRIGGER/DTrigger.h"
using namespace jana;
using namespace std;
// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_ST_ZEff());
}
} // "C"
//------------------
// JEventProcessor_ST_ZEff (Constructor)
//------------------
JEventProcessor_ST_ZEff::JEventProcessor_ST_ZEff()
{

}
//------------------
// ~JEventProcessor_ST_ZEff (Destructor)
//------------------
JEventProcessor_ST_ZEff::~JEventProcessor_ST_ZEff()
{

}
//------------------
// init
//------------------
jerror_t JEventProcessor_ST_ZEff::init(void)
{
	// ROOT mutex like this:
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();
  //****** Define Some Constants*********************************
  int NoBins_z = 240;
  double z_lower_limit = 38.5;
  double z_upper_limit = 98.5;
  // Do not reconstruct tracks with start counter time
  gPARMS->SetParameter("TRKFIT:USE_SC_TIME",false);
  int USE_SC_TIME = 0;
  if(gPARMS->Exists("TRKFIT:USE_SC_TIME"))
    gPARMS->GetParameter("TRKFIT:USE_SC_TIME", USE_SC_TIME);
    
  //cout << "USE_SC_TIME = " << USE_SC_TIME << endl;
  // Warning message if sc time is used in track reconstruction
  if (USE_SC_TIME == 0)
    {
      cout << "=========================================================================="<< endl;
      cout << "TRKFIT: USE_SC_TIME = 0; WARNING SC TIME WILL NOT BE USED IN TRACK FITTING"<< endl;
      cout << "This is required in ST_ZEff plugin which calculate SC efficiency          "<< endl;
      cout << "=========================================================================="<< endl;
    }
  // Create root folder for ST and cd to it, store main dir
  TDirectory *main = gDirectory;
  gDirectory->mkdir("st_efficiency")->cd();
  //eff histos
  h_fom =new TH1D("h_fom", " FOM ; FOM ; Counts", 1000,0.0, 1.0);
  h_N_Hit_in_track = new TH1D("h_N_Hit_in_track", " Ndof + 5; Number of Hits per track ; Counts", 70,0.0, 70.0);
  h1_qVectorSize=new TH1D(" h1_qVectorSize", " number of q tracks; number of q tracks; Counts", 50,0.0, 50.0);
  h1_qVectorSize_ACuts=new TH1D(" h1_qVectorSize_ACuts", " number of q tracks; number of q tracks; Counts", 50,0.0, 50.0);
  h1_tDiff= new TH1D("h1_tDiff", " SC_time - RF; t(ns) ; Counts", 600,-30.0, 30.0);
  h1_RFtime= new TH1D("h1_RFtime", " t0; t0(ns) ; Counts", 600,-30.0, 30.0);
  h1_Centered_RFtime= new TH1D("h1_Centered_RFtime", " RF_C; RF_C(ns) ; Counts", 600,-30.0, 30.0);
  h1_SC_ShiftedTime= new TH1D("h1_SC_ShiftedTime", " SC_shifted Time; Time(ns) ; Counts", 600,-30.0, 30.0);
  
  h_z_v = new TH1D("h_z_v", " z_v; Z (cm); Counts", 1000,0.0, 100.0);
  h1_st_pred_id = new TH1D("h1_st_prd_id", "h1_st_prd_id; Sector; Predicted Hit Counts", 31, -0.5, 30.5);

  h2_z_vs_r    = new TH2I("h2_z_vs_r", "Z vs R; Z (cm); R (cm)", NoBins_z, z_lower_limit, z_upper_limit, 100, 0.0, 10.0);
  h2_x_vs_y = new TH2I ("h2_x_vs_y", "X vs Y; X (cm); Y (cm)", 1000, -20., 20., 500, -10., 10.);
 
 
  for (unsigned int i = 0; i < Nof_ss_intervals; i++)
    {
      h_N_trck_prd_z_ss[i] = new TH1D(Form("h_N_trck_prd_z_ss_%i",i+1), Form("z #in interval %i; Sector; Predicted Hit Counts", i+1), 31, -0.5, 30.5);  
      h_N_recd_hit_z_ss[i] = new TH1D(Form("h_N_recd_hit_z_ss_%i",i+1), Form("z #in interval %i; Sector; Recorded Hit Counts", i+1), 31, -0.5, 30.5); 
      h_N_recd_hit_z_ss_ACC[i] = new TH1D(Form("h_N_recd_hit_z_ss_ACC_%i",i+1), Form("z #in interval %i; Sector; ACC Hit Counts", i+1), 31, -0.5, 30.5); 
      h_N_trck_prd_z_ss_eff[i] = new TH1D(Form("h_N_trck_prd_z_ss_eff_%i",i+1), Form("z #in interval %i; Sector; SS Effeiciency ", i+1), 31, -0.5, 30.5); 
    }
  for (unsigned int i = 0; i < Nof_bs_intervals; i++)
    {
      h_N_trck_prd_z_bs[i] = new TH1D(Form("h_N_trck_prd_z_bs_%i",i+1), Form("z #in interval %i; Sector; Predicted Hit Counts", i+1), 31, -0.5, 30.5);  
      h_N_recd_hit_z_bs[i] = new TH1D(Form("h_N_recd_hit_z_bs_%i",i+1), Form("z #in interval %i; Sector; Recorded Hit Counts", i+1), 31, -0.5, 30.5);  
      h_N_recd_hit_z_bs_ACC[i] = new TH1D(Form("h_N_recd_hit_z_bs_ACC_%i",i+1), Form("z #in interval %i; Sector; ACC Hit Counts", i+1), 31, -0.5, 30.5); 
      h_N_trck_prd_z_bs_eff[i] = new TH1D(Form("h_N_trck_prd_z_bs_eff_%i",i+1), Form("z #in interval %i; Sector; BS Effeiciency ", i+1), 31, -0.5, 30.5);
    }
  for (unsigned int i = 0; i < Nof_ns_intervals; i++)
    {
      h_N_trck_prd_z_ns[i] = new TH1D(Form("h_N_trck_prd_z_ns_%i",i+1), Form("z #in interval %i; Sector; Predicted Hit Counts", i+1), 31, -0.5, 30.5);  
      h_N_recd_hit_z_ns[i] = new TH1D(Form("h_N_recd_hit_z_ns_%i",i+1), Form("z #in interval %i; Sector; Recorded Hit Counts", i+1), 31, -0.5, 30.5); 
      h_N_recd_hit_z_ns_ACC[i] = new TH1D(Form("h_N_recd_hit_z_ns_ACC_%i",i+1), Form("z #in interval %i; Sector; ACC Hit Counts", i+1), 31, -0.5, 30.5); 
      h_N_trck_prd_z_ns_eff[i] = new TH1D(Form("h_N_trck_prd_z_ns_eff_%i",i+1), Form("z #in interval %i; Sector; NS Effeiciency ", i+1), 31, -0.5, 30.5);
    }
  
  // cd back to main directory
  gDirectory->cd("../");
  main->cd();
  return NOERROR;
}
//------------------
// brun
//------------------
jerror_t JEventProcessor_ST_ZEff::brun(JEventLoop *eventLoop, int32_t runnumber)
{
  // This is called whenever the run number changes
  // Get the particleID object for each run
  vector<const DParticleID* > dParticleID_algos;
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
  if (locGeometry->GetStartCounterGeom(sc_pos, sc_norm))
    {
      double theta = sc_norm[0][sc_norm[0].size()-2].Theta(); 
      sc_angle_corr = 1./cos(M_PI_2 - theta);
    }  

  // Propagation Time constant
  if(eventLoop->GetCalib("START_COUNTER/propagation_time_corr", propagation_time_corr))
      jout << "Error loading /START_COUNTER/propagation_time_corr !" << endl;
  // Propagation Time fit Boundaries
  if(eventLoop->GetCalib("START_COUNTER/PTC_Boundary", PTC_Boundary))
    jout << "Error loading /START_COUNTER/PTC_Boundary !" << endl;

  return NOERROR;
}
//------------------
// evnt
//------------------
jerror_t JEventProcessor_ST_ZEff::evnt(JEventLoop *eventLoop, uint64_t eventnumber)
{
  //cout << USE_SC_TIME << endl;
  //cout << " *********************** Event number" << eventnumber << "**************"<<endl;
	// Here's an example:
	//
	// vector<const MyDataClass*> mydataclasses;
	// loop->Get(mydataclasses);
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();

  vector<const DSCHit*>           st_hits;
  vector<const DChargedTrack*> chargedTrackVector;
  
  eventLoop->Get(st_hits);
  eventLoop->Get(chargedTrackVector);
  vector<DVector3> sc_track_position;
  // select events with physics events, i.e., not LED and other front panel triggers
  const DTrigger* locTrigger = NULL; 
  eventLoop->GetSingle(locTrigger); 
  if(locTrigger->Get_L1FrontPanelTriggerBits() != 0) 
    return NOERROR;
  double speed_light = 29.9792458;
   const DEventRFBunch* locEventRFBunch = NULL;
   eventLoop->GetSingle(locEventRFBunch);
  if(locEventRFBunch->dNumParticleVotes <= 1)
   return NOERROR; //don't trust PID: beta-dependence
 
  // Grab the associated detector matches object
  const DDetectorMatches* locDetectorMatches = NULL;
  eventLoop->GetSingle(locDetectorMatches);

  japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
  sc_track_position.clear();
  for (uint32_t i = 0; i < chargedTrackVector.size(); i++)
    {
      // Grab the charged track
      const DChargedTrack *thisChargedTrack = chargedTrackVector[i];
      // Grab associated time based track object by selecting charged track with best FOM
      const DTrackTimeBased *tb_track = thisChargedTrack->Get_BestTrackingFOM()->Get_TrackTimeBased();

      //const DReferenceTrajectory* rt = tb_track->rt;
      // Implement quality cuts for the time based tracks 
      float trackingFOMCut = 2.699796E-3; // +/- 3 sigma: CL
      float NHitsTrack = 14.0; // Number of hits per track 
  
      h_fom->Fill(tb_track->FOM);
      h_N_Hit_in_track->Fill(tb_track->Ndof + 5); // at lest 10 hits per track
      h1_qVectorSize->Fill(chargedTrackVector.size());

      if(tb_track->Ndof + 5 < NHitsTrack) continue;
      if(tb_track->FOM  < trackingFOMCut)  continue;
      // Grab the ST hit match params object and cut on only tracks matched to the ST
      shared_ptr<const DSCHitMatchParams> locBestSCHitMatchParams;
      bool foundSC = dParticleID->Get_BestSCMatchParams(tb_track, locDetectorMatches, locBestSCHitMatchParams);
      if (!foundSC) continue;

      // Define vertex vector
      DVector3 vertex,Momentum;
      // Vertex info
      vertex = tb_track->position();
      Momentum = tb_track->momentum();
      // Cartesian Coordinates
      double z_v = vertex.z();
      double r_v = vertex.Perp();
      bool z_vertex_cut = fabs(z_target_center - z_v) <= 15.0;
      bool r_vertex_cut = r_v < 1.0;
      
       // applied vertex cut
       if (!z_vertex_cut) continue;
       if (!r_vertex_cut) continue;

       // Grab the TOF hit match params object and cut on tracks matched to the TOF
       // Want to use TOF/RF time for t0 for SC hit time
       shared_ptr<const DTOFHitMatchParams> locTOFHitMatchParams;
       bool foundTOF = dParticleID->Get_BestTOFMatchParams(tb_track, locDetectorMatches, locTOFHitMatchParams);

       //BCAL Match Param
       shared_ptr<const DBCALShowerMatchParams> locBCALHitMatchParams;        
       bool foundBCAL = dParticleID->Get_BestBCALMatchParams(tb_track, locDetectorMatches, locBCALHitMatchParams);
       //FCAL Match Param
       shared_ptr<const DFCALShowerMatchParams> locFCALHitMatchParams;
       bool foundFCAL = dParticleID->Get_BestFCALMatchParams(tb_track, locDetectorMatches, locFCALHitMatchParams);
       // Hit must be matched to TOF, BCAL, or FCAL
       if (!(foundBCAL || (foundFCAL && foundTOF))) continue;
 

       // Grab the paramteres associated to a track matched to the ST
       vector<shared_ptr<const DSCHitMatchParams>> st_params;
       bool st_match = locDetectorMatches->Get_SCMatchParams(tb_track, st_params); 
       // If st_match = true, there is a match between this track and the ST
       if (!st_match) continue;

      DVector3 IntersectionPoint, IntersectionMomentum,locProjMom, locPaddleNorm;
      shared_ptr<DSCHitMatchParams> locSCHitMatchParams;
      vector<DTrackFitter::Extrapolation_t>extrapolations=tb_track->extrapolations.at(SYS_START);
      bool st_match_pid = dParticleID->Cut_MatchDistance(extrapolations, st_params[0]->dSCHit, st_params[0]->dSCHit->t, locSCHitMatchParams, true, &IntersectionPoint, &IntersectionMomentum);

      if(!st_match_pid) continue;
 
      int st_pred_id=locBestSCHitMatchParams->dSCHit->sector;  
      double t0 = locEventRFBunch->dTime;
     
      h1_qVectorSize_ACuts->Fill(chargedTrackVector.size());
      h1_st_pred_id->Fill(st_pred_id);
      int st_pred_id_index = st_pred_id - 1;
      // Z intersection of charged track and SC 
      // Acquire the intersection point
      sc_track_position.push_back(IntersectionPoint);
     
      h1_RFtime->Fill(t0); 
       
      double locCenteredRFTime       = t0;
      h1_Centered_RFtime->Fill(locCenteredRFTime);
      // RF time at center of target
      double locCenterToVertexRFTime = (tb_track->z() - z_target_center)*(1.0/speed_light);  // Time correction for photon from target center to vertex of track
      double locVertexRFTime         = locCenteredRFTime + locCenterToVertexRFTime;
      // Define sector array index
      int sc_index =  st_pred_id_index;
      // Start Counter geometry in hall coordinates 
      double sc_pos_soss, sc_pos_eoss, sc_pos_eobs, sc_pos_eons;
      
      sc_pos_soss = sc_pos[sc_index][0].z();   // Start of straight section
      sc_pos_eoss = sc_pos[sc_index][1].z();   // End of straight section
      sc_pos_eobs = sc_pos[sc_index][11].z();  // End of bend section
      sc_pos_eons = sc_pos[sc_index][12].z();  // End of nose section

      // Read the constants from CCDB
      incpt_ss   = propagation_time_corr[sc_index][0];
      slope_ss   = propagation_time_corr[sc_index][1];
      incpt_bs   = propagation_time_corr[sc_index][2];
      slope_bs   = propagation_time_corr[sc_index][3];
      incpt_ns   = propagation_time_corr[sc_index][4];
      slope_ns   = propagation_time_corr[sc_index][5];
      //Read fit boundary from CCDB
      Bound1 = PTC_Boundary[0][0];
      Bound2 = PTC_Boundary[1][0];
      for (uint32_t i = 0; i < sc_track_position.size(); i++)
	{
	  //double locSCzIntersection = IntersectionPoint.z();
	  h_z_v->Fill(sc_track_position[i].z());
	  h2_z_vs_r->Fill(sc_track_position[i].z(),sc_track_position[i].Perp());
	  h2_x_vs_y->Fill(sc_track_position[i].x(), sc_track_position[i].y());
	  double locSCzIntersection = sc_track_position[i].z();
	  // double locSCrIntersection =sc_track_position[i].Perp();
	  	  
	  double ss_interval = (sc_pos_eoss - sc_pos_soss)/Nof_ss_intervals;
	  double bs_interval = (sc_pos_eobs - sc_pos_eoss)/Nof_bs_intervals;
	  double ns_interval = (sc_pos_eons - sc_pos_eobs)/Nof_ns_intervals;
	  // SS intervals
	  for (uint32_t k = 0; k < Nof_ss_intervals; k++)
	    { 
	      z_ss[k] = sc_pos_soss + ss_interval * k  ;
	      //cout << z_ss[k]<<endl;
	      if ((z_ss[k] <= locSCzIntersection)  && (locSCzIntersection < (z_ss[k]+ss_interval)))
		{
		  //predicted tracks to hit SS interval
		  h_N_trck_prd_z_ss[k]->Fill(st_pred_id);
		  
		  for (uint32_t j = 0; j < st_hits.size(); j++)
		    {
		      int hit_sector       = st_hits[j]->sector;
		      if ((st_pred_id == hit_sector)||(st_pred_id == hit_sector-1) || (st_pred_id == hit_sector+1) ||(st_pred_id == hit_sector+29)||(st_pred_id == hit_sector-29))
			{
			  // Get the Flight time 
			  double FlightTime = locBestSCHitMatchParams->dFlightTime; 
			  //St time corrected for the flight time
			  double st_corr_FlightTime =  st_hits[j]->t - FlightTime;
			  			   
			  double pathlength = locSCzIntersection - sc_pos_soss;
			  double corr_time = st_corr_FlightTime  - (incpt_ss + (slope_ss *  pathlength));
			  SC_RFShiftedTime = dRFTimeFactory->Step_TimeToNearInputTime(locVertexRFTime,  corr_time);
			  //  double t_diff  =SC_RFShiftedTime - t0;
			  double t_diff  =corr_time - t0; 
			  h1_tDiff->Fill(t_diff);
			  h1_SC_ShiftedTime->Fill(SC_RFShiftedTime);
			  //SS accidentals
			  if (!((-10 < t_diff) && (t_diff <10)) && (-20 < t_diff) && (t_diff < 20))
			    {  
			      h_N_recd_hit_z_ss_ACC[k]->Fill(st_pred_id);
			    }
			  //SS recorded hits
			  if ((-10 < t_diff) && (t_diff <10))
			    { 
			      h_N_recd_hit_z_ss[k]->Fill(st_pred_id);
			      break;
			    }
			}// sector cretria 
		    }// recorded hits loop
		}// z verification in an interval
	    }// loop over ss intervals
	  //BS intervals
	  for (uint32_t p = 0; p < Nof_bs_intervals; p++)
	    { 
	      z_bs[p] = sc_pos_eoss + bs_interval * p  ;
	      //cout << z_bs[p]<<endl;
	      if ((z_bs[p] <= locSCzIntersection)  && (locSCzIntersection < (z_bs[p]+bs_interval)))
		{
		  //predicted tracks to hit BS interval
		  h_N_trck_prd_z_bs[p]->Fill(st_pred_id);
		  for (uint32_t j = 0; j < st_hits.size(); j++)
		    {
		      int hit_sector       = st_hits[j]->sector;
		      if ((st_pred_id == hit_sector)||(st_pred_id == hit_sector-1) || (st_pred_id == hit_sector+1) ||(st_pred_id == hit_sector+29)||(st_pred_id == hit_sector-29))
			{
			  //BS efficiency
			  // Get the Flight time 
			  double FlightTime = locBestSCHitMatchParams->dFlightTime; 
			  //St time corrected for the flight time
			  double st_corr_FlightTime =  st_hits[j]->t - FlightTime;
			  double SS_Length = sc_pos_eoss - sc_pos_soss;// same for along z or along the paddle
			  
			  double path_bs = SS_Length +  (locSCzIntersection - sc_pos_eoss)*sc_angle_corr;
			  double corr_time_bs=  st_corr_FlightTime  - (incpt_bs + (slope_bs *  path_bs));
			  
			  //  double t_diff  =SC_RFShiftedTime - t0;
			  double t_diff  = corr_time_bs - t0; 
			  h1_tDiff->Fill(t_diff);
			  h1_SC_ShiftedTime->Fill(SC_RFShiftedTime);
			  
			  SC_RFShiftedTime = dRFTimeFactory->Step_TimeToNearInputTime(locVertexRFTime,  Corr_Time_bs);
			       			       
			  //BS accidentals
			  if (!((-10 < t_diff) && (t_diff <10)) && (-20 < t_diff) && (t_diff < 20))
			    {  
			      h_N_recd_hit_z_bs_ACC[p]->Fill(st_pred_id);
			    }
			  //BS recorded hits
			  if ((-10 < t_diff) && (t_diff <10))
			    { 
			      h_N_recd_hit_z_bs[p]->Fill(st_pred_id);
			      break;
			    }
			} // sector cretria 
		    } // recorded hits loop
		} // z verification in an interval
	    }// loop over ss intervals
	  for (uint32_t k = 0; k < Nof_ns_intervals; k++)
	    { 
	      z_ns[k] = sc_pos_eobs + ns_interval * k  ;
	      //cout << z_ns[k]<<endl;
	      if ((z_ns[k] <= locSCzIntersection)  && (locSCzIntersection <= (z_ns[k]+ns_interval)))
		{
		  h_N_trck_prd_z_ns[k]->Fill(st_pred_id);
		  for (uint32_t j = 0; j < st_hits.size(); j++)
		    {
		      int hit_sector       = st_hits[j]->sector;
		      if ((st_pred_id == hit_sector)||(st_pred_id == hit_sector-1) || (st_pred_id == hit_sector+1) ||(st_pred_id == hit_sector+29)||(st_pred_id == hit_sector-29))
			{
			  // Get the Flight time 
			  double FlightTime = locBestSCHitMatchParams->dFlightTime; 
			  double SS_Length = sc_pos_eoss - sc_pos_soss;// same for along z or along the paddle
			  //St time corrected for the flight time
			  double st_corr_FlightTime =  st_hits[j]->t - FlightTime;
			  //NS efficiency
			  double path_ns = SS_Length +  (locSCzIntersection - sc_pos_eoss)*sc_angle_corr;
			  if (path_ns <= Bound1)
			    {       
			      Corr_Time_ns=  st_corr_FlightTime  - (incpt_ss + (slope_ss *  path_ns));
			      SC_RFShiftedTime = dRFTimeFactory->Step_TimeToNearInputTime(locVertexRFTime,  Corr_Time_ns);
			    }
			  else if ((Bound1 < path_ns)&&(path_ns <= Bound2))
			    {
			      Corr_Time_ns = st_corr_FlightTime - (incpt_bs + (slope_bs *  path_ns));
			      SC_RFShiftedTime = dRFTimeFactory->Step_TimeToNearInputTime(locVertexRFTime,  Corr_Time_ns);
			    }
			  else
			    {
			      Corr_Time_ns = st_corr_FlightTime - (incpt_ns + (slope_ns *  path_ns));
			      SC_RFShiftedTime = dRFTimeFactory->Step_TimeToNearInputTime(locVertexRFTime,  Corr_Time_ns);
			    }
			  double t_diff  =Corr_Time_ns - t0; 
			  //NS Accidentals
			  if (!((-10 < t_diff) && (t_diff <10)) && (-20 < t_diff) && (t_diff < 20))
			    {  
			      h_N_recd_hit_z_ns_ACC[k]->Fill(st_pred_id);
			    }
			  // NS recorded hits
			  if ((-10 < t_diff) && (t_diff <10))
			    { 
			      h_N_recd_hit_z_ns[k]->Fill(st_pred_id);
			      break;
			    }
			}// sector cretria 
		    } // recorded hits loop
		} // z verification in an interval
	    }// loop over ss intervals
	}//end of trcak position loop
      
    }// end of charged track loop	
  japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK
  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_ST_ZEff::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_ST_ZEff::fini(void)
{

	// Called before program exit after event processing is finished.
  japp->RootUnLock();
  return NOERROR;
}

