// $Id$
//
//    File: JEventProcessor_ST_online_efficiency.cc
// Created: Wed Jan 20 10:35:58 EST 2016
// Creator: mkamel (on Linux ifarm1102 2.6.32-431.el6.x86_64 x86_64)
//

#include "JEventProcessor_ST_online_efficiency.h"
#include "TRIGGER/DTrigger.h"
using namespace jana;
using namespace std;

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_ST_online_efficiency());
}
} // "C"


//------------------
// JEventProcessor_ST_online_efficiency (Constructor)
//------------------
JEventProcessor_ST_online_efficiency::JEventProcessor_ST_online_efficiency()
{

}

//------------------
// ~JEventProcessor_ST_online_efficiency (Destructor)
//------------------
JEventProcessor_ST_online_efficiency::~JEventProcessor_ST_online_efficiency()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_ST_online_efficiency::init(void)
{
	// This is called once at program startup. If you are creating
	// and filling historgrams in this plugin, you should lock the
	// ROOT mutex like this:
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();
	//
  // USE_SC_TIME = 0 in order to Not reconstruct tracks with start counter time
  gPARMS->SetParameter("TRKFIT:USE_SC_TIME",false);
  
  // Create root folder for ST and cd to it, store main dir
  TDirectory *main = gDirectory;
  gDirectory->mkdir("st_efficiency")->cd();
  //eff histos
  h_N_trck_prd_All = new TH1D("h_N_trck_prd_All", "h_N_trck_prd_All; Sector; Predicted Hit Counts", 31, -0.5, 30.5);
  h_N_recd_hit_All = new TH1D("h_N_recd_hit_All", "h_N_recd_hit_All; Sector; Recorded Hit Counts", 31, -0.5, 30.5);
  
  h_N_trck_prd_ss = new TH1D("h_N_trck_prd_ss", "h_N_trck_prd_ss; Sector; Predicted Hit Counts", 31, -0.5, 30.5);
  h_N_recd_hit_ss = new TH1D("h_N_recd_hit_ss", "h_N_recd_hit_ss; Sector; Recorded Hit Counts", 31, -0.5, 30.5);

  h_N_trck_prd_bs = new TH1D("h_N_trck_prd_bs", "h_N_trck_prd_bs; Sector; Predicted Hit Counts", 31, -0.5, 30.5);
  h_N_recd_hit_bs = new TH1D("h_N_recd_hit_bs", "h_N_recd_hit_bs; Sector; Recorded Hit Counts", 31, -0.5, 30.5);

  h_N_trck_prd_ns = new TH1D("h_N_trck_prd_ns", "h_N_trck_prd_ns; Sector; Predicted Hit Counts", 31, -0.5, 30.5);
  h_N_recd_hit_ns = new TH1D("h_N_recd_hit_ns", "h_N_recd_hit_ns; Sector; Recorded Hit Counts", 31, -0.5, 30.5);

  h_ST_Eff_All= new TH1D("h_ST_Eff_All", " Efficiency; Sector; N_{RECD}/N_{TRCK}", 31, -0.5, 30.5);
  h_ST_Eff_ss = new TH1D("h_ST_Eff_ss", " SS Efficiency; Sector; N_{RECD}/N_{TRCK}", 31, -0.5, 30.5);
  h_ST_Eff_bs = new TH1D("h_ST_Eff_bs", " BS Efficiency; Sector; N_{RECD}/N_{TRCK}", 31, -0.5, 30.5);
  h_ST_Eff_ns = new TH1D("h_ST_Eff_ns", " NS Efficiency; Sector; N_{RECD}/N_{TRCK}", 31, -0.5, 30.5);
 
  gDirectory->cd("../");
  main->cd();

  // Initialize counters
  memset(N_trck_prd_All, 0, sizeof(N_trck_prd_All));
  memset(N_recd_hit_All, 0, sizeof(N_recd_hit_All));
  memset(N_trck_prd_ss, 0, sizeof(N_trck_prd_ss));
  memset(N_recd_hit_ss, 0, sizeof(N_recd_hit_ss));
  memset(N_trck_prd_bs, 0, sizeof(N_trck_prd_bs));
  memset(N_recd_hit_bs, 0, sizeof(N_recd_hit_bs));
  memset(N_trck_prd_ns, 0, sizeof(N_trck_prd_ns));
  memset(N_recd_hit_ns, 0, sizeof(N_recd_hit_ns));
  
  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_ST_online_efficiency::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes
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
  locGeometry->GetStartCounterGeom(sc_pos, sc_norm);
  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_ST_online_efficiency::evnt(JEventLoop *eventLoop, uint64_t eventnumber)
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

  const DTrigger* locTrigger = NULL; 
  eventLoop->GetSingle(locTrigger); 
  if(locTrigger->Get_L1FrontPanelTriggerBits() != 0)
    return NOERROR;

  vector<const DSCDigiHit*>       st_adc_digi_hits;
  vector<const DParticleID*>      pid_algorithm;
  vector<const DSCHit*>           st_hits;
  vector<const DChargedTrack*> chargedTrackVector;
  eventLoop->Get(st_adc_digi_hits);
  eventLoop->Get(pid_algorithm);
  eventLoop->Get(st_hits);
  eventLoop->Get(chargedTrackVector);
  // Grab the associated detector matches object
  const DDetectorMatches* locDetectorMatches = NULL;
  eventLoop->GetSingle(locDetectorMatches);
  
	// FILL HISTOGRAMS
	// Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
  
  // Loop over charged tracks
  for (uint32_t i = 0; i < chargedTrackVector.size(); i++)
    {
      // Grab the charged track
      const DChargedTrack *thisChargedTrack = chargedTrackVector[i];
      // Declare the time based track object
      const DTrackTimeBased *timeBasedTrack;
      // Grab associated time based track object by selecting charged track with best FOM
      thisChargedTrack->Get_BestTrackingFOM()->GetSingle(timeBasedTrack);
      float trackingFOMCut = 0.0027;  // 3 sigma cut
      if(timeBasedTrack->FOM  < trackingFOMCut)  continue;
      // Define vertex vector
      DVector3 vertex;
      // Vertex info
      vertex = timeBasedTrack->position();
      // Cartesian Coordinates
      double z_v = vertex.z();
      double r_v = vertex.Perp();
      bool z_vertex_cut = fabs(z_target_center - z_v) <= 15.0;
      bool r_vertex_cut = r_v < 0.3;
      // applied vertex cut
      if (!z_vertex_cut) continue;
      if (!r_vertex_cut) continue;
      int st_pred_id = pid_algorithm[0]->PredictSCSector(timeBasedTrack->rt,0.05235,&locProjPos,&Barrel);
      int st_pred_id_index = st_pred_id - 1;
      // Z intersection of charged track and SC 
      locSCzIntersection = locProjPos.z();
      locSCrIntersection = locProjPos.Perp();
      // Get the direction of the track momentum
      // Momentum of the track
      DVector3 Momentum;
      Momentum = timeBasedTrack->momentum();
      if (st_pred_id != 0) 
	{ 
	  // Define sector array index
	  sc_index =  st_pred_id_index;
	  // Start Counter geometry in hall coordinates 
	  sc_pos_soss = sc_pos[sc_index][0].z();   // Start of straight section
	  sc_pos_eoss = sc_pos[sc_index][1].z();   // End of straight section
	  sc_pos_eobs = sc_pos[sc_index][11].z();  // End of bend section
	  sc_pos_eons = sc_pos[sc_index][12].z();  // End of nose section
	  ss_interval = (sc_pos_eoss - sc_pos_soss)/Nof_ss_intervals;
	  bs_interval = (sc_pos_eobs - sc_pos_eoss)/Nof_bs_intervals;
	  ns_interval = (sc_pos_eons - sc_pos_eobs)/Nof_ns_intervals;
	  
	  //****************************
	  // Efficiency for the whole ST
	  //****************************
	  N_trck_prd_All[st_pred_id_index] += 1;
	  //loop over the Hit object and get the real sector hit at ST
	  for (uint32_t j = 0; j < st_hits.size(); j++)
	    {
	      int phi_sec_hit_sector       = st_hits[j]->sector;
	      //********************************************************
	      //total efficiency of the ST
	      //********************************************************
	      if (st_pred_id == phi_sec_hit_sector){
		N_recd_hit_All[st_pred_id_index] += 1;
		break;
	      }
	    }
	  // ********************************
	  // Efficiency for Straight Section
	  // ********************************
	  if ( sc_pos_soss <= locSCzIntersection  && locSCzIntersection <= sc_pos_eoss)
	    {
	      N_trck_prd_ss[st_pred_id_index] += 1;
	      //loop over the Hit object and get the real sector hit at ST
	      for (uint32_t j = 0; j < st_hits.size(); j++)
		{
		  int phi_sec_hit_sector       = st_hits[j]->sector;
		  if (st_pred_id == phi_sec_hit_sector)
		    {
		      N_recd_hit_ss[st_pred_id_index] += 1;
		      break;
		    }
		}
	    }// end of straight section loop
	  //********************************
	  // Efficiency for Bend Section
	  //********************************
	  if ( sc_pos_eoss < locSCzIntersection  && locSCzIntersection <= sc_pos_eobs)
	    {
	      N_trck_prd_bs[st_pred_id_index] += 1;
	      //loop over the Hit object and get the real sector hit at ST
	      for (uint32_t j = 0; j < st_hits.size(); j++)
		{
		  int phi_sec_hit_sector       = st_hits[j]->sector;
		  if (st_pred_id == phi_sec_hit_sector)
		    {
		      N_recd_hit_bs[st_pred_id_index] += 1;
		      break;
		    }
		}
	    }// end of bend section loop
	  //********************************
	  // Efficiency for Nose Section
	  //********************************
	  if ( sc_pos_eobs < locSCzIntersection  && locSCzIntersection <= sc_pos_eons)
	    {
	      N_trck_prd_ns[st_pred_id_index] += 1;
	      
	      //loop over the Hit object and get the real sector hit at ST
	      for (uint32_t j = 0; j < st_hits.size(); j++)
		{
		  int phi_sec_hit_sector       = st_hits[j]->sector;
		  if (st_pred_id == phi_sec_hit_sector)
		    {
		      N_recd_hit_ns[st_pred_id_index] += 1;
		      break;
		    }
		}
	    }// end of nose section loop
	} // end if (st_pred_id != 0) 
    }// end of charged track loop

	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_ST_online_efficiency::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_ST_online_efficiency::fini(void)
{
	// Called before program exit after event processing is finished.
  for (uint32_t i = 0; i < NCHANNELS; i++)
    {
      //hit object
      h_N_trck_prd_All->Fill(i+1,double(N_trck_prd_All[i]));
      h_N_recd_hit_All->Fill(i+1,double(N_recd_hit_All[i]));
      h_N_trck_prd_ss->Fill(i+1,double(N_trck_prd_ss[i]));
      h_N_recd_hit_ss->Fill(i+1,double(N_recd_hit_ss[i]));
      h_N_trck_prd_bs->Fill(i+1,double(N_trck_prd_bs[i]));
      h_N_recd_hit_bs->Fill(i+1,double(N_recd_hit_bs[i]));
      h_N_trck_prd_ns->Fill(i+1,double(N_trck_prd_ns[i]));
      h_N_recd_hit_ns->Fill(i+1,double(N_recd_hit_ns[i]));
    }
  return NOERROR;
}

