// $Id$
//
//    File: JEventProcessor_ST_online_tracking.cc
// Created: Fri Jun 19 13:22:21 EDT 2015
// Creator: mkamel (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//
// $Id$
#include "JEventProcessor_ST_online_tracking.h"
using namespace jana;
using namespace std;
// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
// ROOT header files
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
// ST header files
#include "START_COUNTER/DSCHit.h"
#include "START_COUNTER/DSCDigiHit.h"
// C++ header files
#include <stdint.h>
#include <vector>
#include <stdio.h>
// Include some JANA libraries
// PID libraries
#include "PID/DDetectorMatches.h"
#include "PID/DParticleID.h"
#include "PID/DChargedTrack.h"
// Tracking libraries
#include "TRACKING/DTrackTimeBased.h"

// Define some constants
const uint32_t  RAD2DEG       = 57.29577951;      // Convert radians to degrees
const uint32_t  NCHANNELS     = 30;              // number of scintillator paddles

// Declare 2D tracking histos
static TH2I *h2_r_vs_z; 
static TH2I *h2_y_vs_x;
static TH2I *h2_phi_vs_sector;
static TH2I *h2_dphi_sector;

static TH2I *h2_dedx_P_mag;
static TH2I *h2_dedx_P_mag_postv;
static TH2I *h2_dedx_P_mag_negtv;
// Declare 1D tracking histos
static TH1D *h_phi_sec_hit_cntr;
static TH1D *h_phi_sec_pred_hit_cntr;
static TH1D *pEff;
static TH1D *MacropEff;
static TH1D *h_phi_sec_adc_cntr;
static TH1D *pEff_adc;
static TH1D *MacropEff_adc;
// Declare detection efficency counters
uint32_t phi_sec_cntr[NCHANNELS];
uint32_t phi_sec_pred_hit_cntr[NCHANNELS];
uint32_t phi_sec_adc_hit_cntr[NCHANNELS];
uint32_t phi_sec_hit_cntr[NCHANNELS];
double   phi_sec[NCHANNELS][2]; 


extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_ST_online_tracking());
}
} // "C"


//------------------
// JEventProcessor_ST_online_tracking (Constructor)
//------------------
JEventProcessor_ST_online_tracking::JEventProcessor_ST_online_tracking()
{

}

//------------------
// ~JEventProcessor_ST_online_tracking (Destructor)
//------------------
JEventProcessor_ST_online_tracking::~JEventProcessor_ST_online_tracking()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_ST_online_tracking::init(void)
{
	// This is called once at program startup. If you are creating
	// and filling historgrams in this plugin, you should lock the
	// ROOT mutex like this:
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();
	//
  japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
  // Create root folder for ST and cd to it, store main dir
  TDirectory *main = gDirectory;
  gDirectory->mkdir("st_tracking")->cd();
  // Book some 2D histos
  h2_r_vs_z    = new TH2I("r_vs_z", "Charged Track & ST Intersection; Z (cm); R (cm)", 280, 30, 100, 500, 0, 10);
  h2_y_vs_x    = new TH2I("y_vs_x", "Charged Track & ST Intersection; X (cm); Y (cm)", 2000, -20, 20, 1000, -10, 10);
  
  h2_phi_vs_sector = new TH2I("phi_vs_sector", "Charged Track & ST Intersection; Sector Number; #phi (deg)", NCHANNELS, 0.5, NCHANNELS + 0.5, 360, 0, 360);
  h2_dphi_sector  = new TH2I("h2_dphi_sector", "#delta #phi vs Sector; Sector Number; #delta #phi (deg)", NCHANNELS, 0.5, NCHANNELS + 0.5, 600, -15, 15);
 
  h2_dedx_P_mag      = new TH2I("h2_dedx_P_mag", "#frac{dE}{dx} vs Momentum; Momentum (Gev/c); dEdx (au)", 1000, 0, 2, 150, 0, .015);
  h2_dedx_P_mag_postv= new TH2I("h2_dedx_P_mag_postv", "#frac{dE}{dx} vs Momentum (q > 0); Momentum (Gev/c); #frac{dE}{dx} (au)", 1000, 0, 2, 150, 0, .015);
  h2_dedx_P_mag_negtv= new TH2I("h2_dedx_P_mag_negtv", "#frac{dE}{dx} vs Momentum (q < 0); Momentum (Gev/c); #frac{dE}{dx} (au)", 1000, 0, 2, 150, 0, .015);

  //eff histos
  h_phi_sec_pred_hit_cntr = new TH1D("h_phi_sec_pred_hit_cntr", "phi_sec_pred_hit_cntr; Sector; Predicted Hit Counts", 31, -0.5, 30.5);
  h_phi_sec_hit_cntr      = new TH1D("h_phi_sec_hit_cntr", "phi_sec_hit_cntr; Sector; Hit Counts", 31, -0.5, 30.5);
  pEff= new TH1D("pEff", "Hit Efficiency; Sector; N_{HIT}/N_{TRK}", 31, -0.5, 30.5);
  MacropEff= new TH1D("MacropEff", "Hit Efficiency; Sector; N_{HIT}/N_{TRK}", 31, -0.5, 30.5);
  h_phi_sec_adc_cntr = new TH1D("h_phi_sec_adc_cntr", "phi_sec_adc_cntr; Sector; adc Counts", 31, -0.5, 30.5);
  pEff_adc = new TH1D("pEff_adc", "ADC Efficiency; Sector; N_{ADC}/N_{TRK}", 31, -0.5, 30.5);
  MacropEff_adc = new TH1D("MacropEff_adc", "ADC Efficiency; Sector; N_{ADC}/N_{TRK}", 31, -0.5, 30.5);
  // cd back to main directory
  gDirectory->cd("../");
  main->cd();
  japp->RootUnLock(); //RELEASE ROOT LOCK!!
  // Initialize counters
  memset(phi_sec_pred_hit_cntr, 0, sizeof(phi_sec_pred_hit_cntr));
  memset(phi_sec_hit_cntr, 0, sizeof(phi_sec_hit_cntr));
  memset(phi_sec_adc_hit_cntr, 0, sizeof(phi_sec_adc_hit_cntr));
  for (uint32_t i = 0; i < NCHANNELS; i++)
    {
      phi_sec[i][0] = (6.0 + 12.0*double(i)) - 4.0;
      phi_sec[i][1] = (6.0 + 12.0*double(i)) + 4.0;
    }
  return NOERROR;
}
//------------------
// brun
//------------------
jerror_t JEventProcessor_ST_online_tracking::brun(JEventLoop *eventLoop, int runnumber)
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
  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_ST_online_tracking::evnt(JEventLoop *eventLoop, int eventnumber)
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
  vector<const DSCDigiHit*>       st_adc_digi_hits;
  vector<const DParticleID*>      pid_algorithm;
  vector<const DSCHit*>           st_hits;
  eventLoop->Get(st_adc_digi_hits);
  eventLoop->Get(pid_algorithm);
  eventLoop->Get(st_hits);
  vector<DVector3> sc_track_position;
  vector<const DChargedTrack*> chargedTrackVector;
  eventLoop->Get(chargedTrackVector);
  // Grab the associated detector matches object
  const DDetectorMatches* locDetectorMatches = NULL;
  eventLoop->GetSingle(locDetectorMatches);
  // Lock ROOT mutex so other threads won't interfere 
  japp->RootWriteLock();

  sc_track_position.clear();
 
// Loop over charged tracks
  for (uint32_t i = 0; i < chargedTrackVector.size(); i++)
    {
      // Grab the charged track
      const DChargedTrack *thisChargedTrack = chargedTrackVector[i];
      
      // Declare the time based track object
      const DTrackTimeBased *timeBasedTrack;

      // Grab associated time based track object by selecting charged track with best FOM
      thisChargedTrack->Get_BestTrackingFOM()->GetSingle(timeBasedTrack);
      // Implement quality cuts for the time based tracks 
      float trackingFOMCut = 0.0027;  // 3 sigma cut
      if(timeBasedTrack->FOM  < trackingFOMCut)  continue;
      // Get the charge of the track and cut on charged tracks
      int q = timeBasedTrack->charge();
      // Grab the ST hit match params object and cut on only tracks matched to the ST
      DSCHitMatchParams locSCHitMatchParams;
      bool foundSC = dParticleID->Get_BestSCMatchParams(timeBasedTrack, locDetectorMatches, locSCHitMatchParams);
      if (!foundSC) continue;
      // Define vertex vector
      DVector3 vertex;
      // Vertex info
      vertex = timeBasedTrack->position();
      // Cartesian Coordinates
      double z_v = vertex.z();
      double r_v = vertex.Perp();
      // Liquid hydrogen target cuts 
      // Target length 30 cm, centered at z = 65.0 cm
      // Upstream Diameter of target cell = 2.42 cm
      // ID of scattering chamber = 7.49 cm
      bool z_vertex_cut = 50.0 <= z_v && z_v <= 80.0; 
      bool r_vertex_cut = r_v < 3.745;
      // applied vertex cut
      if (!z_vertex_cut) continue;
      if (!r_vertex_cut) continue;
      // Declare a vector which quantizes the point of the intersection of a charged particle 
      //   with a plane in the middle of the scintillator 
      DVector3 IntersectionPoint;
      // Declare a vector which quantizes the unit vector of the charged particle track traversing
      //   through the scintillator with its origin at the intersection point
      DVector3 IntersectionDir;
      // Grab the paramteres associated to a track matched to the ST
      vector<DSCHitMatchParams> st_params;
      bool st_match = locDetectorMatches->Get_SCMatchParams(timeBasedTrack, st_params); 
      // If st_match = true, there is a match between this track and the ST
      if (!st_match) continue;

      bool st_match_pid = dParticleID->MatchToSC(timeBasedTrack, 
						 timeBasedTrack->rt, 
						 st_params[0].dSCHit, 
						 st_params[0].dSCHit->t, 
						 locSCHitMatchParams, 
						 true,
						 &IntersectionPoint, &IntersectionDir);      
      if(!st_match_pid) continue;  

      DVector3 momentum_vec;
      momentum_vec = timeBasedTrack->momentum();
      // applied vertex cut
      DLorentzVector Lor_Mom = timeBasedTrack->lorentzMomentum();
      double P_mag = Lor_Mom.P(); 
      double phi_mom = momentum_vec.Phi()*RAD2DEG;
      if (phi_mom < 0.0) phi_mom += 360.0;
      for (uint32_t j = 0; j < NCHANNELS; j++)
	{
	  if (phi_mom >= phi_sec[j][0] && phi_mom <= phi_sec[j][1])
	    {
	      phi_sec_cntr[j] += 1;
	    }
	}
		  // Grab the ST sector predicted to be hit by the charged time based track
		  // 0.069 rad = 4 deg, 0.087 rad = 5 deg, 0.105 rad = 6 deg
		  uint32_t st_pred_id = pid_algorithm[0]->PredictSCSector(timeBasedTrack->rt, 0.069);
		  uint32_t st_pred_id_index = st_pred_id - 1;
		  
		  if (st_pred_id != 0) 
		    {
		      phi_sec_pred_hit_cntr[st_pred_id_index] += 1;
		      //loop over the Hit object and get the real sector hit at ST
		      for (uint32_t j = 0; j < st_hits.size(); j++)
			{
			  uint32_t phi_sec_hit_sector       = st_hits[j]->sector;
			  uint32_t phi_sec_hit_sector_index = phi_sec_hit_sector - 1;
			  if (st_pred_id == phi_sec_hit_sector)
			    phi_sec_hit_cntr[phi_sec_hit_sector_index] += 1;
			}
		      //loop over the ADC object and get the real sector hit of ADC
		      for (uint32_t j = 0; j < st_adc_digi_hits.size(); j++)
			{
			  uint32_t phi_sec_adc_hit_sector       = st_adc_digi_hits[j]->sector;
			  uint32_t phi_sec_adc_hit_sector_index = phi_sec_adc_hit_sector - 1;
			  if (st_pred_id == phi_sec_adc_hit_sector)
			    phi_sec_adc_hit_cntr[phi_sec_adc_hit_sector_index] += 1;
			}
		    } // end if (st_pred_id != 0) 	
		  
		  if (st_match_pid)  // Get the intersection point which can not be obtained from st_match
		    {
		      // Grab the sector
		      Int_t sector_m = st_params[0].dSCHit->sector;
		      //Acquire the energy loss per unit length in the ST (arbitrary units)
		      double dEdx = st_params[0].dEdx;
		      double dphi = st_params[0].dDeltaPhiToHit*RAD2DEG;
		      // Fill dEdx vs Momentum 
		      h2_dedx_P_mag->Fill(P_mag,dEdx);
		      // Fill dEdx vs Momentum with cut on positive charges  
		      if (q > 0)
			{
			  h2_dedx_P_mag_postv->Fill(P_mag,dEdx);
			}
		      // Fill dEdx vs Momentum with cut on negative charges  
		      if (q < 0)
			{
			  h2_dedx_P_mag_negtv->Fill(P_mag,dEdx);
			}
		      // Obtain the intersection point with the ST		  
		      double phi_ip   = IntersectionPoint.Phi()*RAD2DEG;
		      // Correct phi calculation
		      if (phi_ip < 0.0) phi_ip += 360.0;
		      // Acquire the intersection point
		      sc_track_position.push_back(IntersectionPoint);
		      //Fill 2D histos
		      h2_phi_vs_sector->Fill(sector_m,phi_ip);
		      h2_dphi_sector->Fill(sector_m,dphi);
		    }  // PID match cut
    }                  // Charged track loop
  // Fill 2D histo
  for (uint32_t i = 0; i < sc_track_position.size(); i++)
    {
      h2_r_vs_z->Fill(sc_track_position[i].z(),sc_track_position[i].Perp());
      h2_y_vs_x->Fill(sc_track_position[i].x(),sc_track_position[i].y());
    }
  japp->RootUnLock(); 
  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_ST_online_tracking::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_ST_online_tracking::fini(void)
{
  // Called before program exit after event processing is finished.
  for (uint32_t i = 0; i < NCHANNELS; i++)
    {
      //hit object
      h_phi_sec_hit_cntr->Fill(i+1,double(phi_sec_hit_cntr[i]));
      h_phi_sec_pred_hit_cntr->Fill(i+1,double(phi_sec_pred_hit_cntr[i]));

      h_phi_sec_adc_cntr->Fill(i+1,double(phi_sec_adc_hit_cntr[i]));
     

    }
  // //How to get Binomial errors in an efficiency plot
  // // hit object efficiency
  // h_phi_sec_hit_cntr->Sumw2();
  // h_phi_sec_pred_hit_cntr->Sumw2();
  // pEff->Sumw2();
  // pEff->Divide(h_phi_sec_hit_cntr,h_phi_sec_pred_hit_cntr,1,1,"B");
  // // adc efficiency
  // h_phi_sec_adc_cntr->Sumw2();
  // pEff_adc->Sumw2();
  // pEff_adc->Divide(h_phi_sec_adc_cntr,h_phi_sec_pred_hit_cntr,1,1,"B");
  return NOERROR;
}





