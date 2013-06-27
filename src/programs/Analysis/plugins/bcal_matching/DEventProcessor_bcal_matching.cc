#include "DEventProcessor_bcal_matching.h"

#include "DANA/DApplication.h"
#include "BCAL/DBCALShower.h"
#include "BCAL/DBCALGeometry.h"
#include "TRACKING/DMCThrown.h"
#include "PID/DParticleID.h"
#include "units.h"

#define DPRINT(x) cout << #x "=" << x << endl

// The executable should define the ROOTfile global variable. It will
// be automatically linked when dlopen is called.
extern TFile *ROOTfile;

// Routine used to create our DEventProcessor
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new DEventProcessor_bcal_matching());
  }
} // "C"


//------------------
// init
//------------------
jerror_t DEventProcessor_bcal_matching::init(void)
{

  swim_steps = new DReferenceTrajectory::swim_step_t[max_swim_steps];

  //set up the tree
  bcal_matching_tree = new TTree("bcal_matching_tree","bcal_matching_tree");

  bcal_matching_tree->Branch("p_pi_thrown",&m_p_pi_thrown,"p_pi_thrown/F");
  bcal_matching_tree->Branch("phi_pi_thrown",&m_phi_pi_thrown,"phi_pi_thrown/F");
  bcal_matching_tree->Branch("theta_pi_thrown",&m_theta_pi_thrown,"theta_pi_thrown/F");
  bcal_matching_tree->Branch("BCAL_intersect_z_thrown",&m_BCAL_intersect_z_thrown,"BCAL_intersect_z_thrown/F");
  bcal_matching_tree->Branch("BCAL_intersect_thrown",&m_BCAL_intersect_thrown,"BCAL_intersect_thrown/O");

  bcal_matching_tree->Branch("n_matched_showers",&m_n_matched_showers,"n_matched_showers/I");

  bcal_matching_tree->Branch("n_showers",&m_n_showers,"n_showers/I");

  bcal_matching_tree->Branch("E_shower",m_E_shower,"E_shower[n_showers]/F");
  bcal_matching_tree->Branch("z_shower",m_z_shower,"z_shower[n_showers]/F");
  bcal_matching_tree->Branch("phi_shower",m_phi_shower,"phi_shower[n_showers]/F");
  bcal_matching_tree->Branch("r_shower",m_r_shower,"r_shower[n_showers]/F");
  bcal_matching_tree->Branch("r_poca_thrown",m_r_poca_thrown,"r_poca_thrown[n_showers]/F");
  bcal_matching_tree->Branch("d_shower_track_thrown",m_d_shower_track_thrown,"d_shower_track_thrown[n_showers]/F");
  bcal_matching_tree->Branch("dphi_thrown",m_dphi_thrown,"dphi_thrown[n_showers]/F");
  bcal_matching_tree->Branch("dz_thrown",m_dz_thrown,"dz_thrown[n_showers]/F");
  bcal_matching_tree->Branch("dt_thrown",m_dt_thrown,"dt_thrown[n_showers]/F");


  bcal_matching_tree->Branch("track_recon",&m_track_recon,"track_recon/O");
  bcal_matching_tree->Branch("r_poca_recon",m_r_poca_recon,"r_poca_recon[n_showers]/F");
  bcal_matching_tree->Branch("d_shower_track_recon",m_d_shower_track_recon,"d_shower_track_recon[n_showers]/F");
  bcal_matching_tree->Branch("dphi_recon",m_dphi_recon,"dphi_recon[n_showers]/F");
  bcal_matching_tree->Branch("dz_recon",m_dz_recon,"dz_recon[n_showers]/F");
  bcal_matching_tree->Branch("dt_recon",m_dt_recon,"dt_recon[n_showers]/F");

  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_bcal_matching::brun(JEventLoop *eventLoop, int runnumber)
{
  DApplication* dapp = dynamic_cast<DApplication*>(eventLoop->GetJApplication());    
  if (!dapp) {
    throw RESOURCE_UNAVAILABLE;
  }

  bfield = dapp->GetBfield();

  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_bcal_matching::evnt(JEventLoop *loop, int eventnumber)
{
  vector<const DBCALShower*> bcal_showers;
  vector<const DMCThrown*> mcthrowns;
  vector<const DParticleID*> pid_vect;
  vector<const DChargedTrackHypothesis*> charged_tracks;
  loop->Get(bcal_showers);
  loop->Get(mcthrowns);
  loop->Get(pid_vect);
  loop->Get(charged_tracks);

  if (pid_vect.size() == 0) return NOERROR;
  const DParticleID* pid = pid_vect[0];

  //the assumption here is that we are analyzing single-charged-pion events
  //however there may be >1 mcthrowns due to in-flight decays or something
  if (mcthrowns.size() == 0) return NOERROR;
  if (mcthrowns[0]->PID() != PiPlus && mcthrowns[0]->PID() != PiMinus) return NOERROR;
  DVector3 thrown_pos = mcthrowns[0]->position();
  DVector3 thrown_mom = mcthrowns[0]->momentum();
  double q = mcthrowns[0]->charge();

  //check if there is a "good_track" i.e. a reconstructed track that matches the thrown PID and also the momentum, within some tolerance
  //this assumes that there is only one "good" track, which I hope is the case, at least most of the time
  const DChargedTrackHypothesis* good_track = NULL;
  for (int i=0; i<(int)charged_tracks.size(); i++) {
    if (charged_tracks[i]->PID() != mcthrowns[0]->PID()) continue;

    DVector3 recon_mom = charged_tracks[i]->momentum();
    //I should check if this matching criteria is optimal
    if (fabs(recon_mom.Mag() - thrown_mom.Mag())/thrown_mom.Mag() > .2) continue;
    if (fabs(recon_mom.Theta() - thrown_mom.Theta()) > .09) continue;
    if (fabs(recon_mom.Phi() - thrown_mom.Phi()) > .09) continue;

    good_track = charged_tracks[i];
  }
  if (good_track!=NULL) {
    m_track_recon = true;
  } else {
    m_track_recon = false;
  }

  LockState();

  //fill branches with information of thrown track
  m_p_pi_thrown = thrown_mom.Mag();
  m_phi_pi_thrown = thrown_mom.Phi();
  m_theta_pi_thrown = thrown_mom.Theta();

  //swim the track using the thrown parameters
  //doing this inside this lock will slow things down if we run multi-threaded
  DReferenceTrajectory rt_thrown(bfield, 1.0, swim_steps, max_swim_steps);
  rt_thrown.Swim(thrown_pos,thrown_mom,q);

  //get point of intersection with inner radius of BCAL
  DVector3 intersect_pos;
  jerror_t result = rt_thrown.GetIntersectionWithRadius(DBCALGeometry::BCALINNERRAD, intersect_pos);
  if (result == VALUE_OUT_OF_RANGE || intersect_pos.Z() > 407. || intersect_pos.Z() < 17.) {
    m_BCAL_intersect_thrown = false;
    m_BCAL_intersect_z_thrown = 0.0;
  } else {
    m_BCAL_intersect_thrown = true;
    m_BCAL_intersect_z_thrown = intersect_pos.Z();
  }

  //use the default BCAL matching algorithm (whatever that may be)
  //to see if any (and how many) cluster get matched to the track
  deque<const DBCALShower*> matched_bcal_showers;
  double projected_time=0.0, path_length=0.0, flight_time=0.0;
  pid->MatchToBCAL(&rt_thrown, bcal_showers, matched_bcal_showers, projected_time, path_length, flight_time);
  m_n_matched_showers = matched_bcal_showers.size();

  m_n_showers = 0;

  //loop thru all BCAL showers and compare each one with both the thrown track and the "good" reconstructed track, if one exists
  for (int i = 0; i < (int) bcal_showers.size() && i < max_n_showers; i++) {
    const DBCALShower* shower = bcal_showers[i];

    //first fill branches with shower info
    m_n_showers++;
    m_E_shower[i] = shower->E;
    m_phi_shower[i] = atan2(shower->y,shower->x);
    m_z_shower[i] = shower->z;
    m_r_shower[i] = sqrt(shower->y*shower->y + shower->x*shower->x);

    DVector3 bcal_pos(shower->x, shower->y, shower->z);

    //thrown track POCA
    double path_length=0.0;
    double flight_time=0.0;
    double d = rt_thrown.DistToRTwithTime(bcal_pos,&path_length,&flight_time);
    DVector3 track_pos = rt_thrown.GetLastDOCAPoint();
    m_r_poca_thrown[i] = track_pos.Perp(); //r-coordinate of POCA

    double dz = track_pos.z() - bcal_pos.z();
    double dphi = track_pos.Phi() - bcal_pos.Phi();
    double dphi_alt = (track_pos.Phi() > bcal_pos.Phi()) ? 
                      track_pos.Phi() - bcal_pos.Phi() - 2*M_PI :
                      bcal_pos.Phi() - track_pos.Phi() - 2*M_PI;

    dphi = min(fabs(dphi),fabs(dphi_alt));

    m_d_shower_track_thrown[i] = d;
    m_dphi_thrown[i] = dphi;
    m_dz_thrown[i] = dz;
    m_dt_thrown[i] = shower->t - flight_time;

    //reconstructed track POCA
    if (good_track!=NULL) {
      const DReferenceTrajectory* rt_recon = good_track->dRT; //should probably null-check this

      d = rt_recon->DistToRTwithTime(bcal_pos,&path_length,&flight_time);
      track_pos = rt_recon->GetLastDOCAPoint();
      m_r_poca_recon[i] = track_pos.Perp();

      dz = track_pos.z() - bcal_pos.z();
      dphi = track_pos.Phi() - bcal_pos.Phi();
      dphi_alt = (track_pos.Phi() > bcal_pos.Phi()) ? 
                      track_pos.Phi() - bcal_pos.Phi() - 2*M_PI :
                      bcal_pos.Phi() - track_pos.Phi() - 2*M_PI;

      dphi = min(fabs(dphi),fabs(dphi_alt));

      m_d_shower_track_recon[i] = d;
      m_dphi_recon[i] = dphi;
      m_dz_recon[i] = dz;
      m_dt_recon[i] = shower->t - flight_time;
    } else {
      m_r_poca_recon[i] = 0.0;
      m_d_shower_track_recon[i] = 0.0;
      m_dphi_recon[i] = 0.0;
      m_dz_recon[i] = 0.0;
      m_dt_recon[i] = 0.0;
    }
  }

  bcal_matching_tree->Fill();

  UnlockState();

  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_bcal_matching::erun(void)
{
  // Any final calculations on histograms (like dividing them)
  // should be done here. This may get called more than once.
  return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_bcal_matching::fini(void)
{
  delete[] swim_steps;
  return NOERROR;
}

