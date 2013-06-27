#ifndef _DEventProcessor_bcal_matching_
#define _DEventProcessor_bcal_matching_

#include <JANA/JEventProcessor.h>
using namespace jana;

#include <TFile.h>
#include <TTree.h>

#include "units.h"
#include <HDGEOMETRY/DMagneticFieldMap.h>
#include <TRACKING/DReferenceTrajectory.h>

class DEventProcessor_bcal_matching:public JEventProcessor{
 public:
  DEventProcessor_bcal_matching(){}
  ~DEventProcessor_bcal_matching(){}
  const char* className(void){return "DEventProcessor_bcal_matching";}

 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.

  const DMagneticFieldMap* bfield;

  static const int max_swim_steps = 10000;
  DReferenceTrajectory::swim_step_t *swim_steps;

  TTree *bcal_matching_tree;

  //each member below corresponds to a branch on the tree
  float m_p_pi_thrown;
  float m_phi_pi_thrown;
  float m_theta_pi_thrown;
  //the z-position at which the track (reference trajectory) intersects the inner radius of the BCAL
  bool m_BCAL_intersect_thrown; //does it cross the BCAL at all?
  float m_BCAL_intersect_z_thrown; //if so, where?

  int m_n_matched_showers;

  static const int max_n_showers=20;
  int m_n_showers;

  //shower properties
  float m_E_shower[max_n_showers];
  float m_z_shower[max_n_showers];
  float m_phi_shower[max_n_showers];
  float m_r_shower[max_n_showers];

  //things related to thrown track POCA to shower
  float m_r_poca_thrown[max_n_showers]; //POCA r-coordinate
  float m_d_shower_track_thrown[max_n_showers]; //distance from track to cluster
  float m_dphi_thrown[max_n_showers]; //azimuthal separation between track and cluster
  float m_dz_thrown[max_n_showers]; //z separation
  float m_dt_thrown[max_n_showers]; //(shower time) - (track flight time)

  //reconstructed track POCA to shower
  bool m_track_recon; //Is there a reconstructed track that matches the thrown one?
  float m_r_poca_recon[max_n_showers]; //POCA r-coordinate
  float m_d_shower_track_recon[max_n_showers]; //distance from track to cluster
  float m_dphi_recon[max_n_showers]; //separation between track and cluster
  float m_dz_recon[max_n_showers];
  float m_dt_recon[max_n_showers];


};

#endif // _DEventProcessor_bcal_matching_
