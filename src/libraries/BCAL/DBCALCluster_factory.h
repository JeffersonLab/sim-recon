#ifndef _DBCALCluster_factory_
#define _DBCALCluster_factory_

/*
 *  DBCALCluster_factory.h
 *
 *  Created by Matthew Shepherd on 3/12/11.
 *
 */

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>

using namespace jana;

#include "BCAL/DBCALHit.h"
#include "BCAL/DBCALCluster.h"
#include "BCAL/DBCALUnifiedHit.h"
#include "TRACKING/DTrackWireBased.h"

//#include "TTree.h"
//#include "TFile.h"
#include "TF1.h"
#include "TH2F.h"

class DBCALCluster_factory : public JFactory< DBCALCluster > {
  
public:
  
  DBCALCluster_factory();
  ~DBCALCluster_factory(){}
  
private:

  jerror_t evnt(JEventLoop *loop, uint64_t eventnumber);	
  jerror_t brun(JEventLoop *loop, int32_t runnumber);
  
  void clearPoints();
  
  // these routines combine points and clusters together

  vector<DBCALCluster*> clusterize( vector< const DBCALPoint* > points, vector< const DBCALPoint* > usedPoints,  vector< const DBCALUnifiedHit* > hits, vector< const DTrackWireBased* > tracks ) const;
  void merge( vector<DBCALCluster*>& clusters ) const;

  // This routine removes a point from its original cluster and adds it to its closest cluster if applicable.
  void recycle_points( vector<const DBCALPoint*> usedPoints, vector<DBCALCluster*>& clusters ) const; 
 
  // these are the routines used for testing whether things should be
  // combined -- right now very basic, but can be fine tuned in the future

  bool overlap( const DBCALCluster& highEClust,
                const DBCALCluster& lowEClust ) const;
  
  bool overlap( const DBCALCluster& clust,
                const DBCALPoint* point ) const;
 
  bool overlap_charged( const DBCALCluster& clust, 
			const DBCALPoint* point, DVector3 track_pos ) const;
 
  bool overlap( const DBCALCluster& clust, 
                const DBCALUnifiedHit* hit ) const; 
  
  uint32_t BCALCLUSTERVERBOSE;
  float m_mergeSig;
  float m_moliereRadius;
  float m_clust_hit_timecut;
  float m_timeCut;
  double m_z_target_center;
  vector<double> effective_velocities;
  vector< vector<double > > attenuation_parameters;
  TH2F* charged_dist;
  TF1* charged_fit;

  /*
  TF1* sep_inclusion_curve;
  TF1* dtheta_inclusion_curve;
  TF1* dphi_inclusion_curve;
  TF1* C1_parm;
  TF1* C2_parm;
  */
  jerror_t init();
  jerror_t fini();
  
};

#endif 

