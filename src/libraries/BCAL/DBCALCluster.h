#ifndef _DBCALCluster_
#define _DBCALCluster_

/*
 *  DBCALCluster.h
 *
 *  Created by Matthew Shepherd on 3/12/11.
 *
 */

#include "BCAL/DBCALGeometry.h"
#include "BCAL/DBCALPoint.h"

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

#include <vector>

using namespace jana;
using namespace std;

class DBCALCluster : public JObject {

public:
  
  JOBJECT_PUBLIC( DBCALCluster );
  
  DBCALCluster(double z_target_center, const DBCALGeometry *locGeom);
  DBCALCluster(const DBCALPoint* point, double z_target_center, double q, const DBCALGeometry *locGeom);

  vector< const DBCALPoint* > points() const { return m_points; }
  // Returns a vector of the single-ended hits used in the cluster.
  // This returns a pair because we need to store the attenuated-corrected energy
  // and the DBCALUnifiedHit cannot hold this information.
  vector< pair<const DBCALUnifiedHit*,double> > hits() const { return m_single_ended_hits; }

  int nCells() const { return m_points.size(); }
  
  // the total energy in the cluster
  
  float E() const { return m_E; }
  float E_preshower() const { return m_E_preshower; }
  float E_L2() const { return m_E_L2; }
  float E_L3() const { return m_E_L3; }
  float E_L4() const { return m_E_L4; }
  
  // this is the time at the inner radius of BCAL assuming shower
  // particles propagte into module at the speed of light
  float t() const { return m_t; }
  float sigT() const { return m_sig_t; }
  float rmsTime() const { return m_t_rms; } 
 
  // assuming a photon leaving the target, this the estimate of t0
  // (could be helpful for photon/pion discrimination)
  float t0() const; 
  
  // location of cluster in spherical coordinates with the origin
  // at the center of the target -- WARNING: errors are not rigorously derived!
  float rho() const { return m_rho; }
  float sigRho() const { return m_sig_rho; }
  
  float theta() const { return m_theta; }
  float sigTheta() const { return m_sig_theta; }
  
  float phi() const { return m_phi; }
  float sigPhi() const { return m_sig_phi; }
 
  int Q() const { return charge; }
 
  // these functions modify the cluster
  void addPoint( const DBCALPoint* point , int q );
  void addHit ( const DBCALUnifiedHit* hit, double hit_E_unattenuated );
  void mergeClust( const DBCALCluster& clust, double point_reatten_E );
  void removePoint( const DBCALPoint* point ); 
 
  // this prints out info
  void toStrings( vector< pair < string, string > > &items ) const;

private:
  
  void makeFromPoints();
  void clear();
  
  vector< const DBCALPoint* > m_points;
  vector< const DBCALPoint* > m_points_remove;
  vector< pair<const DBCALUnifiedHit*,double> > m_single_ended_hits; //Store single-ended hits together with their unattenuated energies

  float m_hit_E_unattenuated_sum; //attenuation-corrected sum of energies from single-ended hits
  float m_point_reatten_E_sum;
  float m_E_points;
  float m_E;
  float m_E_preshower;  
  float m_E_L2;
  float m_E_L3;
  float m_E_L4;

  float m_t;
  float m_sig_t;
  float m_t_rms;  

  float m_rho;
  float m_sig_rho;
  float m_theta;
  float m_sig_theta;
  float m_phi;
  float m_sig_phi;

  float m_z_target_center;
  int new_point_q; 
  int charge;
  const DBCALGeometry *m_BCALGeom;
 
};

#endif
