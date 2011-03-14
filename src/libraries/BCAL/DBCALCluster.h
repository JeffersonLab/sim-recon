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
  
  DBCALCluster(){}
  DBCALCluster( const DBCALPoint* point );

  vector< const DBCALPoint* > points() const { return m_points; }

  int nCells() const { return m_points.size(); }
  
  // the total energy in the cluster
  
  float E() const { return m_E; }
  
  // this is the time at the inner radius of BCAL assuming shower
  // particles propagte into module at the speed of light
  float t() const { return m_t; }
  float sigT() const { return m_sig_t; }
  
  // assuming a photon leaving the target, this the estimate of t0
  // (could be helpful for photon/pion discrimination)
  float t0() const; 
  
  // location of cluster in spherical coordinates with the origin
  // at the center of the target -- errors are not rigrously derived!
  float rho() const { return m_rho; }
  float sigRho() const { return m_sig_rho; }
  
  float theta() const { return m_theta; }
  float sigTheta() const { return m_sig_theta; }
  
  float phi() const { return m_phi; }
  float sigPhi() const { return m_sig_phi; }
  
  // these functions modify the cluster
  void addPoint( const DBCALPoint* point );
  void mergeClust( const DBCALCluster& clust );
  
  // this prints out info
  void toStrings( vector< pair < string, string > > &items ) const;

private:
  
  void makeFromPoints();
  void clear();
  
  vector< const DBCALPoint* > m_points;
  
  float m_E;
  
  float m_t;
  float m_sig_t;
  
  float m_rho;
  float m_sig_rho;
  float m_theta;
  float m_sig_theta;
  float m_phi;
  float m_sig_phi;
  
};

#endif
