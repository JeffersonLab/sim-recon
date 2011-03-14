/*
 *  DBCALCluster.cc
 *
 *  Created by Matthew Shepherd on 3/13/11.
 *
 */

#include "BCAL/DBCALCluster.h"
#include "BCAL/DBCALPoint.h"
#include "units.h"

#include <math.h>

DBCALCluster::DBCALCluster( const DBCALPoint* point ){

  m_points.push_back( point );
  makeFromPoints();
}

float
DBCALCluster::t0() const {
  
  float path = DBCALGeometry::BCALINNERRAD / sin( m_theta );
  
  return m_t - ( path / ( 30 * k_cm / k_nsec ) );
}

void
DBCALCluster::addPoint( const DBCALPoint* point ){
  
  // offset phi of the point by +- 2PI to match the cluster
  if( phi() > point->phi() ){
    
    if( fabs( phi() - point->phi() - 2*PI ) < PI ) point->add2Pi();
  }
  else{
    
    if( fabs( point->phi() - phi() - 2*PI ) < PI ) point->sub2Pi();
  }
  
  m_points.push_back( point );
  
  makeFromPoints();
}

void
DBCALCluster::mergeClust( const DBCALCluster& clust ){

  vector< const DBCALPoint* > otherPoints = clust.points();
  
  for( vector< const DBCALPoint* >::const_iterator pt = otherPoints.begin();
      pt != otherPoints.end();
      ++pt ){

    // offset phi of the point by +- 2PI to match the cluster if needed
    if( phi() > (**pt).phi() ){
      
      if( fabs( phi() - (**pt).phi() - 2*PI ) < PI ) (**pt).add2Pi();
    }
    else{
      
      if( fabs( (**pt).phi() - phi() - 2*PI ) < PI ) (**pt).sub2Pi();
    }
    
    m_points.push_back( *pt );
  }
  
  makeFromPoints();
}

void
DBCALCluster::makeFromPoints(){
  
  clear();

  float wt;
  float sum_wt = 0;
  
  for( vector< const DBCALPoint* >::const_iterator pt = m_points.begin();
       pt != m_points.end();
      ++pt ){
   
    m_E += (**pt).E();

    // the choice of weight should be considered and optimized 
    // start with energy weighting for now

    wt = (**pt).E();
    sum_wt += wt;
    
    m_t += (**pt).tInnerRadius() * wt;
    m_sig_t += (**pt).tInnerRadius() * (**pt).tInnerRadius() * wt;

    m_theta += (**pt).theta() * wt;
    m_sig_theta += (**pt).theta() * (**pt).theta() * wt;

    m_phi += (**pt).phi() * wt;
    m_sig_phi += (**pt).phi() * (**pt).phi() * wt;

    m_rho += (**pt).rho() * wt;
    m_sig_rho += (**pt).rho() * (**pt).rho() * wt;
  }
  
  // now adjust the standard deviations and averages
  
  m_t /= sum_wt;
  m_sig_t /= sum_wt;
  m_sig_t -= ( m_t * m_t );
  m_sig_t = sqrt( m_sig_t );
  
  m_theta /= sum_wt;
  m_sig_theta /= sum_wt;
  m_sig_theta -= ( m_theta * m_theta );
  m_sig_theta = sqrt( fabs( m_sig_theta ) );
  
  m_phi /= sum_wt;
  m_sig_phi /= sum_wt;
  m_sig_phi -= ( m_phi * m_phi );
  m_sig_phi = sqrt( fabs( m_sig_phi ) );

  m_rho /= sum_wt;
  m_sig_rho /= sum_wt;
  m_sig_rho -= ( m_rho * m_rho );
  m_sig_rho = sqrt( fabs( m_sig_rho ) );
  
  // if we only have one point, then set sizes based on the size
  // of the point
  
  if( m_points.size() == 1 ){
    
    m_sig_theta = m_points[0]->sigTheta();
    m_sig_phi = m_points[0]->sigPhi();
    m_sig_rho = m_points[0]->sigRho();
  }
  
  // it is also possible to have sigPhi be zero when all points
  // lie in the same radial "tower" and sigTheta can be zero
  // if two points are at the end of the module (same z) in the
  // same layer in these cases set sigmas to the characteristic 
  // size of the cell
  
  if( m_sig_phi   < 1E-6 ) m_sig_phi   = m_points.at(0)->sigPhi();
  if( m_sig_theta < 1E-6 ) m_sig_theta = m_points.at(0)->sigPhi();
  
  // correct phi of the cluster if we've drifted outside of 0 - 2PI
  if( m_phi > 2*PI ) m_phi -= 2*PI;
  if( m_phi < 0 ) m_phi += 2*PI;
}

void 
DBCALCluster::toStrings( vector< pair < string, string > > &items) const {

  AddString(items, "r", "%5.2f", m_rho * sin( m_theta ) );
  AddString(items, "phi", "%5.2f", m_phi );
//  AddString(items, "dphi", "%5.2f", m_sig_phi );
  AddString(items, "z", "%5.2f", m_rho * cos( m_theta ) + 65 );
  AddString(items, "theta", "%5.2f", m_theta);
//  AddString(items, "dtheta", "%5.2f", m_sig_theta);
  AddString(items, "t", "%5.2f", m_t );
  AddString(items, "E", "%5.2f", m_E );
  AddString(items, "N_cell", "%d", m_points.size() );
}


void
DBCALCluster::clear(){
 
  m_E = 0;
  
  m_t = 0;
  m_sig_t = 0;
  
  m_theta = 0;
  m_sig_theta = 0;

  m_phi = 0;
  m_sig_phi = 0;
  
  m_rho = 0;
  m_sig_rho = 0;
  
}

