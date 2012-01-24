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

DBCALCluster::DBCALCluster( const DBCALPoint* point ) : m_points ( 0 ) {

  m_points.push_back( point );
  AddAssociatedObject( point );
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
  AddAssociatedObject( point );
  
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
    AddAssociatedObject( *pt );
  }
  
  makeFromPoints();
}

void
DBCALCluster::makeFromPoints(){
  
  clear();

  //In this function we take a weighted average of the phis/thetas/etc
  //of the individual DBCALPoint's to get the phi/theta/etc of the cluster. The
  //average of phi is weighted by the energy of the hit, while the other
  //quantities are weighted by energy squared. This sort of makes sense because
  //the error of the phi measurement is independent of the energy of the hit,
  //so it is only weighted because higher energy hits make up "more" of the
  //cluster. In the case of theta/timing measurements, the error decreases with
  //hit energy, so it makes sense to weight these with an extra factor of E, at
  //least to some low level of rigor. In any case, the weighting scheme below
  //produces far better angular resolutions for clusters than other
  //weightings tested.
  double wt1;
  double wt2;
  double sum_wt1 = 0;
  double sum_wt1_sq = 0;
  double sum_wt2 = 0;
  double sum_wt2_sq = 0;
  //to do an average of phi correctly for cases with angles near 0, we need to average cos(phi) and sin(phi) instead of phi itself
  double sum_sin_phi=0;
  double sum_cos_phi=0;

  for( vector< const DBCALPoint* >::const_iterator pt = m_points.begin();
       pt != m_points.end();
      ++pt ){
   
    double E = (**pt).E();

    m_E += E;

    wt1 = E;
    wt2 = E*E;

    sum_wt1 += wt1;
    sum_wt1_sq += wt1*wt1;
    sum_wt2 += wt2;
    sum_wt2_sq += wt2*wt2;
    
    m_t += (**pt).tInnerRadius() * wt2;
    m_sig_t += (**pt).tInnerRadius() * (**pt).tInnerRadius() * wt2;

    m_theta += (**pt).theta() * wt2;
    m_sig_theta += (**pt).theta() * (**pt).theta() * wt2;

    sum_sin_phi += sin((**pt).phi()) * wt1;
    sum_cos_phi += cos((**pt).phi()) * wt1;

    m_rho += (**pt).rho() * wt2;
    m_sig_rho += (**pt).rho() * (**pt).rho() * wt2;
  }
  
  // now adjust the standard deviations and averages
  // the variance of the mean of a weighted distribution is s^2/n_eff, where s^2 is the variance of the sample and n_eff is as calculated below
  
  //double n_eff1 = sum_wt1*sum_wt1/sum_wt1_sq;
  double n_eff2 = sum_wt2*sum_wt2/sum_wt2_sq;
  double n = m_points.size();

  m_t /= sum_wt2;
  m_sig_t /= sum_wt2;
  m_sig_t -= ( m_t * m_t );
  m_sig_t = sqrt( m_sig_t );
  m_sig_t /= sqrt(n_eff2);
  
  m_theta /= sum_wt2;
  /*m_sig_theta /= sum_wt2;
  m_sig_theta -= ( m_theta * m_theta );
  m_sig_theta = sqrt( fabs( m_sig_theta ) );
  m_sig_theta /= sqrt(n_eff2);*/

  // The method below for determining sig_theta works better than the one
  // above. sigma_z is determined by fitting. sigma_z is supposed to be of
  // the form a/sqrt(E), but this form fits the data better. This discrepancy
  // is probably due to contributions from dark hits and incomplete
  // gathering of hits to form a cluster.
  double sigma_z = 1.17/sqrt(m_E)+.17/m_E;
  m_sig_theta = sigma_z*sin(m_theta)*sin(m_theta)/DBCALGeometry::BCALINNERRAD;
  
  m_phi = atan2(sum_sin_phi,sum_cos_phi);
  if( m_phi < 0 ) m_phi += 2*PI;
  // calculate the RMS of phi
  m_sig_phi=0;
  for( vector< const DBCALPoint* >::const_iterator pt = m_points.begin();
       pt != m_points.end();
       ++pt ){

    double E = (**pt).E();

    wt1 = E;

    double deltaPhi = (**pt).phi() - m_phi;
    double deltaPhiAlt = ( (**pt).phi() > m_phi ?
			   (**pt).phi() - m_phi - 2*PI :
			   m_phi - (**pt).phi() - 2*PI );
    deltaPhi = min( fabs( deltaPhi ), fabs( deltaPhiAlt ) );
    m_sig_phi += deltaPhi * deltaPhi * wt1;

  }
  m_sig_phi /= sum_wt1;
  m_sig_phi = sqrt( fabs(m_sig_phi) );
  //this should be division, by sqrt(n_eff1), but this works better
  m_sig_phi /= sqrt(n);

  m_rho /= sum_wt2;
  m_sig_rho /= sum_wt2;
  m_sig_rho -= ( m_rho * m_rho );
  m_sig_rho = sqrt( fabs( m_sig_rho ) );
  m_sig_rho /= sqrt(n_eff2);
  
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
  
  if( m_sig_phi   < 1E-6 ) m_sig_phi   = m_points.at(0)->sigPhi()/sqrt(n);
  if( m_sig_theta < 1E-6 ) m_sig_theta = m_points.at(0)->sigPhi()/sqrt(n);
  
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

