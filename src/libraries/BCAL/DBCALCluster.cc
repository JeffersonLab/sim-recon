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

DBCALCluster::DBCALCluster( const DBCALPoint* point, double z_target_center )
  : m_points ( 0 ),  m_hit_E_unattenuated_sum(0.0),  m_z_target_center(z_target_center) {

  m_points.push_back( point );
  AddAssociatedObject( point );
  makeFromPoints();
}

DBCALCluster::DBCALCluster(double z_target_center) : m_z_target_center(z_target_center) {
  clear(); //initialize all values to zero
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
DBCALCluster::addHit( const DBCALUnifiedHit* hit, double hit_E_unattenuated ){
 
  m_hit_E_unattenuated_sum += hit_E_unattenuated; // add the energy of all hits in a cluster
  m_single_ended_hits.push_back(make_pair(hit,hit_E_unattenuated));
  AddAssociatedObject( hit );
 
  makeFromPoints();  // call makeFromPoints so we can include hit energy in the clusterization, but don't use any of the hit positions or times to include in the cluster averaging.

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

  vector<pair<const DBCALUnifiedHit*,double> > otherHits = clust.hits();

  for( vector<pair<const DBCALUnifiedHit*,double> >::const_iterator hit = otherHits.begin();
      hit != otherHits.end();
      ++hit ){

    m_hit_E_unattenuated_sum += hit->second; //add the unattenuated energy
    m_single_ended_hits.push_back( *hit );
    AddAssociatedObject( hit->first );
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

  //Also, don't include hits from the 4th layer in the averages.
  //Three reasons for doing this:
  //1) Averaging the quantities theta, rho, and phi (rather than z,r,phi) is
  //based on the assumption that the shower develops in a linear fashion inside
  //the BCAL. This assumption breaks down
  //deep in the BCAL (see GlueX doc 2229 slides 21-22).
  //(this also applies to layer 3, so at some point, we may want to
  //consider if we want treat layer 3 differently as well)
  //2) The errors in the outer layers for measurements of t and z are
  //relatively large due to the lack
  //of TDC readout so we don't gain much by including them.
  //3) Contamination by noise. (Noise is greatest in the 4th layer)
  //These outer layer hits will of course still contribute to the energy sum.
  //It can happen that a cluster is made up entirely of fourth layer hits.
  //In this case, we must use all hits in the average.

  //Add single-ended hit energies to the cluster energy, but don't use the single-ended hits 
  //to calculate the cluster centroid or time.

  int n = m_points.size();
  int n4 = 0; //number of 4th layer points in the cluster
  for( vector< const DBCALPoint* >::const_iterator pt = m_points.begin();
       pt != m_points.end();
      ++pt ){
    if ((**pt).layer() == 4) n4++;
  }

  bool average_layer4;
  int n_avg; //number of points involved in the average
  if (n == n4) { //if all points are in the 4th layer
    average_layer4 = true; //include 4th layer in average
    n_avg = n4;
  } else {
    average_layer4 = false;
    n_avg = n - n4; //all points except 4th layer involved in average
  }

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

    m_E_points += E;
    m_E = m_E_points + m_hit_E_unattenuated_sum;  // add the energy sum from points to the energy sum from single ended hits
    double wt1, wt2;
    if ((**pt).layer() != 4 || average_layer4) {
      wt1 = E;
      wt2 = E*E;
    } else {
      wt1 = 0;
      wt2 = 0;
    }

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
  // above. parameters of sigma_z are determined using errors when reconstructing MC data.
  // Using m_E_points because the cluster z position only depends on the point z positions
  // and point energies.
  double sigma_z = sqrt(1.394*1.394/m_E_points + 0.859*0.859);
  m_sig_theta = sigma_z*sin(m_theta)*sin(m_theta)/DBCALGeometry::BCALINNERRAD;
  
  m_phi = atan2(sum_sin_phi,sum_cos_phi);
  if( m_phi < 0 ) m_phi += 2*PI;
  // calculate the RMS of phi
  m_sig_phi=0;
  for( vector< const DBCALPoint* >::const_iterator pt = m_points.begin();
       pt != m_points.end();
       ++pt ){

    double E = (**pt).E();

    double wt1;
    if ((**pt).layer() != 4) {
      wt1 = E;
    } else {
      wt1 = 0;
    }

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
  m_sig_phi /= sqrt(n_avg);

  m_rho /= sum_wt2;
  m_sig_rho /= sum_wt2;
  m_sig_rho -= ( m_rho * m_rho );
  m_sig_rho = sqrt( fabs( m_sig_rho ) );
  m_sig_rho /= sqrt(n_eff2);
  
  //Since the z-positions of the DBCALPoint's are not constrained to be
  //physical (be inside the BCAL), we can end up with clusters outside of the
  //BCAL.
  //Also, in extreme cases, the averaging procedure itself can lead to positions
  //outside the BCAL.
  //So at this point, convert our cluster position to cylindrical coordinates,
  //constrain r and z to be inside the BCAL and then convert back to
  //spherical coordinates.
  double z = m_rho*cos(m_theta) + m_z_target_center;
  double r = m_rho*sin(m_theta); 

  double bcal_down = DBCALGeometry::GLOBAL_CENTER + DBCALGeometry::BCALFIBERLENGTH/2.0;
  double bcal_up = DBCALGeometry::GLOBAL_CENTER - DBCALGeometry::BCALFIBERLENGTH/2.0;
  if (z > bcal_down) z = bcal_down;
  if (z < bcal_up) z = bcal_up;

  double r_min = DBCALGeometry::r(DBCALGeometry::cellId(1,1,1)); //This is the center of the first layer of the BCAL. This is the minimum value of r that a DBCALPoint will report.
  double r_max = DBCALGeometry::r(DBCALGeometry::cellId(1,4,1)); //Maximum r that a DBCALPoint will report
  if (r > r_max) r = r_max;
  if (r < r_min) r = r_min;

  m_theta = atan2(r, (z-m_z_target_center));
  m_rho = sqrt(r*r + (z-m_z_target_center)*(z-m_z_target_center));

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
  // (this condition can also occur if we have one or more 4th layer points
  // in addition to a single inner layer hit
  
  //need to think about what to do with these sqrt(n)'s
  //also should add sig_phi and sig_theta here
  if( m_sig_phi   < 1E-6 ) m_sig_phi   = m_points.at(0)->sigPhi()/sqrt(n_avg);
  if( m_sig_theta < 1E-6 ) m_sig_theta = m_points.at(0)->sigPhi()/sqrt(n_avg);
  
  // correct phi of the cluster if we've drifted outside of 0 - 2PI
  if( m_phi > 2*PI ) m_phi -= 2*PI;
  if( m_phi < 0 ) m_phi += 2*PI;
}

void 
DBCALCluster::toStrings( vector< pair < string, string > > &items) const {

  AddString(items, "r", "%5.2f", m_rho * sin( m_theta ) );
  AddString(items, "phi", "%5.2f", m_phi );
//  AddString(items, "dphi", "%5.2f", m_sig_phi );
  AddString(items, "z", "%5.2f", m_rho * cos( m_theta ) + m_z_target_center );
  AddString(items, "theta", "%5.2f", m_theta);
//  AddString(items, "dtheta", "%5.2f", m_sig_theta);
  AddString(items, "t", "%5.2f", m_t );
  AddString(items, "E", "%5.2f", m_E );
  AddString(items, "N_cell", "%i", m_points.size() );
}


void
DBCALCluster::clear(){
 
  m_E = 0;
  m_E_points = 0; 
  m_t = 0;
  m_sig_t = 0;
  
  m_theta = 0;
  m_sig_theta = 0;

  m_phi = 0;
  m_sig_phi = 0;
  
  m_rho = 0;
  m_sig_rho = 0;
  
}

