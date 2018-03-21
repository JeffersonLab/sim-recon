/*
 *  DBCALPoint.cc
 *
 *  Created by Matthew Shepherd on 3/13/11.
 *
 */

#include <iostream>

using namespace std;

#include "BCAL/DBCALHit.h"
#include "BCAL/DBCALPoint.h"
#include "BCAL/DBCALGeometry.h"

#include "units.h"

DBCALPoint::DBCALPoint(const DBCALUnifiedHit& hit1, const DBCALUnifiedHit& hit2, double z_target_center, 
		double attenuation_length, double c_effective, double track_p0, double track_p1, 
		double track_p2, const DBCALGeometry *locGeom) : m_BCALGeom(locGeom)
{
  
  // this is a problem -- both hits are on the same end...
  // this constructor is not equipped for this!

  assert( hit1.end != hit2.end );

  // check to be sure both hits are in the same cell
  
  int cellId = m_BCALGeom->cellId( hit1.module, hit1.layer, hit1.sector );
  assert( cellId == 
             m_BCALGeom->cellId( hit2.module, hit2.layer, hit2.sector ) );

  // save typing
  
  float fibLen = m_BCALGeom->GetBCAL_length();

  // figure out which hit is upstream and which is downstream
  // (downstream means farthest from the target)

  const DBCALUnifiedHit& upHit = 
     ( hit1.end == m_BCALGeom->kUpstream ? hit1 : hit2 );
  const DBCALUnifiedHit& downHit = 
     ( hit1.end == m_BCALGeom->kDownstream ? hit1 : hit2 );
  
  double tUp = upHit.t;
  double tDown = downHit.t;
  m_t_US = tUp;
  m_t_DS = tDown;

  // get the position with respect to the beginning of the module
  // the parameters were extracted from quadratic fits in histograms of z_track = f(tUp - tDown) 
  // we are thus defining z_point so that it matches the z-coordinate of the track (at the middle of each layer)
  // the p0, p1 and p2 tags are the usual names for the parameters of a quadratic fit in ROOT
  m_zGlobal = track_p0 + track_p1 * ( tUp - tDown ) + track_p2 * ( tUp - tDown ) * ( tUp - tDown ); 
  
  // get the position with respect to the center of the module -- positive
  // z in the downstream direction
  m_zLocal = m_zGlobal - m_BCALGeom->GetBCAL_center(); 

  // set the z position relative to the center of the target
  m_z = m_zGlobal - z_target_center;

/* the old method -- kept now for reference. Will be deleted
 *
  // get the position with respect to the center of the module -- positive
  // z in the downstream direction
  m_zLocal = 0.5 * c_effective * ( tUp - tDown ); 

  // set the z position relative to the center of the target
  m_z = m_zLocal + m_BCALGeom->GetBCAL_center() - z_target_center;
*/

  //At this point m_z may be unphysical, i.e. it may be outside the BCAL.
  //For the time being, this is okay. Forcing the z-position inside the
  //BCAL at this point will bias the clustering procedure:
  //1) Consider a cluster with a true position of z=400 cm (where the
  //downstream end of the BCAL is at 407 cm). Say we have three points, one each
  //at z=380 cm, 400 cm, 420 cm. If we force the latter point inside the BCAL
  //the average of the z-positions will be biased toward the upstream direction.
  //2) A point reconstructed at "410 cm" should not be treated equivalently to
  //a point reconstructed at "430 cm", though these would both be truncated to
  //407 cm. The latter should be more likely to
  //match a cluster at 390 cm.
  //The z-position will be forced inside the BCAL *after* clustering
  //and averaging.
  
  // compute the arrival time of the energy at the cell
  m_t = 0.5 * ( tUp + tDown - fibLen / c_effective );
  
  // now compute attentuation factors for each end based on distance
  // the light must travel
  // Make sure not to correct the energy by a distance longer than the length
  // of the BCAL.

  float dUp = 0.5 * fibLen + m_zLocal;
  float dDown = 0.5 * fibLen - m_zLocal;
  if (dUp>fibLen)   dUp=fibLen;
  if (dUp<0)        dUp=0;
  if (dDown>fibLen) dDown=fibLen;
  if (dDown<0)      dDown=0;
  float attUp = exp( -dUp / attenuation_length );
  float attDown = exp( -dDown / attenuation_length );
 
  // use these to correct the energy
  m_E_US =  ( upHit.E / attUp );
  m_E_DS =  ( downHit.E / attDown );
  m_E =  ( m_E_US + m_E_DS ) / 2;
  
  m_r = m_BCALGeom->r( cellId );
  //for a uniform distribution of width a, the RMS is a/sqrt(12)
  m_sig_r = m_BCALGeom->rSize( cellId )/sqrt(12.0);
  
  m_phi = m_BCALGeom->phi( cellId );
  m_sig_phi = m_BCALGeom->phiSize( cellId )/sqrt(12.0);
  
  //make a rough guess of sigma_z for now
  //if we have TDC info
  if (upHit.has_TDC_hit && downHit.has_TDC_hit) {
    //We have TDC hits at both ends, timing is more precise
    //Set z resolution to 1.1 cm/sqrt(E) as in NIM paper

    m_sig_z = 1.1 / sqrt(m_E);
  } else {
    //If we don't have TDC info at both ends then timing is less precise.
    //A reasonable value might be 4 ns/sqrt(12)/sqrt(2)*c_eff=14 cm
    //Although the ADC timing resolution will actually be better than
    //4 ns/sqrt(12) due to FPGA algorithm.
    //For now just set the value as large as needed.

    m_sig_z = 10.0;
  }

  m_module = hit1.module;
  m_layer = hit1.layer;
  m_sector = hit1.sector;
  
  // recast in terms of spherical coordinates
  convertCylindricalToSpherical();
}

float
DBCALPoint::tInnerRadius() const {
 
  // the path length in the module
  
  float modulePath = m_rho - m_BCALGeom->GetBCAL_inner_rad() / sin( m_theta );
  
  // retard the time by that distance divided by the speed of light

  return m_t - modulePath / ( 30 * k_cm / k_nsec );
}

void
DBCALPoint::convertCylindricalToSpherical(){
  
  m_rho = sqrt( m_r * m_r + m_z * m_z );
  m_theta = fabs( atan2( m_r, m_z ) );
  
  float d_rho_d_r = m_r / sqrt( m_r * m_r + m_z * m_z );
  float d_rho_d_z = m_z / sqrt( m_r * m_r + m_z * m_z );
  
  m_sig_rho = sqrt( m_sig_r * m_sig_r * d_rho_d_r * d_rho_d_r +
                    m_sig_z * m_sig_z * d_rho_d_z * d_rho_d_z );

  float d_theta_d_r =  m_z / ( m_r * m_r + m_z * m_z );
  float d_theta_d_z = -m_r / ( m_r * m_r + m_z * m_z );
  
  m_sig_theta = sqrt( m_sig_r * m_sig_r * d_theta_d_r * d_theta_d_r +
                      m_sig_z * m_sig_z * d_theta_d_z * d_theta_d_z );
  
}


// Using const_cast is generally not a good idea, but we
// really want these to be const functions since, in practice,
// they don't change the location of the point.
// Clearly we can't make m_phi mutable -- that would be worse.

void
DBCALPoint::add2Pi() const {

  const_cast< DBCALPoint* >( this )->m_phi += 2*M_PI;
}

void
DBCALPoint::sub2Pi() const {
  
  const_cast< DBCALPoint* >( this )->m_phi -= 2*M_PI;
}
