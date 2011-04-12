/*
 *  DBCALShower_factory.cc
 *
 *  Created by Matthew Shepherd on 3/24/11.
 *
 */

#include "DBCALShower_factory.h"
#include "DBCALCluster.h"

#include "units.h"

DBCALShower_factory::DBCALShower_factory(){
  
  m_zTarget = 65*k_cm;

  if( ! DBCALGeometry::summingOn() ) {
  
    // these are energy calibration parameters -- no summing of cells
    
    m_scaleZ_p0 =  0.876376;
    m_scaleZ_p1 =  0.00150273;
    m_scaleZ_p2 =  -6.57424e-06;
    m_scaleZ_p3 =  7.79246e-09;
    
    m_nonlinZ_p0 =  0.0305786;
    m_nonlinZ_p1 =  -0.00014641;
    m_nonlinZ_p2 =  3.26065e-07;    
    m_nonlinZ_p3 =  0;
  }
  else{
    
    // these are energy calibration parameters -- 2.2.4.2 summing
    
    m_scaleZ_p0 =  0.785604;
    m_scaleZ_p1 =  0.0023703;
    m_scaleZ_p2 =  -9.47468e-06;
    m_scaleZ_p3 =  1.06317e-08;
    
    m_nonlinZ_p0 =  0.0468635;
    m_nonlinZ_p1 =  -9.69989e-05;
    m_nonlinZ_p2 =  0;    
    m_nonlinZ_p3 =  0;
  }
}

jerror_t
DBCALShower_factory::evnt( JEventLoop *loop, int eventnumber ){
 
  vector< const DBCALCluster* > clusters;
  loop->Get( clusters );
  
  // loop through and fill the shower structure from the cluster
  // right now just a simple 1 to 1 correspondence with 
  // an overall energy correction
  
  for( vector< const DBCALCluster* >::const_iterator clItr = clusters.begin();
       clItr != clusters.end();
      ++clItr ){
   
    float cosTh = cos( (**clItr).theta() );
    float sinTh = sin( (**clItr).theta() );
    float cosPhi = cos( (**clItr).phi() );
    float sinPhi = sin( (**clItr).phi() );
    float rho = (**clItr).rho();
    
    DBCALShower* shower = new DBCALShower();
    
    shower->E_raw = (**clItr).E();
    shower->x = rho * sinTh * cosPhi;
    shower->y = rho * sinTh * sinPhi;
    shower->z = rho * cosTh + m_zTarget;
    shower->t = (**clItr).t();
    shower->N_cell = (**clItr).nCells();
    
    float dx_drho = sinTh * cosPhi;
    float dy_drho = sinTh * sinPhi;
    float dz_drho = cosTh;
    float dx_dth  = rho * cosTh * cosPhi;
    float dy_dth  = rho * cosTh * sinPhi;
    float dz_dth  = -rho * sinTh;
    float dx_dphi = -rho * sinTh * sinPhi;
    float dy_dphi = rho * sinTh * cosPhi;
    float dz_dphi = 0;
    
    float drho = (**clItr).sigRho();
    float dphi = (**clItr).sigPhi();
    float dth = (**clItr).sigTheta();
    
    shower->xErr = sqrt( drho * drho * dx_drho * dx_drho + 
                         dphi * dphi * dx_dphi * dx_dphi +
                         dth * dth * dx_dth * dx_dth );

    shower->yErr = sqrt( drho * drho * dy_drho * dy_drho + 
                         dphi * dphi * dy_dphi * dy_dphi +
                         dth * dth * dy_dth * dy_dth );

    shower->zErr = sqrt( drho * drho * dz_drho * dz_drho + 
                         dphi * dphi * dz_dphi * dz_dphi +
                         dth * dth * dz_dth * dz_dth );
    
    shower->tErr = (**clItr).sigT();
    
    // calibrate energy:
    // Energy calibration has a z dependence -- the
    // calibration comes from fitting E_rec / E_gen to scale * E_gen^nonlin
    // for slices of z.  These fit parameters (scale and nonlin) are then plotted 
    // as a function of z and fit.
    
    // center of target should be a geometry lookup
    float zTarget = 65*k_cm;
    float r = sqrt( shower->x * shower->x + shower->y * shower->y );
    
    float zEntry = ( shower->z - zTarget ) * ( DBCALGeometry::BCALINNERRAD / r );
    
    float scale = m_scaleZ_p0  + m_scaleZ_p1*zEntry + 
    m_scaleZ_p2*(zEntry*zEntry) + m_scaleZ_p3*(zEntry*zEntry*zEntry);
    float nonlin = m_nonlinZ_p0  + m_nonlinZ_p1*zEntry + 
    m_nonlinZ_p2*(zEntry*zEntry) + m_nonlinZ_p3*(zEntry*zEntry*zEntry);
    
    shower->E = pow( (shower->E_raw ) / scale, 1 / ( 1 + nonlin ) );
    
    _data.push_back( shower );
  }
  
  return NOERROR;
}
