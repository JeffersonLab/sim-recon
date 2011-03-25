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
    
    // set these the same for now -- need calibration correction
    shower->E = shower->E_raw;

    _data.push_back( shower );
  }
  
  return NOERROR;
}
