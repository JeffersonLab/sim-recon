/*
 *  DBCALShower_factory_IU.cc
 *  (formerly DBCALShower_factory.cc)
 *
 *  Created by Matthew Shepherd on 3/24/11.
 *
 */

#include "DBCALShower_factory_IU.h"
#include "DBCALCluster.h"

#include "DANA/DApplication.h"

#include "units.h"

DBCALShower_factory_IU::DBCALShower_factory_IU(){
}

jerror_t DBCALShower_factory_IU::brun(JEventLoop *loop, int32_t runnumber) {
  DApplication* app = dynamic_cast<DApplication*>(loop->GetJApplication());
  DGeometry* geom = app->GetDGeometry(runnumber);
  geom->GetTargetZ(m_zTarget);

  return NOERROR;
}

jerror_t
DBCALShower_factory_IU::evnt( JEventLoop *loop, uint64_t eventnumber ){
 
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
    shower->E_preshower = (**clItr).E_preshower();
    shower->x = rho * sinTh * cosPhi;
    shower->y = rho * sinTh * sinPhi;
    shower->z = rho * cosTh + m_zTarget;

    //DBCALCluster::t() returns the time at the inner radius
    //so we need to make an adjustment so that the shower t is the time at
    //the shower location (x,y,z)
    double t = (**clItr).t();
    double inner_rad = DBCALGeometry::GetBCAL_inner_rad();
    double dist_in_BCAL = rho - inner_rad/sinTh;
    t = t + dist_in_BCAL/(30*k_cm/k_nsec);
    shower->t = t;

    shower->N_cell = (**clItr).nCells();
    
    //create matrices to rotate errors from cylindrical coordinates to Cartesian coordinates
    float dx_drho = sinTh * cosPhi;
    float dy_drho = sinTh * sinPhi;
    float dz_drho = cosTh;
    float dx_dth  = rho * cosTh * cosPhi;
    float dy_dth  = rho * cosTh * sinPhi;
    float dz_dth  = -rho * sinTh;
    float dx_dphi = -rho * sinTh * sinPhi;
    float dy_dphi = rho * sinTh * cosPhi;
    float dz_dphi = 0;

    DMatrix rotation(3,3);
    DMatrix rotationT(3,3);

    rotation[0][0] = dx_drho;
    rotation[0][1] = dy_drho;
    rotation[0][2] = dz_drho;
    rotation[1][0] = dx_dth;
    rotation[1][1] = dy_dth;
    rotation[1][2] = dz_dth;
    rotation[2][0] = dx_dphi;
    rotation[2][1] = dy_dphi;
    rotation[2][2] = dz_dphi;

    rotationT.Transpose(rotation);
    
    //create covariance matrix in cylindrical coordinates
    //for now assume that these measurements are independent (uncorrelated)
    //will need to think harder to add correlations
    float drho = (**clItr).sigRho();
    float dphi = (**clItr).sigPhi();
    float dth = (**clItr).sigTheta();
    DMatrix errors(3,3);
    errors[0][0] = drho*drho;
    errors[1][1] = dth*dth;
    errors[2][2] = dphi*dphi;

    //do the rotation
    shower->xyzCovariance.ResizeTo(3,3);
    shower->xyzCovariance = rotationT*errors*rotation;
    
    //fill (redundant) x/y/zErr members
    shower->xErr = sqrt(shower->xyzCovariance[0][0]);
    shower->yErr = sqrt(shower->xyzCovariance[1][1]);
    shower->zErr = sqrt(shower->xyzCovariance[2][2]);
    
    shower->tErr = (**clItr).sigT();
    
    shower->E = shower->E_raw;

    shower->AddAssociatedObject(*clItr);
    
    _data.push_back( shower );
  }
  
  return NOERROR;
}
