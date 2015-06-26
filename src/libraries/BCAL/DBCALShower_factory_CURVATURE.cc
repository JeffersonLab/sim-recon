// $Id$
//
//    File: DBCALShower_factory_CURVATURE.cc
// Created: Fri Mar 27 10:57:45 CST 2015
// Creator: beattite (on Linux eos.phys.uregina.ca 2.6.32-504.12.2.el6.x86_64 x86_64)
//

// This factory takes showers from the IU shower factory and alters them based on shower
// curvature tables produced by AS.  These can be found in the CCDB under
// BCAL/curvature_central and BCAL/curvature_side.  They give central z-positions of
// energy depositions in the BCAL for each layer based on the reconstructed IU shower's
// energy and theta values.  The two tables are for 'central' sectors and 'side' sectors.
// A central sector (the one or two sectors closest to the centroid of the shower) exhibit
// different deposition ranges in z than a side sector.

// We first grab the IU showers and check for overlap in z.  We want to merge showers that
// are close in z, t, and phi that the IU code may have reconstructed as two separate showers.
// Next, we extract the energies and theta values for each of the IU (and merged) showers.
// Using the curvature tables, we find the appropriate z-bin position for each layer, then
// loop through the points in the event and put any point that falls within the z-bins for a
// shower into that shower.  The z-bin width is determined by the error on the z-positions of
// the depositions and on a set width factor.

// These new colletions of points are then averaged (energy-squared-weighted averages in x, y,
// z, and t) to get the shower values.

#include <iostream>
#include <iomanip>
using namespace std;

#include "DBCALShower_factory_CURVATURE.h"
#include "DBCALPoint.h"
#include "DBCALShower_factory_IU.h"

#include "DANA/DApplication.h"

#include "units.h"
using namespace jana;
#include "TMath.h"

//------------------
// init
//------------------
jerror_t DBCALShower_factory_CURVATURE::init(void)
{
  if( ! DBCALGeometry::summingOn() ) {
    // in libraries/PID/DNeutralShowerCandidate.h, there are some constants used to calculate the energy uncertainty. If you are updating these constants, you might want to update that also...

    // these are energy calibration parameters -- no summing of cells
    
    m_scaleZ_p0 =  0.950774;
    m_scaleZ_p1 =  0.000483979;
    m_scaleZ_p2 =  -2.08086e-06;
    m_scaleZ_p3 =  8.08534e-10;
    
    m_nonlinZ_p0 =  0.0152548;
    m_nonlinZ_p1 =  0;
    m_nonlinZ_p2 =  0;    
    m_nonlinZ_p3 =  0;
  }
  else{
    
    // these are energy calibration parameters -- 1.2.3.4 summing
    
    //last updated for svn revision 9233 
    m_scaleZ_p0 =  0.992437;
    m_scaleZ_p1 =  0.00039242;
    m_scaleZ_p2 =  -2.23135e-06;
    m_scaleZ_p3 =  1.40158e-09;
    
    m_nonlinZ_p0 =  -0.0147086;
    m_nonlinZ_p1 =  9.69207e-05;
    m_nonlinZ_p2 =  0;    
    m_nonlinZ_p3 =  0;

  }
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DBCALShower_factory_CURVATURE::brun(jana::JEventLoop *loop, int runnumber)
{
  DApplication* app = dynamic_cast<DApplication*>(loop->GetJApplication());
  DGeometry* geom = app->GetDGeometry(runnumber);
  geom->GetTargetZ(m_zTarget);

  vector< vector<double> > curvature_parameters;

  // Read in curvature parameters from two tables in the database.
  // We'll use two four-dimentional arrays (position and sigma).
  // The indices are: [sector type (central or side)][layer (from 0 to 4)][angle bin (from 0 to 11)][energy bin (from 0 to 31)]
  for (int ii = 0; ii < 2; ii++){
    if (ii == 0){
      if (loop->GetCalib("/BCAL/curvature_central", curvature_parameters)) jout << "Error loading /BCAL/curvature_central !" << endl;
    }
    else {
      if (loop->GetCalib("/BCAL/curvature_side", curvature_parameters)) jout << "Error loading /BCAL/curvature_side !" << endl;
    }
    for (int line = 0; line < 1536; line++){
      layer = curvature_parameters[line][2];
      angle = curvature_parameters[line][4];
      energy = curvature_parameters[line][1];
      dataposition = curvature_parameters[line][0];
      datasigma = curvature_parameters[line][3];
      if (angle == 115){ // angle bins are 5, 10, 20, 30, ..., 100, 110, 115.  [angle/10 - 1] won't work for the 115 bin, so we explicitly set these.
        position[ii][layer - 1][11][energy/50 - 1] = dataposition - 48; //table values measure from the end of the BCAL.  Points measure from the center of the target.
        sigma[ii][layer - 1][11][energy/50 - 1] = datasigma;
      }
      else {
        position[ii][layer - 1][angle/10 - 1][energy/50 - 1] = dataposition - 48; //table values measure from the end of the BCAL.  Points measure from the center of the target.
        sigma[ii][layer - 1][angle/10 - 1][energy/50 - 1] = datasigma;
      }
    }
    curvature_parameters.clear();
  }

  // Thresholds for merging showers and including points.
  PHITHRESHOLD = 4*0.06544792; // 4*2 column widths, in radians.
  ZTHRESHOLD = 35.*k_cm; // Range in z, in cm.
  TTHRESHOLD = 8.0*k_nsec; // Loose time range, in ns.
  ETHRESHOLD = 1.0*k_keV; // Points minimal energy, in keV.

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DBCALShower_factory_CURVATURE::evnt(JEventLoop *loop, int eventnumber)
{
  recon_showers_phi.clear();
  recon_showers_theta.clear();
  recon_showers_E.clear();
  recon_showers_t.clear();
  vector<const DBCALShower*> bcal_showers;
  vector<const DBCALPoint*> Points;
  vector<const DBCALPoint*> CurvaturePoints;

  loop->Get(bcal_showers,"IU");
  loop->Get(Points);
  
// Merge overlapping showers and calculate new shower parameters (phi, theta, E)
  while (bcal_showers.size() != 0){
    double recon_x = bcal_showers[0]->x;
    double recon_y = bcal_showers[0]->y;
    double recon_phi = atan2(recon_y,recon_x);
    double recon_theta = atan2(sqrt(bcal_showers[0]->x*bcal_showers[0]->x + bcal_showers[0]->y*bcal_showers[0]->y),bcal_showers[0]->z-m_zTarget);
    double recon_E = bcal_showers[0]->E;
    double recon_t = bcal_showers[0]->t/* - mcgentime*/;
    double total_E = bcal_showers[0]->E;
    for (int ii = 1; ii < (int)bcal_showers.size(); ii++){ // compare each shower in bcal_showers to the first to check for overlap.
      double delPhi = atan2(bcal_showers[0]->y,bcal_showers[0]->x) - atan2(bcal_showers[ii]->y,bcal_showers[ii]->x);
      if (delPhi > TMath::Pi()) delPhi -= 2*TMath::Pi();
      if (delPhi < -1*TMath::Pi()) delPhi += 2*TMath::Pi();
      double delZ = bcal_showers[0]->z - bcal_showers[ii]->z;
      double delT = bcal_showers[0]->t - bcal_showers[ii]->t;
      if (fabs(delPhi) < PHITHRESHOLD && fabs(delZ) < ZTHRESHOLD && fabs(delT) < TTHRESHOLD) overlap.push_back(ii); // Showers overlap!
    }
    if (overlap.size() != 0){
      recon_x *= recon_E;
      recon_y *= recon_E;
      recon_theta *= recon_E;
      recon_E *= recon_E;
      recon_t *= recon_E;
    }
    while (overlap.size() != 0){ // If we had an overlap, average the overlapping showers together.
      recon_x += bcal_showers[overlap.back()]->x*bcal_showers[overlap.back()]->E;
      recon_y += bcal_showers[overlap.back()]->y*bcal_showers[overlap.back()]->E; // Average over x and y rather than over phi to avoid boundary errors.
      recon_theta += atan2(sqrt(bcal_showers[overlap.back()]->x*bcal_showers[overlap.back()]->x + bcal_showers[overlap.back()]->y*bcal_showers[overlap.back()]->y),bcal_showers[overlap.back()]->z-m_zTarget)*bcal_showers[overlap.back()]->E;
      recon_E += bcal_showers[overlap.back()]->E*bcal_showers[overlap.back()]->E;
      recon_t += (bcal_showers[overlap.back()]->t/* - mcgentime*/)*bcal_showers[overlap.back()]->E;
      total_E += bcal_showers[overlap.back()]->E;

      bcal_showers.erase(bcal_showers.begin()+overlap.back()); // Erase the extra overlapping showers.
      overlap.pop_back();
      if (overlap.size() == 0){
        recon_x /= total_E;
        recon_y /= total_E;
        recon_phi = atan2(recon_y,recon_x);
        recon_theta /= total_E;
        recon_E /= total_E;
        recon_t /= total_E;
      }
    }
    // Output new shower parameters for use in creating new showers.
    recon_showers_phi.push_back(recon_phi);
    recon_showers_theta.push_back(recon_theta*180.0/TMath::Pi());
    recon_showers_E.push_back(recon_E);
    recon_showers_t.push_back(recon_t);
    bcal_showers.erase(bcal_showers.begin());
  }

//Create new showers.  Grab curvature parameters based on theta and energy of the old showers.
  for (int m = 0; m < (int)recon_showers_phi.size(); m++){ // Loop over our newly output shower parameters.
    // First, figure out which theta and energy bins we need from the tables.
    // We have to be a little careful with the first and last bins in both
    // theta and energy, since we want to interpolate between two adjacent
    // bin values.  The following carefully sets the two adjacent bins based
    // on the theta and energy values of the IU shower.
    // We do an interpolation between values in adjacent theta bins, then
    // a similar interpolation in adjacent energy bins so that our final z-position
    // will be a more-or-less smooth function of both theta and energy.
    temptheta = recon_showers_theta[m];
    tempenergy = recon_showers_E[m];
    k = l = 1;
    bin = 0;
    while (temptheta > 15){
      temptheta -= 10;
      bin++;
    }
    k = bin;
    if (k > 11) k = 11;
    bin = 0;
    while (tempenergy > 75){
      tempenergy -= 50;
      bin++;
    }
    l = bin;
    if (l > 31) l = 31;

    if (k == 11){
      k2 = k;
      k--;
    }
    else if ((recon_showers_theta[m] - 10*(k+1)) > 0 || k == 0) k2 = k + 1;
    else {
      k2 = k;
      k--;
    }
    if (l == 31){
      l2 = l;
      l--;
    }
    else if ((recon_showers_E[m] - 50*(l+1)) > 0 || l == 0) l2 = l + 1;
    else {
      l2 = l;
      l--;
    }

    // Now that we have the bins, we can use the position and sigma arrays to extrapolate or interpolate from the nearest two bins.
    // We have to again be a bit careful with the first and last angle bins, since they don't follow the regular 10 degree
    // interval pattern, sitting at 5 and 115 degrees respectively.  We find two interpolation offsets for the position and two
    // for the sigma for each layer.  These are then added to the positions and widths to give us our final z-bin position with error
    // for each layer and both central and side sectors (8 pairs in all).
    for (int ii = 0; ii < 2; ii++){
      for (int jj = 0; jj < 4; jj++){
        if (k == 0){ // Treat the k = 0 bin differently.  It lies at 15 degrees rather than the expected (from the pattern) 10 degrees.
          posoffset1 = (position[ii][jj][k][l] - position[ii][jj][k2][l])/(-5)*(recon_showers_theta[m] - 15);
          sigoffset1 = (sigma[ii][jj][k][l] - sigma[ii][jj][k2][l])/(-5)*(recon_showers_theta[m] - 15);
        }
        else if (k2 == 11){ // Treat the k = 11 bin differently.  It lies at 115 degrees rather than the expected (from the pattern) 120 degrees.
          posoffset1 = (position[ii][jj][k][l] - position[ii][jj][k2][l])/(-5)*(recon_showers_theta[m] - 110);
          sigoffset1 = (sigma[ii][jj][k][l] - sigma[ii][jj][k2][l])/(-5)*(recon_showers_theta[m] - 110);
        }
        else {
          posoffset1 = (position[ii][jj][k][l] - position[ii][jj][k2][l])/(-10)*(recon_showers_theta[m] - 10*((double)k+1));
          sigoffset1 = (sigma[ii][jj][k][l] - sigma[ii][jj][k2][l])/(-10)*(recon_showers_theta[m] - 10*((double)k+1));
        }
        posoffset2 = (position[ii][jj][k][l] - position[ii][jj][k][l2])/(-50)*(recon_showers_E[m] - 50*((double)l+1));
        sigoffset2 = (sigma[ii][jj][k][l] - sigma[ii][jj][k][l2])/(-50)*(recon_showers_E[m] - 50*((double)l+1));
        zbinposition[ii][jj] = position[ii][jj][k][l] + posoffset1 + posoffset2;
        zbinwidth[ii][jj] = (sigma[ii][jj][k][l] + sigoffset1 + sigoffset2);
      }
    }

    // We can now create the new shower.
    DBCALShower *shower = new DBCALShower;

    for (int ii = 0; ii < (int)Points.size(); ii++){
      double delPhi;
      if (Points[ii]->phi() > TMath::Pi()) delPhi = recon_showers_phi[m] - (Points[ii]->phi() - 2*TMath::Pi()); // Points->Phi() is measured from 0 to 2pi rather than from -pi to pi.
      else delPhi = recon_showers_phi[m] - Points[ii]->phi();
      if (delPhi > TMath::Pi()) delPhi -= 2*TMath::Pi();
      if (delPhi < -1*TMath::Pi()) delPhi += 2*TMath::Pi();
      if (fabs(delPhi) > PHITHRESHOLD) continue; // Don't include the point in this shower if it is far away in phi.
      if (fabs(delPhi) < 0.020452475) i = 0; // 5/8 of a column width in phi.  This is a 'central' shower (i = 1).
      else i = 1;
      j = Points[ii]->layer() - 1;

      // Below, we compare the point to the shower to see if it should be included in the shower.
      // This uses the thresholds defined at the top of this document.  However, we have to be careful
      // around the ends of the calorimeter.  Points that are reconstructed outside of the calorimeter
      // need to still be considered.  So, for theta >= 119 and for theta <= 13, we must include all
      // points past the end of the BCAL.  Also, for theta <= 17, we increase the z-bin width by a factor
      // of 1.7.  These angles and widths might have to be played with a bit to find the best values,
      // but for now, this seems to give good results (meaning, we catch most of the points in the
      // event and put them in reasonable showers).
      if (recon_showers_theta[m] >= 119.){ // Get all points past the edge of the calorimeter for showers near the edge.
        if (((fabs(Points[ii]->z() - zbinposition[i][j]) < (zbinwidth[i][j] + ZTHRESHOLD)) || (Points[ii]->z() < -48)) && (fabs(Points[ii]->t() - recon_showers_t[m]) < TTHRESHOLD) && (Points[ii]->E() > ETHRESHOLD)){
          CurvaturePoints.push_back(Points[ii]);
          overlap.push_back(ii);
        }
      }
      else if (recon_showers_theta[m] <= 13.){ // Get all points past the edge of the calorimeter for showers near the edge.
        if (((fabs(Points[ii]->z() - zbinposition[i][j]) < 1.7*(zbinwidth[i][j] + ZTHRESHOLD)) || (Points[ii]->z() > 342)) && (fabs(Points[ii]->t() - recon_showers_t[m]) < TTHRESHOLD) && (Points[ii]->E() > ETHRESHOLD)){
          CurvaturePoints.push_back(Points[ii]);
          overlap.push_back(ii);
        }
      }
      else if (recon_showers_theta[m] <= 17.){ // Use a larger z-bin-width for more forward-directed showers.
        if ((fabs(Points[ii]->z() - zbinposition[i][j]) < 1.7*(zbinwidth[i][j] + ZTHRESHOLD)) && (fabs(Points[ii]->t() - recon_showers_t[m]) < TTHRESHOLD) && (Points[ii]->E() > ETHRESHOLD)){
          CurvaturePoints.push_back(Points[ii]);
          overlap.push_back(ii);
        }
      }
      else if ((fabs(Points[ii]->z() - zbinposition[i][j]) < (zbinwidth[i][j] + ZTHRESHOLD)) && (fabs(Points[ii]->t() - recon_showers_t[m]) < TTHRESHOLD) && (Points[ii]->E() > ETHRESHOLD)){
        CurvaturePoints.push_back(Points[ii]);
        overlap.push_back(ii);
      }
    }

    // Average out the showers.  This is done in much the same way as the KLOE code was.
    // x, y, z, and t are averaged weighted by the energy of the point squared, and
    // sig_x, sig_y, sig_z, and sig_t are calculated alongside them.
    double x=0,y=0,z=0,t=0;
    int N_cell=0;
    double sig_x=0,sig_y=0,sig_z=0,sig_t=0;
    double total_E = 0; // Total Energy.
    double total_E2 = 0; // Total Energy squared.  Weight averages by E^2.
    double total_E2_squared = 0; // Used for calculating n_eff.
    double wt = 0; // Save writing: wt = CurvaturePoints[ii]->E().
    bool average_layer4 = false; // Keep track of showers with only layer four points.
    int n = CurvaturePoints.size();
    int n4 = 0; // Number of layer four points in the shower.

    for (int ii = 0; ii < (int)CurvaturePoints.size(); ii++){
      if (CurvaturePoints[ii]->layer() == 4) n4++;
    }
    if (n == n4) average_layer4 = true; // If all points are layer four points, include layer four points in the averages below.

    for (int ii = 0; ii < (int)CurvaturePoints.size(); ii++){
      if (CurvaturePoints[ii]->layer() != 4 || average_layer4) wt = CurvaturePoints[ii]->E();
      else wt = 0;
      total_E += CurvaturePoints[ii]->E(); // We still count energy from points in layer four.
      total_E2 += wt*wt;
      total_E2_squared += wt*wt*wt*wt;
      x += CurvaturePoints[ii]->r()*cos(CurvaturePoints[ii]->phi())*wt*wt;
      y += CurvaturePoints[ii]->r()*sin(CurvaturePoints[ii]->phi())*wt*wt;
      sig_x += CurvaturePoints[ii]->r()*CurvaturePoints[ii]->r()*cos(CurvaturePoints[ii]->phi())*cos(CurvaturePoints[ii]->phi())*wt*wt;
      sig_y += CurvaturePoints[ii]->r()*CurvaturePoints[ii]->r()*sin(CurvaturePoints[ii]->phi())*sin(CurvaturePoints[ii]->phi())*wt*wt;
      z += CurvaturePoints[ii]->z()*wt*wt;
      sig_z += CurvaturePoints[ii]->z()*CurvaturePoints[ii]->z()*wt*wt;
      t += CurvaturePoints[ii]->t()*wt*wt;
      sig_t += CurvaturePoints[ii]->t()*CurvaturePoints[ii]->t()*wt*wt;
      N_cell++;
      shower->AddAssociatedObject(CurvaturePoints[ii]);
    }
    while (overlap.size() != 0){
      Points.erase(Points.begin()+overlap.back());
      overlap.pop_back();
    }
    CurvaturePoints.clear();

    // the variance of the mean of a weighted distribution is s^2/n_eff, where s^2 is the variance of the sample and n_eff is as calculated below.
    double n_eff = total_E2*total_E2/total_E2_squared;

    x /= total_E2;
    sig_x /= total_E2;
    sig_x = sqrt(sig_x - x*x)/sqrt(n_eff);
    y /= total_E2;
    sig_y /= total_E2;
    sig_y = sqrt(sig_y - y*y)/sqrt(n_eff);
    z /= total_E2;
    sig_z /= total_E2;
    sig_z = sqrt(sig_z - z*z)/sqrt(n_eff);
    t /= total_E2;
    sig_t /= total_E2;
    sig_t = sqrt(sig_t - t*t)/sqrt(n_eff);

    // Force showers to be inside the BCal's z-coordinate range.
    double bcal_down = DBCALGeometry::GLOBAL_CENTER + DBCALGeometry::BCALFIBERLENGTH/2.0 - m_zTarget;
    double bcal_up = DBCALGeometry::GLOBAL_CENTER - DBCALGeometry::BCALFIBERLENGTH/2.0 - m_zTarget;
    if (z > bcal_down) z = bcal_down;
    if (z < bcal_up) z = bcal_up;

    shower->E_raw = total_E;
    shower->x = x;
    shower->y = y;
    shower->z = z + m_zTarget;
    shower->t = t;
    shower->N_cell = N_cell;
    shower->xErr = sig_x;
    shower->yErr = sig_y;
    shower->zErr = sig_z;
    shower->tErr = sig_t;
      
    // calibrate energy:
    // Energy calibration has a z dependence -- the
    // calibration comes from fitting E_rec / E_gen to scale * E_gen^nonlin
    // for slices of z.  These fit parameters (scale and nonlin) are then plotted 
    // as a function of z and fit.
    float r = sqrt( shower->x * shower->x + shower->y * shower->y );
    float zEntry = ( shower->z - m_zTarget ) * ( DBCALGeometry::BCALINNERRAD / r );
    float scale = m_scaleZ_p0  + m_scaleZ_p1*zEntry + m_scaleZ_p2*(zEntry*zEntry) + m_scaleZ_p3*(zEntry*zEntry*zEntry);
    float nonlin = m_nonlinZ_p0  + m_nonlinZ_p1*zEntry + m_nonlinZ_p2*(zEntry*zEntry) + m_nonlinZ_p3*(zEntry*zEntry*zEntry);

    shower->E = pow( (shower->E_raw ) / scale, 1 / ( 1 + nonlin ) );

    //copy xyz errors into covariance matrix
    shower->xyzCovariance.ResizeTo(3,3);
    shower->xyzCovariance[0][0] = shower->xErr*shower->xErr;
    shower->xyzCovariance[1][1] = shower->yErr*shower->yErr;
    shower->xyzCovariance[2][2] = shower->zErr*shower->zErr;

    _data.push_back(shower); 
  }

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DBCALShower_factory_CURVATURE::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DBCALShower_factory_CURVATURE::fini(void)
{
	return NOERROR;
}

