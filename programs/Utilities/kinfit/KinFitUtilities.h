// KinFit class header file. -*- C++ -*-
/** @file kinfit/KinFit.h
 * @brief KinFit class defintion file.
 */
#ifndef _KinFitUtilities_H
#define _KinFitUtilities_H
// System Headers:
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
// ROOT Headers:
#include "TMatrixD.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector3.h"
// Local Headers
using namespace std;
//_____________________________________________________________________________

int GetSectorFromP4(const TLorentzVector &__p4);
double AlphaTrack(const TLorentzVector &__p4);
double LambdaTrack(const TLorentzVector &__p4);
double PhiTrack(const TLorentzVector &__p4) ;
double *momentum2tracking(const TLorentzVector &__p4, int __flag);
double *tracking2momentum(const double *__currentTrackingParameters, const double *__origTrackingParameters, int __flag);
double **constraintDerivs(const double *__currentTrackingParameters, const double *__origTrackingParameters, int __flag);

/*
// Convert from 4-momenta to tracking parameters
double *momentum2tracking(const TLorentzVector &__p4, int __flag)
{
  double trackingParameters[3];
  if(__flag == 1) // CLAS
  {
    trackingParameters[0] = __p4.P(); // |p|
    trackingParameters[1] = LambdaTrack(__p4); // tracking angle lambda
    trackingParameters[2] = PhiTrack(__p4); // tracking angle phi
  }
  return trackingParameters;
  //alpha[i] = AlphaTrack(_p4in[i]); // this angle doesn't change
}

/// Get CLAS sector number from 4-momentum
int GetSectorFromP4(const TLorentzVector &__p4){

  int sector = 0;
  double pi = 3.14159;
  double phi_lab = __p4.Phi();
  double phi = (180./pi)*phi_lab;

  if(std::abs(phi) <= 30.) sector = 1;
  else if(phi > 0.){
    if(phi <= 90.) sector = 2;
    else if(phi <= 150) sector = 3;
    else sector = 4;
  }
  else {
    // phi < 0
    if(std::abs(phi) <= 90.) sector = 6;
    else if(std::abs(phi) <= 150.) sector = 5;
    else sector = 4;
  }
  return sector;
}

/// Calculates the tracking quantity \f$ \alpha = \frac{\pi}{3}(sector - 1) \f$
double AlphaTrack(const TLorentzVector &__p4){

  int sector = 0;
  double pi = 3.14159;
  double phi_lab = __p4.Phi();
  double phi = (180./pi)*phi_lab;

  if(std::abs(phi) <= 30.) sector = 1;
  else if(phi > 0.){
    if(phi <= 90.) sector = 2;
    else if(phi <= 150) sector = 3;
    else sector = 4;
  }
  else {
    // phi < 0
    if(std::abs(phi) <= 90.) sector = 6;
    else if(std::abs(phi) <= 150.) sector = 5;
    else sector = 4;
  }
  return (pi/3.)*(sector - 1);
}
//_____________________________________________________________________________

/// Calculates the tracking angle \f$ \lambda \f$.
double LambdaTrack(const TLorentzVector &__p4){

  double lambda;
  double p_mag = __p4.P();
  double x = __p4.X()/p_mag,y = __p4.Y()/p_mag;

  double alpha = AlphaTrack(__p4);
  lambda = asin(cos(alpha)*y - sin(alpha)*x);
  return lambda;
}
//_____________________________________________________________________________

/// Calculates the tracking angle \f$ \phi \f$.
double PhiTrack(const TLorentzVector &__p4) {

  double phi;
  double z = __p4.Pz()/__p4.P(); // normalized z_lab
  double lambda = LambdaTrack(__p4);

  phi = acos(z/cos(lambda));
  return phi;
}
//_____________________________________________________________________________
*/
#endif /* _KinFitUtilities_H */


