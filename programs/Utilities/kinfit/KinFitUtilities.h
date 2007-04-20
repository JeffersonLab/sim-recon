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
double *momentum2tracking(const TLorentzVector &__p4, const string __reconstruction, const string __trackingConversion);
double *tracking2momentum(const double *__currentTrackingParameters, const double *__origTrackingParameters, 
    const string __reconstruction, const string __trackingConversion);
double **constraintDerivs(const double *__currentTrackingParameters, const double *__origTrackingParameters, 
    const string __reconstruction, const string __trackingConversion);

#endif /* _KinFitUtilities_H */


