#include "KinFitUtilities.h"
using namespace std;
//_____________________________________________________________________________

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
// Convert from 4-momenta to tracking parameters
double *momentum2tracking(const TLorentzVector &__p4, const string __reconstruction, const string __trackingConversion)
{
  double *trackingParameters = new double[6]; 
  for(int i = 0; i < 6; ++i) 
  {
    trackingParameters[i] = 0.0;
  }

  if(__reconstruction=="DC" && __trackingConversion=="EASY") // momentum = tracking parameters
  {
    trackingParameters[0] = __p4.X(); // x
    trackingParameters[1] = __p4.Y(); // y
    trackingParameters[2] = __p4.Z(); // z
  }
  else if(__reconstruction=="DC" && __trackingConversion=="CLAS") // CLAS Drift Chambers
  {
    trackingParameters[0] = __p4.P(); // |p|
    trackingParameters[1] = LambdaTrack(__p4); // tracking angle lambda
    trackingParameters[2] = PhiTrack(__p4); // tracking angle phi
    trackingParameters[4] = AlphaTrack(__p4); // this angle doesn't change
  }
  else
  {
    cerr << "Not set up for this reconstruction with this tracking converstion!!!!" << endl;
    cerr << "Fit values are meaningless!!!" << endl;
    exit(-1);
  }
  trackingParameters[3] = __p4.M(); // Store the mass for later use
  return trackingParameters;

}

//_____________________________________________________________________________
// Convert from tracking parameters to 4-momenta
double *tracking2momentum(const double *__currentTrackingParameters, const double *__origTrackingParameters, const string __reconstruction, const string __trackingConversion)
{
  double *p3return = new double[3];

  if(__reconstruction=="DC" && __trackingConversion=="EASY") // momentum = tracking parameters
  {
    p3return[0] = __currentTrackingParameters[0]; // x
    p3return[1] = __currentTrackingParameters[1]; // y
    p3return[2] = __currentTrackingParameters[2]; // z
  }
  else if(__reconstruction=="DC" && __trackingConversion=="CLAS") // CLAS Drift Chambers
  {
    double p = __currentTrackingParameters[0];
    double lambda = __currentTrackingParameters[1];
    double phi = __currentTrackingParameters[2];

    double alpha = __origTrackingParameters[4]; // Calculated from the beginning

    p3return[0] = p*(cos(lambda)*sin(phi)*cos(alpha) - sin(lambda)*sin(alpha));
    p3return[1] = p*(cos(lambda)*sin(phi)*sin(alpha) + sin(lambda)*cos(alpha));
    p3return[2] = p*cos(lambda)*cos(phi);
  }
  else
  {
    cerr << "Not set up for this reconstruction with this tracking converstion!!!!" << endl;
    cerr << "Fit values are meaningless!!!" << endl;
    exit(-1);
  }
  return p3return;
}

//_____________________________________________________________________________
// Get the derivs of the constraint eqn wrt tracking paramaters
double **constraintDerivs(const double *__currentTrackingParameters, const double *__origTrackingParameters, const string __reconstruction, const string __trackingConversion)
{
  //double derivs[4][3];
  double **derivs;
  derivs = new double*[4]; 
  for(int i = 0; i < 4; ++i) 
           derivs[i] = new double[3]; 

  if(__reconstruction=="DC" && __trackingConversion=="EASY") // momentum = tracking parameters
  {
    double px = __currentTrackingParameters[0];
    double py = __currentTrackingParameters[1];
    double pz = __currentTrackingParameters[2];
    double mass = __origTrackingParameters[3]; // Calculated from the beginning
    double erg = sqrt( px*px + py*py + pz*pz + mass*mass);

    /* Energy constraint derivatives */
    derivs[0][0] = -1.0*px/erg;
    derivs[0][1] = -1.0*py/erg;
    derivs[0][2] = -1.0*pz/erg;

    /* Momentum constraint derivatives */
    derivs[1][0]=1;
    derivs[1][1]=0;
    derivs[1][2]=0;

    derivs[2][0]=0;
    derivs[2][1]=1;
    derivs[2][2]=0;

    derivs[3][0]=0;
    derivs[3][1]=0;
    derivs[3][2]=1;

  }
  else if(__reconstruction=="DC" && __trackingConversion=="CLAS") // CLAS Drift Chambers
  {
    double p = __currentTrackingParameters[0];
    double lambda = __currentTrackingParameters[1];
    double phi = __currentTrackingParameters[2];

    double mass  = __origTrackingParameters[3]; // Calculated from the beginning
    double alpha = __origTrackingParameters[4]; // Calculated from the beginning
    double erg = sqrt(p*p + mass*mass);
    //cerr << "in util: " << p << " " << erg << endl;

    /* derivatives of constraint eq's w/r to p */
    derivs[0][0] = -p/erg;
    derivs[1][0] = cos(lambda)*sin(phi)*cos(alpha) - sin(lambda)*sin(alpha);
    derivs[2][0] = cos(lambda)*sin(phi)*sin(alpha) + sin(lambda)*cos(alpha);
    derivs[3][0] = cos(lambda)*cos(phi);
    /* derivatives of constraint eq's w/r to lambda */
    derivs[0][1] = 0.;
    derivs[1][1] = p*(-sin(lambda)*sin(phi)*cos(alpha) - cos(lambda)*sin(alpha));
    derivs[2][1] = p*(-sin(lambda)*sin(phi)*sin(alpha) + cos(lambda)*cos(alpha));
    derivs[3][1] = p*(-sin(lambda)*cos(phi));
    /* derivatives of constraint eq's w/r to phi */
    derivs[0][2] = 0.;
    derivs[1][2] = p*cos(lambda)*cos(phi)*cos(alpha);
    derivs[2][2] = p*cos(lambda)*cos(phi)*sin(alpha);
    derivs[3][2] = p*cos(lambda)*(-sin(phi));

  }
  else
  {
    cerr << "Not set up for this reconstruction with this tracking converstion!!!!" << endl;
    cerr << "Fit values are meaningless!!!" << endl;
    exit(-1);
  }
  return derivs;
}


