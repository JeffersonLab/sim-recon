// Base class for dealing with the Lorentz effect
#ifndef _DLorentzDeflections_
#define _DLorentzDeflections_

#include <JANA/jerror.h>
/* The folowing are for interpreting grid of Lorentz deflection data */
#define PACKAGE_Z_POINTS 10
#define LORENTZ_X_POINTS 21
#define LORENTZ_Z_POINTS (4*PACKAGE_Z_POINTS)

class DLorentzDeflections{
 public:
  
  DLorentzDeflections(){};
  virtual ~DLorentzDeflections(){};
  jerror_t GetLorentzCorrectionParameters(double x,double y,double z,
				double &tanz, double &tanr) const;
  double GetLorentzCorrection(double x,double y,double z,double alpha,
			      double dx) const;
  
 protected:
  // Variables for implementing lorentz effect (deflections of avalanche 
  // position due to the magnetic field).
  double lorentz_x[LORENTZ_X_POINTS];
  double lorentz_z[LORENTZ_Z_POINTS];
  double lorentz_nx[LORENTZ_X_POINTS][LORENTZ_Z_POINTS];
  double lorentz_nz[LORENTZ_X_POINTS][LORENTZ_Z_POINTS];
  
};

#endif // _DLorentzDeflections_
